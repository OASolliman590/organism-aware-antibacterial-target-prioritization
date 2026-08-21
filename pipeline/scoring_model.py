"""Leakage-checked regularised logistic combiner with Platt calibration.

Only the Platt output is named a probability. Rows with missing feature evidence
are excluded and counted; feature values are never imputed.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    from pipeline.config import ProjectConfig, load_config
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig, load_config


@dataclass(frozen=True)
class CalibratedCombiner:
    features: tuple[str, ...]
    base_model: Pipeline
    platt_model: LogisticRegression
    label_column: str
    train_query_ids: tuple[str, ...]
    calibration_query_ids: tuple[str, ...]


def _complete_rows(
    frame: pd.DataFrame, features: list[str], label_column: str
) -> tuple[pd.DataFrame, int]:
    required = {"query_id", label_column, *features}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Model frame is missing columns: {missing}")
    numeric = frame.copy()
    for column in [*features, label_column]:
        numeric[column] = pd.to_numeric(numeric[column], errors="coerce")
    complete = numeric[[*features, label_column]].apply(
        lambda column: np.isfinite(column)
    ).all(axis=1)
    complete &= numeric[label_column].isin([0, 1])
    return numeric.loc[complete].copy(), int((~complete).sum())


def _assert_disjoint_query_ids(*frames: pd.DataFrame) -> None:
    sets = [set(frame["query_id"].astype(str)) for frame in frames]
    for left in range(len(sets)):
        for right in range(left + 1, len(sets)):
            overlap = sets[left] & sets[right]
            if overlap:
                raise ValueError(
                    f"Query leakage across train/calibration/test: {sorted(overlap)[:5]}"
                )


def fit_calibrated_combiner(
    train: pd.DataFrame,
    calibration: pd.DataFrame,
    *,
    config: ProjectConfig,
    label_column: str = "is_correct_target",
) -> tuple[CalibratedCombiner, dict[str, Any]]:
    """Fit training-only scaling/logistic model and disjoint Platt calibrator."""

    _assert_disjoint_query_ids(train, calibration)
    features = list(config.value("model.features"))
    train_complete, train_excluded = _complete_rows(train, features, label_column)
    calibration_complete, calibration_excluded = _complete_rows(
        calibration, features, label_column
    )
    minimum = int(config.value("model.minimum_class_count"))
    for label, frame in [
        ("training", train_complete),
        ("calibration", calibration_complete),
    ]:
        counts = frame[label_column].value_counts().to_dict()
        if counts.get(0, 0) < minimum or counts.get(1, 0) < minimum:
            raise ValueError(
                f"{label} requires at least {minimum} rows from each class; got {counts}"
            )

    base = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "logistic",
                LogisticRegression(
                    C=float(config.value("model.regularization_c")),
                    class_weight="balanced",
                    random_state=int(config.value("seeds.model")),
                    solver="lbfgs",
                    max_iter=int(config.value("model.max_iterations")),
                ),
            ),
        ]
    )
    base.fit(train_complete[features], train_complete[label_column].astype(int))
    calibration_logits = base.decision_function(calibration_complete[features]).reshape(
        -1, 1
    )
    platt = LogisticRegression(
        C=float(config.value("model.platt_regularization_c")),
        random_state=int(config.value("seeds.model")),
        solver="lbfgs",
        max_iter=int(config.value("model.max_iterations")),
    )
    platt.fit(calibration_logits, calibration_complete[label_column].astype(int))
    model = CalibratedCombiner(
        tuple(features),
        base,
        platt,
        label_column,
        tuple(sorted(train_complete.query_id.astype(str).unique())),
        tuple(sorted(calibration_complete.query_id.astype(str).unique())),
    )
    manifest = {
        "model_type": "standardized_regularized_logistic_plus_platt",
        "features": features,
        "label_semantics": "correct target-class candidate, not molecular inactivity",
        "n_train_rows": len(train_complete),
        "n_calibration_rows": len(calibration_complete),
        "n_train_rows_excluded_missing": train_excluded,
        "n_calibration_rows_excluded_missing": calibration_excluded,
        "train_query_ids": list(model.train_query_ids),
        "calibration_query_ids": list(model.calibration_query_ids),
        "random_seed": int(config.value("seeds.model")),
    }
    return model, manifest


def predict_calibrated(
    model: CalibratedCombiner, test: pd.DataFrame
) -> pd.DataFrame:
    """Predict without imputation; only the Platt output is a probability."""

    _assert_disjoint_query_ids(
        pd.DataFrame({"query_id": model.train_query_ids}),
        pd.DataFrame({"query_id": model.calibration_query_ids}),
        test,
    )
    result = test.copy()
    for feature in model.features:
        if feature not in result:
            raise ValueError(f"Test frame is missing model feature: {feature}")
        result[feature] = pd.to_numeric(result[feature], errors="coerce")
    complete = result[list(model.features)].apply(
        lambda column: np.isfinite(column)
    ).all(axis=1)
    result["model_evidence_status"] = np.where(
        complete, "complete", "excluded_missing_feature_no_imputation"
    )
    result["uncalibrated_logistic_score"] = np.nan
    result["uncalibrated_logistic_score_is_probability"] = False
    result["calibrated_probability"] = np.nan
    result["calibrated_probability_method"] = "Platt"
    if complete.any():
        logits = model.base_model.decision_function(
            result.loc[complete, list(model.features)]
        )
        result.loc[complete, "uncalibrated_logistic_score"] = 1.0 / (
            1.0 + np.exp(-logits)
        )
        result.loc[complete, "calibrated_probability"] = model.platt_model.predict_proba(
            np.asarray(logits).reshape(-1, 1)
        )[:, 1]
    return result


def reliability_curve(
    predictions: pd.DataFrame,
    *,
    label_column: str = "is_correct_target",
    n_bins: int = 5,
) -> pd.DataFrame:
    """Return fixed-width reliability bins, retaining empty bins as gaps."""

    labels = pd.to_numeric(predictions[label_column], errors="coerce")
    probabilities = pd.to_numeric(
        predictions["calibrated_probability"], errors="coerce"
    )
    valid = labels.isin([0, 1]) & np.isfinite(probabilities)
    labels = labels[valid]
    probabilities = probabilities[valid]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins = np.minimum(np.searchsorted(edges, probabilities, side="right") - 1, n_bins - 1)
    rows = []
    for index in range(n_bins):
        selected = bins == index
        rows.append(
            {
                "bin_index": index,
                "bin_lower": float(edges[index]),
                "bin_upper": float(edges[index + 1]),
                "n": int(selected.sum()),
                "mean_calibrated_probability": (
                    float(probabilities[selected].mean()) if selected.any() else np.nan
                ),
                "observed_fraction_correct": (
                    float(labels[selected].mean()) if selected.any() else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def calibrated_metrics_with_bootstrap(
    predictions: pd.DataFrame,
    *,
    config: ProjectConfig,
    label_column: str = "is_correct_target",
) -> pd.DataFrame:
    """Held-out AUROC and Brier score with deterministic query bootstrap CIs."""

    required = {"query_id", label_column, "calibrated_probability"}
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise ValueError(f"Calibrated predictions are missing: {missing}")
    frame = predictions.copy()
    frame[label_column] = pd.to_numeric(frame[label_column], errors="coerce")
    frame["calibrated_probability"] = pd.to_numeric(
        frame["calibrated_probability"], errors="coerce"
    )
    valid = frame[label_column].isin([0, 1]) & np.isfinite(
        frame["calibrated_probability"]
    )
    frame = frame[valid]

    def calculate(data: pd.DataFrame, metric: str) -> float:
        labels = data[label_column].astype(int)
        probability = data["calibrated_probability"].astype(float)
        if data.empty:
            return np.nan
        if metric == "brier_score":
            return float(brier_score_loss(labels, probability))
        return (
            float(roc_auc_score(labels, probability))
            if labels.nunique() == 2
            else np.nan
        )

    query_ids = sorted(frame.query_id.astype(str).unique())
    rng = np.random.default_rng(int(config.value("seeds.bootstrap")))
    samples = [
        rng.choice(query_ids, size=len(query_ids), replace=True).tolist()
        for _ in range(int(config.value("benchmark.bootstrap_n")))
    ] if query_ids else []
    rows = []
    for metric in ["auroc", "brier_score"]:
        estimate = calculate(frame, metric)
        bootstrapped = []
        for selection in samples:
            sampled = pd.concat(
                [frame[frame.query_id.astype(str) == query_id] for query_id in selection],
                ignore_index=True,
            )
            value = calculate(sampled, metric)
            if math.isfinite(value):
                bootstrapped.append(value)
        lower, upper = (
            np.quantile(bootstrapped, [0.025, 0.975]).tolist()
            if bootstrapped
            else (np.nan, np.nan)
        )
        rows.append(
            {
                "metric": metric,
                "estimate": estimate,
                "ci_lower_95": lower,
                "ci_upper_95": upper,
                "n_queries": len(query_ids),
                "n_rows": len(frame),
                "bootstrap_unit": "query_id",
                "bootstrap_n": int(config.value("benchmark.bootstrap_n")),
                "status": "available" if math.isfinite(estimate) else "unavailable",
            }
        )
    return pd.DataFrame(rows)


def assess_training_readiness(config: ProjectConfig) -> dict[str, Any]:
    """Require actual split roles and all configured benchmark splits before fitting."""

    scores_path = config.path_for("results") / "benchmark_target_scores_by_split_v3.csv"
    provenance_path = config.path_for("results") / "benchmark_split_provenance_v3.csv"
    reasons = []
    if not scores_path.is_file():
        reasons.append("benchmark_scores_missing")
    if not provenance_path.is_file():
        reasons.append("split_provenance_missing")
    if scores_path.is_file():
        columns = pd.read_csv(scores_path, nrows=0).columns
        if "model_split_role" not in columns:
            reasons.append("disjoint_train_calibration_test_roles_missing")
    if provenance_path.is_file():
        provenance = pd.read_csv(provenance_path)
        for split_type in config.value("benchmark.splits"):
            available = provenance[
                (provenance.split_type == split_type)
                & (provenance.status == "available")
            ]
            if available.empty:
                reasons.append(f"{split_type}_split_unavailable")
    return {
        "status": "ready" if not reasons else "pending_sparse_or_incomplete_labels",
        "status_reason": ";".join(reasons) if reasons else "all requirements satisfied",
        "scores_path": str(scores_path),
        "split_provenance_path": str(provenance_path),
        "calibration_method": str(config.value("model.calibration_method")),
        "only_calibrated_output_named_probability": True,
    }


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    results_dir.mkdir(parents=True, exist_ok=True)
    status = assess_training_readiness(config)
    pd.DataFrame([status]).to_csv(results_dir / "scoring_model_status_v3.csv", index=False)
    pd.DataFrame(
        columns=["metric", "estimate", "ci_lower_95", "ci_upper_95", "status"]
    ).to_csv(results_dir / "scoring_model_heldout_metrics_v3.csv", index=False)
    pd.DataFrame(
        columns=[
            "bin_index",
            "bin_lower",
            "bin_upper",
            "n",
            "mean_calibrated_probability",
            "observed_fraction_correct",
        ]
    ).to_csv(results_dir / "scoring_model_reliability_v3.csv", index=False)
    print(pd.DataFrame([status]).to_string(index=False))


if __name__ == "__main__":
    main()
