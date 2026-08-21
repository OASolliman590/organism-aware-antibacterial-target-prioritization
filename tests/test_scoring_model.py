from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd

from pipeline.config import ProjectConfig, load_config
from pipeline.scoring_model import (
    calibrated_metrics_with_bootstrap,
    fit_calibrated_combiner,
    predict_calibrated,
    reliability_curve,
)


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["model"]["features"] = ["feature_2d", "feature_3d"]
    data["benchmark"]["bootstrap_n"] = 50
    return ProjectConfig(path=base.path, root=base.root, data=data)


def _frame(prefix: str, offset: float = 0.0) -> pd.DataFrame:
    rows = []
    for index in range(8):
        label = index % 2
        rows.append(
            {
                "query_id": f"{prefix}{index // 2}",
                "target_class": f"t{index}",
                "feature_2d": 0.8 + offset if label else 0.2 + offset,
                "feature_3d": 0.7 + offset if label else 0.3 + offset,
                "is_correct_target": label,
            }
        )
    return pd.DataFrame(rows)


def test_calibrated_combiner_uses_disjoint_sets_and_names_only_platt_probability() -> None:
    config = _config()
    model, manifest = fit_calibrated_combiner(
        _frame("train"), _frame("cal"), config=config
    )
    test = _frame("test")
    predictions = predict_calibrated(model, test)
    curve = reliability_curve(predictions, n_bins=5)
    metrics = calibrated_metrics_with_bootstrap(predictions, config=config)

    assert predictions.calibrated_probability.between(0, 1).all()
    assert not predictions.uncalibrated_logistic_score_is_probability.any()
    assert (predictions.calibrated_probability_method == "Platt").all()
    assert manifest["label_semantics"] == (
        "correct target-class candidate, not molecular inactivity"
    )
    assert len(curve) == 5
    assert set(metrics.metric) == {"auroc", "brier_score"}
    assert metrics[["ci_lower_95", "ci_upper_95"]].notna().all().all()


def test_missing_features_are_excluded_not_imputed() -> None:
    config = _config()
    model, _ = fit_calibrated_combiner(
        _frame("train"), _frame("cal"), config=config
    )
    test = _frame("test")
    test.loc[0, "feature_3d"] = np.nan

    predictions = predict_calibrated(model, test)

    assert np.isnan(predictions.loc[0, "calibrated_probability"])
    assert predictions.loc[0, "model_evidence_status"] == (
        "excluded_missing_feature_no_imputation"
    )


def test_query_overlap_is_rejected_as_leakage() -> None:
    config = _config()
    try:
        fit_calibrated_combiner(_frame("same"), _frame("same"), config=config)
    except ValueError as error:
        assert "leakage" in str(error).lower()
    else:
        raise AssertionError("overlapping query IDs must be rejected")
