"""Weight/layer perturbation and reference-bootstrap rank stability for v3."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

try:
    from pipeline.config import ProjectConfig, load_config
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig, load_config


REFERENCE_COMPONENT_COLUMNS = {
    "ecfp4_max": "ecfp4_similarity",
    "maccs_max": "maccs_similarity",
    "usrcat_max": "usrcat_similarity",
    "o3a_shape_tanimoto_max": "o3a_shape_tanimoto",
    "o3a_color_max": "o3a_color",
    "pharmacophore_2d_gobbi_sim_max": "pharmacophore_2d_gobbi_similarity",
    "pharmacophore_3d_sim_max": "pharmacophore_3d_similarity",
}

BIOLOGY_COMPONENT_COLUMNS = {
    "organism_scope": "organism_scope_score",
    "clinical": "clinical_priority_score",
    "essentiality": "essentiality_score",
    "accessibility": "cellular_access_score",
    "resistance": "resistance_relevance_score",
    "card_context": "card_resistance_context_score",
}

OVERALL_AFFINE_LAYERS = {
    "transfer": "species_transfer_score",
    "pocket": "pocket_evidence_score",
    "biology": "biological_priority_score",
}


def perturbation_scenarios(components: list[str], factors: list[float]):
    baseline = {component: 1.0 for component in components}
    scenarios: list[tuple[str, str, dict[str, float]]] = []
    for component in components:
        for factor in factors:
            weights = baseline.copy()
            weights[component] = float(factor)
            scenarios.append(
                (
                    f"{component}_x{factor:g}",
                    "single_component_weight_perturbation",
                    weights,
                )
            )
        weights = baseline.copy()
        weights[component] = 0.0
        scenarios.append(
            (f"leave_out_{component}", "leave_one_component_out", weights)
        )
    layer_groups = {
        "2d": ["ecfp4_max", "maccs_max"],
        "shape": ["usrcat_max", "o3a_shape_tanimoto_max", "o3a_color_max"],
        "pharmacophore": [
            "pharmacophore_2d_gobbi_sim_max",
            "pharmacophore_3d_sim_max",
        ],
    }
    for layer, members in layer_groups.items():
        present = [component for component in members if component in baseline]
        if not present:
            continue
        weights = baseline.copy()
        for component in present:
            weights[component] = 0.0
        scenarios.append(
            (f"leave_out_{layer}_layer", "leave_one_layer_out", weights)
        )
    return baseline, scenarios


def weighted_rrf_scores(
    evidence: pd.DataFrame,
    weights: dict[str, float],
    *,
    reciprocal_rank_constant: float,
) -> pd.Series:
    """Weighted RRF with fixed missing-evidence denominator; never an imputation."""

    missing = sorted(set(weights) - set(evidence.columns))
    if missing:
        raise ValueError(f"Sensitivity evidence is missing components: {missing}")
    positive_weight = sum(float(weight) for weight in weights.values() if weight > 0)
    if positive_weight <= 0:
        raise ValueError("At least one sensitivity weight must be positive")
    contributions = []
    available = pd.Series(False, index=evidence.index)
    for component, weight in weights.items():
        if weight <= 0:
            continue
        values = pd.to_numeric(evidence[component], errors="coerce").where(
            lambda series: np.isfinite(series)
        )
        ranks = values.groupby(evidence["query_id"], sort=False).rank(
            method="average", ascending=False, na_option="keep"
        )
        contributions.append(float(weight) / (reciprocal_rank_constant + ranks))
        available |= values.notna()
    raw = pd.concat(contributions, axis=1).sum(axis=1, min_count=1)
    ideal = positive_weight / (reciprocal_rank_constant + 1.0)
    return (raw / ideal).clip(0.0, 1.0).where(available)


def _rank_order(frame: pd.DataFrame, score_column: str, top_k: int) -> list[str]:
    available = frame[np.isfinite(pd.to_numeric(frame[score_column], errors="coerce"))]
    return (
        available.sort_values(
            [score_column, "target_class"],
            ascending=[False, True],
            kind="mergesort",
        )["target_class"]
        .astype(str)
        .head(top_k)
        .tolist()
    )


def rank_biased_overlap(left: list[str], right: list[str], persistence: float) -> float:
    """Finite extrapolated RBO for two deterministic ranked lists."""

    depth = min(len(left), len(right))
    if depth == 0:
        return np.nan
    overlap = 0
    weighted = 0.0
    left_seen: set[str] = set()
    right_seen: set[str] = set()
    for index in range(depth):
        left_seen.add(left[index])
        right_seen.add(right[index])
        overlap = len(left_seen & right_seen)
        weighted += (overlap / (index + 1)) * persistence**index
    return float(
        (1.0 - persistence) * weighted
        + (overlap / depth) * persistence**depth
    )


def scenario_rank_stability(
    evidence: pd.DataFrame,
    baseline_weights: dict[str, float],
    scenario_weights: dict[str, float],
    *,
    reciprocal_rank_constant: float,
    top_k: int,
    rbo_persistence: float,
) -> tuple[float, float, int]:
    compared = evidence[["query_id", "target_class", *baseline_weights]].copy()
    compared["baseline"] = weighted_rrf_scores(
        compared,
        baseline_weights,
        reciprocal_rank_constant=reciprocal_rank_constant,
    )
    compared["variant"] = weighted_rrf_scores(
        compared,
        scenario_weights,
        reciprocal_rank_constant=reciprocal_rank_constant,
    )
    taus = []
    rbos = []
    for _, group in compared.groupby("query_id", sort=True):
        complete = group.dropna(subset=["baseline", "variant"]).copy()
        if len(complete) < 2:
            continue
        baseline_order = _rank_order(complete, "baseline", len(complete))
        variant_order = _rank_order(complete, "variant", len(complete))
        baseline_position = {target: index for index, target in enumerate(baseline_order)}
        variant_position = {target: index for index, target in enumerate(variant_order)}
        common = sorted(set(baseline_position) & set(variant_position))
        tau = kendalltau(
            [baseline_position[target] for target in common],
            [variant_position[target] for target in common],
        ).statistic
        if np.isfinite(tau):
            taus.append(float(tau))
        rbos.append(
            rank_biased_overlap(
                baseline_order[:top_k], variant_order[:top_k], rbo_persistence
            )
        )
    rbos = [value for value in rbos if np.isfinite(value)]
    return (
        float(np.mean(taus)) if taus else np.nan,
        float(np.mean(rbos)) if rbos else np.nan,
        max(len(taus), len(rbos)),
    )


def final_ranking_scenarios(config: ProjectConfig):
    """Enumerate every published coefficient that affects the final heuristic rank."""

    factors = [float(value) for value in config.value("sensitivity.weight_factors")]
    baseline = {
        "specificity": {
            key: float(value)
            for key, value in config.value("v2_scoring.specificity").items()
        },
        "specificity_enabled": True,
        "reference_quality": {
            key: float(value)
            for key, value in config.value("v2_scoring.reference_quality").items()
        },
        "reference_quality_enabled": True,
        "biology_weights": {
            key: float(value)
            for key, value in config.value("v2_scoring.biology_weights").items()
        },
        "overall_factors": {
            layer: {
                "base": float(config.value(f"v2_scoring.overall_factors.{layer}.base")),
                "weight": float(
                    config.value(f"v2_scoring.overall_factors.{layer}.weight")
                ),
            }
            for layer in OVERALL_AFFINE_LAYERS
        },
        "anti_target_penalty": float(
            config.value("v2_scoring.overall_factors.anti_target_penalty")
        ),
    }
    scenarios: list[tuple[str, str, dict[str, Any]]] = []

    for component in BIOLOGY_COMPONENT_COLUMNS:
        for factor in factors:
            variant = json.loads(json.dumps(baseline))
            variant["biology_weights"][component] *= factor
            scenarios.append(
                (
                    f"biology_weight.{component}_x{factor:g}",
                    "final_biology_weight_perturbation",
                    variant,
                )
            )
        variant = json.loads(json.dumps(baseline))
        variant["biology_weights"][component] = 0.0
        scenarios.append(
            (
                f"leave_out_biology_component.{component}",
                "leave_one_final_biology_component_out",
                variant,
            )
        )

    for layer in OVERALL_AFFINE_LAYERS:
        for coefficient in ("base", "weight"):
            for factor in factors:
                variant = json.loads(json.dumps(baseline))
                variant["overall_factors"][layer][coefficient] *= factor
                scenarios.append(
                    (
                        f"overall.{layer}.{coefficient}_x{factor:g}",
                        "final_overall_coefficient_perturbation",
                        variant,
                    )
                )
        variant = json.loads(json.dumps(baseline))
        variant["overall_factors"][layer] = {"base": 1.0, "weight": 0.0}
        scenarios.append(
            (
                f"leave_out_overall_layer.{layer}",
                "leave_one_final_overall_layer_out",
                variant,
            )
        )

    for factor in factors:
        variant = json.loads(json.dumps(baseline))
        variant["anti_target_penalty"] *= factor
        scenarios.append(
            (
                f"overall.anti_target_penalty_x{factor:g}",
                "final_overall_coefficient_perturbation",
                variant,
            )
        )
    variant = json.loads(json.dumps(baseline))
    variant["anti_target_penalty"] = 0.0
    scenarios.append(
        (
            "leave_out_overall_layer.anti_target",
            "leave_one_final_overall_layer_out",
            variant,
        )
    )

    for coefficient in ("baseline", "margin_weight", "margin_scale"):
        for factor in factors:
            variant = json.loads(json.dumps(baseline))
            variant["specificity"][coefficient] *= factor
            scenarios.append(
                (
                    f"specificity.{coefficient}_x{factor:g}",
                    "final_chemical_coefficient_perturbation",
                    variant,
                )
            )
    variant = json.loads(json.dumps(baseline))
    variant["specificity_enabled"] = False
    scenarios.append(
        (
            "leave_out_final_chemical_layer.specificity",
            "leave_one_final_chemical_layer_out",
            variant,
        )
    )

    for grade in sorted(baseline["reference_quality"]):
        for factor in factors:
            variant = json.loads(json.dumps(baseline))
            variant["reference_quality"][grade] *= factor
            scenarios.append(
                (
                    f"reference_quality.{grade}_x{factor:g}",
                    "final_chemical_coefficient_perturbation",
                    variant,
                )
            )
    variant = json.loads(json.dumps(baseline))
    variant["reference_quality_enabled"] = False
    scenarios.append(
        (
            "leave_out_final_chemical_layer.reference_quality",
            "leave_one_final_chemical_layer_out",
            variant,
        )
    )
    return baseline, scenarios


def recompute_final_priority(
    predictions: pd.DataFrame,
    parameters: dict[str, Any],
) -> pd.DataFrame:
    """Recompute the published heuristic from separate evidence fields."""

    required = {
        "organism",
        "query_id",
        "target_class",
        "chemical_evidence_score_v3",
        "target_specificity_margin",
        "reference_quality_grade",
        "species_transfer_score",
        "pocket_evidence_score",
        "anti_target_risk_score",
        *BIOLOGY_COMPONENT_COLUMNS.values(),
    }
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise ValueError(f"Final-ranking sensitivity is missing fields: {missing}")

    result = predictions.copy()
    numeric_columns = sorted(
        required
        - {"organism", "query_id", "target_class", "reference_quality_grade"}
    )
    numeric = result[numeric_columns].apply(pd.to_numeric, errors="coerce")
    evaluable = np.isfinite(numeric.to_numpy(dtype=float)).all(axis=1)
    if not evaluable.any():
        raise ValueError(
            "Final-ranking sensitivity has no rows with measured finite component "
            "fields; missing values are not imputed"
        )
    result = result.loc[evaluable].copy()
    numeric = numeric.loc[evaluable]
    for column in numeric_columns:
        result[column] = numeric[column]

    biology = pd.Series(0.0, index=result.index, dtype=float)
    for component, column in BIOLOGY_COMPONENT_COLUMNS.items():
        biology += float(parameters["biology_weights"][component]) * result[column]
    result["sensitivity_biological_priority_score"] = biology

    if parameters["specificity_enabled"]:
        specificity_parameters = parameters["specificity"]
        margin_scale = float(specificity_parameters["margin_scale"])
        if margin_scale <= 0:
            raise ValueError("Specificity margin_scale must remain positive")
        specificity = (
            float(specificity_parameters["baseline"])
            + float(specificity_parameters["margin_weight"])
            * result["target_specificity_margin"]
            / margin_scale
        ).clip(0.0, 1.0)
    else:
        specificity = pd.Series(1.0, index=result.index, dtype=float)
    if parameters["reference_quality_enabled"]:
        reference_quality = result["reference_quality_grade"].map(
            parameters["reference_quality"]
        )
        if reference_quality.isna().any():
            unknown = sorted(
                result.loc[reference_quality.isna(), "reference_quality_grade"]
                .astype(str)
                .unique()
            )
            raise ValueError(f"Unconfigured reference-quality grades: {unknown}")
        reference_quality = reference_quality.astype(float)
    else:
        reference_quality = pd.Series(1.0, index=result.index, dtype=float)
    chemical_quality = (
        result["chemical_evidence_score_v3"].astype(float)
        * specificity
        * reference_quality
    )
    result["sensitivity_target_specificity_score"] = specificity
    result["sensitivity_reference_quality_score"] = reference_quality
    result["sensitivity_chemical_hypothesis_score"] = chemical_quality

    overall = chemical_quality.copy()
    overall_parameters = parameters["overall_factors"]
    layer_values = {
        "transfer": result["species_transfer_score"],
        "pocket": result["pocket_evidence_score"],
        "biology": biology,
    }
    for layer, values in layer_values.items():
        coefficients = overall_parameters[layer]
        overall *= float(coefficients["base"]) + float(
            coefficients["weight"]
        ) * values
    overall *= 1.0 - float(parameters["anti_target_penalty"]) * result[
        "anti_target_risk_score"
    ]
    result["sensitivity_overall_priority_score"] = overall
    result["sensitivity_overall_priority_score_is_probability"] = False
    return result


def _group_rank_stability(
    baseline: pd.DataFrame,
    variant: pd.DataFrame,
    *,
    top_k: int,
    rbo_persistence: float,
) -> pd.DataFrame:
    joined = baseline[
        ["organism", "query_id", "target_class", "sensitivity_overall_priority_score"]
    ].merge(
        variant[
            [
                "organism",
                "query_id",
                "target_class",
                "sensitivity_overall_priority_score",
            ]
        ],
        on=["organism", "query_id", "target_class"],
        validate="one_to_one",
        suffixes=("_baseline", "_variant"),
    )
    rows = []
    for (organism, query_id), group in joined.groupby(
        ["organism", "query_id"], sort=True
    ):
        if len(group) < 2:
            continue
        baseline_order = _rank_order(
            group.rename(
                columns={
                    "sensitivity_overall_priority_score_baseline": "baseline"
                }
            ),
            "baseline",
            len(group),
        )
        variant_order = _rank_order(
            group.rename(
                columns={"sensitivity_overall_priority_score_variant": "variant"}
            ),
            "variant",
            len(group),
        )
        baseline_position = {
            target: index for index, target in enumerate(baseline_order)
        }
        variant_position = {target: index for index, target in enumerate(variant_order)}
        common = sorted(set(baseline_position) & set(variant_position))
        tau = kendalltau(
            [baseline_position[target] for target in common],
            [variant_position[target] for target in common],
        ).statistic
        rows.append(
            {
                "organism": organism,
                "query_id": query_id,
                "kendall_tau": float(tau) if np.isfinite(tau) else np.nan,
                "rbo": rank_biased_overlap(
                    baseline_order[:top_k], variant_order[:top_k], rbo_persistence
                ),
            }
        )
    return pd.DataFrame(rows)


def run_final_ranking_sensitivity(
    predictions: pd.DataFrame,
    config: ProjectConfig,
) -> pd.DataFrame:
    """Quantify sensitivity of every heuristic coefficient used in final ranks."""

    baseline_parameters, scenarios = final_ranking_scenarios(config)
    baseline = recompute_final_priority(predictions, baseline_parameters)
    if "chemical_hypothesis_score" in predictions:
        recorded_chemical = pd.to_numeric(
            predictions.loc[baseline.index, "chemical_hypothesis_score"],
            errors="coerce",
        )
        recalculated_chemical = baseline["sensitivity_chemical_hypothesis_score"]
        if not np.allclose(
            recorded_chemical, recalculated_chemical, rtol=1e-10, atol=1e-12
        ):
            raise ValueError(
                "Configured chemical-quality formula does not reproduce recorded "
                "baseline scores"
            )
    if "overall_priority_score" in predictions:
        recorded = pd.to_numeric(
            predictions.loc[baseline.index, "overall_priority_score"],
            errors="coerce",
        )
        recalculated = baseline["sensitivity_overall_priority_score"]
        if not np.allclose(recorded, recalculated, rtol=1e-10, atol=1e-12):
            raise ValueError(
                "Configured final-ranking formula does not reproduce recorded baseline scores"
            )

    top_k = int(config.value("sensitivity.top_k"))
    persistence = float(config.value("sensitivity.rbo_persistence"))
    bootstrap_n = int(config.value("sensitivity.bootstrap_n"))
    seed = int(config.value("seeds.bootstrap"))
    rows = []
    for scenario_name, scenario_type, parameters in scenarios:
        variant = recompute_final_priority(predictions, parameters)
        grouped = _group_rank_stability(
            baseline,
            variant,
            top_k=top_k,
            rbo_persistence=persistence,
        )
        tau_values = grouped["kendall_tau"].dropna().to_numpy(dtype=float)
        rbo_values = grouped["rbo"].dropna().to_numpy(dtype=float)
        tau_bootstrap = []
        rbo_bootstrap = []
        if len(grouped):
            rng = np.random.default_rng(seed)
            for _ in range(bootstrap_n):
                sampled = grouped.iloc[rng.integers(0, len(grouped), size=len(grouped))]
                sampled_tau = sampled["kendall_tau"].dropna()
                sampled_rbo = sampled["rbo"].dropna()
                if len(sampled_tau):
                    tau_bootstrap.append(float(sampled_tau.mean()))
                if len(sampled_rbo):
                    rbo_bootstrap.append(float(sampled_rbo.mean()))
        tau_ci = (
            np.quantile(tau_bootstrap, [0.025, 0.975]).tolist()
            if tau_bootstrap
            else [np.nan, np.nan]
        )
        rbo_ci = (
            np.quantile(rbo_bootstrap, [0.025, 0.975]).tolist()
            if rbo_bootstrap
            else [np.nan, np.nan]
        )
        rows.append(
            {
                "scenario": scenario_name,
                "scenario_type": scenario_type,
                "parameters_json": json.dumps(parameters, sort_keys=True),
                "mean_kendall_tau": (
                    float(tau_values.mean()) if len(tau_values) else np.nan
                ),
                "kendall_tau_ci_lower_95": tau_ci[0],
                "kendall_tau_ci_upper_95": tau_ci[1],
                "mean_rbo": float(rbo_values.mean()) if len(rbo_values) else np.nan,
                "rbo_ci_lower_95": rbo_ci[0],
                "rbo_ci_upper_95": rbo_ci[1],
                "n_ranked_lists_evaluable": int(len(grouped)),
                "n_prediction_rows_input": int(len(predictions)),
                "n_prediction_rows_evaluable": int(len(baseline)),
                "n_prediction_rows_excluded_missing": int(
                    len(predictions) - len(baseline)
                ),
                "bootstrap_n": bootstrap_n,
                "bootstrap_seed": seed,
                "bootstrap_unit": "organism_query_ranked_list",
                "score_is_probability": False,
            }
        )
    return pd.DataFrame(rows)


def _bootstrap_reference_aggregates(
    reference_evidence: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    for (query_id, target_class), group in reference_evidence.groupby(
        ["query_id", "target_class"], sort=True
    ):
        indices = rng.integers(0, len(group), size=len(group))
        sampled = group.iloc[indices]
        row = {"query_id": query_id, "target_class": target_class}
        for component, reference_column in REFERENCE_COMPONENT_COLUMNS.items():
            values = pd.to_numeric(sampled[reference_column], errors="coerce")
            values = values[np.isfinite(values)]
            row[component] = float(values.max()) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def run_sensitivity_analysis(
    scores: pd.DataFrame,
    reference_evidence: pd.DataFrame,
    config: ProjectConfig,
) -> pd.DataFrame:
    """Perturb all fusion weights/layers and bootstrap reference ligands."""

    components = list(config.value("fusion.components"))
    factors = [float(value) for value in config.value("sensitivity.weight_factors")]
    baseline, scenarios = perturbation_scenarios(components, factors)
    rrf_constant = float(config.value("fusion.reciprocal_rank_constant"))
    top_k = int(config.value("sensitivity.top_k"))
    persistence = float(config.value("sensitivity.rbo_persistence"))
    bootstrap_n = int(config.value("sensitivity.bootstrap_n"))
    seed = int(config.value("seeds.bootstrap"))
    rows = []
    for split_type in config.value("benchmark.splits"):
        split_scores = scores[scores.split_type == split_type]
        split_references = reference_evidence[
            reference_evidence.split_type == split_type
        ] if not reference_evidence.empty else reference_evidence
        point_values: dict[str, tuple[float, float, int]] = {}
        for scenario_name, _, weights in scenarios:
            point_values[scenario_name] = scenario_rank_stability(
                split_scores,
                baseline,
                weights,
                reciprocal_rank_constant=rrf_constant,
                top_k=top_k,
                rbo_persistence=persistence,
            ) if not split_scores.empty else (np.nan, np.nan, 0)
        bootstrap_values = {
            scenario_name: {"kendall": [], "rbo": []}
            for scenario_name, _, _ in scenarios
        }
        if not split_references.empty:
            rng = np.random.default_rng(seed)
            for _ in range(bootstrap_n):
                aggregated = _bootstrap_reference_aggregates(split_references, rng)
                for scenario_name, _, weights in scenarios:
                    tau, rbo, _ = scenario_rank_stability(
                        aggregated,
                        baseline,
                        weights,
                        reciprocal_rank_constant=rrf_constant,
                        top_k=top_k,
                        rbo_persistence=persistence,
                    )
                    if np.isfinite(tau):
                        bootstrap_values[scenario_name]["kendall"].append(tau)
                    if np.isfinite(rbo):
                        bootstrap_values[scenario_name]["rbo"].append(rbo)
        for scenario_name, scenario_type, weights in scenarios:
            tau, rbo, n_queries = point_values[scenario_name]
            tau_samples = bootstrap_values[scenario_name]["kendall"]
            rbo_samples = bootstrap_values[scenario_name]["rbo"]
            tau_ci = (
                np.quantile(tau_samples, [0.025, 0.975]).tolist()
                if tau_samples
                else [np.nan, np.nan]
            )
            rbo_ci = (
                np.quantile(rbo_samples, [0.025, 0.975]).tolist()
                if rbo_samples
                else [np.nan, np.nan]
            )
            rows.append(
                {
                    "split_type": split_type,
                    "scenario": scenario_name,
                    "scenario_type": scenario_type,
                    "weights_json": json.dumps(weights, sort_keys=True),
                    "mean_kendall_tau": tau,
                    "kendall_tau_ci_lower_95": tau_ci[0],
                    "kendall_tau_ci_upper_95": tau_ci[1],
                    "mean_rbo": rbo,
                    "rbo_ci_lower_95": rbo_ci[0],
                    "rbo_ci_upper_95": rbo_ci[1],
                    "n_queries_evaluable": n_queries,
                    "reference_bootstrap_n": bootstrap_n,
                    "reference_bootstrap_seed": seed,
                    "bootstrap_unit": "reference_ligand_within_query_target",
                    "bootstrap_status": (
                        "available"
                        if tau_samples or rbo_samples
                        else "unavailable_missing_reference_evidence"
                    ),
                    "score_is_probability": False,
                }
            )
    return pd.DataFrame(rows)


def plot_sensitivity_report(report: pd.DataFrame, output_path: Path) -> None:
    """Plot layer-ablation RBO with its reference-bootstrap interval."""

    import matplotlib.pyplot as plt

    selected = report[
        (report.scenario_type == "leave_one_layer_out")
        & report.mean_rbo.notna()
    ].copy()
    if selected.empty:
        return
    selected = selected.sort_values(["split_type", "scenario"], kind="mergesort")
    labels = (
        selected["split_type"].astype(str)
        + " | "
        + selected["scenario"].str.replace("leave_out_", "", regex=False)
    )
    lower = selected["mean_rbo"] - selected["rbo_ci_lower_95"]
    upper = selected["rbo_ci_upper_95"] - selected["mean_rbo"]
    positions = np.arange(len(selected))
    figure, axis = plt.subplots(figsize=(10, max(4, 0.45 * len(selected))))
    axis.errorbar(
        selected["mean_rbo"],
        positions,
        xerr=np.vstack([lower.clip(lower=0), upper.clip(lower=0)]),
        fmt="o",
        capsize=3,
    )
    axis.set_yticks(positions, labels)
    axis.set_xlim(0, 1.02)
    axis.set_xlabel("Rank-biased overlap with equal-weight fusion")
    axis.set_title("V3 layer-ablation rank stability (95% reference bootstrap CI)")
    axis.grid(axis="x", alpha=0.25)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def plot_final_ranking_sensitivity_report(
    report: pd.DataFrame, output_path: Path
) -> None:
    """Plot final-ranking layer ablations with grouped-bootstrap intervals."""

    import matplotlib.pyplot as plt

    selected = report[
        report.scenario_type.str.startswith("leave_one_final_", na=False)
        & report.mean_rbo.notna()
    ].copy()
    if selected.empty:
        return
    selected = selected.sort_values("scenario", kind="mergesort")
    lower = selected["mean_rbo"] - selected["rbo_ci_lower_95"]
    upper = selected["rbo_ci_upper_95"] - selected["mean_rbo"]
    positions = np.arange(len(selected))
    figure, axis = plt.subplots(figsize=(10, max(5, 0.34 * len(selected))))
    axis.errorbar(
        selected["mean_rbo"],
        positions,
        xerr=np.vstack([lower.clip(lower=0), upper.clip(lower=0)]),
        fmt="o",
        capsize=3,
    )
    axis.set_yticks(positions, selected["scenario"])
    axis.set_xlim(0, 1.02)
    axis.set_xlabel("Rank-biased overlap with the configured final ranking")
    axis.set_title("Final heuristic rank sensitivity (95% grouped-bootstrap CI)")
    axis.grid(axis="x", alpha=0.25)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    score_path = results_dir / "benchmark_target_scores_by_split_v3.csv"
    reference_path = results_dir / "benchmark_reference_evidence_by_split_v3.csv"
    if score_path.is_file():
        scores = pd.read_csv(score_path)
        references = pd.read_csv(reference_path) if reference_path.is_file() else pd.DataFrame()
        report = run_sensitivity_analysis(scores, references, config)
        plot_sensitivity_report(
            report, results_dir / "figures" / "sensitivity_rank_stability_v3.png"
        )
        fusion_status = "available"
        fusion_reason = (
            "point perturbations and reference bootstrap complete"
            if not references.empty
            else "point perturbations complete; reference bootstrap unavailable"
        )
    else:
        report = pd.DataFrame()
        fusion_status = "pending_missing_v3_benchmark_scores"
        fusion_reason = str(score_path)
    report.to_csv(results_dir / "sensitivity_analysis_v3.csv", index=False)

    prediction_path = results_dir / "open_target_predictions_by_organism_v3.csv"
    if prediction_path.is_file():
        predictions = pd.read_csv(prediction_path)
        final_report = run_final_ranking_sensitivity(predictions, config)
        plot_final_ranking_sensitivity_report(
            final_report,
            results_dir / "figures" / "final_ranking_sensitivity_v3.png",
        )
        final_status = "available"
        final_reason = (
            "all configured final-ranking coefficients perturbed with "
            "organism-query grouped bootstrap"
        )
    else:
        final_report = pd.DataFrame()
        final_status = "pending_missing_v3_organism_predictions"
        final_reason = str(prediction_path)
    final_report.to_csv(
        results_dir / "final_ranking_sensitivity_v3.csv", index=False
    )

    status = {
        "status": (
            "available"
            if fusion_status == "available" and final_status == "available"
            else "partial_or_pending"
        ),
        "fusion_sensitivity_status": fusion_status,
        "fusion_sensitivity_reason": fusion_reason,
        "fusion_scenarios": len(report),
        "final_ranking_sensitivity_status": final_status,
        "final_ranking_sensitivity_reason": final_reason,
        "final_ranking_scenarios": len(final_report),
    }
    pd.DataFrame([status]).to_csv(
        results_dir / "sensitivity_analysis_status_v3.csv", index=False
    )
    print(pd.DataFrame([status]).to_string(index=False))


if __name__ == "__main__":
    main()
