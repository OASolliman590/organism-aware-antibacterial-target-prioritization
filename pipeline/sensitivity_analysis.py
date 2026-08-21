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
    "pharmacophore_sim_max": "pharmacophore_similarity",
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
        "pharmacophore": ["pharmacophore_sim_max"],
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
        status = {
            "status": "available",
            "status_reason": (
                "point perturbations and reference bootstrap complete"
                if not references.empty
                else "point perturbations complete; reference bootstrap unavailable"
            ),
            "n_scenarios": len(report),
        }
    else:
        report = pd.DataFrame()
        status = {
            "status": "pending_missing_v3_benchmark_scores",
            "status_reason": str(score_path),
            "n_scenarios": 0,
        }
    report.to_csv(results_dir / "sensitivity_analysis_v3.csv", index=False)
    pd.DataFrame([status]).to_csv(
        results_dir / "sensitivity_analysis_status_v3.csv", index=False
    )
    print(pd.DataFrame([status]).to_string(index=False))


if __name__ == "__main__":
    main()
