"""Render v3 figures only from completed run tables; never invent empty panels."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

try:
    from pipeline.config import load_config
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import load_config


def _benchmark_figure(frame: pd.DataFrame, path: Path) -> bool:
    import matplotlib.pyplot as plt

    selected_metrics = ["bedroc_alpha_80_5", "ef_1pct", "mrr", "coverage"]
    frame = frame[
        frame.metric.isin(selected_metrics) & frame.estimate.notna()
    ].copy()
    if frame.empty:
        return False
    splits = list(dict.fromkeys(frame.split_type.astype(str)))
    figure, axes = plt.subplots(
        len(splits), 1, figsize=(10, max(4, 3.5 * len(splits))), squeeze=False
    )
    modes = ["2d_only", "3d_only", "fusion"]
    width = 0.24
    for axis, split in zip(axes[:, 0], splits):
        subset = frame[frame.split_type == split]
        positions = np.arange(len(selected_metrics))
        for mode_index, mode in enumerate(modes):
            mode_rows = subset[subset.score_mode == mode].set_index("metric")
            estimates = [mode_rows.estimate.get(metric, np.nan) for metric in selected_metrics]
            lower = [
                max(0.0, estimate - mode_rows.ci_lower_95.get(metric, estimate))
                if np.isfinite(estimate)
                else 0.0
                for metric, estimate in zip(selected_metrics, estimates)
            ]
            upper = [
                max(0.0, mode_rows.ci_upper_95.get(metric, estimate) - estimate)
                if np.isfinite(estimate)
                else 0.0
                for metric, estimate in zip(selected_metrics, estimates)
            ]
            axis.bar(
                positions + (mode_index - 1) * width,
                estimates,
                width,
                yerr=np.vstack([lower, upper]),
                capsize=2,
                label=mode,
            )
        axis.set_xticks(positions, selected_metrics)
        axis.set_ylabel("metric value")
        axis.set_title(f"{split} split")
        axis.legend()
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle("V3 retrieval comparison (95% query-bootstrap CI)")
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return True


def _ad_figure(frame: pd.DataFrame, path: Path) -> bool:
    import matplotlib.pyplot as plt

    if frame.empty or "applicability_domain_flag" not in frame:
        return False
    counts = frame["applicability_domain_flag"].value_counts().sort_index()
    if counts.empty:
        return False
    figure, axis = plt.subplots(figsize=(8, 4.5))
    counts.plot.bar(ax=axis)
    axis.set_ylabel("prediction rows")
    axis.set_title("V3 applicability-domain coverage")
    axis.tick_params(axis="x", rotation=25)
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return True


def _disagreement_figure(frame: pd.DataFrame, path: Path) -> bool:
    import matplotlib.pyplot as plt

    if frame.empty:
        return False
    selected = frame.nlargest(20, "absolute_rank_shift").copy()
    labels = selected.query_id.astype(str) + " | " + selected.target_class.astype(str)
    colors = np.where(selected.rank_shift_2d_to_v3 > 0, "#2a9d8f", "#e76f51")
    figure, axis = plt.subplots(figsize=(9, max(4, 0.35 * len(selected))))
    positions = np.arange(len(selected))
    axis.barh(positions, selected.rank_shift_2d_to_v3, color=colors)
    axis.set_yticks(positions, labels)
    axis.axvline(0, color="black", linewidth=0.8)
    axis.set_xlabel("rank shift (positive = promoted by v3)")
    axis.set_title("Material 2D-vs-v3 target-rank disagreements")
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return True


def _reliability_figure(frame: pd.DataFrame, path: Path) -> bool:
    import matplotlib.pyplot as plt

    observed = frame[(frame.n > 0) & frame.mean_calibrated_probability.notna()]
    if observed.empty:
        return False
    figure, axis = plt.subplots(figsize=(5.5, 5.5))
    axis.plot([0, 1], [0, 1], linestyle="--", color="gray", label="ideal")
    axis.plot(
        observed.mean_calibrated_probability,
        observed.observed_fraction_correct,
        marker="o",
        label="held out",
    )
    axis.set(xlim=(0, 1), ylim=(0, 1), xlabel="calibrated probability", ylabel="observed fraction correct")
    axis.set_title("V3 reliability diagram")
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return True


def generate_v3_figures(results_dir: Path) -> pd.DataFrame:
    figure_dir = results_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    specifications = [
        (
            "benchmark_mode_comparison",
            results_dir / "benchmark_mode_comparison_v3.csv",
            figure_dir / "benchmark_mode_comparison_v3.png",
            _benchmark_figure,
        ),
        (
            "applicability_domain",
            results_dir / "benchmark_target_scores_by_split_v3.csv",
            figure_dir / "applicability_domain_v3.png",
            _ad_figure,
        ),
        (
            "rank_disagreement",
            results_dir / "chemical_evidence_disagreements_v3.csv",
            figure_dir / "chemical_evidence_disagreements_v3.png",
            _disagreement_figure,
        ),
        (
            "calibration_reliability",
            results_dir / "scoring_model_reliability_v3.csv",
            figure_dir / "scoring_model_reliability_v3.png",
            _reliability_figure,
        ),
    ]
    statuses = []
    for figure_name, source, output, renderer in specifications:
        if not source.is_file():
            status = "unavailable_source_missing"
        else:
            frame = pd.read_csv(source)
            status = "created" if renderer(frame, output) else "unavailable_no_evaluable_rows"
        statuses.append(
            {
                "figure": figure_name,
                "status": status,
                "source": str(source),
                "output": str(output) if status == "created" else None,
            }
        )
    return pd.DataFrame(statuses)


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    status = generate_v3_figures(results_dir)
    status.to_csv(results_dir / "v3_figure_status.csv", index=False)
    print(status.to_string(index=False))


if __name__ == "__main__":
    main()
