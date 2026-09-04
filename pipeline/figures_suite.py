"""Compound-facing visualization suite for the organism-aware pipeline.

Every panel is rendered only from a completed run table. A panel that has no
source file, no evaluable rows, or an unavailable optional dependency reports
that reason in ``results/figure_suite_status.csv`` instead of emitting an empty
or invented figure. The status table is the audit trail: a missing PNG must
always be explainable from it.

This module covers the *per-compound* view of a run. It complements
``pipeline/v3_figures.py`` (benchmark/calibration diagnostics) and
``pipeline/v2_figures.py`` (aggregate organism views) rather than replacing
them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

try:  # pragma: no cover - import shim for direct module execution
    from pipeline.config import load_config
    from pipeline import figure_style as style
    from pipeline.compound_manifest import CompoundManifest, load_from_config
except ModuleNotFoundError:  # pragma: no cover
    from config import load_config
    import figure_style as style
    from compound_manifest import CompoundManifest, load_from_config


STATUS_CREATED = "created"
STATUS_SOURCE_MISSING = "unavailable_source_missing"
STATUS_NO_ROWS = "unavailable_no_evaluable_rows"
STATUS_DEPENDENCY = "unavailable_dependency_unavailable"
STATUS_NOT_APPLICABLE = "not_applicable_for_single_organism"

PREDICTIONS_BY_ORGANISM = "v2_open_target_predictions_by_organism.csv"
PREDICTIONS_V3 = "open_target_predictions_by_organism_v3.csv"
DISAGREEMENTS_V3 = "chemical_evidence_disagreements_v3.csv"
SENSITIVITY_V3 = "final_ranking_sensitivity_v3.csv"
UNCERTAINTY_V2 = "v2_uncertainty_private.csv"
BENCHMARK_SUMMARY_V2 = "benchmark_v2_summary.csv"


@dataclass(frozen=True)
class SuiteContext:
    """Inputs a renderer may need beyond the run tables themselves."""

    private_compounds: Path | None = None
    benchmark_compounds: Path | None = None
    umap_seed: int = 17
    figure_sample_seed: int = 7
    top_target_classes: int = 18
    max_disagreements_per_compound: int = 6
    # Design intent from the private manifest: which organism each compound was
    # prepared against. Used to annotate and focus figures, never to score.
    manifest: CompoundManifest | None = None
    # Set when this suite is rendered for a single organism.
    organism: str | None = None
    # When true, a single-organism suite covers only the compounds the manifest
    # assigned to that organism, rather than the whole series scored against it.
    restrict_to_assigned: bool = False

    def scoped_compounds(self) -> list[str] | None:
        """Compounds this suite is limited to, or None for no compound limit."""

        if not self.restrict_to_assigned or self.organism is None:
            return None
        if self.manifest is None:
            return None
        return self.manifest.compounds_for(self.organism)

    def label(self, compound: str) -> str:
        """Compound label, marked when it was prepared against this organism.

        The marker distinguishes assigned compounds from the rest, so it is
        redundant once the suite is already scoped to the assigned set.
        """

        if self.scoped_compounds() is not None:
            return compound
        if self.manifest is not None and self.manifest.is_assigned(
            compound, self.organism
        ):
            return f"{compound}*"
        return compound

    def labels(self, compounds) -> list[str]:
        return [self.label(str(compound)) for compound in compounds]

    def scope_label(self) -> str:
        """Title suffix naming the organism a suite is filtered to, if any."""

        return f" — {self.organism}" if self.organism else ""

    def assignment_note(self) -> str:
        """Footnote stating which compounds the figure covers, and why."""

        if self.manifest is None or self.organism is None:
            return ""
        assigned = self.manifest.compounds_for(self.organism)
        if self.scoped_compounds() is not None:
            if not assigned:
                return "no compound in this series was prepared against this organism"
            return (
                "Limited to the compounds prepared against this organism: "
                + ", ".join(assigned)
            )
        if not assigned:
            return "* no compound in this series was prepared against this organism"
        return "* prepared against this organism: " + ", ".join(assigned)


@dataclass(frozen=True)
class PanelSpec:
    """One figure: the tables it needs, where it lands, and how it draws."""

    name: str
    sources: dict[str, str]
    output: str
    renderer: Callable[[dict[str, pd.DataFrame], Path, SuiteContext], str]
    description: str = ""
    optional_sources: dict[str, str] = field(default_factory=dict)
    # Panels that compare organisms against each other carry no meaning once the
    # run is filtered to one organism, and would render a misleading figure from
    # the surviving rows. They are skipped, and recorded as skipped, in the
    # per-organism suites.
    cross_organism_only: bool = False


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    """Coerce ``column`` to float, treating unparseable entries as missing."""

    return pd.to_numeric(frame.get(column), errors="coerce")


def _compound_order(frame: pd.DataFrame) -> list[str]:
    """Compounds ordered by their strongest hypothesis, best first."""

    if "query_id" not in frame:
        return []
    score = (
        "overall_priority_score"
        if "overall_priority_score" in frame
        else "chemical_evidence_score"
    )
    if score not in frame:
        return sorted(frame["query_id"].dropna().astype(str).unique())
    ranked = (
        frame.assign(_score=_numeric(frame, score))
        .groupby(frame["query_id"].astype(str))["_score"]
        .max()
        .sort_values(ascending=False)
    )
    return list(ranked.index)


# --------------------------------------------------------------------------
# Panel renderers
# --------------------------------------------------------------------------


def _compound_target_priority(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Small-multiple heatmaps: one organism x target grid per compound."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    frame = frames["predictions"].copy()
    required = {"query_id", "organism", "target_class", "overall_priority_score"}
    if not required.issubset(frame.columns):
        return STATUS_NO_ROWS
    frame["overall_priority_score"] = _numeric(frame, "overall_priority_score")
    frame = frame.dropna(subset=["overall_priority_score"])
    if frame.empty:
        return STATUS_NO_ROWS

    compounds = _compound_order(frame)
    targets = style.top_categories(
        frame, "target_class", "overall_priority_score", context.top_target_classes
    )
    organisms = sorted(frame["organism"].dropna().astype(str).unique())
    if not compounds or not targets or not organisms:
        return STATUS_NO_ROWS

    if context.organism is not None:
        # Filtered to one organism, the per-compound grids collapse to single
        # rows. A compound x target heatmap carries the same numbers legibly.
        matrix = (
            frame.pivot_table(
                index="query_id",
                columns="target_class",
                values="overall_priority_score",
                aggfunc="max",
            )
            .reindex(index=compounds, columns=targets)
        )
        figure, axis = plt.subplots(
            figsize=(
                max(9.0, 0.52 * len(targets)),
                max(4.0, 0.42 * len(compounds)) + 2.4,
            )
        )
        sns.heatmap(
            matrix,
            cmap=style.SEQUENTIAL_CMAP,
            vmin=0.0,
            linewidths=0.3,
            linecolor="white",
            cbar_kws={"label": "Overall priority score"},
            ax=axis,
        )
        axis.set_yticks(
            axis.get_yticks(), context.labels(matrix.index), rotation=0, fontsize=8
        )
        axis.set_xticks(
            axis.get_xticks(),
            style.wrap_labels(targets, width=18),
            rotation=90,
            fontsize=7,
        )
        axis.set_xlabel("Open target class")
        axis.set_ylabel("Private compound")
        axis.set_title(f"Target priority in {context.organism}")
        note = context.assignment_note()
        if note:
            figure.text(0.01, 0.005, note, fontsize=7.5, color="#4a5058")
        style.save_figure(figure, path)
        return STATUS_CREATED

    columns = min(3, len(compounds))
    rows = int(np.ceil(len(compounds) / columns))
    vmax = float(frame["overall_priority_score"].max()) or 1.0
    cmap = sns.color_palette(style.SEQUENTIAL_CMAP, as_cmap=True)

    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(6.4 * columns, 1.05 * len(organisms) * rows + 2.2),
        squeeze=False,
        layout="constrained",
    )
    image = None
    for index, axis in enumerate(axes.flat):
        if index >= len(compounds):
            axis.axis("off")
            continue
        compound = compounds[index]
        grid = (
            frame[frame["query_id"].astype(str) == compound]
            .pivot_table(
                index="organism",
                columns="target_class",
                values="overall_priority_score",
                aggfunc="max",
            )
            .reindex(index=organisms, columns=targets)
        )
        image = axis.imshow(
            grid.to_numpy(dtype=float),
            aspect="auto",
            cmap=cmap,
            vmin=0.0,
            vmax=vmax,
        )
        axis.set_title(context.label(compound))
        # Tick labels only on the outer edges: repeating six organism names and
        # eighteen target names in every cell of the grid would overrun the
        # neighbouring panels.
        if index % columns == 0:
            axis.set_yticks(range(len(organisms)), organisms, fontsize=7.5)
        else:
            axis.set_yticks([])
        if index >= len(compounds) - columns:
            axis.set_xticks(
                range(len(targets)),
                style.wrap_labels(targets, width=16),
                rotation=90,
                fontsize=6,
            )
        else:
            axis.set_xticks([])
        axis.grid(False)

    if image is not None:
        colorbar = figure.colorbar(
            image, ax=axes.ravel().tolist(), fraction=0.02, pad=0.015
        )
        colorbar.set_label("Overall priority score")
    figure.suptitle(
        "Organism-aware target priority per private compound "
        f"(top {len(targets)} target classes)"
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _organism_target_atlas(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Single atlas heatmap of mean priority, pooled across all compounds."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    frame = frames["predictions"].copy()
    if not {"organism", "target_class", "overall_priority_score"}.issubset(frame.columns):
        return STATUS_NO_ROWS
    frame["overall_priority_score"] = _numeric(frame, "overall_priority_score")
    frame = frame.dropna(subset=["overall_priority_score"])
    if frame.empty:
        return STATUS_NO_ROWS

    matrix = frame.pivot_table(
        index="organism",
        columns="target_class",
        values="overall_priority_score",
        aggfunc="mean",
    )
    ordered = matrix.max().sort_values(ascending=False).index
    matrix = matrix[ordered]

    figure, axis = plt.subplots(
        figsize=(max(10.0, 0.34 * matrix.shape[1]), max(4.0, 0.55 * matrix.shape[0]) + 2.0)
    )
    sns.heatmap(
        matrix,
        cmap=style.SEQUENTIAL_CMAP,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Mean overall priority"},
        ax=axis,
    )
    n_compounds = frame["query_id"].nunique() if "query_id" in frame else 0
    axis.set_title(
        f"Target-class atlas across {n_compounds} private compounds"
    )
    axis.set_xlabel("Open target class")
    axis.set_ylabel("Organism")
    axis.set_xticks(
        axis.get_xticks(),
        style.wrap_labels([tick.get_text() for tick in axis.get_xticklabels()], width=20),
        rotation=90,
        fontsize=6.5,
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _evidence_decomposition(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Which evidence layer carries each compound's best hypothesis."""

    import matplotlib.pyplot as plt

    frame = frames["predictions"].copy()
    columns = [column for column, _ in style.EVIDENCE_COLUMNS]
    if "query_id" not in frame or not set(columns).issubset(frame.columns):
        return STATUS_NO_ROWS
    for column in columns:
        frame[column] = _numeric(frame, column)
    frame = frame.dropna(subset=["overall_priority_score"])
    if frame.empty:
        return STATUS_NO_ROWS

    best = (
        frame.sort_values("overall_priority_score", ascending=False)
        .groupby(frame["query_id"].astype(str), as_index=False)
        .head(1)
    )
    order = _compound_order(frame)
    best = best.set_index(best["query_id"].astype(str)).reindex(order).dropna(how="all")
    if best.empty:
        return STATUS_NO_ROWS

    labels = [label for _, label in style.EVIDENCE_COLUMNS]
    positions = np.arange(len(best))
    width = 0.16
    figure, axis = plt.subplots(figsize=(max(10.0, 1.15 * len(best)), 6.0))
    for offset, (column, label, color) in enumerate(
        zip(columns, labels, style.EVIDENCE_COLORS)
    ):
        axis.bar(
            positions + (offset - (len(columns) - 1) / 2) * width,
            best[column].to_numpy(dtype=float),
            width,
            label=label,
            color=color,
        )
    annotations = [
        f"{context.label(compound)}\n"
        f"{str(row.get('target_class', ''))} | {str(row.get('organism', ''))}"
        for compound, row in best.iterrows()
    ]
    axis.set_xticks(positions, annotations, rotation=40, ha="right", fontsize=7)
    axis.set_ylim(0, 1)
    axis.set_ylabel("Score")
    axis.grid(axis="x", visible=False)
    axis.legend(
        ncol=len(labels), loc="upper center", bbox_to_anchor=(0.5, -0.28)
    )
    figure.suptitle(
        "Evidence decomposition of each compound's highest-priority hypothesis"
        + context.scope_label()
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _confidence_profile(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Confidence-class composition, per compound and per organism."""

    import matplotlib.pyplot as plt

    frame = frames["predictions"]
    classes = style.confidence_columns(frame)
    if not classes or "query_id" not in frame:
        return STATUS_NO_ROWS

    facets: list[tuple[str, str]] = [("query_id", "Private compound")]
    # An organism facet is a single bar once the run is filtered to one organism.
    if "organism" in frame and context.organism is None:
        facets.append(("organism", "Organism"))

    figure, axes = plt.subplots(
        1, len(facets), figsize=(7.5 * len(facets), 5.4), squeeze=False
    )
    for axis, (column, label) in zip(axes[0], facets):
        counts = (
            pd.crosstab(frame[column].astype(str), frame["confidence_class"].astype(str))
            .reindex(columns=classes, fill_value=0)
        )
        if column == "query_id":
            counts = counts.reindex(_compound_order(frame)).dropna(how="all")
        counts = counts.fillna(0)
        bottom = np.zeros(len(counts))
        values = counts.to_numpy(dtype=float)
        for index, name in enumerate(classes):
            axis.bar(
                np.arange(len(counts)),
                values[:, index],
                bottom=bottom,
                label=name,
                color=style.CONFIDENCE_COLORS.get(name, "#9aa2ab"),
            )
            bottom += values[:, index]
        tick_labels = (
            context.labels(counts.index)
            if column == "query_id"
            else list(counts.index)
        )
        axis.set_xticks(
            np.arange(len(counts)), tick_labels, rotation=40, ha="right", fontsize=8
        )
        axis.set_xlabel(label)
        axis.set_ylabel("Compound-target hypotheses")
        axis.grid(axis="x", visible=False)
    axes[0][-1].legend(
        title="Confidence", loc="upper left", bbox_to_anchor=(1.01, 1.0)
    )
    figure.suptitle(
        "Confidence composition of organism-aware hypotheses" + context.scope_label()
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _fusion_contribution(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Mean reciprocal-rank contribution of each fusion component, by compound."""

    import matplotlib.pyplot as plt

    frame = frames["predictions_v3"]
    contribution_columns = sorted(
        column for column in frame.columns if column.endswith("_fusion_contribution")
    )
    if "query_id" not in frame or not contribution_columns:
        return STATUS_NO_ROWS

    numeric = frame[contribution_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.notna().to_numpy().sum() == 0:
        return STATUS_NO_ROWS
    means = numeric.groupby(frame["query_id"].astype(str)).mean()
    means = means.loc[[c for c in _compound_order(frame) if c in means.index]]
    if means.empty:
        return STATUS_NO_ROWS

    labels = [
        column.removesuffix("_fusion_contribution").replace("_", " ")
        for column in contribution_columns
    ]
    figure, axis = plt.subplots(figsize=(max(10.0, 1.0 * len(means)), 6.0))
    bottom = np.zeros(len(means))
    values = means.to_numpy(dtype=float)
    for index, label in enumerate(labels):
        column_values = np.nan_to_num(values[:, index])
        axis.bar(
            np.arange(len(means)),
            column_values,
            bottom=bottom,
            label=label,
            color=style.COMPOUND_COLORS[index % len(style.COMPOUND_COLORS)],
        )
        bottom += column_values
    axis.set_xticks(
        np.arange(len(means)), context.labels(means.index), rotation=40, ha="right"
    )
    axis.set_ylabel("Mean reciprocal-rank contribution")
    axis.set_xlabel("Private compound")
    axis.set_title(
        "Fusion component contributions to the v3 chemical evidence score"
        + context.scope_label()
    )
    axis.legend(ncol=2, loc="upper right")
    axis.grid(axis="x", visible=False)
    style.save_figure(figure, path)
    return STATUS_CREATED


def _rank_shift_disagreements(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Where 3D fusion moved a target relative to the 2D-only ranking."""

    import matplotlib.pyplot as plt

    frame = frames["disagreements"].copy()
    required = {"query_id", "target_class", "rank_shift_2d_to_v3"}
    if not required.issubset(frame.columns):
        return STATUS_NO_ROWS
    frame["rank_shift_2d_to_v3"] = _numeric(frame, "rank_shift_2d_to_v3")
    frame["absolute_rank_shift"] = (
        _numeric(frame, "absolute_rank_shift")
        if "absolute_rank_shift" in frame
        else frame["rank_shift_2d_to_v3"].abs()
    )
    frame = frame.dropna(subset=["rank_shift_2d_to_v3", "absolute_rank_shift"])
    if frame.empty:
        return STATUS_NO_ROWS

    selected = (
        frame.sort_values("absolute_rank_shift", ascending=False)
        .groupby(frame["query_id"].astype(str), as_index=False)
        .head(context.max_disagreements_per_compound)
    )
    if selected.empty:
        return STATUS_NO_ROWS
    selected = selected.sort_values(
        ["query_id", "rank_shift_2d_to_v3"], ascending=[True, True]
    )
    labels = (
        selected["query_id"].astype(str) + "  ·  " + selected["target_class"].astype(str)
    )
    colors = np.where(
        selected["rank_shift_2d_to_v3"] > 0,
        style.DIVERGING_POSITIVE,
        style.DIVERGING_NEGATIVE,
    )

    figure, axis = plt.subplots(figsize=(9.5, max(5.0, 0.24 * len(selected))))
    positions = np.arange(len(selected))
    axis.barh(positions, selected["rank_shift_2d_to_v3"].to_numpy(dtype=float), color=colors)
    axis.set_yticks(positions, list(labels), fontsize=6.5)
    axis.axvline(0, color="#22262b", linewidth=0.9)
    axis.set_xlabel("Rank shift (positive = promoted by 3D fusion)")
    axis.set_title(
        "Largest 2D-versus-v3 target-rank disagreements per compound "
        f"(top {context.max_disagreements_per_compound} each)"
    )
    axis.grid(axis="y", visible=False)
    style.save_figure(figure, path)
    return STATUS_CREATED


def _ranking_stability(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """How far the final ranking moves when scoring weights are perturbed."""

    import matplotlib.pyplot as plt

    frame = frames["sensitivity"].copy()
    if "scenario_type" not in frame:
        return STATUS_NO_ROWS
    frame["mean_kendall_tau"] = _numeric(frame, "mean_kendall_tau")
    frame["mean_rbo"] = _numeric(frame, "mean_rbo")
    frame = frame.dropna(subset=["mean_kendall_tau", "mean_rbo"], how="all")
    if frame.empty:
        return STATUS_NO_ROWS

    scenario_types = list(
        frame.groupby("scenario_type")["mean_kendall_tau"].min().sort_values().index
    )
    figure, axes = plt.subplots(
        1, 2, figsize=(13.5, max(4.0, 0.55 * len(scenario_types)) + 1.5), sharey=True
    )
    metrics = [
        ("mean_kendall_tau", "kendall_tau_ci_lower_95", "kendall_tau_ci_upper_95",
         "Kendall tau vs configured ranking", style.EVIDENCE_COLORS[0]),
        ("mean_rbo", "rbo_ci_lower_95", "rbo_ci_upper_95",
         "Rank-biased overlap vs configured ranking", style.EVIDENCE_COLORS[1]),
    ]
    for axis, (value_column, lower_column, upper_column, title, color) in zip(axes, metrics):
        for position, scenario_type in enumerate(scenario_types):
            subset = frame[frame["scenario_type"] == scenario_type]
            values = subset[value_column].to_numpy(dtype=float)
            jitter = np.linspace(-0.18, 0.18, num=max(len(values), 1))
            axis.scatter(
                values,
                np.full(len(values), position) + jitter,
                s=22,
                color=color,
                alpha=0.75,
                edgecolor="white",
                linewidth=0.4,
            )
            if lower_column in subset and upper_column in subset:
                low = _numeric(subset, lower_column).min()
                high = _numeric(subset, upper_column).max()
                if np.isfinite(low) and np.isfinite(high):
                    axis.plot([low, high], [position, position], color=color, alpha=0.35, linewidth=6)
        axis.set_yticks(
            range(len(scenario_types)),
            [name.replace("_", " ") for name in scenario_types],
            fontsize=8,
        )
        axis.set_xlim(min(0.5, float(frame[value_column].min(skipna=True)) - 0.02), 1.005)
        axis.set_xlabel(title)
        axis.grid(axis="y", visible=False)
    figure.suptitle(
        f"Final-ranking stability across {len(frame)} weight-perturbation scenarios"
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _uncertainty_landscape(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Bootstrap stability against the empirical decoy p-value, per compound."""

    import matplotlib.pyplot as plt

    frame = frames["uncertainty"].copy()
    required = {"query_id", "bootstrap_stability_score", "empirical_decoy_p_value"}
    if not required.issubset(frame.columns):
        return STATUS_NO_ROWS
    frame["bootstrap_stability_score"] = _numeric(frame, "bootstrap_stability_score")
    frame["empirical_decoy_p_value"] = _numeric(frame, "empirical_decoy_p_value")
    frame = frame.dropna(subset=["bootstrap_stability_score", "empirical_decoy_p_value"])
    if frame.empty:
        return STATUS_NO_ROWS

    compounds = sorted(frame["query_id"].astype(str).unique())
    palette = style.compound_palette(compounds)
    sizes = _numeric(frame, "bootstrap_top1_probability")
    marker_size = (
        30.0 + 170.0 * sizes.fillna(0.0).clip(0, 1) if sizes.notna().any() else 60.0
    )

    figure, axis = plt.subplots(figsize=(11.0, 7.0))
    for compound in compounds:
        subset = frame[frame["query_id"].astype(str) == compound]
        axis.scatter(
            subset["bootstrap_stability_score"],
            subset["empirical_decoy_p_value"],
            s=marker_size.loc[subset.index] if hasattr(marker_size, "loc") else marker_size,
            color=palette[compound],
            alpha=0.75,
            edgecolor="white",
            linewidth=0.5,
            label=context.label(compound),
        )
    axis.axhline(0.05, color="#22262b", linestyle="--", linewidth=0.9)
    axis.text(
        axis.get_xlim()[0], 0.055, "decoy p = 0.05", fontsize=7.5, color="#22262b"
    )
    axis.set_xlabel("Bootstrap rank-stability score")
    axis.set_ylabel("Empirical decoy p-value (lower is stronger)")
    axis.set_title(
        "Hypothesis uncertainty landscape"
        + context.scope_label()
        + " (marker size = bootstrap top-1 probability)"
    )
    axis.legend(
        title="Compound",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        fontsize=8,
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _benchmark_enrichment(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Retrieval enrichment over a random ranking, by split and scoring mode."""

    import matplotlib.pyplot as plt

    frame = frames["benchmark"].copy()
    metrics = ["top1_enrichment_over_random", "top3_enrichment_over_random"]
    if not {"split", "mode"}.issubset(frame.columns):
        return STATUS_NO_ROWS
    available = [metric for metric in metrics if metric in frame.columns]
    for metric in available:
        frame[metric] = _numeric(frame, metric)
    frame = frame.dropna(subset=available, how="all") if available else frame.iloc[0:0]
    if frame.empty or not available:
        return STATUS_NO_ROWS

    modes = list(dict.fromkeys(frame["mode"].astype(str)))
    splits = list(dict.fromkeys(frame["split"].astype(str)))
    figure, axes = plt.subplots(
        1, len(available), figsize=(6.8 * len(available), 5.4), squeeze=False
    )
    width = 0.8 / max(len(modes), 1)
    for axis, metric in zip(axes[0], available):
        positions = np.arange(len(splits))
        for index, mode in enumerate(modes):
            subset = frame[frame["mode"].astype(str) == mode].set_index(
                frame.loc[frame["mode"].astype(str) == mode, "split"].astype(str)
            )
            values = [float(subset[metric].get(split, np.nan)) for split in splits]
            axis.bar(
                positions + (index - (len(modes) - 1) / 2) * width,
                values,
                width,
                label=mode.replace("_", " "),
                color=style.COMPOUND_COLORS[index % len(style.COMPOUND_COLORS)],
            )
        axis.axhline(1.0, color="#22262b", linestyle="--", linewidth=0.9)
        axis.set_xticks(positions, splits, rotation=15, ha="right")
        axis.set_ylabel("Enrichment over random")
        axis.set_title(metric.replace("_", " "))
        axis.grid(axis="x", visible=False)
    axes[0][0].legend(title="Scoring mode")
    figure.suptitle("Benchmark retrieval enrichment against a random ranking")
    style.save_figure(figure, path)
    return STATUS_CREATED


def _assignment_concordance(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Does each compound score best against the organism it was made for?

    The manifest records design intent. This panel puts that intent next to the
    evidence: the compound's priority in its assigned organism against its best
    score in any other organism. It reports concordance, it does not enforce it.
    """

    import matplotlib.pyplot as plt

    if context.manifest is None or not context.manifest:
        return STATUS_SOURCE_MISSING
    frame = frames["predictions"].copy()
    if not {"query_id", "organism", "overall_priority_score"}.issubset(frame.columns):
        return STATUS_NO_ROWS
    frame["overall_priority_score"] = _numeric(frame, "overall_priority_score")
    frame = frame.dropna(subset=["overall_priority_score"])
    if frame.empty:
        return STATUS_NO_ROWS

    assignments = context.manifest.assignments
    rows = []
    for compound, assigned in sorted(assignments.items()):
        subset = frame[frame["query_id"].astype(str) == compound]
        if subset.empty:
            continue
        own = subset[subset["organism"].astype(str) == assigned]
        others = subset[subset["organism"].astype(str) != assigned]
        if own.empty:
            continue
        best_other = others.loc[others["overall_priority_score"].idxmax()] if not others.empty else None
        rows.append(
            {
                "compound": compound,
                "assigned": assigned,
                "assigned_score": float(own["overall_priority_score"].max()),
                "best_other_score": (
                    float(best_other["overall_priority_score"]) if best_other is not None else np.nan
                ),
                "best_other_organism": (
                    str(best_other["organism"]) if best_other is not None else ""
                ),
            }
        )
    if not rows:
        return STATUS_NO_ROWS
    table = pd.DataFrame(rows).sort_values("assigned_score", ascending=False)

    positions = np.arange(len(table))
    width = 0.38
    figure, axis = plt.subplots(figsize=(max(9.0, 1.05 * len(table)), 6.2))
    axis.bar(
        positions - width / 2,
        table["assigned_score"],
        width,
        label="Assigned organism",
        color=style.EVIDENCE_COLORS[1],
    )
    axis.bar(
        positions + width / 2,
        table["best_other_score"],
        width,
        label="Best other organism",
        color=style.EVIDENCE_COLORS[4],
    )
    labels = [
        f"{row.compound}\n{row.assigned}" for row in table.itertuples()
    ]
    axis.set_xticks(positions, labels, rotation=40, ha="right", fontsize=7.5)
    axis.set_ylabel("Overall priority score")
    axis.grid(axis="x", visible=False)
    axis.legend(loc="upper right")
    concordant = int((table["assigned_score"] >= table["best_other_score"].fillna(-1)).sum())
    figure.suptitle(
        "Design intent versus evidence: assigned organism against best alternative "
        f"({concordant} of {len(table)} compounds score highest in their assigned organism)"
    )
    style.save_figure(figure, path)
    return STATUS_CREATED


def _chemical_space(
    frames: dict[str, pd.DataFrame], path: Path, context: SuiteContext
) -> str:
    """Private compounds against public benchmark drugs in ECFP4 space.

    Requires RDKit and UMAP. Both ship compiled extensions that are unavailable
    on some hosts, so an import failure is reported as a dependency status
    rather than raised: the rest of the suite must still render.
    """

    import matplotlib.pyplot as plt

    if context.private_compounds is None or not context.private_compounds.is_file():
        return STATUS_SOURCE_MISSING
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem
        from umap import UMAP
    except Exception:  # pragma: no cover - host-dependent native import
        return STATUS_DEPENDENCY

    records: list[tuple[str, str, object]] = []
    for molecule in Chem.SDMolSupplier(str(context.private_compounds), removeHs=True):
        if molecule is None:
            continue
        name = molecule.GetProp("_Name") if molecule.HasProp("_Name") else "compound"
        records.append((name, "private compound", molecule))
    if context.benchmark_compounds is not None and context.benchmark_compounds.is_file():
        benchmark = pd.read_csv(context.benchmark_compounds)
        smiles_column = next(
            (column for column in ("canonical_smiles", "smiles") if column in benchmark),
            None,
        )
        name_column = next(
            (column for column in ("drug", "name", "compound") if column in benchmark),
            None,
        )
        if smiles_column and name_column:
            for _, row in benchmark.iterrows():
                molecule = Chem.MolFromSmiles(str(row[smiles_column]))
                if molecule is not None:
                    records.append((str(row[name_column]), "benchmark drug", molecule))
    if len(records) < 5:
        return STATUS_NO_ROWS

    generator = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    matrix = np.zeros((len(records), 2048), dtype=np.float32)
    for index, (_, _, molecule) in enumerate(records):
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(molecule), matrix[index])
    embedding = UMAP(
        n_neighbors=min(10, len(records) - 1),
        min_dist=0.25,
        metric="jaccard",
        random_state=context.umap_seed,
    ).fit_transform(matrix)

    frame = pd.DataFrame(
        {
            "id": [name for name, _, _ in records],
            "source": [source for _, source, _ in records],
            "x": embedding[:, 0],
            "y": embedding[:, 1],
        }
    )
    figure, axis = plt.subplots(figsize=(10.0, 8.0))
    for source, color, marker in (
        ("benchmark drug", "#7a869a", "o"),
        ("private compound", style.DIVERGING_NEGATIVE, "D"),
    ):
        subset = frame[frame["source"] == source]
        axis.scatter(
            subset["x"], subset["y"], s=70, color=color, marker=marker,
            alpha=0.85, edgecolor="white", linewidth=0.6, label=source,
        )
    for _, row in frame[frame["source"] == "private compound"].iterrows():
        axis.annotate(row["id"], (row["x"], row["y"]), fontsize=7.5,
                      xytext=(4, 4), textcoords="offset points")
    axis.set_xlabel("UMAP-1")
    axis.set_ylabel("UMAP-2")
    axis.set_title("ECFP4 chemical space: private compounds versus benchmark drugs")
    axis.legend()
    style.save_figure(figure, path)
    return STATUS_CREATED


PANELS: tuple[PanelSpec, ...] = (
    PanelSpec(
        name="compound_target_priority",
        sources={"predictions": PREDICTIONS_BY_ORGANISM},
        output="compound_target_priority.png",
        renderer=_compound_target_priority,
        description="Per-compound organism x target-class priority heatmaps.",
    ),
    PanelSpec(
        name="organism_target_atlas",
        sources={"predictions": PREDICTIONS_BY_ORGANISM},
        output="organism_target_atlas.png",
        renderer=_organism_target_atlas,
        description="Mean priority per organism and target class across all compounds.",
        cross_organism_only=True,
    ),
    PanelSpec(
        name="evidence_decomposition",
        sources={"predictions": PREDICTIONS_BY_ORGANISM},
        output="evidence_decomposition.png",
        renderer=_evidence_decomposition,
        description="Evidence layers behind each compound's strongest hypothesis.",
    ),
    PanelSpec(
        name="confidence_profile",
        sources={"predictions": PREDICTIONS_BY_ORGANISM},
        output="confidence_profile.png",
        renderer=_confidence_profile,
        description="Confidence-class composition by compound and by organism.",
    ),
    PanelSpec(
        name="fusion_contribution",
        sources={"predictions_v3": PREDICTIONS_V3},
        output="fusion_contribution.png",
        renderer=_fusion_contribution,
        description="Reciprocal-rank contribution of each fusion component.",
    ),
    PanelSpec(
        name="rank_shift_disagreements",
        sources={"disagreements": DISAGREEMENTS_V3},
        output="rank_shift_disagreements.png",
        renderer=_rank_shift_disagreements,
        description="Targets that 3D fusion promoted or demoted against 2D-only.",
    ),
    PanelSpec(
        name="ranking_stability",
        sources={"sensitivity": SENSITIVITY_V3},
        output="ranking_stability.png",
        renderer=_ranking_stability,
        description="Kendall tau and RBO under scoring-weight perturbation.",
    ),
    PanelSpec(
        name="uncertainty_landscape",
        sources={"uncertainty": UNCERTAINTY_V2},
        output="uncertainty_landscape.png",
        renderer=_uncertainty_landscape,
        description="Bootstrap stability against empirical decoy p-values.",
    ),
    PanelSpec(
        name="benchmark_enrichment",
        sources={"benchmark": BENCHMARK_SUMMARY_V2},
        output="benchmark_enrichment.png",
        renderer=_benchmark_enrichment,
        description="Retrieval enrichment over random, by split and scoring mode.",
    ),
    PanelSpec(
        name="assignment_concordance",
        sources={"predictions": PREDICTIONS_BY_ORGANISM},
        output="assignment_concordance.png",
        renderer=_assignment_concordance,
        description="Manifest-assigned organism against the best alternative organism.",
        cross_organism_only=True,
    ),
    PanelSpec(
        name="chemical_space",
        sources={},
        output="chemical_space.png",
        renderer=_chemical_space,
        description="ECFP4 UMAP of private compounds against benchmark drugs.",
    ),
)


def organism_slug(organism: str) -> str:
    """Filesystem-safe directory name for one organism."""

    return "_".join(str(organism).split()).lower()


def _restrict_to_organism(frame: pd.DataFrame, organism: str) -> pd.DataFrame:
    """Keep only rows belonging to ``organism``.

    Tables carrying no organism column are run-level and organism-agnostic
    (ranking stability, benchmark summaries). They are passed through unchanged
    rather than silently emptied: filtering them by organism would be inventing
    a distinction the pipeline did not compute.
    """

    if "organism" not in frame.columns:
        return frame
    return frame[frame["organism"].astype(str) == organism]


def _restrict_to_compounds(frame: pd.DataFrame, compounds: list[str]) -> pd.DataFrame:
    """Keep only rows for ``compounds``, matched on ``query_id``.

    Tables with no ``query_id`` are run-level and pass through unchanged, for
    the same reason as ``_restrict_to_organism``.
    """

    if "query_id" not in frame.columns:
        return frame
    return frame[frame["query_id"].astype(str).isin(compounds)]


def generate_figure_suite(
    results_dir: Path,
    *,
    context: SuiteContext | None = None,
    figure_dirname: str = "figures_suite",
) -> pd.DataFrame:
    """Render every panel it can and return one status row per panel."""

    style.apply_theme()
    context = context or SuiteContext()
    figure_dir = results_dir / figure_dirname
    figure_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for panel in PANELS:
        output = figure_dir / panel.output
        if panel.cross_organism_only and context.organism is not None:
            rows.append(
                {
                    "organism": context.organism,
                    "figure": panel.name,
                    "status": STATUS_NOT_APPLICABLE,
                    "source": "",
                    "output": None,
                    "description": panel.description,
                }
            )
            continue
        source_paths = {alias: results_dir / name for alias, name in panel.sources.items()}
        missing = [str(path) for path in source_paths.values() if not path.is_file()]
        if missing:
            status = STATUS_SOURCE_MISSING
        else:
            frames = {alias: pd.read_csv(path) for alias, path in source_paths.items()}
            if context.organism is not None:
                frames = {
                    alias: _restrict_to_organism(frame, context.organism)
                    for alias, frame in frames.items()
                }
                scope = context.scoped_compounds()
                if scope is not None:
                    frames = {
                        alias: _restrict_to_compounds(frame, scope)
                        for alias, frame in frames.items()
                    }
            if any(frame.empty for frame in frames.values()):
                status = STATUS_NO_ROWS
            else:
                status = panel.renderer(frames, output, context)
        rows.append(
            {
                "organism": context.organism or "",
                "figure": panel.name,
                "status": status,
                "source": ";".join(str(path) for path in source_paths.values()),
                "output": str(output) if status == STATUS_CREATED else None,
                "description": panel.description,
            }
        )
    return pd.DataFrame(rows)


def generate_per_organism_suites(
    results_dir: Path,
    *,
    context: SuiteContext,
    organisms: list[str],
    figure_dirname: str = "figures_suite",
) -> pd.DataFrame:
    """Render one focused suite per organism under ``figures_suite/by_organism``.

    Each suite covers the compounds the manifest prepared against that organism,
    scored against it: the Klebsiella folder is about the Klebsiella compounds.
    The cross-series comparison is not lost, it lives in the overall suite, whose
    per-compound panels already show every compound against every organism.

    With ``restrict_to_assigned`` false, or with no manifest, a suite falls back
    to the whole series scored against that organism, with assigned compounds
    marked in the axis labels.
    """

    from dataclasses import replace

    statuses = []
    for organism in organisms:
        organism_context = replace(
            context, organism=organism, restrict_to_assigned=True
        )
        statuses.append(
            generate_figure_suite(
                results_dir,
                context=organism_context,
                figure_dirname=f"{figure_dirname}/by_organism/{organism_slug(organism)}",
            )
        )
    if not statuses:
        return pd.DataFrame(
            columns=["organism", "figure", "status", "source", "output", "description"]
        )
    return pd.concat(statuses, ignore_index=True)


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    manifest = load_from_config(config)
    context = SuiteContext(
        private_compounds=config.path_for("private_compounds"),
        benchmark_compounds=config.path_for("benchmark"),
        umap_seed=int(config.value("seeds.umap")),
        figure_sample_seed=int(config.value("seeds.figure_sampling")),
        manifest=manifest,
    )
    if manifest.unresolved_groups:
        print(
            "WARNING: manifest microbe groups with no organisms.manifest_aliases "
            f"entry, left unassigned: {', '.join(manifest.unresolved_groups)}"
        )

    status = generate_figure_suite(results_dir, context=context)
    organisms = list(config.value("organisms.names"))
    status = pd.concat(
        [
            status,
            generate_per_organism_suites(
                results_dir, context=context, organisms=organisms
            ),
        ],
        ignore_index=True,
    )
    status.to_csv(results_dir / "figure_suite_status.csv", index=False)

    summary = status.drop(columns=["description", "source"])
    print(summary.to_string(index=False))
    print(
        f"\ncreated {(status.status == STATUS_CREATED).sum()} of {len(status)} panels "
        f"across 1 overall suite and {len(organisms)} per-organism suites"
    )


if __name__ == "__main__":
    main()
