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
except ModuleNotFoundError:  # pragma: no cover
    from config import load_config
    import figure_style as style


STATUS_CREATED = "created"
STATUS_SOURCE_MISSING = "unavailable_source_missing"
STATUS_NO_ROWS = "unavailable_no_evaluable_rows"
STATUS_DEPENDENCY = "unavailable_dependency_unavailable"

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


@dataclass(frozen=True)
class PanelSpec:
    """One figure: the tables it needs, where it lands, and how it draws."""

    name: str
    sources: dict[str, str]
    output: str
    renderer: Callable[[dict[str, pd.DataFrame], Path, SuiteContext], str]
    description: str = ""
    optional_sources: dict[str, str] = field(default_factory=dict)


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
        axis.set_title(compound)
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
        f"{compound}\n{str(row.get('target_class', ''))} | {str(row.get('organism', ''))}"
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
    if "organism" in frame:
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
        axis.set_xticks(
            np.arange(len(counts)), list(counts.index), rotation=40, ha="right", fontsize=8
        )
        axis.set_xlabel(label)
        axis.set_ylabel("Compound-target hypotheses")
        axis.grid(axis="x", visible=False)
    axes[0][-1].legend(
        title="Confidence", loc="upper left", bbox_to_anchor=(1.01, 1.0)
    )
    figure.suptitle("Confidence composition of organism-aware hypotheses")
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
    axis.set_xticks(np.arange(len(means)), list(means.index), rotation=40, ha="right")
    axis.set_ylabel("Mean reciprocal-rank contribution")
    axis.set_xlabel("Private compound")
    axis.set_title("Fusion component contributions to the v3 chemical evidence score")
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
            label=compound,
        )
    axis.axhline(0.05, color="#22262b", linestyle="--", linewidth=0.9)
    axis.text(
        axis.get_xlim()[0], 0.055, "decoy p = 0.05", fontsize=7.5, color="#22262b"
    )
    axis.set_xlabel("Bootstrap rank-stability score")
    axis.set_ylabel("Empirical decoy p-value (lower is stronger)")
    axis.set_title(
        "Hypothesis uncertainty landscape "
        "(marker size = bootstrap top-1 probability)"
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
        name="chemical_space",
        sources={},
        output="chemical_space.png",
        renderer=_chemical_space,
        description="ECFP4 UMAP of private compounds against benchmark drugs.",
    ),
)


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
        source_paths = {alias: results_dir / name for alias, name in panel.sources.items()}
        missing = [str(path) for path in source_paths.values() if not path.is_file()]
        if missing:
            status = STATUS_SOURCE_MISSING
        else:
            frames = {alias: pd.read_csv(path) for alias, path in source_paths.items()}
            if any(frame.empty for frame in frames.values()):
                status = STATUS_NO_ROWS
            else:
                status = panel.renderer(frames, output, context)
        rows.append(
            {
                "figure": panel.name,
                "status": status,
                "source": ";".join(str(path) for path in source_paths.values()),
                "output": str(output) if status == STATUS_CREATED else None,
                "description": panel.description,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    context = SuiteContext(
        private_compounds=config.path_for("private_compounds"),
        benchmark_compounds=config.path_for("benchmark"),
        umap_seed=int(config.value("seeds.umap")),
        figure_sample_seed=int(config.value("seeds.figure_sampling")),
    )
    status = generate_figure_suite(results_dir, context=context)
    status.to_csv(results_dir / "figure_suite_status.csv", index=False)
    print(status.drop(columns=["description"]).to_string(index=False))


if __name__ == "__main__":
    main()
