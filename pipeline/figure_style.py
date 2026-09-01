"""Shared visual language for every generated figure.

One palette, one theme, one save path. Panels import from here instead of
choosing their own colours so that a figure suite rendered across several
modules still reads as a single document.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

# Ordered confidence vocabulary emitted by the v2 scoring layer. The order is
# meaningful: stacked bars and legends must follow it so that "High" is always
# the first band and "Insufficient" always the last.
CONFIDENCE_ORDER: tuple[str, ...] = ("High", "Moderate", "Low", "Insufficient")

CONFIDENCE_COLORS: dict[str, str] = {
    "High": "#1b7f5f",
    "Moderate": "#d9a441",
    "Low": "#d1683a",
    "Insufficient": "#9aa2ab",
}

# Evidence layers of the overall priority score, in the order the pipeline
# composes them (chemistry -> organism transfer -> structure -> biology).
EVIDENCE_COLUMNS: tuple[tuple[str, str], ...] = (
    ("chemical_quality_adjusted_score", "Chemical"),
    ("species_transfer_score", "Species transfer"),
    ("pocket_evidence_score", "Pocket"),
    ("biological_priority_score", "Biology"),
    ("overall_priority_score", "Overall"),
)

EVIDENCE_COLORS: tuple[str, ...] = (
    "#3b6ea5",
    "#2a9d8f",
    "#e9a03c",
    "#8b5fbf",
    "#c1443c",
)

# Qualitative palette for per-compound series. Twelve private compounds is the
# current study size; the sequence repeats safely if that grows.
COMPOUND_COLORS: tuple[str, ...] = (
    "#3b6ea5",
    "#c1443c",
    "#2a9d8f",
    "#e9a03c",
    "#8b5fbf",
    "#5d7a3a",
    "#c2557f",
    "#4aa3c7",
    "#a5673f",
    "#7a7f8a",
    "#d4823a",
    "#5b5fa8",
)

SEQUENTIAL_CMAP = "mako"
DIVERGING_POSITIVE = "#2a7f62"
DIVERGING_NEGATIVE = "#c1443c"
GRID_COLOR = "#c9ced4"

FIGURE_DPI = 300


def apply_theme() -> None:
    """Install the shared matplotlib/seaborn theme. Safe to call repeatedly."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#4a5058",
            "axes.labelcolor": "#22262b",
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "legend.frameon": False,
            "grid.color": GRID_COLOR,
            "grid.linewidth": 0.6,
            "savefig.bbox": "tight",
        }
    )


def save_figure(figure, path: Path) -> None:
    """Write one figure at publication resolution and release its memory."""

    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=FIGURE_DPI)
    plt.close(figure)


def compound_palette(compounds: Sequence[str]) -> dict[str, str]:
    """Map compound identifiers to stable colours, ordered as given."""

    return {
        compound: COMPOUND_COLORS[index % len(COMPOUND_COLORS)]
        for index, compound in enumerate(compounds)
    }


def confidence_columns(frame: pd.DataFrame) -> list[str]:
    """Confidence classes present in ``frame``, in the canonical order.

    Unknown classes are appended rather than dropped so that an unexpected
    scoring vocabulary is visible in the figure instead of silently missing.
    """

    if "confidence_class" not in frame:
        return []
    observed = set(frame["confidence_class"].dropna().astype(str))
    ordered = [name for name in CONFIDENCE_ORDER if name in observed]
    return ordered + sorted(observed - set(CONFIDENCE_ORDER))


def top_categories(
    frame: pd.DataFrame, column: str, value_column: str, limit: int
) -> list[str]:
    """The ``limit`` categories with the highest peak ``value_column``.

    Used to keep wide target-class axes legible without hiding the strongest
    hypotheses: selection is by maximum, never by mean, so a target that is
    outstanding for a single compound still survives the trim.
    """

    if column not in frame or value_column not in frame:
        return []
    ranked = (
        frame.groupby(column)[value_column]
        .max()
        .sort_values(ascending=False)
        .head(limit)
    )
    return list(ranked.index.astype(str))


def wrap_labels(labels: Iterable[str], width: int = 22) -> list[str]:
    """Soft-wrap long ontology labels so dense axes stay readable."""

    import textwrap

    return ["\n".join(textwrap.wrap(str(label), width=width)) or str(label) for label in labels]
