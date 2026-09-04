"""Tie an external docking/MD campaign to the organism-aware pipeline ranking.

The pipeline proposes target hypotheses; a docking campaign tests some of them.
This module joins the two so a run can be reported against what was actually
docked, and answers three questions explicitly:

1. **Coverage** — which docked receptors map onto pipeline target classes, and
   which of the pipeline's strongest hypotheses were never docked at all.
2. **Agreement** — over the shared targets only, does docking affinity rank
   targets the way the pipeline ranks them?
3. **Size bias** — docking scores grow with molecular size, so raw affinity is
   also reported per heavy atom (ligand efficiency) against the native
   co-crystal ligand, which is the control that says whether an affinity is
   good *for a molecule that size*.

Nothing here rescores the pipeline. Docking is evidence about a pose in one
structure; the pipeline ranks target-class plausibility. They are different
quantities, so only their rank agreement over shared targets is compared, and
coverage is reported alongside it: a correlation over five targets means little
without knowing which targets were left out.

The mapping from docking receptor to target class is declared in ``config.yaml``
and never inferred. An unmapped receptor is reported, not dropped silently — a
cross-check that quietly matches zero rows reports agreement it never measured.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

try:  # pragma: no cover - import shim for direct module execution
    from pipeline.config import load_config
    from pipeline import figure_style as style
except ModuleNotFoundError:  # pragma: no cover
    from config import load_config
    import figure_style as style


PREDICTIONS = "v2_open_target_predictions_by_organism.csv"
CROSSCHECK_CSV = "docking_pipeline_crosscheck.csv"
COVERAGE_CSV = "docking_pipeline_coverage.csv"
AGREEMENT_CSV = "docking_pipeline_agreement.csv"


class DockingCrosscheckError(ValueError):
    """Raised when the cross-check cannot be performed honestly."""


@dataclass(frozen=True)
class DockingSpec:
    """How to read one docking campaign's summary table."""

    organism: str
    target_column: str = "target_folder"
    ligand_column: str = "ligand"
    affinity_column: str = "docking_affinity_kcal_mol"
    target_aliases: dict[str, str] | None = None
    excluded_targets: dict[str, str] | None = None
    control_ligand_prefixes: tuple[str, ...] = ()

    @classmethod
    def from_config(cls, config) -> "DockingSpec":
        return cls(
            organism=str(config.value("docking.organism")),
            target_column=str(config.value("docking.target_column")),
            ligand_column=str(config.value("docking.ligand_column")),
            affinity_column=str(config.value("docking.affinity_column")),
            target_aliases=dict(config.value("docking.target_aliases")),
            excluded_targets=dict(config.value("docking.excluded_targets")),
            control_ligand_prefixes=tuple(
                config.value("docking.control_ligand_prefixes")
            ),
        )

    def is_control(self, ligand: str) -> bool:
        return str(ligand).startswith(self.control_ligand_prefixes)


def _best_affinity(frame: pd.DataFrame, spec: DockingSpec) -> pd.DataFrame:
    """Most favourable (most negative) affinity per receptor and ligand."""

    working = frame.copy()
    working[spec.affinity_column] = pd.to_numeric(
        working[spec.affinity_column], errors="coerce"
    )
    working = working.dropna(subset=[spec.affinity_column])
    index = working.groupby([spec.target_column, spec.ligand_column])[
        spec.affinity_column
    ].idxmin()
    return working.loc[index]


def build_coverage(
    docking: pd.DataFrame, spec: DockingSpec, predictions: pd.DataFrame
) -> pd.DataFrame:
    """One row per docked receptor: how it maps, or why it does not."""

    aliases = spec.target_aliases or {}
    excluded = spec.excluded_targets or {}
    known_classes = set(predictions["target_class"].astype(str))

    rows = []
    for receptor in sorted(docking[spec.target_column].astype(str).unique()):
        if receptor in excluded:
            rows.append(
                {
                    "docking_receptor": receptor,
                    "target_class": None,
                    "mapping_status": "excluded_by_config",
                    "detail": excluded[receptor],
                }
            )
            continue
        target_class = aliases.get(receptor)
        if target_class is None:
            rows.append(
                {
                    "docking_receptor": receptor,
                    "target_class": None,
                    "mapping_status": "unmapped_no_alias",
                    "detail": "add an entry to docking.target_aliases in config.yaml",
                }
            )
        elif target_class not in known_classes:
            rows.append(
                {
                    "docking_receptor": receptor,
                    "target_class": target_class,
                    "mapping_status": "unmapped_class_absent_from_run",
                    "detail": f"{target_class!r} is not a target class in this run",
                }
            )
        else:
            rows.append(
                {
                    "docking_receptor": receptor,
                    "target_class": target_class,
                    "mapping_status": "mapped",
                    "detail": "",
                }
            )
    return pd.DataFrame(rows)


def build_crosscheck(
    predictions: pd.DataFrame,
    docking: pd.DataFrame,
    spec: DockingSpec,
    *,
    heavy_atoms: dict[str, int] | None = None,
) -> pd.DataFrame:
    """Join docked affinity to the pipeline ranking, one row per target/compound."""

    coverage = build_coverage(docking, spec, predictions)
    mapped = coverage[coverage.mapping_status == "mapped"]
    mapping = dict(zip(mapped.docking_receptor, mapped.target_class))
    if not mapping:
        # Fail loudly and early. The alternative is an empty join that reads as
        # "no disagreement found" while nothing was ever compared.
        raise DockingCrosscheckError(
            "no docked receptor mapped onto a target class in this run; "
            "check docking.target_aliases in config.yaml"
        )

    best = _best_affinity(docking, spec).copy()
    best["target_class"] = best[spec.target_column].astype(str).map(mapping)
    best = best.dropna(subset=["target_class"])
    best["is_control"] = best[spec.ligand_column].map(spec.is_control)

    organism_rows = predictions[
        predictions["organism"].astype(str) == spec.organism
    ].copy()
    if organism_rows.empty:
        raise DockingCrosscheckError(
            f"no prediction rows for organism {spec.organism!r}"
        )
    organism_rows["overall_priority_score"] = pd.to_numeric(
        organism_rows["overall_priority_score"], errors="coerce"
    )
    organism_rows["pipeline_rank"] = organism_rows.groupby("query_id")[
        "overall_priority_score"
    ].rank(ascending=False, method="min")
    n_classes = organism_rows["target_class"].nunique()

    heavy = heavy_atoms or {}

    # Queries: the compounds the pipeline scored. Controls keep their own rows so
    # the native co-crystal ligand can anchor the ligand-efficiency comparison.
    queries = best[~best.is_control].copy()
    joined = queries.merge(
        organism_rows[
            [
                "query_id",
                "target_class",
                "overall_priority_score",
                "chemical_evidence_score",
                "confidence_class",
                "species_transfer_score",
                "pocket_evidence_score",
                "pipeline_rank",
            ]
        ],
        left_on=[spec.ligand_column, "target_class"],
        right_on=["query_id", "target_class"],
        how="left",
    )
    joined["n_target_classes_in_run"] = n_classes
    joined = joined.rename(columns={spec.affinity_column: "docking_affinity_kcal_mol"})
    joined["heavy_atoms"] = joined[spec.ligand_column].map(heavy)
    joined["ligand_efficiency"] = (
        -joined["docking_affinity_kcal_mol"] / joined["heavy_atoms"]
    )

    # Native control for each receptor, as the size-matched reference point.
    controls = best[best.is_control].copy()
    controls["heavy_atoms"] = controls[spec.ligand_column].map(heavy)
    controls["control_ligand_efficiency"] = (
        -controls[spec.affinity_column] / controls["heavy_atoms"]
    )
    control_best = (
        controls.sort_values("control_ligand_efficiency", ascending=False)
        .groupby("target_class")
        .head(1)
        .set_index("target_class")
    )
    joined["control_ligand"] = joined["target_class"].map(
        control_best[spec.ligand_column]
    )
    joined["control_affinity_kcal_mol"] = joined["target_class"].map(
        control_best[spec.affinity_column]
    )
    joined["control_ligand_efficiency"] = joined["target_class"].map(
        control_best["control_ligand_efficiency"]
    )
    joined["beats_control_on_affinity"] = (
        joined["docking_affinity_kcal_mol"] < joined["control_affinity_kcal_mol"]
    )
    joined["beats_control_on_efficiency"] = (
        joined["ligand_efficiency"] > joined["control_ligand_efficiency"]
    )

    columns = [
        spec.ligand_column,
        "target_class",
        spec.target_column,
        "docking_affinity_kcal_mol",
        "heavy_atoms",
        "ligand_efficiency",
        "control_ligand",
        "control_affinity_kcal_mol",
        "control_ligand_efficiency",
        "beats_control_on_affinity",
        "beats_control_on_efficiency",
        "overall_priority_score",
        "pipeline_rank",
        "n_target_classes_in_run",
        "confidence_class",
        "chemical_evidence_score",
        "species_transfer_score",
        "pocket_evidence_score",
    ]
    columns = [column for column in columns if column in joined.columns]
    return joined[columns].rename(columns={spec.ligand_column: "query_id"}).sort_values(
        ["query_id", "docking_affinity_kcal_mol"]
    )


def collapse_to_target_class(crosscheck: pd.DataFrame) -> pd.DataFrame:
    """One row per compound and target class, keeping the best affinity.

    A campaign may dock several receptor structures for the same target class
    (two KPC-2 crystal forms, say). Those share a single pipeline priority, so
    leaving them as separate rows would pseudo-replicate that value and inflate
    both the sample size and the correlation.
    """

    if crosscheck.empty:
        return crosscheck
    index = crosscheck.groupby(["query_id", "target_class"])[
        "docking_affinity_kcal_mol"
    ].idxmin()
    return crosscheck.loc[index].reset_index(drop=True)


def build_agreement(crosscheck: pd.DataFrame) -> pd.DataFrame:
    """Rank agreement between docking affinity and pipeline priority.

    Computed on one row per compound and target class; see
    ``collapse_to_target_class``.
    """

    from scipy.stats import spearmanr

    crosscheck = collapse_to_target_class(crosscheck)
    rows = []

    def _one(label: str, frame: pd.DataFrame) -> None:
        usable = frame.dropna(
            subset=["docking_affinity_kcal_mol", "overall_priority_score"]
        )
        n = len(usable)
        if n < 3 or usable["overall_priority_score"].nunique() < 2:
            rows.append(
                {
                    "scope": label,
                    "n_shared_targets": n,
                    "spearman_rho": np.nan,
                    "p_value": np.nan,
                    "status": "unavailable_too_few_comparable_targets",
                }
            )
            return
        rho, p = spearmanr(
            -usable["docking_affinity_kcal_mol"], usable["overall_priority_score"]
        )
        rows.append(
            {
                "scope": label,
                "n_shared_targets": n,
                "spearman_rho": float(rho),
                "p_value": float(p),
                "status": "available",
            }
        )

    _one("pooled", crosscheck)
    for compound, frame in crosscheck.groupby("query_id"):
        _one(str(compound), frame)
    return pd.DataFrame(rows)


def undocked_top_targets(
    predictions: pd.DataFrame,
    crosscheck: pd.DataFrame,
    spec: DockingSpec,
    *,
    top_n: int = 6,
) -> pd.DataFrame:
    """The pipeline's strongest hypotheses that the campaign never tested."""

    rows = predictions[predictions["organism"].astype(str) == spec.organism].copy()
    rows["overall_priority_score"] = pd.to_numeric(
        rows["overall_priority_score"], errors="coerce"
    )
    compounds = set(crosscheck["query_id"].astype(str))
    if compounds:
        rows = rows[rows["query_id"].astype(str).isin(compounds)]
    ranked = (
        rows.groupby("target_class")["overall_priority_score"]
        .mean()
        .sort_values(ascending=False)
        .head(top_n)
        .rename("mean_overall_priority_score")
        .reset_index()
    )
    docked = set(crosscheck["target_class"].astype(str))
    ranked["pipeline_rank"] = range(1, len(ranked) + 1)
    ranked["docked"] = ranked["target_class"].astype(str).isin(docked)
    return ranked


def _figure(crosscheck: pd.DataFrame, path: Path, organism: str) -> bool:
    """Docking affinity against pipeline priority, with the size control shown."""

    import matplotlib.pyplot as plt

    usable = collapse_to_target_class(crosscheck).dropna(
        subset=["docking_affinity_kcal_mol", "overall_priority_score"]
    )
    if usable.empty:
        return False

    style.apply_theme()
    compounds = sorted(usable["query_id"].astype(str).unique())
    palette = style.compound_palette(compounds)
    figure, axes = plt.subplots(1, 2, figsize=(13.5, 6.0))

    for compound in compounds:
        sub = usable[usable["query_id"].astype(str) == compound]
        axes[0].scatter(
            -sub["docking_affinity_kcal_mol"],
            sub["overall_priority_score"],
            s=90,
            color=palette[compound],
            edgecolor="white",
            linewidth=0.6,
            label=compound,
        )
        for row in sub.itertuples():
            label = str(row.target_class)
            if len(label) > 14:
                # Long ontology names collide at this density; the axis label
                # and the companion panel carry the full name.
                label = label.split("_")[0]
            axes[0].annotate(
                label,
                (-row.docking_affinity_kcal_mol, row.overall_priority_score),
                fontsize=7.5,
                xytext=(6, 4),
                textcoords="offset points",
            )
    axes[0].set_xlabel("Docking affinity (−kcal/mol, higher is stronger)")
    axes[0].set_ylabel("Pipeline overall priority score")
    axes[0].set_title("Shared targets: docking against pipeline")
    axes[0].legend(title="Compound", fontsize=8)

    efficiency = usable.dropna(subset=["ligand_efficiency"])
    if efficiency.empty:
        axes[1].axis("off")
    else:
        targets = list(dict.fromkeys(efficiency["target_class"].astype(str)))
        positions = np.arange(len(targets))
        width = 0.8 / (len(compounds) + 1)
        for index, compound in enumerate(compounds):
            sub = efficiency[efficiency["query_id"].astype(str) == compound]
            values = [
                float(sub.loc[sub.target_class == t, "ligand_efficiency"].max())
                if (sub.target_class == t).any()
                else np.nan
                for t in targets
            ]
            axes[1].bar(
                positions + index * width,
                values,
                width,
                color=palette[compound],
                label=compound,
            )
        control = [
            float(
                efficiency.loc[
                    efficiency.target_class == t, "control_ligand_efficiency"
                ].max()
            )
            for t in targets
        ]
        axes[1].bar(
            positions + len(compounds) * width,
            control,
            width,
            color="#7a869a",
            label="native control",
        )
        axes[1].set_xticks(
            positions + width * len(compounds) / 2, targets, rotation=30, ha="right",
            fontsize=8,
        )
        axes[1].set_ylabel("Ligand efficiency (−kcal/mol per heavy atom)")
        axes[1].set_title("Size-corrected: is the affinity good for this size?")
        axes[1].legend(fontsize=8)
        axes[1].grid(axis="x", visible=False)

    figure.suptitle(f"Docking campaign against pipeline ranking — {organism}")
    style.save_figure(figure, path)
    return True


def run_crosscheck(
    results_dir: Path,
    docking_summary: Path,
    spec: DockingSpec,
    *,
    heavy_atoms: dict[str, int] | None = None,
    figure_dirname: str = "figures_suite",
) -> dict[str, pd.DataFrame]:
    """Write the cross-check tables and figure; return them for inspection."""

    predictions_path = results_dir / PREDICTIONS
    if not predictions_path.is_file():
        raise DockingCrosscheckError(
            f"{predictions_path} is missing; run the pipeline before cross-checking"
        )
    if not Path(docking_summary).is_file():
        raise DockingCrosscheckError(f"docking summary is missing: {docking_summary}")

    predictions = pd.read_csv(predictions_path)
    docking = pd.read_csv(docking_summary)

    coverage = build_coverage(docking, spec, predictions)
    crosscheck = build_crosscheck(
        predictions, docking, spec, heavy_atoms=heavy_atoms
    )
    if crosscheck.empty:
        raise DockingCrosscheckError(
            "no docked receptor mapped onto a target class in this run; "
            "check docking.target_aliases in config.yaml"
        )
    agreement = build_agreement(crosscheck)
    undocked = undocked_top_targets(predictions, crosscheck, spec)

    coverage.to_csv(results_dir / COVERAGE_CSV, index=False, lineterminator="\n")
    crosscheck.to_csv(results_dir / CROSSCHECK_CSV, index=False, lineterminator="\n")
    agreement.to_csv(results_dir / AGREEMENT_CSV, index=False, lineterminator="\n")

    figure_path = results_dir / figure_dirname / "docking_pipeline_crosscheck.png"
    _figure(crosscheck, figure_path, spec.organism)

    return {
        "coverage": coverage,
        "crosscheck": crosscheck,
        "agreement": agreement,
        "undocked_top_targets": undocked,
    }


def load_heavy_atoms(path: Path) -> dict[str, int]:
    if not Path(path).is_file():
        return {}
    frame = pd.read_csv(path)
    if not {"ligand", "heavy_atoms"}.issubset(frame.columns):
        return {}
    return dict(zip(frame["ligand"].astype(str), frame["heavy_atoms"].astype(int)))


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    spec = DockingSpec.from_config(config)
    tables = run_crosscheck(
        results_dir,
        config.path_for("docking_summary"),
        spec,
        heavy_atoms=load_heavy_atoms(config.path_for("docking_ligand_heavy_atoms")),
    )

    print("=== receptor mapping ===")
    print(tables["coverage"].to_string(index=False))
    print("\n=== rank agreement (docking affinity vs pipeline priority) ===")
    print(tables["agreement"].to_string(index=False))
    print(f"\n=== pipeline top targets for {spec.organism}: docked? ===")
    print(tables["undocked_top_targets"].to_string(index=False))


if __name__ == "__main__":
    main()
