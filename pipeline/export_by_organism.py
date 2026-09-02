"""Split the organism-aware predictions into one CSV set per organism.

Two views per organism, written to ``results/by_organism/<organism>/``:

``targets_predicted.csv``
    Every compound-target hypothesis scored against that organism, best first.

``target_summary.csv``
    One row per target class: how each target fared across the whole compound
    series in that organism.

Nothing is recomputed here. This is a projection of the run tables, so the
numbers match ``v2_open_target_predictions_by_organism.csv`` exactly.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

try:  # pragma: no cover - import shim for direct module execution
    from pipeline.config import load_config
    from pipeline.compound_manifest import CompoundManifest, load_from_config
    from pipeline.figures_suite import organism_slug
except ModuleNotFoundError:  # pragma: no cover
    from config import load_config
    from compound_manifest import CompoundManifest, load_from_config
    from figures_suite import organism_slug


PREDICTIONS = "v2_open_target_predictions_by_organism.csv"

#: Columns carried into the per-organism prediction export, in reading order:
#: identity, then the evidence layers, then the context that justifies them.
EXPORT_COLUMNS: tuple[str, ...] = (
    "query_id",
    "manifest_assigned",
    "organism",
    "target_class",
    "parent_target_class",
    "target_subtype",
    "binding_site_or_mechanism",
    "overall_priority_score",
    "confidence_class",
    "chemical_evidence_score",
    "chemical_quality_adjusted_score",
    "ecfp4_max",
    "maccs_max",
    "target_specificity_score",
    "reference_quality_grade",
    "species_transfer_score",
    "sequence_mapping_status",
    "organism_transfer_source",
    "pocket_evidence_score",
    "rcsb_structure_candidate_count",
    "rcsb_co_crystal_ligand_count",
    "biological_priority_score",
    "clinical_translation_score",
    "organism_scope_score",
    "essentiality_score",
    "cellular_access_score",
    "resistance_relevance_score",
    "card_model_count",
    "card_snp_row_count",
    "organism_specific_snp_row_count",
    "anti_target_risk_score",
    "anti_target_evidence_status",
    "uncertainty_reasons",
    "recommended_validation",
    "best_reference_molecule",
    "best_reference_organism",
    "clinical_status",
    "target_role",
    "organism_scope",
    "cellular_localization",
    "resistance_relevance",
)


def _prepare(frame: pd.DataFrame, manifest: CompoundManifest) -> pd.DataFrame:
    prepared = frame.copy()
    prepared["overall_priority_score"] = pd.to_numeric(
        prepared.get("overall_priority_score"), errors="coerce"
    )
    assignments = manifest.assignments if manifest is not None else {}
    prepared["manifest_assigned"] = [
        assignments.get(str(compound)) == str(organism)
        for compound, organism in zip(prepared["query_id"], prepared["organism"])
    ]
    return prepared


def _summarize(frame: pd.DataFrame) -> pd.DataFrame:
    """Per target class: how the whole series scored against it."""

    best_rows = frame.loc[
        frame.groupby("target_class")["overall_priority_score"].idxmax()
    ].set_index("target_class")

    grouped = frame.groupby("target_class")
    summary = pd.DataFrame(
        {
            "n_compound_hypotheses": grouped["query_id"].size(),
            "max_overall_priority_score": grouped["overall_priority_score"].max(),
            "mean_overall_priority_score": grouped["overall_priority_score"].mean(),
            "best_compound": best_rows["query_id"],
            "best_compound_is_manifest_assigned": best_rows["manifest_assigned"],
            "best_compound_confidence": best_rows["confidence_class"],
        }
    )
    if "confidence_class" in frame:
        counts = (
            pd.crosstab(frame["target_class"], frame["confidence_class"])
            .add_prefix("n_")
            .reindex(summary.index, fill_value=0)
        )
        summary = summary.join(counts)
    for column in (
        "species_transfer_score",
        "pocket_evidence_score",
        "biological_priority_score",
        "clinical_translation_score",
    ):
        if column in best_rows:
            summary[column] = best_rows[column]
    for column in ("sequence_mapping_status", "target_role", "clinical_status"):
        if column in best_rows:
            summary[column] = best_rows[column]
    return summary.sort_values(
        "max_overall_priority_score", ascending=False
    ).reset_index()


def export_by_organism(
    results_dir: Path,
    *,
    manifest: CompoundManifest | None = None,
    dirname: str = "by_organism",
) -> pd.DataFrame:
    """Write one CSV set per organism; return an index of what was written."""

    source = results_dir / PREDICTIONS
    if not source.is_file():
        raise FileNotFoundError(
            f"{source} is missing; run the pipeline before exporting per-organism tables"
        )
    frame = _prepare(pd.read_csv(source), manifest)
    if "organism" not in frame.columns:
        raise ValueError(f"{source} has no organism column")

    root = results_dir / dirname
    root.mkdir(parents=True, exist_ok=True)

    index: list[dict[str, object]] = []
    for organism in sorted(frame["organism"].dropna().astype(str).unique()):
        subset = frame[frame["organism"].astype(str) == organism].copy()
        subset = subset.sort_values(
            ["overall_priority_score", "query_id"], ascending=[False, True]
        )
        columns = [column for column in EXPORT_COLUMNS if column in subset.columns]
        directory = root / organism_slug(organism)
        directory.mkdir(parents=True, exist_ok=True)

        predictions_path = directory / "targets_predicted.csv"
        subset[columns].to_csv(predictions_path, index=False, lineterminator="\n")

        summary = _summarize(subset)
        summary_path = directory / "target_summary.csv"
        summary.to_csv(summary_path, index=False, lineterminator="\n")

        assigned = sorted(
            subset.loc[subset["manifest_assigned"], "query_id"].astype(str).unique()
        )
        index.append(
            {
                "organism": organism,
                "n_rows": len(subset),
                "n_compounds": subset["query_id"].nunique(),
                "n_target_classes": subset["target_class"].nunique(),
                "manifest_assigned_compounds": ";".join(assigned),
                "targets_predicted_csv": str(predictions_path),
                "target_summary_csv": str(summary_path),
            }
        )

    index_frame = pd.DataFrame(index)
    index_frame.to_csv(root / "index.csv", index=False, lineterminator="\n")
    return index_frame


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    index = export_by_organism(results_dir, manifest=load_from_config(config))
    print(
        index[
            [
                "organism",
                "n_rows",
                "n_compounds",
                "n_target_classes",
                "manifest_assigned_compounds",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
