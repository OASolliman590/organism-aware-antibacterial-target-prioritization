"""Strict ingestion for externally sourced property-matched benchmark decoys.

Cross-target reference ligands are deliberately excluded from this interface:
they support the specificity margin but are not experimentally confirmed inactive.
Likewise, DUD-E-style decoys are labelled as presumed, not confirmed, inactive.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

import pandas as pd

_RDKIT_IMPORT_ERROR: ImportError | None = None
try:
    from rdkit import Chem
except ImportError as error:  # missing-artifact status remains available
    Chem = None
    _RDKIT_IMPORT_ERROR = error


REQUIRED_DECOY_COLUMNS = (
    "decoy_id",
    "canonical_smiles",
    "matched_active_query_id",
    "matched_target_class",
    "source_dataset",
    "source_version",
    "source_record_id",
    "property_matching_method",
)

CANDIDATE_DATE_COLUMNS = (
    "activity_date",
    "measurement_date",
    "publication_date",
    "document_date",
    "activity_year",
    "publication_year",
    "document_year",
    "year",
)

CANDIDATE_COLUMNS = (
    "candidate_id",
    "source_candidate_id",
    "query_id",
    "canonical_smiles",
    "canonical_smiles_source",
    "source_target_label",
    "target_class",
    "matched_target_class",
    "matched_active_query_id",
    "retrieval_task_id",
    "target_mapping_status",
    "is_active",
    "candidate_type",
    "label_semantics",
    "decoy_source_dataset",
    "decoy_source_version",
    "decoy_source_record_id",
    "property_matching_method",
    *CANDIDATE_DATE_COLUMNS,
)


@dataclass(frozen=True)
class DecoyLoadResult:
    decoys: pd.DataFrame
    status: dict[str, object]


def _empty_decoys() -> pd.DataFrame:
    return pd.DataFrame(columns=[*REQUIRED_DECOY_COLUMNS, "canonical_smiles_rdkit"])


def _require_rdkit() -> None:
    if _RDKIT_IMPORT_ERROR is not None:
        raise RuntimeError("RDKit is required to validate decoy structures") from (
            _RDKIT_IMPORT_ERROR
        )


def load_property_matched_decoys(path: str | Path) -> DecoyLoadResult:
    """Load a versioned decoy table or return a machine-readable pending gap."""

    selected = Path(path).resolve()
    base_status: dict[str, object] = {
        "dataset_path": str(selected),
        "decoy_requirement": "DUD-E-style property-matched presumed inactives",
        "negative_label_semantics": (
            "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
        ),
        "cross_target_decoy_policy": "specificity_margin_only_not_inactive",
    }
    if not selected.is_file():
        return DecoyLoadResult(
            _empty_decoys(),
            {
                **base_status,
                "status": "pending_missing_property_matched_decoy_dataset",
                "status_reason": (
                    "No versioned property-matched decoy artifact exists in the pinned "
                    "snapshot; cross-target ligands were not relabelled as inactive."
                ),
                "n_decoys": 0,
            },
        )

    decoys = pd.read_csv(selected)
    missing = sorted(set(REQUIRED_DECOY_COLUMNS) - set(decoys.columns))
    if missing:
        raise ValueError(f"Property-matched decoy table is missing columns: {missing}")
    if decoys.empty:
        raise ValueError("Property-matched decoy table is empty")
    for column in REQUIRED_DECOY_COLUMNS:
        if decoys[column].isna().any() or decoys[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Property-matched decoy provenance is blank: {column}")
    if decoys["decoy_id"].duplicated().any():
        raise ValueError("Property-matched decoy_id values must be unique")
    if decoys["source_version"].astype(str).str.lower().isin(
        {"unknown", "unrecorded", "none", "nan"}
    ).any():
        raise ValueError("Every decoy requires a recorded source_version")

    _require_rdkit()
    canonical: list[str] = []
    invalid: list[str] = []
    for _, row in decoys.iterrows():
        molecule = Chem.MolFromSmiles(str(row["canonical_smiles"]))
        if molecule is None:
            invalid.append(str(row["decoy_id"]))
        else:
            canonical.append(Chem.MolToSmiles(molecule, canonical=True))
    if invalid:
        raise ValueError(f"Invalid decoy structures: {invalid[:5]}")
    decoys = decoys.copy()
    decoys["canonical_smiles_rdkit"] = canonical
    duplicate_structures = decoys.duplicated(
        [
            "matched_active_query_id",
            "matched_target_class",
            "canonical_smiles_rdkit",
        ],
        keep=False,
    )
    if duplicate_structures.any():
        examples = decoys.loc[
            duplicate_structures,
            ["decoy_id", "matched_active_query_id", "matched_target_class"],
        ].head(5)
        raise ValueError(
            "Duplicate decoy structures within an active/target retrieval task "
            "are not permitted: "
            f"{examples.to_dict('records')}"
        )
    source_datasets = sorted(decoys["source_dataset"].astype(str).unique())
    source_versions = sorted(decoys["source_version"].astype(str).unique())
    return DecoyLoadResult(
        decoys,
        {
            **base_status,
            "status": "available",
            "status_reason": "versioned property-matched decoys validated",
            "n_decoys": len(decoys),
            "n_matched_active_queries": int(
                decoys["matched_active_query_id"].nunique()
            ),
            "n_matched_target_classes": int(decoys["matched_target_class"].nunique()),
            "artifact_sha256": hashlib.sha256(selected.read_bytes()).hexdigest(),
            "source_datasets": ";".join(source_datasets),
            "source_versions": ";".join(source_versions),
        },
    )


def build_active_decoy_candidates(
    active_queries: pd.DataFrame,
    decoy_result: DecoyLoadResult,
    *,
    classes_by_alias: dict[str, set[str]] | None = None,
    valid_target_classes: set[str] | None = None,
    unmapped_active_sink: list[dict[str, object]] | None = None,
) -> pd.DataFrame:
    """Create target-specific candidates without inventing negative labels."""

    required_active = {"query_id", "canonical_smiles", "target_class"}
    missing = sorted(required_active - set(active_queries.columns))
    if missing:
        raise ValueError(f"Active benchmark table is missing columns: {missing}")
    for column in required_active:
        if (
            active_queries[column].isna().any()
            or active_queries[column].astype(str).str.strip().eq("").any()
        ):
            raise ValueError(f"Active benchmark provenance is blank: {column}")
    if active_queries["query_id"].astype(str).duplicated().any():
        raise ValueError("Active benchmark query_id values must be unique")
    _require_rdkit()
    mappings = classes_by_alias or {}
    active_rows: list[dict[str, object]] = []
    active_supported_targets: dict[str, set[str]] = {}
    unmapped_gaps_by_active: dict[str, dict[str, object]] = {}
    all_active_structures: set[str] = set()
    for _, row in active_queries.iterrows():
        source_query_id = str(row["query_id"])
        label = str(row["target_class"])
        molecule = Chem.MolFromSmiles(str(row["canonical_smiles"]))
        if molecule is None:
            raise ValueError(f"Invalid active benchmark structure: {source_query_id}")
        canonical = Chem.MolToSmiles(molecule, canonical=True)
        all_active_structures.add(canonical)
        matched_classes = sorted(mappings.get(label, {label}))
        if valid_target_classes is not None:
            matched_classes = [
                target for target in matched_classes if target in valid_target_classes
            ]
        if not matched_classes:
            if unmapped_active_sink is not None:
                gap: dict[str, object] = {
                    "source_candidate_id": source_query_id,
                    "source_target_label": label,
                    "mapping_status": "pending_no_pinned_ontology_mapping",
                    "status_reason": (
                        "Active benchmark label has no exact or declared alias "
                        "mapping in the pinned target ontology"
                    ),
                    "n_linked_decoys_excluded": 0,
                    "linked_decoy_ids_sha256": "",
                }
                unmapped_active_sink.append(gap)
                unmapped_gaps_by_active[source_query_id] = gap
            continue
        active_supported_targets[source_query_id] = set(matched_classes)
        for matched_target_class in matched_classes:
            active_rows.append(
                {
                    **row.to_dict(),
                    "source_candidate_id": source_query_id,
                    "candidate_id": f"{source_query_id}::{matched_target_class}",
                    "canonical_smiles_source": str(row["canonical_smiles"]),
                    "canonical_smiles": canonical,
                    "source_target_label": label,
                    "target_class": matched_target_class,
                    "matched_target_class": matched_target_class,
                    "matched_active_query_id": source_query_id,
                    "retrieval_task_id": (
                        f"{source_query_id}::{matched_target_class}"
                    ),
                    "target_mapping_status": "pinned_ontology_exact_or_declared_alias",
                }
            )
    decoys = decoy_result.decoys.copy()
    if decoys.empty:
        return pd.DataFrame(columns=list(CANDIDATE_COLUMNS))
    active_query_ids = set(active_queries["query_id"].astype(str))
    unknown_active_ids = sorted(
        set(decoys["matched_active_query_id"].astype(str)) - active_query_ids
    )
    if unknown_active_ids:
        raise ValueError(
            "Property-matched decoys reference unknown active query_id values: "
            f"{unknown_active_ids[:5]}"
        )
    linked_to_unmapped = ~decoys["matched_active_query_id"].astype(str).isin(
        active_supported_targets
    )
    if linked_to_unmapped.any():
        if unmapped_active_sink is None:
            raise ValueError(
                "Property-matched decoys linked to unmapped active queries require "
                "an unmapped_active_sink so exclusions cannot be silent"
            )
        for active_query_id, group in decoys[linked_to_unmapped].groupby(
            decoys.loc[linked_to_unmapped, "matched_active_query_id"].astype(str),
            sort=True,
        ):
            decoy_ids = sorted(group["decoy_id"].astype(str))
            gap = unmapped_gaps_by_active[str(active_query_id)]
            gap["n_linked_decoys_excluded"] = len(decoy_ids)
            gap["linked_decoy_ids_sha256"] = hashlib.sha256(
                "\n".join(decoy_ids).encode("utf-8")
            ).hexdigest()
        decoys = decoys.loc[~linked_to_unmapped].copy()

    decoy_targets = set(decoys["matched_target_class"].astype(str))
    if valid_target_classes is not None:
        unmapped = sorted(decoy_targets - valid_target_classes)
        if unmapped:
            raise ValueError(
                f"Decoy matched_target_class values are not in the pinned ontology: {unmapped}"
            )
    unsupported_pairs = []
    for _, row in decoys.iterrows():
        active_query_id = str(row["matched_active_query_id"])
        target_class = str(row["matched_target_class"])
        if target_class not in active_supported_targets.get(active_query_id, set()):
            unsupported_pairs.append(
                {
                    "decoy_id": str(row["decoy_id"]),
                    "matched_active_query_id": active_query_id,
                    "matched_target_class": target_class,
                }
            )
    if unsupported_pairs:
        raise ValueError(
            "Property-matched decoy active/target links are unsupported by the "
            f"pinned ontology mapping: {unsupported_pairs[:5]}"
        )
    supported_pairs = set(
        zip(
            decoys["matched_active_query_id"].astype(str),
            decoys["matched_target_class"].astype(str),
        )
    )
    if active_rows:
        active = pd.DataFrame(active_rows)
        active = active[
            active.apply(
                lambda row: (
                    str(row["matched_active_query_id"]),
                    str(row["matched_target_class"]),
                )
                in supported_pairs,
                axis=1,
            )
        ].copy()
        active["query_id"] = active["candidate_id"]
        active["is_active"] = 1
        active["candidate_type"] = "curated_benchmark_active"
        active["label_semantics"] = "curated_known_mechanism_active"
        active["decoy_source_dataset"] = pd.NA
        active["decoy_source_version"] = pd.NA
        active["decoy_source_record_id"] = pd.NA
        active["property_matching_method"] = pd.NA
        for column in CANDIDATE_DATE_COLUMNS:
            if column not in active:
                active[column] = pd.NA
        active = active[list(CANDIDATE_COLUMNS)]
    else:
        active = pd.DataFrame(columns=list(CANDIDATE_COLUMNS))
    overlaps = sorted(
        set(decoys["canonical_smiles_rdkit"].astype(str))
        & all_active_structures
    )
    if overlaps:
        raise ValueError(
            "Property-matched decoys overlap exact active structures; "
            f"n_overlaps={len(overlaps)}"
        )
    decoys["source_candidate_id"] = decoys["decoy_id"].astype(str)
    decoys["candidate_id"] = (
        decoys["decoy_id"].astype(str)
        + "::"
        + decoys["matched_target_class"].astype(str)
    )
    collisions = sorted(set(active["candidate_id"]) & set(decoys["candidate_id"]))
    if collisions:
        raise ValueError(f"Active/decoy candidate IDs collide: {collisions[:5]}")
    decoys["query_id"] = decoys["candidate_id"]
    decoys["target_class"] = decoys["matched_target_class"].astype(str)
    decoys["retrieval_task_id"] = (
        decoys["matched_active_query_id"].astype(str)
        + "::"
        + decoys["matched_target_class"].astype(str)
    )
    decoys["source_target_label"] = decoys["matched_target_class"].astype(str)
    decoys["canonical_smiles_source"] = decoys["canonical_smiles"].astype(str)
    decoys["canonical_smiles"] = decoys["canonical_smiles_rdkit"].astype(str)
    decoys["target_mapping_status"] = "exact_pinned_ontology_target_class"
    decoys["is_active"] = 0
    decoys["candidate_type"] = "property_matched_decoy"
    decoys["label_semantics"] = (
        "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
    )
    decoys["decoy_source_dataset"] = decoys["source_dataset"]
    decoys["decoy_source_version"] = decoys["source_version"]
    decoys["decoy_source_record_id"] = decoys["source_record_id"]
    for column in CANDIDATE_DATE_COLUMNS:
        decoys[column] = pd.NA

    return pd.concat(
        [
            active,
            decoys[list(CANDIDATE_COLUMNS)],
        ],
        ignore_index=True,
    )
