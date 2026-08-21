"""Strict ingestion for externally sourced property-matched benchmark decoys.

Cross-target reference ligands are deliberately excluded from this interface:
they support the specificity margin but are not experimentally confirmed inactive.
Likewise, DUD-E-style decoys are labelled as presumed, not confirmed, inactive.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem


REQUIRED_DECOY_COLUMNS = (
    "decoy_id",
    "canonical_smiles",
    "matched_target_class",
    "source_dataset",
    "source_version",
    "source_record_id",
    "property_matching_method",
)


@dataclass(frozen=True)
class DecoyLoadResult:
    decoys: pd.DataFrame
    status: dict[str, object]


def _empty_decoys() -> pd.DataFrame:
    return pd.DataFrame(columns=[*REQUIRED_DECOY_COLUMNS, "canonical_smiles_rdkit"])


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
    source_datasets = sorted(decoys["source_dataset"].astype(str).unique())
    source_versions = sorted(decoys["source_version"].astype(str).unique())
    return DecoyLoadResult(
        decoys,
        {
            **base_status,
            "status": "available",
            "status_reason": "versioned property-matched decoys validated",
            "n_decoys": len(decoys),
            "source_datasets": ";".join(source_datasets),
            "source_versions": ";".join(source_versions),
        },
    )


def build_active_decoy_candidates(
    active_queries: pd.DataFrame, decoy_result: DecoyLoadResult
) -> pd.DataFrame:
    """Create auditable candidate labels for an active-vs-decoy retrieval table."""

    required_active = {"query_id", "canonical_smiles", "target_class"}
    missing = sorted(required_active - set(active_queries.columns))
    if missing:
        raise ValueError(f"Active benchmark table is missing columns: {missing}")
    active = active_queries.copy()
    active["candidate_id"] = active["query_id"].astype(str)
    active["matched_target_class"] = active["target_class"].astype(str)
    active["is_active"] = 1
    active["candidate_type"] = "curated_benchmark_active"
    active["label_semantics"] = "curated_known_mechanism_active"
    active["decoy_source_dataset"] = pd.NA
    active["decoy_source_version"] = pd.NA
    active["decoy_source_record_id"] = pd.NA
    active["property_matching_method"] = pd.NA

    decoys = decoy_result.decoys.copy()
    if decoys.empty:
        return active
    decoys["candidate_id"] = decoys["decoy_id"].astype(str)
    decoys["query_id"] = decoys["candidate_id"]
    decoys["target_class"] = decoys["matched_target_class"].astype(str)
    decoys["is_active"] = 0
    decoys["candidate_type"] = "property_matched_decoy"
    decoys["label_semantics"] = (
        "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
    )
    decoys["decoy_source_dataset"] = decoys["source_dataset"]
    decoys["decoy_source_version"] = decoys["source_version"]
    decoys["decoy_source_record_id"] = decoys["source_record_id"]

    columns = [
        "candidate_id",
        "query_id",
        "canonical_smiles",
        "target_class",
        "matched_target_class",
        "is_active",
        "candidate_type",
        "label_semantics",
        "decoy_source_dataset",
        "decoy_source_version",
        "decoy_source_record_id",
        "property_matching_method",
    ]
    return pd.concat([active[columns], decoys[columns]], ignore_index=True)
