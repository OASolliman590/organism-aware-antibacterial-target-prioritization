"""Leakage-controlled v3 benchmark splits and provenance.

Missing family or time metadata is never inferred. Affected split/query pairs are
returned with a pending status and zero evaluable references so downstream code
cannot accidentally report them as a valid benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles

try:
    from pipeline.config import ProjectConfig
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig


SUPPORTED_SPLITS = ("target_family", "scaffold", "temporal")
DATE_FIELDS = (
    "activity_date",
    "measurement_date",
    "publication_date",
    "document_date",
    "activity_year",
    "publication_year",
    "document_year",
    "year",
)


@dataclass(frozen=True)
class SplitResult:
    query_id: str
    split_type: str
    references: dict[str, list[dict[str, Any]]]
    provenance: dict[str, Any]


def _molecule(record: dict[str, Any]) -> Chem.Mol | None:
    existing = record.get("_mol") or record.get("mol")
    if existing is not None:
        return Chem.RemoveHs(Chem.Mol(existing))
    smiles = record.get("canonical_smiles_standardized") or record.get(
        "canonical_smiles"
    )
    molecule = Chem.MolFromSmiles(str(smiles)) if smiles else None
    return Chem.RemoveHs(molecule) if molecule is not None else None


def _fingerprint(record: dict[str, Any], config: ProjectConfig):
    existing = record.get("_fp")
    if existing is not None:
        return existing
    molecule = _molecule(record)
    if molecule is None:
        return None
    return AllChem.GetMorganGenerator(
        radius=int(config.value("chem2d.fingerprint_radius")),
        fpSize=int(config.value("chem2d.fingerprint_bits")),
    ).GetFingerprint(molecule)


def _scaffold(record: dict[str, Any]) -> str | None:
    if "_scaffold" in record:
        value = str(record["_scaffold"] or "")
        return value or None
    molecule = _molecule(record)
    if molecule is None:
        return None
    value = MurckoScaffoldSmiles(mol=molecule)
    return value or None


def _record_date(record: dict[str, Any]) -> tuple[pd.Timestamp | None, str | None]:
    for field in DATE_FIELDS:
        value = record.get(field)
        if value is None or (isinstance(value, float) and np.isnan(value)):
            continue
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        if field.endswith("year") or field == "year":
            try:
                year = int(float(text))
            except ValueError:
                continue
            if 1000 <= year <= 9999:
                return pd.Timestamp(year=year, month=1, day=1), field
            continue
        parsed = pd.to_datetime(text, errors="coerce")
        if pd.notna(parsed):
            return pd.Timestamp(parsed), field
    return None, None


def ontology_family_maps(
    ontology: pd.DataFrame,
) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Map target classes and exact declared benchmark aliases to families."""

    required = {"target_class", "target_family"}
    if not required.issubset(ontology.columns):
        raise ValueError(f"Ontology must contain {sorted(required)}")
    family_by_class: dict[str, str] = {}
    classes_by_alias: dict[str, set[str]] = {}
    for _, row in ontology.iterrows():
        target_class = str(row["target_class"])
        family = row["target_family"]
        if pd.notna(family) and str(family).strip():
            family_by_class[target_class] = str(family).strip()
        classes_by_alias.setdefault(target_class, set()).add(target_class)
        aliases = row.get("benchmark_aliases")
        if pd.notna(aliases):
            for alias in str(aliases).split(";"):
                alias = alias.strip()
                if alias:
                    classes_by_alias.setdefault(alias, set()).add(target_class)
        parent = row.get("parent_target_class")
        if pd.notna(parent) and str(parent).strip():
            classes_by_alias.setdefault(str(parent).strip(), set()).add(target_class)
    return family_by_class, classes_by_alias


def _flatten_references(
    references: dict[str, list[dict[str, Any]]]
) -> list[tuple[str, dict[str, Any]]]:
    return [
        (target_class, record)
        for target_class in sorted(references)
        for record in references[target_class]
    ]


def make_split(
    query: dict[str, Any],
    references: dict[str, list[dict[str, Any]]],
    split_type: str,
    *,
    family_by_class: dict[str, str],
    accepted_query_classes: set[str],
    config: ProjectConfig,
) -> SplitResult:
    """Filter one query's references and assert the analogue leakage guard."""

    if split_type not in SUPPORTED_SPLITS:
        raise ValueError(f"Unsupported split: {split_type}")
    query_id = str(query.get("query_id") or query.get("drug") or "")
    if not query_id:
        raise ValueError("query requires query_id or drug")
    query_molecule = _molecule(query)
    query_fp = _fingerprint(query, config)
    if query_molecule is None or query_fp is None:
        raise ValueError(f"Query {query_id} has no valid molecular structure")
    query_scaffold = _scaffold(query)
    threshold = float(config.value("benchmark.analogue_exclusion_threshold"))
    cutoff = pd.Timestamp(str(config.value("benchmark.time_cutoff")))
    query_date, query_date_field = _record_date(query)
    query_families = {
        family_by_class[target_class]
        for target_class in accepted_query_classes
        if target_class in family_by_class
    }

    counts = {
        "n_references_input": 0,
        "n_removed_close_analogue": 0,
        "n_removed_same_scaffold": 0,
        "n_removed_target_family": 0,
        "n_removed_post_cutoff": 0,
        "n_removed_missing_date": 0,
        "n_removed_invalid_structure": 0,
    }
    kept: dict[str, list[dict[str, Any]]] = {}
    for target_class, record in _flatten_references(references):
        counts["n_references_input"] += 1
        reference_fp = _fingerprint(record, config)
        if reference_fp is None:
            counts["n_removed_invalid_structure"] += 1
            continue
        similarity = float(DataStructs.TanimotoSimilarity(query_fp, reference_fp))
        close = similarity >= threshold
        same_scaffold = bool(
            query_scaffold
            and _scaffold(record)
            and query_scaffold == _scaffold(record)
        )
        target_family = family_by_class.get(target_class)
        same_family = bool(target_family and target_family in query_families)
        reference_date, _ = _record_date(record)
        post_cutoff = bool(reference_date is not None and reference_date >= cutoff)
        missing_date = reference_date is None

        remove = close
        if close:
            counts["n_removed_close_analogue"] += 1
        if split_type == "scaffold" and same_scaffold:
            remove = True
            counts["n_removed_same_scaffold"] += 1
        elif split_type == "target_family" and same_family:
            remove = True
            counts["n_removed_target_family"] += 1
        elif split_type == "temporal":
            if missing_date:
                remove = True
                counts["n_removed_missing_date"] += 1
            elif post_cutoff:
                remove = True
                counts["n_removed_post_cutoff"] += 1
        if not remove:
            kept.setdefault(target_class, []).append(record)

    status = "available"
    reason = "split filters and analogue leakage guard applied"
    if split_type == "target_family" and not query_families:
        status = "pending_missing_target_family"
        reason = "query target label has no target-family mapping in the pinned ontology"
        kept = {}
    elif split_type == "temporal" and query_date is None:
        status = "pending_missing_query_date"
        reason = "query has no dated activity/measurement metadata for temporal assignment"
        kept = {}
    elif split_type == "temporal" and query_date < cutoff:
        status = "not_in_post_cutoff_test_set"
        reason = "query date precedes the configured temporal test cutoff"
        kept = {}

    remaining_similarities = []
    for _, record in _flatten_references(kept):
        reference_fp = _fingerprint(record, config)
        if reference_fp is not None:
            remaining_similarities.append(
                float(DataStructs.TanimotoSimilarity(query_fp, reference_fp))
            )
    if any(value >= threshold for value in remaining_similarities):
        raise AssertionError(
            f"Analogue leakage guard failed for {query_id}/{split_type}"
        )

    provenance = {
        "query_id": query_id,
        "split_type": split_type,
        "status": status,
        "status_reason": reason,
        "query_target_classes": ";".join(sorted(accepted_query_classes)),
        "query_target_families": ";".join(sorted(query_families)),
        "query_date_field": query_date_field,
        "query_date_value": query_date.isoformat() if query_date is not None else None,
        "time_cutoff": cutoff.isoformat(),
        "analogue_exclusion_threshold": threshold,
        "query_murcko_scaffold": query_scaffold,
        **counts,
        "n_references_after_split": sum(len(rows) for rows in kept.values()),
        "n_target_classes_after_split": len(kept),
        "max_remaining_query_reference_tanimoto": (
            max(remaining_similarities) if remaining_similarities else np.nan
        ),
        "analogue_leakage_guard_passed": not any(
            value >= threshold for value in remaining_similarities
        ),
    }
    return SplitResult(query_id, split_type, kept, provenance)


def generate_splits(
    queries: list[dict[str, Any]],
    references: dict[str, list[dict[str, Any]]],
    ontology: pd.DataFrame,
    config: ProjectConfig,
) -> tuple[list[SplitResult], pd.DataFrame]:
    """Generate every configured split and its row-level provenance table."""

    family_by_class, classes_by_alias = ontology_family_maps(ontology)
    results: list[SplitResult] = []
    for split_type in config.value("benchmark.splits"):
        for query in queries:
            label = str(query.get("target_class") or query.get("query_target_label") or "")
            accepted = set(classes_by_alias.get(label, {label} if label else set()))
            results.append(
                make_split(
                    query,
                    references,
                    str(split_type),
                    family_by_class=family_by_class,
                    accepted_query_classes=accepted,
                    config=config,
                )
            )
    provenance = pd.DataFrame([result.provenance for result in results])
    return results, provenance
