"""Leakage-controlled v3 benchmark splits and provenance.

Missing family or time metadata is never inferred. Affected split/query pairs are
returned with a pending status and zero evaluable references so downstream code
cannot accidentally report them as a valid benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

_RDKIT_IMPORT_ERROR: ImportError | None = None
try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
except ImportError as error:  # pure metric/fusion helpers remain testable
    Chem = None
    DataStructs = None
    AllChem = None
    MurckoScaffoldSmiles = None
    _RDKIT_IMPORT_ERROR = error

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


def _require_rdkit() -> None:
    if _RDKIT_IMPORT_ERROR is not None:
        raise RuntimeError(
            "RDKit is required for benchmark structure and split operations"
        ) from _RDKIT_IMPORT_ERROR


def _molecule(record: dict[str, Any]) -> Chem.Mol | None:
    _require_rdkit()
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


def _bedroc(labels: np.ndarray, scores: np.ndarray, alpha: float) -> float:
    """Tie-aware BEDROC using average score ranks (rank 1 is best)."""

    n_molecules = len(labels)
    n_actives = int(labels.sum())
    if n_molecules == 0 or n_actives == 0 or n_actives == n_molecules:
        return np.nan
    ranks = pd.Series(scores).rank(method="average", ascending=False).to_numpy()
    denominator = (1.0 / n_molecules) * (
        (1.0 - math.exp(-alpha)) / (math.exp(alpha / n_molecules) - 1.0)
    )
    rie = float(np.exp(-(alpha * ranks[labels == 1]) / n_molecules).sum()) / (
        n_actives * denominator
    )
    ratio = n_actives / n_molecules
    rie_max = (1.0 - math.exp(-alpha * ratio)) / (
        ratio * (1.0 - math.exp(-alpha))
    )
    rie_min = (1.0 - math.exp(alpha * ratio)) / (
        ratio * (1.0 - math.exp(alpha))
    )
    return float((rie - rie_min) / (rie_max - rie_min))


def _tie_aware_enrichment(
    labels: np.ndarray, scores: np.ndarray, fraction: float
) -> float:
    n_molecules = len(labels)
    n_actives = int(labels.sum())
    if n_molecules == 0 or n_actives == 0 or n_actives == n_molecules:
        return np.nan
    top_n = max(1, int(math.ceil(n_molecules * fraction)))
    ordered = np.sort(scores)[::-1]
    boundary = ordered[top_n - 1]
    above = scores > boundary
    tied = scores == boundary
    slots_at_boundary = top_n - int(above.sum())
    expected_boundary_hits = (
        slots_at_boundary * float(labels[tied].mean()) if tied.any() else 0.0
    )
    expected_hits = float(labels[above].sum()) + expected_boundary_hits
    return float((expected_hits / top_n) / (n_actives / n_molecules))


def query_level_metrics(
    scores: pd.DataFrame,
    *,
    score_column: str,
    label_column: str = "is_active",
    bedroc_alphas: tuple[float, ...] = (20.0, 80.5),
    enrichment_fractions: tuple[float, ...] = (0.01, 0.05),
) -> pd.DataFrame:
    """Compute retrieval metrics per query before any cross-query aggregation."""

    required = {"query_id", "target_class", score_column, label_column}
    missing = sorted(required - set(scores.columns))
    if missing:
        raise ValueError(f"Benchmark scores are missing columns: {missing}")
    metric_names = [
        "auroc",
        *[
            f"bedroc_alpha_{str(float(alpha)).replace('.', '_')}"
            for alpha in bedroc_alphas
        ],
        *[
            f"ef_{int(round(float(fraction) * 100))}pct"
            for fraction in enrichment_fractions
        ],
        "mrr",
        "coverage",
    ]
    rows: list[dict[str, Any]] = []
    for query_id, group in scores.groupby("query_id", sort=True):
        numeric_scores = pd.to_numeric(group[score_column], errors="coerce")
        labels_all = pd.to_numeric(group[label_column], errors="coerce")
        valid = np.isfinite(numeric_scores) & labels_all.isin([0, 1])
        labels = labels_all[valid].astype(int).to_numpy()
        values = numeric_scores[valid].astype(float).to_numpy()
        n_positive = int(labels.sum())
        n_negative = int(len(labels) - n_positive)
        covered = int(n_positive > 0)
        if covered:
            ranks = pd.Series(values).rank(
                method="average", ascending=False
            ).to_numpy()
            reciprocal_rank = float(1.0 / ranks[labels == 1].min())
        else:
            reciprocal_rank = 0.0
        row: dict[str, Any] = {
            "query_id": query_id,
            "n_candidates": len(labels),
            "n_positive": n_positive,
            "n_negative": n_negative,
            "coverage": float(covered),
            "mrr": reciprocal_rank,
            "auroc": (
                float(roc_auc_score(labels, values))
                if n_positive > 0 and n_negative > 0
                else np.nan
            ),
        }
        for alpha in bedroc_alphas:
            name = f"bedroc_alpha_{str(alpha).replace('.', '_')}"
            row[name] = _bedroc(labels, values, float(alpha))
        for fraction in enrichment_fractions:
            name = f"ef_{int(round(float(fraction) * 100))}pct"
            row[name] = _tie_aware_enrichment(labels, values, float(fraction))
        rows.append(row)
    return pd.DataFrame(
        rows,
        columns=[
            "query_id",
            "n_candidates",
            "n_positive",
            "n_negative",
            *metric_names,
        ],
    )


def _provenance_summary(provenance: pd.DataFrame) -> dict[str, Any]:
    if provenance.empty:
        return {
            "split_status_counts": "{}",
            "n_references_input": 0,
            "n_references_after_split": 0,
            "n_removed_close_analogue": 0,
            "n_removed_same_scaffold": 0,
            "n_removed_target_family": 0,
            "n_removed_post_cutoff": 0,
            "n_removed_missing_date": 0,
        }
    count_columns = [
        "n_references_input",
        "n_references_after_split",
        "n_removed_close_analogue",
        "n_removed_same_scaffold",
        "n_removed_target_family",
        "n_removed_post_cutoff",
        "n_removed_missing_date",
    ]
    summary: dict[str, Any] = {
        "split_status_counts": json.dumps(
            provenance["status"].value_counts().sort_index().to_dict(),
            sort_keys=True,
        )
    }
    for column in count_columns:
        summary[column] = (
            int(pd.to_numeric(provenance[column], errors="coerce").sum())
            if column in provenance
            else 0
        )
    return summary


def aggregate_metrics_with_bootstrap(
    query_metrics: pd.DataFrame,
    *,
    split_type: str,
    score_mode: str,
    split_provenance: pd.DataFrame,
    config: ProjectConfig,
) -> pd.DataFrame:
    """Aggregate query metrics with deterministic query-level bootstrap CIs."""

    metric_columns = [
        "auroc",
        *[
            f"bedroc_alpha_{str(float(alpha)).replace('.', '_')}"
            for alpha in config.value("benchmark.bedroc_alphas")
        ],
        *[
            f"ef_{int(round(float(fraction) * 100))}pct"
            for fraction in config.value("benchmark.enrichment_fractions")
        ],
        "mrr",
        "coverage",
    ]
    missing = sorted(set(metric_columns) - set(query_metrics.columns))
    if missing:
        raise ValueError(f"Query metrics are missing configured metrics: {missing}")
    bootstrap_n = int(config.value("benchmark.bootstrap_n"))
    bootstrap_seed = int(config.value("seeds.bootstrap"))
    rng = np.random.default_rng(bootstrap_seed)
    n_queries = len(query_metrics)
    bootstrap_indices = (
        rng.integers(0, n_queries, size=(bootstrap_n, n_queries))
        if n_queries
        else np.empty((bootstrap_n, 0), dtype=int)
    )
    provenance_summary = _provenance_summary(split_provenance)
    rows: list[dict[str, Any]] = []
    for metric in metric_columns:
        values = pd.to_numeric(query_metrics[metric], errors="coerce").to_numpy(
            dtype=float
        )
        finite = np.isfinite(values)
        estimate = float(np.mean(values[finite])) if finite.any() else np.nan
        samples: list[float] = []
        for indices in bootstrap_indices:
            sampled = values[indices]
            sampled = sampled[np.isfinite(sampled)]
            if len(sampled):
                samples.append(float(np.mean(sampled)))
        if samples:
            ci_lower, ci_upper = np.quantile(samples, [0.025, 0.975]).tolist()
        else:
            ci_lower, ci_upper = np.nan, np.nan
        status = "available" if finite.any() else "unavailable_no_evaluable_queries"
        rows.append(
            {
                "split_type": split_type,
                "score_mode": score_mode,
                "metric": metric,
                "estimate": estimate,
                "ci_lower_95": ci_lower,
                "ci_upper_95": ci_upper,
                "n_queries": n_queries,
                "n_evaluable_queries": int(finite.sum()),
                "bootstrap_unit": "query_id",
                "bootstrap_n": bootstrap_n,
                "bootstrap_seed": bootstrap_seed,
                "status": status,
                "snapshot_id": str(config.value("snapshots.snapshot_id")),
                "time_cutoff": str(config.value("benchmark.time_cutoff")),
                "analogue_exclusion_threshold": float(
                    config.value("benchmark.analogue_exclusion_threshold")
                ),
                **provenance_summary,
            }
        )
    return pd.DataFrame(rows)


SCORE_MODES = {
    "2d_only": "chemical_evidence_score",
    "3d_only": "chemical_evidence_score_3d_only",
    "fusion": "chemical_evidence_score_v3",
}


def add_3d_only_score(scores: pd.DataFrame, config: ProjectConfig) -> pd.DataFrame:
    """Rank-fuse only the declared non-2D components, split by split."""

    try:
        from pipeline.evidence_fusion import reciprocal_rank_fusion
    except ModuleNotFoundError:  # direct module execution/import compatibility
        from evidence_fusion import reciprocal_rank_fusion

    three_dimensional = list(
        config.value("benchmark.three_dimensional_components")
    )
    if not three_dimensional:
        raise ValueError("No 3D/pharmacophore fusion components are configured")
    required = {"split_type", "query_id", "target_class", *three_dimensional}
    missing = sorted(required - set(scores.columns))
    if missing:
        raise ValueError(f"Scores are missing 3D-only fields: {missing}")
    if scores.empty:
        result = scores.copy()
        result["chemical_evidence_score_3d_only"] = np.nan
        result["chemical_evidence_score_3d_only_is_probability"] = False
        return result

    pieces = []
    for _, group in scores.groupby("split_type", sort=True):
        source = group[["query_id", "target_class", *three_dimensional]].copy()
        fused = reciprocal_rank_fusion(
            source,
            components=three_dimensional,
            reciprocal_rank_constant=float(
                config.value("fusion.reciprocal_rank_constant")
            ),
        )
        addition = fused[
            ["query_id", "target_class", "chemical_evidence_score_v3"]
        ].rename(
            columns={
                "chemical_evidence_score_v3": "chemical_evidence_score_3d_only"
            }
        )
        piece = group.merge(
            addition,
            on=["query_id", "target_class"],
            how="left",
            validate="one_to_one",
        )
        pieces.append(piece)
    result = pd.concat(pieces, ignore_index=True)
    result["chemical_evidence_score_3d_only_is_probability"] = False
    return result


def compare_score_modes(
    scores: pd.DataFrame, split_provenance: pd.DataFrame, config: ProjectConfig
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return query-level metrics and one honest split×mode comparison table."""

    query_tables = []
    aggregate_tables = []
    for split_type in config.value("benchmark.splits"):
        split_scores = scores[scores["split_type"] == split_type]
        provenance = split_provenance[
            split_provenance["split_type"] == split_type
        ]
        for mode, score_column in SCORE_MODES.items():
            query_table = query_level_metrics(
                split_scores,
                score_column=score_column,
                bedroc_alphas=tuple(
                    float(value) for value in config.value("benchmark.bedroc_alphas")
                ),
                enrichment_fractions=tuple(
                    float(value)
                    for value in config.value("benchmark.enrichment_fractions")
                ),
            )
            query_table["split_type"] = split_type
            query_table["score_mode"] = mode
            query_tables.append(query_table)
            aggregate_tables.append(
                aggregate_metrics_with_bootstrap(
                    query_table,
                    split_type=split_type,
                    score_mode=mode,
                    split_provenance=provenance,
                    config=config,
                )
            )
    query_metrics = pd.concat(query_tables, ignore_index=True)
    comparison = pd.concat(aggregate_tables, ignore_index=True)
    reference = comparison[comparison.score_mode == "2d_only"][
        ["split_type", "metric", "estimate"]
    ].rename(columns={"estimate": "estimate_2d_reference"})
    comparison = comparison.merge(
        reference,
        on=["split_type", "metric"],
        how="left",
        validate="many_to_one",
    )
    comparison["delta_vs_2d"] = (
        comparison["estimate"] - comparison["estimate_2d_reference"]
    )
    available = comparison["estimate"].notna() & comparison[
        "estimate_2d_reference"
    ].notna()
    comparison["performance_vs_2d"] = np.select(
        [
            ~available,
            comparison["score_mode"] == "2d_only",
            comparison["delta_vs_2d"] > 0,
            comparison["delta_vs_2d"] < 0,
        ],
        ["unavailable", "reference", "improved", "worse"],
        default="equal",
    )
    return query_metrics, comparison


_SPLIT_WORKER_CONTEXT = None


def _initialize_split_worker(
    split_results,
    queries,
    ontology,
    quality,
    config,
    reference_evidence_enabled,
):
    global _SPLIT_WORKER_CONTEXT
    _, classes_by_alias = ontology_family_maps(ontology)
    _SPLIT_WORKER_CONTEXT = {
        "split_results": split_results,
        "query_by_id": {
            str(query.get("query_id") or query.get("drug")): query
            for query in queries
        },
        "ontology": ontology,
        "quality": quality,
        "config": config,
        "classes_by_alias": classes_by_alias,
        "reference_evidence_enabled": reference_evidence_enabled,
    }


def _score_split_worker(index):
    """Score one split row in an initialized spawned process."""

    if _SPLIT_WORKER_CONTEXT is None:
        raise RuntimeError("benchmark split worker context was not initialized")
    try:
        from pipeline import open_target_discovery_v2 as discovery
    except ModuleNotFoundError:  # direct module execution/import compatibility
        import open_target_discovery_v2 as discovery

    context = _SPLIT_WORKER_CONTEXT
    split = context["split_results"][index]
    if split.provenance["status"] != "available" or not split.references:
        return None
    source = context["query_by_id"][split.query_id]
    molecule = _molecule(source)
    if molecule is None:
        return None
    query = {
        **source,
        "query_id": split.query_id,
        "mol": molecule,
        "fp": discovery.fp(molecule),
        "maccs": discovery.maccs(molecule),
        "source": "public_benchmark_v3",
    }
    scored = discovery.score_query_v3(
        query,
        split.references,
        context["quality"],
        pd.DataFrame(),
        context["ontology"],
        config=context["config"],
        exclude_close=False,
        return_reference_evidence=context["reference_evidence_enabled"],
    )
    reference_evidence = None
    if context["reference_evidence_enabled"]:
        frame, reference_evidence = scored
        reference_evidence = reference_evidence.copy()
        reference_evidence["split_type"] = split.split_type
    else:
        frame = scored
    if frame.empty:
        return None
    label = str(
        source.get("target_class") or source.get("query_target_label") or ""
    )
    accepted = set(
        context["classes_by_alias"].get(label, {label} if label else set())
    )
    frame["is_active"] = frame["target_class"].isin(accepted).astype(int)
    frame["query_target_label"] = label
    frame["split_type"] = split.split_type
    return frame, reference_evidence


def score_split_results(
    split_results: list[SplitResult],
    queries: list[dict[str, Any]],
    ontology: pd.DataFrame,
    quality: pd.DataFrame,
    config: ProjectConfig,
    *,
    reference_evidence_sink: list[pd.DataFrame] | None = None,
) -> pd.DataFrame:
    """Compute additive 2D/3D/fusion target scores for available split rows."""

    try:
        from pipeline.execution import ordered_process_map
    except ModuleNotFoundError:  # direct module execution/import compatibility
        from execution import ordered_process_map

    scored_splits = ordered_process_map(
        _score_split_worker,
        range(len(split_results)),
        workers=int(config.value("chem3d.scoring_workers")),
        initializer=_initialize_split_worker,
        initargs=(
            split_results,
            queries,
            ontology,
            quality,
            config,
            reference_evidence_sink is not None,
        ),
    )
    rows = []
    for result in scored_splits:
        if result is None:
            continue
        frame, reference_evidence = result
        rows.append(frame)
        if reference_evidence_sink is not None and reference_evidence is not None:
            reference_evidence_sink.append(reference_evidence)
    if not rows:
        return pd.DataFrame(
            columns=[
                "query_id",
                "target_class",
                "split_type",
                "is_active",
                *config.value("fusion.components"),
                "chemical_evidence_score",
                "chemical_evidence_score_v3",
            ]
        )
    return add_3d_only_score(pd.concat(rows, ignore_index=True), config)


def main() -> None:
    """Run the pinned v3 benchmark and write auditable split/comparison tables."""

    try:
        from pipeline import benchmark_v2
        from pipeline.benchmark_decoys import load_property_matched_decoys
        from pipeline.config import load_config, set_global_seed
        from pipeline.snapshots import verify_snapshot
    except ModuleNotFoundError:  # direct module execution/import compatibility
        import benchmark_v2
        from benchmark_decoys import load_property_matched_decoys
        from config import load_config, set_global_seed
        from snapshots import verify_snapshot

    config = load_config()
    set_global_seed(config)
    if bool(config.value("snapshots.verify_on_load")):
        verify_snapshot(config)
    results_dir = config.path_for("results")
    results_dir.mkdir(parents=True, exist_ok=True)
    ontology = pd.read_csv(config.path_for("target_ontology"))
    subtype_path = config.path_for("target_subtype_ontology")
    if subtype_path.exists():
        ontology = pd.concat(
            [ontology, pd.read_csv(subtype_path)], ignore_index=True
        ).drop_duplicates("target_class", keep="first")
    queries = benchmark_v2.load_bench()
    references = benchmark_v2.load_refs()
    split_results, provenance = generate_splits(
        queries, references, ontology, config
    )
    provenance.to_csv(results_dir / "benchmark_split_provenance_v3.csv", index=False)
    decoy_result = load_property_matched_decoys(
        config.path_for("property_matched_decoys")
    )
    pd.DataFrame([decoy_result.status]).to_csv(
        results_dir / "benchmark_decoy_status_v3.csv", index=False
    )
    quality_path = config.path_for("reference_quality")
    quality = (
        pd.read_csv(quality_path)
        if quality_path.exists()
        else pd.DataFrame(columns=["target_class"])
    )
    reference_evidence_tables: list[pd.DataFrame] = []
    scores = score_split_results(
        split_results,
        queries,
        ontology,
        quality,
        config,
        reference_evidence_sink=reference_evidence_tables,
    )
    scores.to_csv(results_dir / "benchmark_target_scores_by_split_v3.csv", index=False)
    reference_evidence = (
        pd.concat(reference_evidence_tables, ignore_index=True)
        if reference_evidence_tables
        else pd.DataFrame()
    )
    reference_evidence.to_csv(
        results_dir / "benchmark_reference_evidence_by_split_v3.csv", index=False
    )
    query_metrics, comparison = compare_score_modes(scores, provenance, config)
    query_metrics.to_csv(results_dir / "benchmark_query_metrics_v3.csv", index=False)
    comparison.to_csv(results_dir / "benchmark_mode_comparison_v3.csv", index=False)
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
