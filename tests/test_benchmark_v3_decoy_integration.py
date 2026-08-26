from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd

from pipeline.benchmark_decoys import CANDIDATE_COLUMNS
from pipeline.benchmark_v3 import (
    SCORE_MODES,
    build_property_decoy_target_score_table,
    compare_property_decoy_score_modes,
)
from pipeline.config import ProjectConfig, load_config


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["benchmark"]["bootstrap_n"] = 50
    return ProjectConfig(path=base.path, root=base.root, data=data)


def _candidates() -> pd.DataFrame:
    rows = []
    for target in ["A", "B"]:
        for source_id, is_active in [("active", 1), ("decoy-1", 0), ("decoy-2", 0)]:
            candidate_id = f"{source_id}::{target}"
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "source_candidate_id": source_id,
                    "query_id": candidate_id,
                    "canonical_smiles": "CCO",
                    "canonical_smiles_source": "CCO",
                    "source_target_label": target,
                    "target_class": target,
                    "matched_target_class": target,
                    "matched_active_query_id": "active",
                    "retrieval_task_id": f"active::{target}",
                    "target_mapping_status": "fixture_exact",
                    "is_active": is_active,
                    "candidate_type": (
                        "curated_benchmark_active"
                        if is_active
                        else "property_matched_decoy"
                    ),
                    "label_semantics": (
                        "curated_known_mechanism_active"
                        if is_active
                        else "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
                    ),
                    "decoy_source_dataset": "fixture" if not is_active else pd.NA,
                    "decoy_source_version": "v1" if not is_active else pd.NA,
                    "decoy_source_record_id": source_id if not is_active else pd.NA,
                    "property_matching_method": "fixture" if not is_active else pd.NA,
                }
            )
    return pd.DataFrame(rows, columns=list(CANDIDATE_COLUMNS))


def _scores(config: ProjectConfig, *, active_high: bool) -> pd.DataFrame:
    rows = []
    ad_flags = [
        "in_domain",
        "near_domain",
        "out_of_domain",
        "unassessable_missing_tanimoto",
        "in_domain",
        "near_domain",
    ]
    for index, candidate in _candidates().iterrows():
        high = bool(candidate.is_active) == active_high
        value = 0.9 if high else 0.1
        tanimoto = [0.8, 0.3, 0.1, np.nan, 0.7, 0.35][index]
        row = {
            "query_id": candidate.candidate_id,
            "target_class": candidate.matched_target_class,
            "split_type": "scaffold",
            "is_active": 1,  # target-retrieval label must not replace external label
            **{component: value for component in config.value("fusion.components")},
            "chemical_evidence_score": value,
            "chemical_evidence_score_3d_only": value,
            "chemical_evidence_score_v3": value,
            "ad_nearest_reference_tanimoto": tanimoto,
            "ad_nearest_reference_usrcat_similarity": 0.6,
            "ad_nearest_reference_usrcat_distance": 0.4,
            "ad_tanimoto_flag": ad_flags[index],
            "ad_usrcat_status": "continuous_similarity_available_no_calibrated_cutoff",
            "applicability_domain_flag": ad_flags[index],
            "applicability_domain_flag_basis": "fixture thresholds",
            "ad_tanimoto_in_threshold": 0.4,
            "ad_tanimoto_out_threshold": 0.25,
            "ad_shortlist_eligible": ad_flags[index] in {"in_domain", "near_domain"},
            "ad_shortlist_discount_policy": "ordering_only_no_score_rescaling",
        }
        rows.append(row)
    return pd.DataFrame(rows)


def _provenance() -> pd.DataFrame:
    rows = []
    for candidate_id in _candidates().candidate_id:
        rows.append(
            {
                "query_id": candidate_id,
                "split_type": "scaffold",
                "status": "available",
                "status_reason": "fixture split",
                "analogue_exclusion_threshold": 0.85,
                "query_murcko_scaffold": "",
                "n_references_input": 10,
                "n_references_after_split": 8,
                "n_removed_close_analogue": 1,
                "n_removed_same_scaffold": 1,
                "n_removed_target_family": 0,
                "n_removed_post_cutoff": 0,
                "n_removed_missing_date": 0,
                "max_remaining_query_reference_tanimoto": 0.8,
                "analogue_leakage_guard_passed": True,
            }
        )
    return pd.DataFrame(rows)


def test_missing_decoy_artifact_produces_empty_scores_and_pending_metrics() -> None:
    config = _config()
    candidates = pd.DataFrame(columns=list(CANDIDATE_COLUMNS))
    target_scores = build_property_decoy_target_score_table(
        pd.DataFrame(), candidates, pd.DataFrame(), config
    )
    query_metrics, metrics = compare_property_decoy_score_modes(
        target_scores, pd.DataFrame(), config
    )

    assert target_scores.empty
    assert "applicability_domain_flag" in target_scores
    assert {"2d_only_status", "3d_only_status", "fusion_status"}.issubset(
        target_scores.columns
    )
    assert query_metrics.empty
    assert len(metrics) == 3 * 3 * 7
    assert metrics["estimate"].isna().all()
    assert metrics["status"].eq("unavailable_no_evaluable_queries").all()
    assert metrics["bootstrap_unit"].eq(
        "matched_target_class_cluster"
    ).all()


def test_property_decoys_are_scored_as_molecules_and_change_enrichment() -> None:
    config = _config()
    candidates = _candidates()
    provenance = _provenance()
    good_scores = build_property_decoy_target_score_table(
        _scores(config, active_high=True), candidates, provenance, config
    )
    bad_scores = build_property_decoy_target_score_table(
        _scores(config, active_high=False), candidates, provenance, config
    )

    assert good_scores.loc[
        good_scores.split_type == "scaffold", "is_active"
    ].tolist() == [1, 0, 0, 1, 0, 0]
    assert not good_scores["scores_are_probabilities"].any()
    assert set(
        good_scores.loc[
            good_scores.split_type == "scaffold", "applicability_domain_flag"
        ]
    ) == {
        "in_domain",
        "near_domain",
        "out_of_domain",
        "unassessable_missing_tanimoto",
    }
    assert good_scores.loc[
        good_scores.split_type == "scaffold", "candidate_score_status"
    ].eq("available_all_modes").all()
    good_query, good_metrics = compare_property_decoy_score_modes(
        good_scores, provenance, config
    )
    _, bad_metrics = compare_property_decoy_score_modes(
        bad_scores, provenance, config
    )

    good = good_metrics[
        (good_metrics.split_type == "scaffold")
        & (good_metrics.score_mode == "fusion")
    ].set_index("metric")
    bad = bad_metrics[
        (bad_metrics.split_type == "scaffold")
        & (bad_metrics.score_mode == "fusion")
    ].set_index("metric")
    assert good.loc["auroc", "estimate"] == 1.0
    assert bad.loc["auroc", "estimate"] == 0.0
    assert good.loc["ef_1pct", "estimate"] > bad.loc["ef_1pct", "estimate"]
    assert set(good_query.score_mode) == set(SCORE_MODES)
    temporal_query = good_query[good_query.split_type == "temporal"]
    assert temporal_query["coverage"].isna().all()
    assert temporal_query["coverage_status"].str.startswith(
        "unavailable_candidate_split_status_"
    ).all()
    assert temporal_query["discrimination_status"].str.startswith(
        "unavailable_candidate_split_status_"
    ).all()
    assert temporal_query["status"].eq("unavailable_split_task").all()
    temporal_coverage = good_metrics[
        (good_metrics.split_type == "temporal")
        & (good_metrics.metric == "coverage")
    ]
    assert temporal_coverage["estimate"].isna().all()
    assert temporal_coverage["status"].eq(
        "unavailable_no_evaluable_queries"
    ).all()
    assert good_metrics["bootstrap_unit"].eq(
        "matched_target_class_cluster"
    ).all()
    assert good_metrics["cross_target_decoy_policy"].eq(
        "specificity_margin_only_not_inactive"
    ).all()


def test_mode_specific_missingness_keeps_coverage_and_scored_counts() -> None:
    config = _config()
    candidates = _candidates()
    provenance = _provenance()
    source = _scores(config, active_high=True)
    source.loc[
        source["query_id"].str.startswith("decoy-"),
        "chemical_evidence_score_3d_only",
    ] = np.nan

    target_scores = build_property_decoy_target_score_table(
        source, candidates, provenance, config
    )
    scaffold_decoys = target_scores[
        (target_scores.split_type == "scaffold")
        & (target_scores.is_active == 0)
    ]
    assert scaffold_decoys["2d_only_status"].eq("available").all()
    assert scaffold_decoys["3d_only_status"].eq(
        "unavailable_no_matched_target_score"
    ).all()
    assert scaffold_decoys["fusion_status"].eq("available").all()
    assert scaffold_decoys["candidate_score_status"].eq(
        "partial_mode_coverage"
    ).all()

    query_metrics, metrics = compare_property_decoy_score_modes(
        target_scores, provenance, config
    )
    three_d_query = query_metrics[
        (query_metrics.split_type == "scaffold")
        & (query_metrics.score_mode == "3d_only")
    ]
    assert three_d_query["coverage"].eq(1 / 3).all()
    assert three_d_query["auroc"].isna().all()
    assert three_d_query["discrimination_status"].eq(
        "unavailable_requires_scored_active_and_property_decoy"
    ).all()
    assert three_d_query["coverage_status"].eq("available").all()
    assert three_d_query["status"].eq("available_coverage_only").all()

    three_d_metrics = metrics[
        (metrics.split_type == "scaffold")
        & (metrics.score_mode == "3d_only")
    ].set_index("metric")
    assert three_d_metrics.loc["coverage", "estimate"] == 1 / 3
    assert three_d_metrics.loc["coverage", "status"] == "available"
    assert pd.isna(three_d_metrics.loc["auroc", "estimate"])
    assert three_d_metrics["n_active_candidate_target_pairs_scored"].eq(2).all()
    assert three_d_metrics["n_decoy_candidate_target_pairs_scored"].eq(0).all()
    assert three_d_metrics["n_candidate_target_pairs_scored"].eq(2).all()


def test_scored_candidates_cannot_drop_applicability_domain_fields() -> None:
    config = _config()
    source = _scores(config, active_high=True).drop(
        columns=["applicability_domain_flag"]
    )

    try:
        build_property_decoy_target_score_table(
            source, _candidates(), _provenance(), config
        )
    except ValueError as error:
        assert "applicability-domain" in str(error)
        assert "applicability_domain_flag" in str(error)
    else:
        raise AssertionError("scored candidates must retain their AD fields")


def test_noncanonical_retrieval_task_id_is_rejected() -> None:
    config = _config()
    candidates = _candidates()
    candidates.loc[0, "retrieval_task_id"] = "corrupt-task-id"

    try:
        build_property_decoy_target_score_table(
            pd.DataFrame(), candidates, pd.DataFrame(), config
        )
    except ValueError as error:
        assert "retrieval_task_id" in str(error)
        assert "matched_active_query_id::matched_target_class" in str(error)
    else:
        raise AssertionError("a noncanonical retrieval task ID must be rejected")


def test_active_specific_decoy_tasks_are_not_pooled_within_target_class() -> None:
    config = _config()
    specifications = [
        ("active-a::A", "active-a", 1, 0.9),
        ("decoy-a::A", "active-a", 0, 0.1),
        ("active-b::A", "active-b", 1, 0.1),
        ("decoy-b::A", "active-b", 0, 0.9),
    ]
    candidate_rows = []
    score_rows = []
    provenance_rows = []
    for candidate_id, matched_active_id, is_active, value in specifications:
        candidate = {column: pd.NA for column in CANDIDATE_COLUMNS}
        candidate.update(
            {
                "candidate_id": candidate_id,
                "source_candidate_id": candidate_id.split("::", 1)[0],
                "query_id": candidate_id,
                "canonical_smiles": "CCO",
                "canonical_smiles_source": "CCO",
                "source_target_label": "A",
                "target_class": "A",
                "matched_target_class": "A",
                "matched_active_query_id": matched_active_id,
                "retrieval_task_id": f"{matched_active_id}::A",
                "target_mapping_status": "fixture_exact",
                "is_active": is_active,
                "candidate_type": (
                    "curated_benchmark_active"
                    if is_active
                    else "property_matched_decoy"
                ),
                "label_semantics": (
                    "curated_known_mechanism_active"
                    if is_active
                    else "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
                ),
            }
        )
        candidate_rows.append(candidate)
        score_rows.append(
            {
                "query_id": candidate_id,
                "target_class": "A",
                "split_type": "scaffold",
                **{
                    component: value
                    for component in config.value("fusion.components")
                },
                "chemical_evidence_score": value,
                "chemical_evidence_score_3d_only": value,
                "chemical_evidence_score_v3": value,
                "ad_nearest_reference_tanimoto": 0.8,
                "ad_nearest_reference_usrcat_similarity": 0.6,
                "ad_nearest_reference_usrcat_distance": 0.4,
                "ad_tanimoto_flag": "in_domain",
                "ad_usrcat_status": (
                    "continuous_similarity_available_no_calibrated_cutoff"
                ),
                "applicability_domain_flag": "in_domain",
                "applicability_domain_flag_basis": "fixture thresholds",
                "ad_tanimoto_in_threshold": 0.4,
                "ad_tanimoto_out_threshold": 0.25,
                "ad_shortlist_eligible": True,
                "ad_shortlist_discount_policy": (
                    "ordering_only_no_score_rescaling"
                ),
            }
        )
        provenance_rows.append(
            {
                "query_id": candidate_id,
                "split_type": "scaffold",
                "status": "available",
                "status_reason": "fixture split",
            }
        )
    candidates = pd.DataFrame(candidate_rows, columns=list(CANDIDATE_COLUMNS))
    provenance = pd.DataFrame(provenance_rows)
    target_scores = build_property_decoy_target_score_table(
        pd.DataFrame(score_rows), candidates, provenance, config
    )

    query_metrics, metrics = compare_property_decoy_score_modes(
        target_scores, provenance, config
    )
    paired = query_metrics[
        (query_metrics.split_type == "scaffold")
        & (query_metrics.score_mode == "fusion")
    ].sort_values("retrieval_task_id")
    assert paired["retrieval_task_id"].tolist() == ["active-a::A", "active-b::A"]
    assert paired["n_candidates_total"].tolist() == [2, 2]
    assert paired["auroc"].tolist() == [1.0, 0.0]

    aggregate = metrics[
        (metrics.split_type == "scaffold")
        & (metrics.score_mode == "fusion")
    ].set_index("metric")
    assert aggregate.loc["auroc", "estimate"] == 0.5
    assert aggregate["n_retrieval_tasks"].eq(2).all()
    assert aggregate["n_matched_target_classes"].eq(1).all()
    assert aggregate["bootstrap_unit"].eq(
        "matched_target_class_cluster"
    ).all()
