from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from pipeline.benchmark_decoys import (
    build_active_decoy_candidates,
    load_property_matched_decoys,
)


def test_missing_decoy_artifact_is_pending_and_never_uses_cross_target_ligands(
    tmp_path: Path,
) -> None:
    result = load_property_matched_decoys(tmp_path / "absent.csv")

    assert result.decoys.empty
    assert result.status["status"] == "pending_missing_property_matched_decoy_dataset"
    assert result.status["n_decoys"] == 0
    assert result.status["cross_target_decoy_policy"] == "specificity_margin_only_not_inactive"


def test_versioned_property_matched_decoys_integrate_with_honest_labels(
    tmp_path: Path,
) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "trimethoprim",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            },
        ]
    ).to_csv(path, index=False)
    result = load_property_matched_decoys(path)
    active = pd.DataFrame(
        [
            {
                "query_id": "trimethoprim",
                "canonical_smiles": "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",
                "target_class": "DHFR",
                "activity_date": "2017-05-01",
            }
        ]
    )

    candidates = build_active_decoy_candidates(active, result).set_index(
        "candidate_id"
    )

    assert result.status["status"] == "available"
    assert result.status["source_datasets"] == "DUD-E"
    assert len(result.status["artifact_sha256"]) == 64
    assert candidates.loc["trimethoprim::DHFR", "is_active"] == 1
    assert candidates.loc["d1::DHFR", "is_active"] == 0
    assert candidates.loc["d1::DHFR", "query_id"] == "d1::DHFR"
    assert candidates.loc["d1::DHFR", "matched_active_query_id"] == "trimethoprim"
    assert candidates.loc["d1::DHFR", "label_semantics"] == (
        "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
    )
    assert candidates.loc["d1::DHFR", "decoy_source_record_id"] == "ZINC000001"
    assert candidates.loc["trimethoprim::DHFR", "activity_date"] == "2017-05-01"
    assert pd.isna(candidates.loc["d1::DHFR", "activity_date"])


def test_decoys_without_recorded_version_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "active",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "unrecorded",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            },
        ]
    ).to_csv(path, index=False)

    try:
        load_property_matched_decoys(path)
    except ValueError as error:
        assert "source_version" in str(error)
    else:
        raise AssertionError("unversioned decoys must not be accepted")


def test_exact_active_structure_is_rejected_as_a_decoy(tmp_path: Path) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCO",
                "matched_active_query_id": "active",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            }
        ]
    ).to_csv(path, index=False)
    result = load_property_matched_decoys(path)
    active = pd.DataFrame(
        [{"query_id": "active", "canonical_smiles": "CCO", "target_class": "DHFR"}]
    )

    try:
        build_active_decoy_candidates(
            active, result, valid_target_classes={"DHFR"}
        )
    except ValueError as error:
        assert "overlap exact active structures" in str(error)
    else:
        raise AssertionError("an exact active structure must not be accepted as a decoy")


def test_unmapped_decoy_target_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "active",
                "matched_target_class": "NOT_IN_ONTOLOGY",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            }
        ]
    ).to_csv(path, index=False)
    result = load_property_matched_decoys(path)
    active = pd.DataFrame(
        [{"query_id": "active", "canonical_smiles": "CCO", "target_class": "DHFR"}]
    )

    try:
        build_active_decoy_candidates(
            active, result, valid_target_classes={"DHFR"}
        )
    except ValueError as error:
        assert "pinned ontology" in str(error)
    else:
        raise AssertionError("an unmapped decoy target must not be accepted")


def test_unknown_matched_active_query_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "not-a-benchmark-query",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            }
        ]
    ).to_csv(path, index=False)
    result = load_property_matched_decoys(path)
    active = pd.DataFrame(
        [{"query_id": "active", "canonical_smiles": "CCO", "target_class": "DHFR"}]
    )

    try:
        build_active_decoy_candidates(
            active, result, valid_target_classes={"DHFR"}
        )
    except ValueError as error:
        assert "unknown active query_id" in str(error)
    else:
        raise AssertionError("a decoy linked to an unknown active must not be accepted")


def test_decoy_active_target_pair_must_follow_pinned_ontology_mapping(
    tmp_path: Path,
) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "active",
                "matched_target_class": "GYRA",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            }
        ]
    ).to_csv(path, index=False)
    result = load_property_matched_decoys(path)
    active = pd.DataFrame(
        [{"query_id": "active", "canonical_smiles": "CCO", "target_class": "DHFR"}]
    )

    try:
        build_active_decoy_candidates(
            active, result, valid_target_classes={"DHFR", "GYRA"}
        )
    except ValueError as error:
        assert "active/target links" in str(error)
    else:
        raise AssertionError("an unsupported active/target link must not be accepted")


def test_unmapped_active_is_logged_without_aborting_supported_decoy_task(
    tmp_path: Path,
) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "mapped-active",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            },
            {
                "decoy_id": "d2",
                "canonical_smiles": "CCC",
                "matched_active_query_id": "unmapped-active",
                "matched_target_class": "NOT_IN_ONTOLOGY",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000002",
                "property_matching_method": "DUD-E supplied target decoy",
            },
        ]
    ).to_csv(path, index=False)
    result = load_property_matched_decoys(path)
    active = pd.DataFrame(
        [
            {
                "query_id": "mapped-active",
                "canonical_smiles": "CCO",
                "target_class": "DHFR",
            },
            {
                "query_id": "unmapped-active",
                "canonical_smiles": "CCN",
                "target_class": "lipid A / membrane",
            },
        ]
    )
    gaps: list[dict[str, object]] = []

    candidates = build_active_decoy_candidates(
        active,
        result,
        classes_by_alias={"DHFR": {"DHFR"}},
        valid_target_classes={"DHFR"},
        unmapped_active_sink=gaps,
    )

    assert set(candidates["candidate_id"]) == {"mapped-active::DHFR", "d1::DHFR"}
    assert gaps == [
        {
            "source_candidate_id": "unmapped-active",
            "source_target_label": "lipid A / membrane",
            "mapping_status": "pending_no_pinned_ontology_mapping",
            "status_reason": (
                "Active benchmark label has no exact or declared alias mapping in "
                "the pinned target ontology"
            ),
            "n_linked_decoys_excluded": 1,
            "linked_decoy_ids_sha256": hashlib.sha256(b"d2").hexdigest(),
        }
    ]


def test_same_structure_may_be_provenanced_for_distinct_active_tasks(
    tmp_path: Path,
) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "active-a",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            },
            {
                "decoy_id": "d2",
                "canonical_smiles": "CCOC(=O)N",
                "matched_active_query_id": "active-b",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "download-2026-08-21",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            },
        ]
    ).to_csv(path, index=False)

    result = load_property_matched_decoys(path)

    assert len(result.decoys) == 2
    assert set(result.decoys["matched_active_query_id"]) == {
        "active-a",
        "active-b",
    }
