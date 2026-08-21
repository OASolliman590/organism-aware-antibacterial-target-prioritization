from __future__ import annotations

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
        [
            {
                "query_id": "trimethoprim",
                "canonical_smiles": "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",
                "target_class": "DHFR",
            }
        ]
    )

    candidates = build_active_decoy_candidates(active, result).set_index(
        "candidate_id"
    )

    assert result.status["status"] == "available"
    assert result.status["source_datasets"] == "DUD-E"
    assert candidates.loc["trimethoprim", "is_active"] == 1
    assert candidates.loc["d1", "is_active"] == 0
    assert candidates.loc["d1", "label_semantics"] == (
        "presumed_inactive_property_matched_decoy_not_confirmed_inactive"
    )
    assert candidates.loc["d1", "decoy_source_record_id"] == "ZINC000001"


def test_decoys_without_recorded_version_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "decoys.csv"
    pd.DataFrame(
        [
            {
                "decoy_id": "d1",
                "canonical_smiles": "CCOC(=O)N",
                "matched_target_class": "DHFR",
                "source_dataset": "DUD-E",
                "source_version": "unrecorded",
                "source_record_id": "ZINC000001",
                "property_matching_method": "DUD-E supplied target decoy",
            }
        ]
    ).to_csv(path, index=False)

    try:
        load_property_matched_decoys(path)
    except ValueError as error:
        assert "source_version" in str(error)
    else:
        raise AssertionError("unversioned decoys must not be accepted")
