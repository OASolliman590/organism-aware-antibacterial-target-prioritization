from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

import run_pipeline
from pipeline.config import load_config
from pipeline.snapshots import (
    finalize_refresh_snapshot,
    RefreshSafetyError,
    SnapshotIntegrityError,
    hash_declared_paths,
    require_refresh_output,
    verify_snapshot,
    write_text_exclusive,
)


ROOT = Path(__file__).resolve().parents[1]


def test_repository_snapshot_verifies() -> None:
    verified = verify_snapshot(load_config(ROOT / "config.yaml"))

    assert verified["snapshot_id"] == "v2-public-baseline-2ed4684"
    assert set(verified["datasets"]) == {
        "chembl_reference_ligands",
        "reference_quality",
        "pubchem_benchmark",
        "uniprot_species_targets",
        "card_derived_annotations",
        "rcsb_structure_catalog",
        "curated_ontologies",
    }
    assert all("fetch_" not in step for step in run_pipeline.ANALYSIS_STEPS)


def test_snapshot_verification_detects_content_changes(tmp_path: Path) -> None:
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["snapshots"]["snapshot_id"] = "fixture-snapshot"
    raw["paths"]["snapshot_manifest"] = "SNAPSHOT_VERSIONS.json"
    config_path = tmp_path / "config.yaml"
    data_path = tmp_path / "fixture.txt"
    data_path.write_text("original\n", encoding="utf-8")
    record = hash_declared_paths(tmp_path, ["fixture.txt"])
    manifest = {
            "schema_version": 1,
            "snapshot_id": "fixture-snapshot",
            "analysis_ready": True,
        "hash_algorithm": "sha256-tree-v1",
        "datasets": {"fixture": {"paths": ["fixture.txt"], **record}},
    }
    (tmp_path / "SNAPSHOT_VERSIONS.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    config_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    config = load_config(config_path)

    verify_snapshot(config)
    data_path.write_text("changed\n", encoding="utf-8")

    with pytest.raises(SnapshotIntegrityError, match="mismatch"):
        verify_snapshot(config)


def test_refresh_requires_new_dated_snapshot_and_never_overwrites(
    tmp_path: Path,
) -> None:
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    snapshot_id = "2026-08-21-fixture"
    raw["run"]["refresh_external_data"] = True
    raw["snapshots"]["snapshot_id"] = snapshot_id
    raw["paths"]["snapshot_manifest"] = (
        f"data/snapshots/{snapshot_id}/SNAPSHOT_VERSIONS.json"
    )
    config_path = tmp_path / "refresh.yaml"
    config_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    config = load_config(config_path)
    output = tmp_path / "data" / "snapshots" / snapshot_id / "source.json"

    assert require_refresh_output(config, output) == output.resolve()
    write_text_exclusive(config, output, "{}\n")

    with pytest.raises(RefreshSafetyError, match="overwrite"):
        write_text_exclusive(config, output, "{}\n")


def test_refresh_finalization_freezes_raw_data_but_marks_analysis_pending(
    tmp_path: Path,
) -> None:
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    snapshot_id = "2026-08-21-raw-fixture"
    prefix = f"data/snapshots/{snapshot_id}"
    raw["run"]["refresh_external_data"] = True
    raw["snapshots"]["snapshot_id"] = snapshot_id
    raw["refresh"]["chembl"]["source_release"] = "fixture-chembl"
    raw["refresh"]["uniprot"]["source_release"] = "fixture-uniprot"
    raw["refresh"]["card"]["source_version"] = "fixture-card"
    path_values = {
        "snapshot_manifest": f"{prefix}/SNAPSHOT_VERSIONS.json",
        "reference_ligands": f"{prefix}/reference_ligands",
        "chembl_cache": f"{prefix}/chembl_cache",
        "benchmark": f"{prefix}/benchmark/eskape_benchmark_drugs.csv",
        "species_proteins": f"{prefix}/species/species_target_proteins.csv",
        "species_fasta": f"{prefix}/species/species_target_proteins.fasta",
        "species_metadata": f"{prefix}/species/species_target_fetch_metadata.json",
        "card_raw": f"{prefix}/card",
        "structure_candidates": f"{prefix}/structures/rcsb_structure_candidates_v2.csv",
        "structure_summary": f"{prefix}/structures/rcsb_structure_summary_v2.csv",
    }
    raw["paths"].update(path_values)
    config_path = tmp_path / "refresh.yaml"
    config_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    config = load_config(config_path)

    files = [
        config.path_for("reference_ligands") / "ref_ligands_fixture.json",
        config.path_for("chembl_cache") / "manifest.json",
        config.path_for("benchmark"),
        config.path_for("benchmark").with_name("eskape_benchmark_sources.json"),
        config.path_for("species_proteins"),
        config.path_for("species_fasta"),
        config.path_for("species_metadata"),
        config.path_for("card_raw") / "card-data.tar.bz2",
        config.path_for("structure_candidates"),
        config.path_for("structure_summary"),
    ]
    for path in files:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n", encoding="utf-8")

    manifest = finalize_refresh_snapshot(config)

    assert manifest["analysis_ready"] is False
    assert manifest["status"] == "raw_refresh_complete_derived_layers_pending"
    assert config.path_for("snapshot_manifest").is_file()
    with pytest.raises(SnapshotIntegrityError, match="not marked analysis_ready"):
        verify_snapshot(config)
