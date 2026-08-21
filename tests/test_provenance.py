from __future__ import annotations

import json
from pathlib import Path

import yaml

import run_pipeline
from pipeline.config import load_config
from pipeline.provenance import start_run, utc_now, write_run_manifest


ROOT = Path(__file__).resolve().parents[1]


def _temporary_config(tmp_path: Path):
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["paths"]["results"] = str(tmp_path / "results")
    raw["paths"]["snapshot_manifest"] = str(tmp_path / "missing-snapshot.json")
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    return load_config(path, project_root=ROOT)


def test_manifest_is_written_and_finalized_without_inventing_missing_versions(
    tmp_path: Path,
) -> None:
    config = _temporary_config(tmp_path)
    context = start_run(config)

    running_path = config.path_for("results") / "run_manifest.json"
    running = json.loads(running_path.read_text(encoding="utf-8"))
    assert running["status"] == "running"
    assert running["data_snapshots"]["status"] == "manifest_missing"
    assert running["data_snapshots"]["versions"] is None

    final = write_run_manifest(
        config, context, status="completed", completed_at=utc_now()
    )
    saved = json.loads(running_path.read_text(encoding="utf-8"))
    archived = (
        config.path_for("results")
        / "run_manifests"
        / f"{context.run_id}.json"
    )

    assert final == saved
    assert saved["status"] == "completed"
    assert saved["config"]["sha256"] == config.config_hash
    assert saved["seeds"]["conformer"] == config.value("seeds.conformer")
    assert saved["code"]["commit"]
    assert saved["environment"]["packages"]["rdkit"]
    assert archived.is_file()


def test_configured_runner_emits_a_completed_manifest(
    tmp_path: Path, monkeypatch
) -> None:
    config = _temporary_config(tmp_path)
    monkeypatch.setattr(run_pipeline, "STEPS", [])

    run_pipeline.main(["--config", str(config.path)])

    manifest = json.loads(
        (config.path_for("results") / "run_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["status"] == "completed"
    assert manifest["timestamps"]["completed_at_utc"] is not None
