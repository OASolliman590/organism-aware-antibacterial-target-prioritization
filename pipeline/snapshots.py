"""Verify immutable data snapshots and guard refresh destinations."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import re
from typing import Any, Iterable

try:
    from pipeline.config import ProjectConfig
except ModuleNotFoundError:  # imported by direct ``python pipeline/<script>.py`` execution
    from config import ProjectConfig


TREE_HASH_ALGORITHM = "sha256-tree-v1"
DATED_SNAPSHOT_ID = re.compile(r"^\d{4}-\d{2}-\d{2}(?:[._-][A-Za-z0-9._-]+)?$")


class SnapshotError(RuntimeError):
    """Base class for snapshot validation and refresh errors."""


class SnapshotIntegrityError(SnapshotError):
    """Raised when pinned files do not match their declared snapshot."""


class RefreshSafetyError(SnapshotError):
    """Raised when refresh mode could alter an existing or undated snapshot."""


def _files_for_paths(root: Path, declared_paths: Iterable[str]) -> list[Path]:
    files: list[Path] = []
    for declared in declared_paths:
        path = (root / declared).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError as exc:
            raise SnapshotIntegrityError(
                f"Snapshot path leaves the project root: {declared}"
            ) from exc
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(item for item in path.rglob("*") if item.is_file())
        else:
            raise SnapshotIntegrityError(f"Pinned snapshot path is missing: {declared}")
    return sorted(set(files), key=lambda item: item.relative_to(root).as_posix())


def hash_declared_paths(root: Path, declared_paths: Iterable[str]) -> dict[str, Any]:
    """Hash relative filenames and file digests to bind content and layout."""

    digest = hashlib.sha256()
    files = _files_for_paths(root, declared_paths)
    total_bytes = 0
    for path in files:
        relative = path.relative_to(root).as_posix()
        content = path.read_bytes()
        total_bytes += len(content)
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(hashlib.sha256(content).digest())
        digest.update(b"\n")
    return {
        "hash_algorithm": TREE_HASH_ALGORITHM,
        "sha256": digest.hexdigest(),
        "file_count": len(files),
        "total_bytes": total_bytes,
    }


def load_snapshot_manifest(config: ProjectConfig) -> dict[str, Any]:
    path = config.path_for("snapshot_manifest")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise SnapshotIntegrityError(f"Snapshot manifest is missing: {path}") from exc
    except json.JSONDecodeError as exc:
        raise SnapshotIntegrityError(f"Snapshot manifest is invalid JSON: {path}") from exc
    if not isinstance(manifest, dict):
        raise SnapshotIntegrityError("Snapshot manifest root must be an object")
    return manifest


def verify_snapshot(config: ProjectConfig) -> dict[str, Any]:
    """Fail closed if any pinned input differs from the declared snapshot."""

    manifest = load_snapshot_manifest(config)
    expected_id = str(config.value("snapshots.snapshot_id"))
    if manifest.get("snapshot_id") != expected_id:
        raise SnapshotIntegrityError(
            f"Configured snapshot {expected_id!r} does not match manifest "
            f"{manifest.get('snapshot_id')!r}"
        )
    if manifest.get("hash_algorithm") != TREE_HASH_ALGORITHM:
        raise SnapshotIntegrityError(
            f"Unsupported snapshot hash algorithm: {manifest.get('hash_algorithm')!r}"
        )
    if manifest.get("analysis_ready") is not True:
        raise SnapshotIntegrityError(
            "Snapshot is not marked analysis_ready; incomplete refreshes cannot be "
            "used for reported analysis"
        )
    datasets = manifest.get("datasets")
    if not isinstance(datasets, dict) or not datasets:
        raise SnapshotIntegrityError("Snapshot manifest has no datasets")
    verified: dict[str, Any] = {}
    for name, record in datasets.items():
        if not isinstance(record, dict) or not isinstance(record.get("paths"), list):
            raise SnapshotIntegrityError(f"Dataset {name!r} has no path list")
        actual = hash_declared_paths(config.root, record["paths"])
        for field in ("sha256", "file_count", "total_bytes"):
            if actual[field] != record.get(field):
                raise SnapshotIntegrityError(
                    f"Dataset {name!r} {field} mismatch: "
                    f"expected {record.get(field)!r}, observed {actual[field]!r}"
                )
        verified[name] = actual
    return {
        "snapshot_id": expected_id,
        "manifest_path": str(config.path_for("snapshot_manifest")),
        "datasets": verified,
    }


def refresh_snapshot_root(config: ProjectConfig) -> Path:
    """Return a safe new refresh root or fail before external access."""

    if not bool(config.value("run.refresh_external_data")):
        raise RefreshSafetyError(
            "Live data access is refresh-only; set run.refresh_external_data: true "
            "in a dedicated refresh config"
        )
    snapshot_id = str(config.value("snapshots.snapshot_id"))
    if not DATED_SNAPSHOT_ID.fullmatch(snapshot_id):
        raise RefreshSafetyError(
            "A refresh snapshot_id must start with an ISO date (YYYY-MM-DD)"
        )
    root = (config.root / "data" / "snapshots" / snapshot_id).resolve()
    manifest_path = config.path_for("snapshot_manifest")
    try:
        manifest_path.relative_to(root)
    except ValueError as exc:
        raise RefreshSafetyError(
            "Refresh snapshot_manifest must be inside the new dated snapshot root"
        ) from exc
    if manifest_path.exists():
        raise RefreshSafetyError(
            f"Refusing to overwrite existing snapshot manifest: {manifest_path}"
        )
    return root


def require_refresh_output(config: ProjectConfig, path: Path) -> Path:
    """Validate that one refresh output is new and contained in its snapshot."""

    root = refresh_snapshot_root(config)
    resolved = path.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise RefreshSafetyError(
            f"Refresh output must be inside {root}: {resolved}"
        ) from exc
    if resolved.exists():
        raise RefreshSafetyError(f"Refusing to overwrite existing snapshot file: {resolved}")
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def write_bytes_exclusive(config: ProjectConfig, path: Path, payload: bytes) -> None:
    target = require_refresh_output(config, path)
    with target.open("xb") as handle:
        handle.write(payload)


def write_text_exclusive(
    config: ProjectConfig, path: Path, payload: str, *, encoding: str = "utf-8"
) -> None:
    write_bytes_exclusive(config, path, payload.encode(encoding))


def _relative_refresh_path(config: ProjectConfig, path: Path) -> str:
    refresh_root = refresh_snapshot_root(config)
    resolved = path.resolve()
    try:
        resolved.relative_to(refresh_root)
        return resolved.relative_to(config.root).as_posix()
    except ValueError as exc:
        raise RefreshSafetyError(
            f"Refresh dataset path must be inside {refresh_root}: {resolved}"
        ) from exc


def finalize_refresh_snapshot(config: ProjectConfig) -> dict[str, Any]:
    """Freeze fetched raw layers while marking derivation/analysis work pending."""

    refresh_snapshot_root(config)
    dataset_specs = {
        "chembl_reference_ligands": {
            "source": "ChEMBL",
            "source_version": config.value("refresh.chembl.source_release"),
            "paths": [
                _relative_refresh_path(config, config.path_for("reference_ligands")),
                _relative_refresh_path(config, config.path_for("chembl_cache")),
            ],
        },
        "pubchem_benchmark": {
            "source": "PubChem PUG REST plus curated mechanism labels",
            "source_version": None,
            "paths": [
                _relative_refresh_path(config, config.path_for("benchmark")),
                _relative_refresh_path(
                    config,
                    config.path_for("benchmark").with_name(
                        "eskape_benchmark_sources.json"
                    ),
                ),
            ],
        },
        "uniprot_species_targets": {
            "source": "UniProt REST API",
            "source_version": config.value("refresh.uniprot.source_release"),
            "paths": [
                _relative_refresh_path(config, config.path_for("species_proteins")),
                _relative_refresh_path(config, config.path_for("species_fasta")),
                _relative_refresh_path(config, config.path_for("species_metadata")),
            ],
        },
        "card_raw": {
            "source": "CARD data archive",
            "source_version": config.value("refresh.card.source_version"),
            "paths": [_relative_refresh_path(config, config.path_for("card_raw"))],
        },
        "rcsb_structure_catalog": {
            "source": "RCSB PDB search and data APIs",
            "source_version": None,
            "paths": [
                _relative_refresh_path(config, config.path_for("structure_candidates")),
                _relative_refresh_path(config, config.path_for("structure_summary")),
            ],
        },
    }
    datasets: dict[str, Any] = {}
    for name, spec in dataset_specs.items():
        record = hash_declared_paths(config.root, spec["paths"])
        datasets[name] = {**spec, **record}
    recorded_at = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    manifest = {
        "schema_version": 1,
        "snapshot_id": config.value("snapshots.snapshot_id"),
        "analysis_ready": False,
        "snapshot_recorded_at": recorded_at,
        "hash_algorithm": TREE_HASH_ALGORITHM,
        "status": "raw_refresh_complete_derived_layers_pending",
        "provenance_warning": (
            "Raw API refresh is frozen, but derived quality, compatibility, and "
            "resistance layers have not been built. This snapshot must not be used "
            "for reported analysis."
        ),
        "datasets": datasets,
    }
    write_text_exclusive(
        config,
        config.path_for("snapshot_manifest"),
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return manifest
