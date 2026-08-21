"""Emit an auditable provenance manifest for every configured pipeline run."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any

try:
    from pipeline.config import ProjectConfig
except ModuleNotFoundError:  # imported by direct ``python pipeline/<script>.py`` execution
    from config import ProjectConfig


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_output(root: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def git_provenance(root: Path) -> dict[str, Any]:
    commit = _git_output(root, "rev-parse", "HEAD")
    branch = _git_output(root, "branch", "--show-current")
    status = _git_output(root, "status", "--porcelain", "--untracked-files=no")
    return {
        "commit": commit,
        "branch": branch,
        "dirty_tracked_files": None if status is None else bool(status),
        "status": "available" if commit else "unavailable",
    }


def package_versions() -> dict[str, str]:
    """Record the complete installed distribution set, not only direct imports."""

    versions = {
        str(distribution.metadata["Name"]): distribution.version
        for distribution in metadata.distributions()
        if distribution.metadata["Name"]
    }
    return dict(sorted(versions.items(), key=lambda item: item[0].casefold()))


def snapshot_provenance(config: ProjectConfig) -> dict[str, Any]:
    manifest_path = config.path_for("snapshot_manifest")
    result: dict[str, Any] = {
        "configured_snapshot_id": config.value("snapshots.snapshot_id"),
        "manifest_path": str(manifest_path),
    }
    if not manifest_path.is_file():
        result.update(
            {
                "status": "manifest_missing",
                "manifest_sha256": None,
                "versions": None,
            }
        )
        return result
    try:
        versions = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        result.update(
            {
                "status": "manifest_unreadable",
                "manifest_sha256": sha256_file(manifest_path),
                "versions": None,
                "error": f"{type(exc).__name__}: {exc}",
            }
        )
        return result
    result.update(
        {
            "status": "available",
            "manifest_sha256": sha256_file(manifest_path),
            "versions": versions,
        }
    )
    return result


@dataclass(frozen=True)
class RunContext:
    run_id: str
    started_at: str


def start_run(config: ProjectConfig) -> RunContext:
    started_at = utc_now()
    compact_time = started_at.replace("-", "").replace(":", "").replace(".", "")
    run_id = f"{compact_time}-{config.config_hash[:12]}"
    context = RunContext(run_id=run_id, started_at=started_at)
    write_run_manifest(config, context, status="running")
    return context


def build_run_manifest(
    config: ProjectConfig,
    context: RunContext,
    *,
    status: str,
    completed_at: str | None = None,
    error: BaseException | None = None,
) -> dict[str, Any]:
    if status not in {"running", "completed", "failed"}:
        raise ValueError(f"Unsupported run status: {status}")
    return {
        "schema_version": 1,
        "run_id": context.run_id,
        "status": status,
        "timestamps": {
            "started_at_utc": context.started_at,
            "completed_at_utc": completed_at,
        },
        "config": {
            "path": str(config.path),
            "sha256": config.config_hash,
            "schema_version": config.value("schema_version"),
        },
        "seeds": {
            "global": config.value("run.seed"),
            **config.value("seeds"),
        },
        "data_snapshots": snapshot_provenance(config),
        "code": git_provenance(config.root),
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "packages": package_versions(),
        },
        "command": [str(item) for item in sys.argv],
        "error": (
            None
            if error is None
            else {"type": type(error).__name__, "message": str(error)}
        ),
    }


def _atomic_json_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def write_run_manifest(
    config: ProjectConfig,
    context: RunContext,
    *,
    status: str,
    completed_at: str | None = None,
    error: BaseException | None = None,
) -> dict[str, Any]:
    """Write both the required latest manifest and an immutable per-run copy."""

    manifest = build_run_manifest(
        config,
        context,
        status=status,
        completed_at=completed_at,
        error=error,
    )
    results = config.path_for("results")
    _atomic_json_write(results / "run_manifest.json", manifest)
    _atomic_json_write(results / "run_manifests" / f"{context.run_id}.json", manifest)
    return manifest
