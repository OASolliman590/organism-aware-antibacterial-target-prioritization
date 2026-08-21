"""Reproducible PIDGINv4 external-baseline adapter with honest gap handling."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any

import pandas as pd

try:
    from pipeline.config import ProjectConfig, load_config
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig, load_config


HEAD_TO_HEAD_COLUMNS = [
    "query_id",
    "baseline_target_id",
    "target_class",
    "baseline_score",
    "baseline_provider",
    "baseline_code_commit",
    "baseline_model_release",
    "baseline_score_is_calibrated_probability",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit(path: Path) -> str | None:
    if not (path / ".git").exists():
        return None
    completed = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def probe_pidginv4(config: ProjectConfig) -> dict[str, Any]:
    """Report every prerequisite; never auto-install or silently change versions."""

    code_dir = config.path_for("pidgin_code")
    model_dir = config.path_for("pidgin_models")
    target_map = config.path_for("external_baseline_target_map")
    python2 = shutil.which(str(config.value("external_baseline.python_executable")))
    expected_commit = str(config.value("external_baseline.code_commit"))
    actual_commit = _git_commit(code_dir) if code_dir.is_dir() else None
    checks = {
        "python2_runtime_available": bool(python2),
        "code_checkout_available": (code_dir / "predict.py").is_file(),
        "code_commit_matches": actual_commit == expected_commit,
        "model_directory_available": model_dir.is_dir(),
        "model_pickles_available": (model_dir / "pkls").is_dir(),
        "model_ad_data_available": (model_dir / "ad_data").is_dir(),
        "model_metadata_available": (model_dir / "uniprot_information.txt").is_file(),
        "target_mapping_available": target_map.is_file(),
    }
    missing = [name for name, available in checks.items() if not available]
    return {
        "provider": "PIDGINv4",
        "status": "available" if not missing else "pending_unavailable_prerequisites",
        "status_reason": (
            "all pinned PIDGINv4 prerequisites are available"
            if not missing
            else ";".join(missing)
        ),
        "repository": str(config.value("external_baseline.repository")),
        "expected_code_commit": expected_commit,
        "actual_code_commit": actual_commit,
        "model_release": str(config.value("external_baseline.model_release")),
        "python_executable": python2,
        "code_dir": str(code_dir),
        "model_dir": str(model_dir),
        "target_mapping_path": str(target_map),
        "score_semantics": (
            "PIDGINv4 random-forest mean class score; official documentation says "
            "v4 outputs are not Platt-scaled"
        ),
        "score_is_calibrated_probability": False,
        **checks,
    }


def _validate_target_map(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path)
    required = {
        "baseline_target_id",
        "target_class",
        "mapping_source",
        "mapping_version",
    }
    missing = sorted(required - set(mapping.columns))
    if missing:
        raise ValueError(f"External-baseline target map is missing: {missing}")
    if mapping[list(required)].isna().any().any():
        raise ValueError("External-baseline target map contains blank provenance")
    if mapping["baseline_target_id"].duplicated().any():
        raise ValueError("baseline_target_id mappings must be unique")
    return mapping


def parse_pidginv4_predictions(
    path: str | Path,
    target_map_path: str | Path,
    *,
    code_commit: str,
    model_release: str,
) -> pd.DataFrame:
    """Map a PIDGIN target×compound matrix through a versioned explicit map."""

    prediction_path = Path(path)
    predictions = pd.read_csv(prediction_path, sep="\t")
    if predictions.empty or len(predictions.columns) < 2:
        raise ValueError("PIDGINv4 prediction matrix is empty or malformed")
    target_id_column = predictions.columns[0]
    long = predictions.melt(
        id_vars=[target_id_column],
        var_name="query_id",
        value_name="baseline_score",
    ).rename(columns={target_id_column: "baseline_target_id"})
    long["baseline_score"] = pd.to_numeric(long["baseline_score"], errors="coerce")
    mapping = _validate_target_map(Path(target_map_path))
    mapped = long.merge(
        mapping[["baseline_target_id", "target_class"]],
        on="baseline_target_id",
        how="inner",
        validate="many_to_one",
    )
    mapped["baseline_provider"] = "PIDGINv4"
    mapped["baseline_code_commit"] = code_commit
    mapped["baseline_model_release"] = model_release
    mapped["baseline_score_is_calibrated_probability"] = False
    return mapped[HEAD_TO_HEAD_COLUMNS]


def run_or_load_pidginv4(
    queries: pd.DataFrame, config: ProjectConfig
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run the exact pinned CLI or load a hash-verified cached raw matrix."""

    status = probe_pidginv4(config)
    if status["status"] != "available":
        return pd.DataFrame(columns=HEAD_TO_HEAD_COLUMNS), status
    required = {"query_id", "canonical_smiles"}
    missing = sorted(required - set(queries.columns))
    if missing:
        raise ValueError(f"External-baseline queries are missing: {missing}")
    cache_dir = config.path_for("external_baseline_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    input_path = cache_dir / "benchmark_queries.smi"
    output_path = cache_dir / "pidginv4_predictions.tsv"
    manifest_path = cache_dir / "pidginv4_predictions.manifest.json"
    input_text = "".join(
        f"{row.canonical_smiles} {row.query_id}\n" for row in queries.itertuples()
    )
    input_path.write_text(input_text, encoding="utf-8", newline="\n")
    input_hash = sha256_file(input_path)
    expected_manifest = {
        "provider": "PIDGINv4",
        "code_commit": status["expected_code_commit"],
        "model_release": status["model_release"],
        "query_sha256": input_hash,
    }
    cache_hit = False
    if output_path.is_file() and manifest_path.is_file():
        saved = json.loads(manifest_path.read_text(encoding="utf-8"))
        cache_hit = all(saved.get(key) == value for key, value in expected_manifest.items())
        cache_hit = cache_hit and saved.get("output_sha256") == sha256_file(output_path)
    if not cache_hit:
        command = [
            str(status["python_executable"]),
            str(Path(str(status["code_dir"])) / "predict.py"),
            "-f",
            str(input_path),
            "-o",
            str(output_path),
            "-n",
            "1",
            "--ad",
            str(config.value("external_baseline.applicability_domain_percentile")),
            "--model_dir",
            str(status["model_dir"]),
        ]
        subprocess.run(command, check=True, cwd=status["code_dir"])
        expected_manifest["output_sha256"] = sha256_file(output_path)
        manifest_path.write_text(
            json.dumps(expected_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    parsed = parse_pidginv4_predictions(
        output_path,
        config.path_for("external_baseline_target_map"),
        code_commit=str(status["expected_code_commit"]),
        model_release=str(status["model_release"]),
    )
    return parsed, {**status, "cache_hit": cache_hit, "n_mapped_scores": len(parsed)}


def main() -> None:
    config = load_config()
    results_dir = config.path_for("results")
    results_dir.mkdir(parents=True, exist_ok=True)
    queries = pd.read_csv(config.path_for("benchmark")).rename(
        columns={"drug": "query_id"}
    )
    scores, status = run_or_load_pidginv4(queries, config)
    scores.to_csv(results_dir / "external_baseline_head_to_head_v3.csv", index=False)
    pd.DataFrame([status]).to_csv(
        results_dir / "external_baseline_status_v3.csv", index=False
    )
    print(pd.DataFrame([status]).to_string(index=False))


if __name__ == "__main__":
    main()
