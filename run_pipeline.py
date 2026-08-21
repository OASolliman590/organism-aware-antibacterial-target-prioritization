"""Run the configured open-target-discovery workflow.

Private structures are optional and are read only from local ignored paths.
Public benchmark/reference/annotation modules use paths declared in config.yaml.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import subprocess
import sys

from pipeline.config import ProjectConfig, load_config, set_global_seed
from pipeline.provenance import start_run, utc_now, write_run_manifest


STEPS = [
    'pipeline/fetch_benchmark_structures.py',
    'pipeline/fetch_chembl_reference_subtypes_v21.py',
    'pipeline/build_reference_quality.py',
    'pipeline/fetch_card_data.py',
    'pipeline/fetch_species_targets.py',
    'pipeline/sequence_compatibility.py',
    'pipeline/build_card_resistance_annotations.py',
    'pipeline/parse_card_snps_v2.py',
    'pipeline/fetch_structure_catalog_v2.py',
    'pipeline/open_target_discovery_v2.py',
    'pipeline/benchmark_v2.py',
    'pipeline/calibrate_uncertainty_v2.py',
    'pipeline/build_validation_plan_v2.py',
    'pipeline/v2_figures.py',
    'pipeline/summarize_v2.py',
]


def child_environment(config: ProjectConfig) -> dict[str, str]:
    """Pass one config path; legacy variables carry paths only, never parameters."""

    env = os.environ.copy()
    env["OATP_CONFIG"] = str(config.path)
    env["PROJECT_ROOT"] = str(config.root)
    env["INPUT_DIR"] = str(config.path_for("inputs"))
    return env


def run_step(relative_path: str, *, config: ProjectConfig) -> None:
    step = config.root / relative_path
    print(f"\n=== {step.name} ===", flush=True)
    subprocess.run(
        [sys.executable, str(step)],
        check=True,
        cwd=config.root,
        env=child_environment(config),
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Single YAML run configuration (default: config.yaml)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = load_config(args.config)
    set_global_seed(config)
    run_context = start_run(config)
    try:
        for relative_path in STEPS:
            if (config.root / relative_path).exists():
                run_step(relative_path, config=config)
    except Exception as exc:
        write_run_manifest(
            config,
            run_context,
            status="failed",
            completed_at=utc_now(),
            error=exc,
        )
        raise
    write_run_manifest(
        config, run_context, status="completed", completed_at=utc_now()
    )
    print(
        "\nConfigured open-target-discovery pipeline completed. "
        f"Snapshot: {config.value('snapshots.snapshot_id')}. "
        f"See {config.path_for('results')} for outputs."
    )


if __name__ == "__main__":
    main()
