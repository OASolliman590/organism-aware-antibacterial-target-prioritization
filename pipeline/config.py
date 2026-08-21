"""Load and validate the single declarative configuration for v2/v3 runs."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import random
from typing import Any

import yaml


DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[1] / "config.yaml"
CONFIG_PATH_ENV = "OATP_CONFIG"


class ConfigError(ValueError):
    """Raised when a run configuration is missing or scientifically invalid."""


@dataclass(frozen=True)
class ProjectConfig:
    """Validated configuration plus deterministic path/hash helpers."""

    path: Path
    root: Path
    data: dict[str, Any]

    def value(self, dotted_key: str) -> Any:
        value: Any = self.data
        for key in dotted_key.split("."):
            if not isinstance(value, dict) or key not in value:
                raise ConfigError(f"Missing required config key: {dotted_key}")
            value = value[key]
        return value

    def path_for(self, key: str) -> Path:
        raw = Path(str(self.value(f"paths.{key}")))
        return raw if raw.is_absolute() else (self.root / raw).resolve()

    @property
    def config_hash(self) -> str:
        canonical = json.dumps(
            self.data, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        return hashlib.sha256(canonical).hexdigest()


def _require_number(
    config: ProjectConfig,
    key: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> float:
    value = config.value(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{key} must be numeric")
    number = float(value)
    if minimum is not None and number < minimum:
        raise ConfigError(f"{key} must be >= {minimum}")
    if maximum is not None and number > maximum:
        raise ConfigError(f"{key} must be <= {maximum}")
    return number


def _validate(config: ProjectConfig) -> None:
    required = [
        "schema_version",
        "run.seed",
        "run.fusion_mode",
        "run.combiner",
        "run.refresh_external_data",
        "snapshots.snapshot_id",
        "paths.snapshot_manifest",
        "paths.results",
        "paths.reference_ligands",
        "paths.benchmark",
        "chem2d.fingerprint_radius",
        "chem2d.fingerprint_bits",
        "chem2d.close_analogue_cutoff",
        "chem2d.top_k",
        "chem3d.n_confs",
        "chem3d.prewarm_workers",
        "chem3d.scoring_workers",
        "chem3d.prune_rms",
        "chem3d.energy_window_kcal",
        "chem3d.o3a_shortlist_top",
        "chem3d.aggregate_top_k",
        "chem3d.o3a_atom_typing",
        "chem3d.o3a_max_iterations",
        "chem3d.max_iterations",
        "chem3d.num_threads",
        "fusion.reciprocal_rank_constant",
        "fusion.disagreement_min_absolute_rank_shift",
        "fusion.components",
        "benchmark.bootstrap_n",
        "benchmark.bedroc_alphas",
        "benchmark.enrichment_fractions",
        "benchmark.splits",
        "benchmark.analogue_exclusion_threshold",
        "benchmark.time_cutoff",
        "applicability_domain.tanimoto_in",
        "applicability_domain.tanimoto_out",
        "external_baseline.provider",
        "external_baseline.repository",
        "external_baseline.code_commit",
        "external_baseline.model_release",
        "external_baseline.python_executable",
        "external_baseline.applicability_domain_percentile",
        "model.features",
        "model.regularization_c",
        "model.platt_regularization_c",
        "model.calibration_method",
        "model.minimum_class_count",
        "model.max_iterations",
        "model.reliability_bins",
        "sensitivity.weight_factors",
        "sensitivity.top_k",
        "sensitivity.rbo_persistence",
        "sensitivity.bootstrap_n",
    ]
    for key in required:
        config.value(key)

    if int(config.value("schema_version")) != 1:
        raise ConfigError("Only config schema_version 1 is supported")
    if config.value("run.fusion_mode") not in {"rank_fusion", "score_fusion"}:
        raise ConfigError("run.fusion_mode must be rank_fusion or score_fusion")
    if config.value("run.combiner") not in {"heuristic", "learned"}:
        raise ConfigError("run.combiner must be heuristic or learned")
    if not isinstance(config.value("run.refresh_external_data"), bool):
        raise ConfigError("run.refresh_external_data must be true or false")

    _require_number(config, "run.seed", minimum=0)
    _require_number(config, "chem2d.fingerprint_radius", minimum=1)
    _require_number(config, "chem2d.fingerprint_bits", minimum=128)
    _require_number(config, "chem2d.close_analogue_cutoff", minimum=0, maximum=1)
    _require_number(config, "chem2d.top_k", minimum=1)
    _require_number(config, "chem3d.n_confs", minimum=1)
    _require_number(config, "chem3d.prewarm_workers", minimum=1)
    _require_number(config, "chem3d.scoring_workers", minimum=1)
    _require_number(config, "chem3d.prune_rms", minimum=0)
    _require_number(config, "chem3d.energy_window_kcal", minimum=0)
    _require_number(config, "chem3d.o3a_shortlist_top", minimum=1)
    _require_number(config, "chem3d.aggregate_top_k", minimum=1)
    if config.value("chem3d.o3a_atom_typing") not in {"mmff94", "crippen"}:
        raise ConfigError("chem3d.o3a_atom_typing must be mmff94 or crippen")
    _require_number(config, "chem3d.o3a_max_iterations", minimum=1)
    _require_number(config, "chem3d.max_iterations", minimum=1)
    if int(_require_number(config, "chem3d.num_threads", minimum=1)) != 1:
        raise ConfigError("chem3d.num_threads must be 1 for deterministic runs")
    _require_number(config, "fusion.reciprocal_rank_constant", minimum=1e-12)
    _require_number(
        config, "fusion.disagreement_min_absolute_rank_shift", minimum=1e-12
    )
    fusion_components = config.value("fusion.components")
    if (
        not isinstance(fusion_components, list)
        or not fusion_components
        or not all(isinstance(component, str) and component for component in fusion_components)
        or len(fusion_components) != len(set(fusion_components))
    ):
        raise ConfigError("fusion.components must be a non-empty unique string list")
    _require_number(config, "benchmark.bootstrap_n", minimum=1)
    bedroc_alphas = config.value("benchmark.bedroc_alphas")
    if not isinstance(bedroc_alphas, list) or {
        float(value) for value in bedroc_alphas
    } != {20.0, 80.5}:
        raise ConfigError("benchmark.bedroc_alphas must contain 20.0 and 80.5")
    enrichment_fractions = config.value("benchmark.enrichment_fractions")
    if not isinstance(enrichment_fractions, list) or {
        float(value) for value in enrichment_fractions
    } != {0.01, 0.05}:
        raise ConfigError(
            "benchmark.enrichment_fractions must contain 0.01 and 0.05"
        )
    benchmark_splits = config.value("benchmark.splits")
    if (
        not isinstance(benchmark_splits, list)
        or set(benchmark_splits) != {"target_family", "scaffold", "temporal"}
        or len(benchmark_splits) != 3
    ):
        raise ConfigError(
            "benchmark.splits must contain target_family, scaffold, and temporal exactly once"
        )
    _require_number(
        config, "benchmark.analogue_exclusion_threshold", minimum=0, maximum=1
    )
    tanimoto_in = _require_number(
        config, "applicability_domain.tanimoto_in", minimum=0, maximum=1
    )
    tanimoto_out = _require_number(
        config, "applicability_domain.tanimoto_out", minimum=0, maximum=1
    )
    if tanimoto_out >= tanimoto_in:
        raise ConfigError(
            "applicability_domain.tanimoto_out must be below tanimoto_in"
        )
    if config.value("external_baseline.provider") != "PIDGINv4":
        raise ConfigError("external_baseline.provider must be PIDGINv4")
    _require_number(
        config,
        "external_baseline.applicability_domain_percentile",
        minimum=0,
        maximum=100,
    )
    model_features = config.value("model.features")
    if (
        not isinstance(model_features, list)
        or not model_features
        or len(model_features) != len(set(model_features))
    ):
        raise ConfigError("model.features must be a non-empty unique list")
    _require_number(config, "model.regularization_c", minimum=1e-12)
    _require_number(config, "model.platt_regularization_c", minimum=1e-12)
    if config.value("model.calibration_method") != "Platt":
        raise ConfigError("model.calibration_method must be Platt")
    _require_number(config, "model.minimum_class_count", minimum=2)
    _require_number(config, "model.max_iterations", minimum=1)
    _require_number(config, "model.reliability_bins", minimum=2)
    sensitivity_factors = config.value("sensitivity.weight_factors")
    if not isinstance(sensitivity_factors, list) or {
        float(value) for value in sensitivity_factors
    } != {0.5, 0.75, 1.25, 1.5}:
        raise ConfigError(
            "sensitivity.weight_factors must contain 0.5, 0.75, 1.25, and 1.5"
        )
    _require_number(config, "sensitivity.top_k", minimum=1)
    _require_number(
        config, "sensitivity.rbo_persistence", minimum=1e-12, maximum=0.999999
    )
    _require_number(config, "sensitivity.bootstrap_n", minimum=1)

    weights = config.value("v2_scoring.chemical.weights")
    if not isinstance(weights, dict) or any(float(v) < 0 for v in weights.values()):
        raise ConfigError("v2_scoring.chemical.weights must be non-negative")
    if sum(float(v) for v in weights.values()) > 1.0 + 1e-12:
        raise ConfigError("v2_scoring.chemical.weights must sum to at most 1")


def load_config(
    path: str | Path | None = None, *, project_root: str | Path | None = None
) -> ProjectConfig:
    """Load a YAML config; only the config path may be selected via environment."""

    selected = path or os.environ.get(CONFIG_PATH_ENV) or DEFAULT_CONFIG_PATH
    config_path = Path(selected).expanduser().resolve()
    if not config_path.is_file():
        raise ConfigError(f"Config file does not exist: {config_path}")
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ConfigError("The config root must be a YAML mapping")
    root = (
        Path(project_root).expanduser().resolve()
        if project_root is not None
        else config_path.parent
    )
    config = ProjectConfig(path=config_path, root=root, data=raw)
    _validate(config)
    return config


def set_global_seed(config: ProjectConfig) -> int:
    """Seed Python and NumPy from the declared global seed."""

    seed = int(config.value("run.seed"))
    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass
    return seed
