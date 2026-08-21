from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from pipeline.config import ConfigError, load_config, set_global_seed


ROOT = Path(__file__).resolve().parents[1]


def test_default_config_resolves_paths_and_declares_scientific_parameters() -> None:
    config = load_config(ROOT / "config.yaml")

    assert config.path_for("benchmark") == (
        ROOT / "data" / "benchmark" / "eskape_benchmark_drugs.csv"
    )
    assert config.value("chem2d.close_analogue_cutoff") == 0.85
    assert config.value("seeds.conformer") == 20240601
    assert config.value("benchmark.bedroc_alphas") == [20.0, 80.5]
    assert len(config.config_hash) == 64
    assert set_global_seed(config) == config.value("run.seed")


def test_config_rejects_overlapping_applicability_domain_thresholds(
    tmp_path: Path,
) -> None:
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["applicability_domain"]["tanimoto_out"] = raw[
        "applicability_domain"
    ]["tanimoto_in"]
    invalid = tmp_path / "invalid.yaml"
    invalid.write_text(yaml.safe_dump(raw), encoding="utf-8")

    with pytest.raises(ConfigError, match="tanimoto_out"):
        load_config(invalid, project_root=ROOT)
