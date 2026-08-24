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


@pytest.mark.parametrize(
    ("key", "value", "message"),
    [
        ("pharmacophore_3d_feature_definition", "custom.fdef", "BaseFeatures"),
        ("pharmacophore_3d_score_mode", "all", "score_mode"),
        ("pharmacophore_3d_profile", "triangle", "profile"),
        ("pharmacophore_3d_radius", 0, "radius"),
        ("pharmacophore_3d_width", 0, "width"),
        ("pharmacophore_3d_direction_mode", "dot", "direction_mode"),
        ("pharmacophore_3d_normalization", "asymmetric", "normalization"),
    ],
)
def test_config_rejects_unpinned_3d_pharmacophore_options(
    tmp_path: Path, key: str, value, message: str
) -> None:
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["chem3d"][key] = value
    invalid = tmp_path / f"invalid-{key}.yaml"
    invalid.write_text(yaml.safe_dump(raw), encoding="utf-8")

    with pytest.raises(ConfigError, match=message):
        load_config(invalid, project_root=ROOT)


def test_config_rejects_3d_benchmark_components_outside_fusion(
    tmp_path: Path,
) -> None:
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["benchmark"]["three_dimensional_components"].append("unknown_3d_score")
    invalid = tmp_path / "invalid-3d-components.yaml"
    invalid.write_text(yaml.safe_dump(raw), encoding="utf-8")

    with pytest.raises(ConfigError, match="must be fusion components"):
        load_config(invalid, project_root=ROOT)
