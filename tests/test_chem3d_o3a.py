from __future__ import annotations

from pathlib import Path

import yaml
from rdkit import Chem

from pipeline.chem3d_matching import score_o3a_by_target
from pipeline.config import load_config


ROOT = Path(__file__).resolve().parents[1]


def _fast_config(tmp_path: Path):
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["chem3d"]["n_confs"] = 3
    raw["chem3d"]["o3a_shortlist_top"] = 1
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    return load_config(path, project_root=ROOT)


def test_o3a_shortlist_shape_and_color_scores_are_bounded(tmp_path: Path) -> None:
    config = _fast_config(tmp_path)
    query = Chem.MolFromSmiles("COC1=CC=CC=C1")
    other = Chem.MolFromSmiles("CC(=O)N")
    assert query is not None and other is not None
    references = {
        "same": [{"_mol": query}, {"_mol": other}],
        "other": [{"_mol": other}],
    }

    scores = score_o3a_by_target(
        "fixture-query", query, references, config, cache_dir=tmp_path / "cache"
    ).set_index("target_class")

    assert scores["o3a_shape_tanimoto_max"].between(0, 1).all()
    assert scores["o3a_color_max"].between(0, 1).all()
    assert (scores["n_o3a_references_shortlisted"] == 1).all()
    assert (scores["o3a_status"] == "ok").all()
    assert scores.loc["same", "o3a_shape_tanimoto_max"] > 0.99
    assert scores.loc["same", "o3a_color_max"] > 0.99
