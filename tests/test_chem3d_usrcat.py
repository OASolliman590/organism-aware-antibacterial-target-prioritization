from __future__ import annotations

from pathlib import Path

import yaml
from rdkit import Chem

from pipeline.chem3d_matching import score_usrcat_by_target
from pipeline.config import load_config


ROOT = Path(__file__).resolve().parents[1]


def _fast_config(tmp_path: Path):
    raw = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    raw["chem3d"]["n_confs"] = 4
    raw["chem3d"]["aggregate_top_k"] = 2
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    return load_config(path, project_root=ROOT)


def test_usrcat_scores_and_aggregates_every_query_target_pair(
    tmp_path: Path,
) -> None:
    config = _fast_config(tmp_path)
    query = Chem.MolFromSmiles("COC1=CC=CC=C1")
    other = Chem.MolFromSmiles("CC(=O)N")
    assert query is not None and other is not None
    references = {
        "same_shape": [{"_mol": query}, {"_mol": other}],
        "other_shape": [{"_mol": other}],
    }

    scores = score_usrcat_by_target(
        "fixture-query", query, references, config, cache_dir=tmp_path / "cache"
    )

    assert set(scores["target_class"]) == {"same_shape", "other_shape"}
    assert scores["usrcat_max"].between(0, 1).all()
    assert scores["usrcat_top5_mean"].between(0, 1).all()
    same = scores.set_index("target_class").loc["same_shape"]
    assert same["usrcat_max"] == 1.0
    assert same["n_usrcat_references_scored"] == 2
    assert same["query_conformer_status"] == "ok"
