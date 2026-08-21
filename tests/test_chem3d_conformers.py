from __future__ import annotations

from pathlib import Path

import numpy as np
from rdkit import Chem

from pipeline import chem3d_matching
from pipeline.chem3d_matching import generate_conformer_ensemble
from pipeline.config import load_config


ROOT = Path(__file__).resolve().parents[1]


def _coordinates(molecule: Chem.Mol) -> list[np.ndarray]:
    return [
        np.asarray(conformer.GetPositions(), dtype=float)
        for conformer in molecule.GetConformers()
    ]


def test_etkdgv3_conformers_are_deterministic_and_cached(tmp_path: Path) -> None:
    config = load_config(ROOT / "config.yaml")
    molecule = Chem.MolFromSmiles("COC1=CC=CC=C1")
    assert molecule is not None

    generated = generate_conformer_ensemble(molecule, config, cache_dir=tmp_path)
    cached = generate_conformer_ensemble(molecule, config, cache_dir=tmp_path)

    assert generated.status == "ok"
    assert generated.cache_hit is False
    assert cached.cache_hit is True
    assert generated.cache_key == cached.cache_key
    assert 1 <= generated.n_conformers <= config.value("chem3d.n_confs")
    assert generated.relative_energies_kcal == cached.relative_energies_kcal
    assert generated.optimization_statuses == cached.optimization_statuses
    assert generated.molecule is not None and cached.molecule is not None
    for first, second in zip(
        _coordinates(generated.molecule), _coordinates(cached.molecule)
    ):
        assert np.array_equal(first, second)
    assert generated.cache_path is not None
    assert generated.cache_path.is_file()
    assert generated.cache_path.with_name(f"{generated.cache_key}.json").is_file()


def test_embedding_failure_is_cached_as_explicit_missingness(
    tmp_path: Path, monkeypatch
) -> None:
    config = load_config(ROOT / "config.yaml")
    molecule = Chem.MolFromSmiles("CCO")
    assert molecule is not None
    calls = 0

    def fail_embedding(*args, **kwargs):
        nonlocal calls
        calls += 1
        return []

    monkeypatch.setattr(
        chem3d_matching.AllChem, "EmbedMultipleConfs", fail_embedding
    )
    first = generate_conformer_ensemble(molecule, config, cache_dir=tmp_path)
    second = generate_conformer_ensemble(molecule, config, cache_dir=tmp_path)

    assert calls == 1
    assert first.status == second.status == "embedding_failed"
    assert first.molecule is second.molecule is None
    assert first.cache_hit is False
    assert second.cache_hit is True
    assert len(list(tmp_path.glob("*.json"))) == 1
    assert not list(tmp_path.glob("*.rdkit.bin"))
