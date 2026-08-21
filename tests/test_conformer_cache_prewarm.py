import json

from rdkit import Chem

from pipeline.config import load_config
from pipeline.prewarm_conformer_cache import prewarm_molecules


def test_parallel_prewarm_is_deterministic_and_reuses_cache(tmp_path):
    config = load_config()
    molecules = [
        Chem.MolFromSmiles("CCO"),
        Chem.MolFromSmiles("c1ccccc1"),
        Chem.MolFromSmiles("CCO"),
    ]
    first = prewarm_molecules(
        molecules, config, cache_dir=tmp_path, workers=2
    ).iloc[0]
    manifests_before = {
        path.name: json.loads(path.read_text(encoding="utf-8"))
        for path in tmp_path.glob("*.json")
    }
    second = prewarm_molecules(
        list(reversed(molecules)), config, cache_dir=tmp_path, workers=2
    ).iloc[0]
    manifests_after = {
        path.name: json.loads(path.read_text(encoding="utf-8"))
        for path in tmp_path.glob("*.json")
    }

    assert first["n_unique_valid_structures"] == 2
    assert first["n_input_records"] == 3
    assert first["n_cache_misses"] == 2
    assert first["n_with_conformers"] == 2
    assert second["n_cache_hits"] == 2
    assert manifests_before == manifests_after
    assert all(
        manifest["parameters"]["num_threads"] == 1
        for manifest in manifests_after.values()
    )
