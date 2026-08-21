from __future__ import annotations

from copy import deepcopy

import pandas as pd
from rdkit import Chem

from pipeline.config import ProjectConfig
from pipeline import open_target_discovery_v2 as discovery


def _fast_config(tmp_path) -> ProjectConfig:
    data = deepcopy(discovery.CONFIG.data)
    data["chem3d"].update(
        {
            "n_confs": 2,
            "o3a_shortlist_top": 1,
            "aggregate_top_k": 1,
            "max_iterations": 100,
            "cache_dir": str(tmp_path / "conformers"),
        }
    )
    return ProjectConfig(
        path=discovery.CONFIG.path,
        root=discovery.CONFIG.root,
        data=data,
    )


def _record(smiles: str, identifier: str):
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return {
        "molecule_chembl_id": identifier,
        "canonical_smiles": smiles,
        "_mol": molecule,
        "_fp": discovery.fp(molecule),
        "_maccs": discovery.maccs(molecule),
        "_smi": Chem.MolToSmiles(molecule),
    }


def _inputs():
    molecule = discovery.mol("COC1=CC(=CC(=C1OC)OC)CC2=CN=C(N=C2N)N")
    assert molecule is not None
    query = {
        "query_id": "public-query",
        "query_name": "public-query",
        "mol": molecule,
        "fp": discovery.fp(molecule),
        "maccs": discovery.maccs(molecule),
        "source": "test_public_fixture",
    }
    references = {
        "class_a": [_record("CC(=O)NC1=CC=CC=C1", "ref-a")],
        "class_b": [_record("CC(=O)OC1=CC=CC=C1C(=O)O", "ref-b")],
    }
    quality = pd.DataFrame(columns=["target_class"])
    empty = pd.DataFrame()
    ontology = pd.DataFrame(columns=["target_class"])
    return query, references, quality, empty, ontology


def test_v3_query_scoring_is_additive_and_deterministic(tmp_path) -> None:
    query, references, quality, empty, ontology = _inputs()
    config = _fast_config(tmp_path)

    first = discovery.score_query_v3(
        query,
        references,
        quality,
        empty,
        ontology,
        config=config,
        cache_dir=tmp_path / "conformers",
    )
    second = discovery.score_query_v3(
        query,
        references,
        quality,
        empty,
        ontology,
        config=config,
        cache_dir=tmp_path / "conformers",
    )
    third, reference_evidence = discovery.score_query_v3(
        query,
        references,
        quality,
        empty,
        ontology,
        config=config,
        cache_dir=tmp_path / "conformers",
        return_reference_evidence=True,
    )

    pd.testing.assert_frame_equal(first, second, check_exact=True)
    pd.testing.assert_frame_equal(first, third, check_exact=True)
    expected = {
        "chemical_evidence_score",
        "usrcat_max",
        "usrcat_top5_mean",
        "o3a_shape_tanimoto_max",
        "o3a_color_max",
        "pharmacophore_sim_max",
        "chemical_evidence_score_v3",
        "chemical_quality_adjusted_score_v3",
    }
    assert expected.issubset(first.columns)
    assert len(first) == 2
    assert first["chemical_evidence_score_v3"].between(0, 1).all()
    assert len(reference_evidence) == 2
    assert {
        "reference_id",
        "ecfp4_similarity",
        "usrcat_similarity",
        "o3a_was_shortlisted",
        "pharmacophore_similarity",
    }.issubset(reference_evidence.columns)


def test_emit_v3_outputs_writes_suffix_without_touching_v2_files(tmp_path) -> None:
    query, references, quality, empty, ontology = _inputs()
    config = _fast_config(tmp_path)
    result_dir = tmp_path / "results"

    written = discovery.emit_v3_outputs(
        [],
        [query],
        references,
        quality,
        ontology,
        empty,
        config=config,
        result_dir=result_dir,
        cache_dir=tmp_path / "conformers",
    )

    assert {path.name for path in written} == {
        "benchmark_open_target_scores_v3.csv",
        "chemical_evidence_disagreements_v3.csv",
        "open_target_reference_coverage_v3.csv",
    }
    assert all(path.is_file() and path.name.endswith("_v3.csv") for path in written)
    assert not list(result_dir.glob("v2_*.csv"))
