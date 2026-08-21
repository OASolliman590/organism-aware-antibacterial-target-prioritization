from __future__ import annotations

from rdkit import Chem

from pipeline.chem3d_matching import score_pharmacophore_by_target


def test_gobbi_pharmacophore_similarity_is_auditable_and_bounded() -> None:
    query = Chem.MolFromSmiles("COC1=CC(=CC(=C1OC)OC)CC2=CN=C(N=C2N)N")
    other = Chem.MolFromSmiles("CC(=O)N")
    assert query is not None and other is not None
    references = {
        "same": [{"_mol": query}, {"_mol": other}],
        "other": [{"_mol": other}],
    }

    scores = score_pharmacophore_by_target(
        "fixture-query", query, references
    ).set_index("target_class")

    assert scores["pharmacophore_sim_max"].between(0, 1).all()
    assert scores.loc["same", "pharmacophore_sim_max"] == 1.0
    assert scores.loc["same", "n_pharmacophore_references_scored"] == 2
    assert set(scores["pharmacophore_method"]) == {"Gobbi_Pharm2D"}
    assert set(scores["pharmacophore_status"]) == {"ok"}
