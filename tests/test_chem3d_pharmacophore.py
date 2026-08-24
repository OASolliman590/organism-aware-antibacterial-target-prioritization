from __future__ import annotations

from copy import deepcopy

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms

import pipeline.chem3d_matching as chem3d_matching

from pipeline.chem3d_matching import (
    ConformerEnsemble,
    Pharmacophore3DScore,
    gobbi_pharmacophore_fingerprint,
    o3a_shape_color_pharmacophore_similarity,
    pharmacophore_3d_similarity,
    score_pharmacophore_by_target,
    score_reference_evidence_by_target,
)
from pipeline.config import ProjectConfig, load_config


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["chem3d"]["n_confs"] = 2
    data["chem3d"]["o3a_shortlist_top"] = 1
    return ProjectConfig(path=base.path, root=base.root, data=data)


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

    assert scores["pharmacophore_2d_gobbi_sim_max"].between(0, 1).all()
    assert scores["pharmacophore_sim_max"].between(0, 1).all()
    assert scores["pharmacophore_sim_max"].equals(
        scores["pharmacophore_2d_gobbi_sim_max"]
    )
    assert scores.loc["same", "pharmacophore_sim_max"] == 1.0
    assert scores.loc["same", "n_pharmacophore_references_scored"] == 2
    assert set(scores["pharmacophore_method"]) == {"Gobbi_Pharm2D"}
    assert set(scores["pharmacophore_status"]) == {"ok"}


def _embedded_ethanol() -> Chem.Mol:
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 41
    assert AllChem.EmbedMolecule(molecule, parameters) == 0
    return molecule


def _embedded_flexible_molecule() -> Chem.Mol:
    molecule = Chem.AddHs(Chem.MolFromSmiles("NCCCCO"))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 41
    assert AllChem.EmbedMolecule(molecule, parameters) == 0
    return molecule


def _single_conformer_ensemble(molecule: Chem.Mol) -> ConformerEnsemble:
    return ConformerEnsemble(
        molecule=molecule,
        relative_energies_kcal=(0.0,),
        optimization_statuses=(0,),
        cache_key="fixture",
        cache_path=None,
        cache_hit=False,
        status="ok",
        n_embedded=1,
    )


def test_3d_pharmacophore_changes_with_internal_geometry_while_gobbi_does_not() -> None:
    config = _config()
    query = _embedded_flexible_molecule()
    distorted = Chem.Mol(query)
    conformer = distorted.GetConformer(0)
    initial_dihedral = rdMolTransforms.GetDihedralDeg(conformer, 0, 1, 2, 3)
    rdMolTransforms.SetDihedralDeg(
        conformer, 0, 1, 2, 3, initial_dihedral + 120.0
    )

    identical = pharmacophore_3d_similarity(query, 0, query, 0, config)
    distorted_score = pharmacophore_3d_similarity(
        query, 0, distorted, 0, config
    )
    query_gobbi = gobbi_pharmacophore_fingerprint(query)
    distorted_gobbi = gobbi_pharmacophore_fingerprint(distorted)

    assert identical.status == "ok"
    assert identical.similarity == 1.0
    assert distorted_score.status == "ok"
    assert distorted_score.similarity is not None
    assert 0 <= distorted_score.similarity < identical.similarity
    assert query_gobbi == distorted_gobbi


def test_o3a_3d_pharmacophore_is_invariant_to_initial_rigid_translation() -> None:
    config = _config()
    query = _embedded_ethanol()
    translated = Chem.Mol(query)
    conformer = translated.GetConformer(0)
    for atom_index in range(translated.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom_index)
        conformer.SetAtomPosition(
            atom_index,
            Geometry.Point3D(point.x + 10.0, point.y - 4.0, point.z + 2.0),
        )

    evidence = o3a_shape_color_pharmacophore_similarity(
        _single_conformer_ensemble(query),
        _single_conformer_ensemble(translated),
        config,
    )

    assert evidence.shape_similarity is not None
    assert evidence.shape_similarity > 0.99
    assert evidence.pharmacophore_3d_similarity is not None
    assert evidence.pharmacophore_3d_similarity > 0.99


def test_color_failure_does_not_suppress_shape_or_3d_pharmacophore(
    monkeypatch,
) -> None:
    config = _config()
    molecule = _embedded_ethanol()

    def fail_color(*args, **kwargs):
        raise RuntimeError("fixture color failure")

    monkeypatch.setattr(chem3d_matching.rdShapeAlign, "ScoreMol", fail_color)
    evidence = o3a_shape_color_pharmacophore_similarity(
        _single_conformer_ensemble(molecule),
        _single_conformer_ensemble(Chem.Mol(molecule)),
        config,
    )

    assert evidence.shape_similarity is not None
    assert evidence.color_similarity is None
    assert evidence.pharmacophore_3d_similarity is not None
    assert evidence.overlay_successes == 0
    assert evidence.overlay_failures == 1
    assert evidence.pharmacophore_3d_successes == 1
    assert evidence.pharmacophore_3d_failures == 0


def test_feature_failure_does_not_suppress_shape_or_color(monkeypatch) -> None:
    config = _config()
    molecule = _embedded_ethanol()

    def fail_features(*args, **kwargs):
        return Pharmacophore3DScore(None, "feature_extraction_failed", 0, 0)

    monkeypatch.setattr(
        chem3d_matching, "pharmacophore_3d_similarity", fail_features
    )
    evidence = o3a_shape_color_pharmacophore_similarity(
        _single_conformer_ensemble(molecule),
        _single_conformer_ensemble(Chem.Mol(molecule)),
        config,
    )

    assert evidence.shape_similarity is not None
    assert evidence.color_similarity is not None
    assert evidence.pharmacophore_3d_similarity is None
    assert evidence.overlay_successes == 1
    assert evidence.overlay_failures == 0
    assert evidence.pharmacophore_3d_successes == 0
    assert evidence.pharmacophore_3d_failures == 1
    assert evidence.pharmacophore_3d_status == "feature_extraction_failed"


def test_3d_pharmacophore_no_features_is_missing_not_zero() -> None:
    config = _config()
    query = _embedded_ethanol()
    helium = Chem.MolFromSmiles("[He]")
    assert helium is not None
    helium.AddConformer(Chem.Conformer(helium.GetNumAtoms()), assignId=True)

    result = pharmacophore_3d_similarity(query, 0, helium, 0, config)

    assert result.similarity is None
    assert result.status == "unavailable_no_features"
    assert result.reference_feature_count == 0


def test_reference_evidence_keeps_2d_and_o3a_aligned_3d_pharmacophores(
    tmp_path,
) -> None:
    config = _config()
    query = Chem.MolFromSmiles("CCO")
    other = Chem.MolFromSmiles("CC(=O)N")
    assert query is not None and other is not None
    evidence = score_reference_evidence_by_target(
        "fixture-query",
        query,
        {"fixture": [{"_mol": query}, {"_mol": other}]},
        config,
        cache_dir=tmp_path / "conformers",
    )

    assert evidence["pharmacophore_2d_gobbi_similarity"].notna().all()
    assert evidence["pharmacophore_similarity"].equals(
        evidence["pharmacophore_2d_gobbi_similarity"]
    )
    assert int(evidence["o3a_was_shortlisted"].sum()) == 1
    shortlisted = evidence[evidence["o3a_was_shortlisted"]]
    not_shortlisted = evidence[~evidence["o3a_was_shortlisted"]]
    assert shortlisted["pharmacophore_3d_similarity"].between(0, 1).all()
    assert shortlisted["pharmacophore_3d_status"].eq("ok").all()
    assert not_shortlisted["pharmacophore_3d_similarity"].isna().all()
    assert not_shortlisted["pharmacophore_3d_status"].eq("not_shortlisted").all()
    assert shortlisted["pharmacophore_3d_feature_definition_sha256"].str.len().eq(64).all()
