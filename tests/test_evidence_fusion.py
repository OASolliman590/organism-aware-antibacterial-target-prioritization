from __future__ import annotations

import math

import numpy as np
import pandas as pd

from pipeline.config import load_config
from pipeline.evidence_fusion import fuse_evidence, reciprocal_rank_fusion


def test_rank_fusion_retains_components_and_exposes_contributions() -> None:
    chem2d = pd.DataFrame(
        [
            {"query_id": "q1", "target_class": "A", "ecfp4_max": 0.9, "maccs_max": 0.7},
            {"query_id": "q1", "target_class": "B", "ecfp4_max": 0.4, "maccs_max": 0.8},
        ]
    )
    chem3d = pd.DataFrame(
        [
            {
                "query_id": "q1",
                "target_class": "A",
                "usrcat_max": 0.5,
                "o3a_shape_tanimoto_max": 0.6,
                "o3a_color_max": 0.3,
                "pharmacophore_sim_max": 0.2,
            },
            {
                "query_id": "q1",
                "target_class": "B",
                "usrcat_max": 0.8,
                "o3a_shape_tanimoto_max": 0.7,
                "o3a_color_max": 0.9,
                "pharmacophore_sim_max": 0.6,
            },
        ]
    )

    fused = fuse_evidence(chem2d, chem3d, load_config())

    assert set(chem2d.columns).issubset(fused.columns)
    assert set(chem3d.columns).issubset(fused.columns)
    assert fused["chemical_evidence_score_v3"].between(0, 1).all()
    assert not fused["chemical_evidence_score_v3_is_probability"].any()
    assert (fused["fusion_component_count"] == 6).all()
    assert (fused["fusion_missing_components"] == "").all()
    assert fused.loc[fused.target_class == "B", "chemical_evidence_score_v3"].iloc[0] > fused.loc[
        fused.target_class == "A", "chemical_evidence_score_v3"
    ].iloc[0]
    for component in load_config().value("fusion.components"):
        assert f"{component}_fusion_rank" in fused
        assert f"{component}_fusion_contribution" in fused


def test_rank_fusion_preserves_missingness_and_uses_fixed_denominator() -> None:
    evidence = pd.DataFrame(
        [
            {"query_id": "q", "target_class": "complete", "a": 1.0, "b": 1.0},
            {"query_id": "q", "target_class": "partial", "a": 1.0, "b": np.nan},
            {"query_id": "q", "target_class": "missing", "a": np.nan, "b": np.nan},
        ]
    )
    fused = reciprocal_rank_fusion(
        evidence, components=["a", "b"], reciprocal_rank_constant=60.0
    ).set_index("target_class")

    assert math.isnan(fused.loc["partial", "b_fusion_contribution"])
    assert fused.loc["partial", "fusion_component_count"] == 1
    assert fused.loc["partial", "fusion_missing_components"] == "b"
    assert fused.loc["partial", "chemical_evidence_score_v3"] < fused.loc[
        "complete", "chemical_evidence_score_v3"
    ]
    assert math.isnan(fused.loc["missing", "chemical_evidence_score_v3"])


def test_rank_fusion_is_independent_between_queries() -> None:
    evidence = pd.DataFrame(
        [
            {"query_id": "q1", "target_class": "A", "x": 0.9},
            {"query_id": "q1", "target_class": "B", "x": 0.1},
            {"query_id": "q2", "target_class": "A", "x": 0.2},
            {"query_id": "q2", "target_class": "B", "x": 0.1},
        ]
    )
    fused = reciprocal_rank_fusion(
        evidence, components=["x"], reciprocal_rank_constant=60.0
    )
    top = fused[fused.target_class == "A"]
    assert (top["x_fusion_rank"] == 1.0).all()
    assert top["chemical_evidence_score_v3"].nunique() == 1
