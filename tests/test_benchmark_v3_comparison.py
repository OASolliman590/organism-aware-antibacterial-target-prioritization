from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd

from pipeline.benchmark_v3 import add_3d_only_score, compare_score_modes
from pipeline.config import ProjectConfig, load_config


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["benchmark"]["bootstrap_n"] = 50
    return ProjectConfig(path=base.path, root=base.root, data=data)


def _scores() -> pd.DataFrame:
    rows = []
    for split in ["target_family", "scaffold"]:
        rows.extend(
            [
                {
                    "split_type": split,
                    "query_id": "q",
                    "target_class": "active",
                    "is_active": 1,
                    "chemical_evidence_score": 0.9,
                    "chemical_evidence_score_v3": 0.1,
                    "ecfp4_max": 0.9,
                    "maccs_max": 0.9,
                    "usrcat_max": 0.1,
                    "o3a_shape_tanimoto_max": 0.1,
                    "o3a_color_max": 0.1,
                    "pharmacophore_2d_gobbi_sim_max": 0.9,
                    "pharmacophore_3d_sim_max": 0.1,
                },
                {
                    "split_type": split,
                    "query_id": "q",
                    "target_class": "decoy",
                    "is_active": 0,
                    "chemical_evidence_score": 0.1,
                    "chemical_evidence_score_v3": 0.9,
                    "ecfp4_max": 0.1,
                    "maccs_max": 0.1,
                    "usrcat_max": 0.9,
                    "o3a_shape_tanimoto_max": 0.9,
                    "o3a_color_max": 0.9,
                    "pharmacophore_2d_gobbi_sim_max": 0.1,
                    "pharmacophore_3d_sim_max": 0.9,
                },
            ]
        )
    return pd.DataFrame(rows)


def _provenance() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "split_type": split,
                "status": "available" if split != "temporal" else "pending_missing_query_date",
                "n_references_input": 2,
                "n_references_after_split": 2 if split != "temporal" else 0,
                "n_removed_close_analogue": 0,
                "n_removed_same_scaffold": 0,
                "n_removed_target_family": 0,
                "n_removed_post_cutoff": 0,
                "n_removed_missing_date": 2 if split == "temporal" else 0,
            }
            for split in ["target_family", "scaffold", "temporal"]
        ]
    )


def test_single_table_compares_all_modes_and_reports_worse_3d_honestly() -> None:
    config = _config()
    scored = add_3d_only_score(_scores(), config)
    query_metrics, comparison = compare_score_modes(
        scored, _provenance(), config
    )

    assert scored["chemical_evidence_score_3d_only"].between(0, 1).all()
    assert len(comparison) == 3 * 3 * 7
    scaffold_auroc = comparison[
        (comparison.split_type == "scaffold") & (comparison.metric == "auroc")
    ].set_index("score_mode")
    assert scaffold_auroc.loc["2d_only", "estimate"] == 1.0
    assert scaffold_auroc.loc["3d_only", "estimate"] == 0.0
    assert scaffold_auroc.loc["3d_only", "performance_vs_2d"] == "worse"
    assert scaffold_auroc.loc["fusion", "performance_vs_2d"] == "worse"
    temporal = comparison[comparison.split_type == "temporal"]
    assert temporal.estimate.isna().all()
    assert (temporal.performance_vs_2d == "unavailable").all()
    assert set(query_metrics.score_mode) == {"2d_only", "3d_only", "fusion"}


def test_3d_only_score_excludes_alignment_free_2d_gobbi_component() -> None:
    config = _config()
    original = _scores()
    reversed_gobbi = original.copy()
    reversed_gobbi["pharmacophore_2d_gobbi_sim_max"] = 1.0 - reversed_gobbi[
        "pharmacophore_2d_gobbi_sim_max"
    ]

    original_score = add_3d_only_score(original, config)[
        "chemical_evidence_score_3d_only"
    ]
    reversed_score = add_3d_only_score(reversed_gobbi, config)[
        "chemical_evidence_score_3d_only"
    ]

    pd.testing.assert_series_equal(original_score, reversed_score)


def test_3d_only_score_uses_aligned_3d_pharmacophore_component() -> None:
    config = _config()
    forward = _scores()
    for component in config.value("benchmark.three_dimensional_components"):
        forward[component] = 0.5
    forward.loc[forward.target_class == "active", "pharmacophore_3d_sim_max"] = 0.9
    forward.loc[forward.target_class == "decoy", "pharmacophore_3d_sim_max"] = 0.1
    reversed_3d = forward.copy()
    reversed_3d["pharmacophore_3d_sim_max"] = 1.0 - reversed_3d[
        "pharmacophore_3d_sim_max"
    ]

    forward_scored = add_3d_only_score(forward, config)
    reversed_scored = add_3d_only_score(reversed_3d, config)

    assert (
        forward_scored.loc[
            forward_scored.target_class == "active",
            "chemical_evidence_score_3d_only",
        ].iloc[0]
        > forward_scored.loc[
            forward_scored.target_class == "decoy",
            "chemical_evidence_score_3d_only",
        ].iloc[0]
    )
    assert (
        reversed_scored.loc[
            reversed_scored.target_class == "active",
            "chemical_evidence_score_3d_only",
        ].iloc[0]
        < reversed_scored.loc[
            reversed_scored.target_class == "decoy",
            "chemical_evidence_score_3d_only",
        ].iloc[0]
    )
