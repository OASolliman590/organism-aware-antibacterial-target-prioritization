from __future__ import annotations

from copy import deepcopy

import pandas as pd

from pipeline.config import ProjectConfig, load_config
from pipeline.sensitivity_analysis import (
    rank_biased_overlap,
    run_sensitivity_analysis,
)


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["sensitivity"]["bootstrap_n"] = 5
    data["sensitivity"]["top_k"] = 3
    return ProjectConfig(path=base.path, root=base.root, data=data)


def _scores(config: ProjectConfig) -> pd.DataFrame:
    rows = []
    for split in ["target_family", "scaffold"]:
        for query in ["q1", "q2"]:
            for index, target in enumerate(["A", "B", "C"]):
                row = {
                    "split_type": split,
                    "query_id": query,
                    "target_class": target,
                }
                for component_index, component in enumerate(
                    config.value("fusion.components")
                ):
                    row[component] = 1.0 - 0.1 * index - 0.01 * component_index
                rows.append(row)
    return pd.DataFrame(rows)


def _references(config: ProjectConfig) -> pd.DataFrame:
    rows = []
    for split in ["target_family", "scaffold"]:
        for query in ["q1", "q2"]:
            for target_index, target in enumerate(["A", "B", "C"]):
                for reference_index in range(3):
                    base = 1.0 - 0.1 * target_index - 0.01 * reference_index
                    rows.append(
                        {
                            "split_type": split,
                            "query_id": query,
                            "target_class": target,
                            "reference_id": f"{target}-{reference_index}",
                            "ecfp4_similarity": base,
                            "maccs_similarity": base - 0.01,
                            "usrcat_similarity": base - 0.02,
                            "o3a_shape_tanimoto": base - 0.03,
                            "o3a_color": base - 0.04,
                            "pharmacophore_similarity": base - 0.05,
                        }
                    )
    return pd.DataFrame(rows)


def test_weight_and_layer_perturbations_have_reference_bootstrap_cis() -> None:
    config = _config()
    first = run_sensitivity_analysis(_scores(config), _references(config), config)
    second = run_sensitivity_analysis(_scores(config), _references(config), config)

    pd.testing.assert_frame_equal(first, second, check_exact=True)
    assert "leave_out_2d_layer" in set(first.scenario)
    assert "leave_out_shape_layer" in set(first.scenario)
    assert "leave_out_pharmacophore_layer" in set(first.scenario)
    available = first[first.split_type != "temporal"]
    temporal = first[first.split_type == "temporal"]
    assert (available.bootstrap_status == "available").all()
    assert available[["kendall_tau_ci_lower_95", "rbo_ci_lower_95"]].notna().all().all()
    assert (temporal.bootstrap_status == "unavailable_missing_reference_evidence").all()
    assert not first.score_is_probability.any()


def test_rbo_rewards_identical_rankings() -> None:
    ranking = ["A", "B", "C"]
    assert rank_biased_overlap(ranking, ranking, 0.9) == 1.0
    assert rank_biased_overlap(ranking, list(reversed(ranking)), 0.9) < 1.0
