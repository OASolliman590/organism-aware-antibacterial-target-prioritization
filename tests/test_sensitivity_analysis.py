from __future__ import annotations

from copy import deepcopy

import pandas as pd

from pipeline.config import ProjectConfig, load_config
from pipeline.sensitivity_analysis import (
    final_ranking_scenarios,
    rank_biased_overlap,
    recompute_final_priority,
    run_final_ranking_sensitivity,
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
                            "pharmacophore_2d_gobbi_similarity": base - 0.05,
                            "pharmacophore_3d_similarity": base - 0.06,
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


def _organism_predictions(config: ProjectConfig) -> pd.DataFrame:
    rows = []
    for organism_index, organism in enumerate(["org-a", "org-b"]):
        for query_index, query in enumerate(["q1", "q2"]):
            for target_index, target in enumerate(["A", "B", "C"]):
                rows.append(
                    {
                        "organism": organism,
                        "query_id": query,
                        "target_class": target,
                        "chemical_evidence_score_v3": 0.82 - 0.12 * target_index,
                        "target_specificity_margin": 0.05 + 0.08 * target_index,
                        "reference_quality_grade": ["A", "B", "C"][target_index],
                        "species_transfer_score": 0.45
                        + 0.20 * ((target_index + query_index) % 3),
                        "pocket_evidence_score": 0.30
                        + 0.25 * ((target_index + organism_index) % 3),
                        "anti_target_risk_score": 0.10 + 0.20 * target_index,
                        "organism_scope_score": 0.95 - 0.20 * target_index,
                        "clinical_priority_score": 0.35 + 0.25 * target_index,
                        "essentiality_score": 0.85 - 0.15 * target_index,
                        "cellular_access_score": 0.40 + 0.15 * target_index,
                        "resistance_relevance_score": 0.25 + 0.25 * target_index,
                        "card_resistance_context_score": float(target_index == 2),
                    }
                )
    frame = pd.DataFrame(rows)
    baseline, _ = final_ranking_scenarios(config)
    recomputed = recompute_final_priority(frame, baseline)
    frame["chemical_hypothesis_score"] = recomputed[
        "sensitivity_chemical_hypothesis_score"
    ]
    frame["overall_priority_score"] = recomputed[
        "sensitivity_overall_priority_score"
    ]
    return frame


def test_final_ranking_sensitivity_covers_all_published_ranking_coefficients() -> None:
    config = _config()
    predictions = _organism_predictions(config)
    first = run_final_ranking_sensitivity(predictions, config)
    second = run_final_ranking_sensitivity(predictions, config)

    pd.testing.assert_frame_equal(first, second, check_exact=True)
    scenarios = set(first.scenario)
    assert "leave_out_overall_layer.transfer" in scenarios
    assert "leave_out_overall_layer.pocket" in scenarios
    assert "leave_out_overall_layer.biology" in scenarios
    assert "leave_out_overall_layer.anti_target" in scenarios
    assert "leave_out_biology_component.organism_scope" in scenarios
    assert "overall.transfer.base_x0.5" in scenarios
    assert "overall.transfer.weight_x1.5" in scenarios
    assert "overall.anti_target_penalty_x1.25" in scenarios
    assert "specificity.margin_weight_x1.5" in scenarios
    assert "reference_quality.A_x0.5" in scenarios
    assert "leave_out_final_chemical_layer.specificity" in scenarios
    assert "leave_out_final_chemical_layer.reference_quality" in scenarios
    assert first[["kendall_tau_ci_lower_95", "rbo_ci_lower_95"]].notna().all().all()
    assert (first.n_ranked_lists_evaluable == 4).all()
    assert (first.bootstrap_unit == "organism_query_ranked_list").all()
    assert not first.score_is_probability.any()


def test_final_ranking_sensitivity_rejects_inconsistent_recorded_baseline() -> None:
    config = _config()
    predictions = _organism_predictions(config)
    predictions.loc[0, "overall_priority_score"] += 0.01

    try:
        run_final_ranking_sensitivity(predictions, config)
    except ValueError as error:
        assert "does not reproduce recorded baseline" in str(error)
    else:
        raise AssertionError("Expected an inconsistent-baseline error")


def test_final_ranking_sensitivity_excludes_missing_rows_without_imputation() -> None:
    config = _config()
    predictions = _organism_predictions(config)
    predictions.loc[0, "species_transfer_score"] = float("nan")
    predictions.loc[0, "overall_priority_score"] = float("nan")

    report = run_final_ranking_sensitivity(predictions, config)

    assert (report.n_prediction_rows_input == len(predictions)).all()
    assert (report.n_prediction_rows_evaluable == len(predictions) - 1).all()
    assert (report.n_prediction_rows_excluded_missing == 1).all()
