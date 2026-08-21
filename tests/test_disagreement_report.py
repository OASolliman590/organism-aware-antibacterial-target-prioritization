from __future__ import annotations

import pandas as pd

from pipeline.disagreement_report import build_disagreement_report


COMPONENTS = ["ecfp4_max", "maccs_max", "usrcat_max"]


def _scores() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "query_id": "q",
                "target_class": "A",
                "chemical_evidence_score": 0.9,
                "chemical_evidence_score_v3": 0.1,
                "ecfp4_max": 0.9,
                "maccs_max": 0.8,
                "usrcat_max": 0.2,
                "fusion_component_count": 3,
                "fusion_missing_components": "",
            },
            {
                "query_id": "q",
                "target_class": "B",
                "chemical_evidence_score": 0.5,
                "chemical_evidence_score_v3": 0.5,
                "ecfp4_max": 0.5,
                "maccs_max": 0.5,
                "usrcat_max": 0.5,
                "fusion_component_count": 3,
                "fusion_missing_components": "",
            },
            {
                "query_id": "q",
                "target_class": "C",
                "chemical_evidence_score": 0.1,
                "chemical_evidence_score_v3": 0.9,
                "ecfp4_max": 0.1,
                "maccs_max": 0.2,
                "usrcat_max": 0.9,
                "fusion_component_count": 3,
                "fusion_missing_components": "",
            },
        ]
    )


def test_report_lists_all_material_promotions_and_demotions_with_components() -> None:
    report = build_disagreement_report(
        _scores(), component_columns=COMPONENTS, minimum_absolute_rank_shift=2
    ).set_index("target_class")

    assert set(report.index) == {"A", "C"}
    assert report.loc["A", "rank_change_direction"] == "demoted_by_v3"
    assert report.loc["C", "rank_change_direction"] == "promoted_by_v3"
    assert report.loc["A", "rank_shift_2d_to_v3"] == -2
    assert report.loc["C", "rank_shift_2d_to_v3"] == 2
    assert set(COMPONENTS).issubset(report.columns)


def test_report_is_honestly_empty_when_no_shift_meets_threshold() -> None:
    report = build_disagreement_report(
        _scores(), component_columns=COMPONENTS, minimum_absolute_rank_shift=3
    )

    assert report.empty
    assert set(COMPONENTS).issubset(report.columns)
    assert "rank_shift_2d_to_v3" in report
