from __future__ import annotations

from pathlib import Path

import pandas as pd

from pipeline.v3_figures import generate_v3_figures


def test_figures_are_created_only_for_evaluable_run_tables(tmp_path: Path) -> None:
    pd.DataFrame(
        [
            {
                "split_type": "scaffold",
                "score_mode": "2d_only",
                "metric": "mrr",
                "estimate": 0.7,
                "ci_lower_95": 0.5,
                "ci_upper_95": 0.9,
            },
            {
                "split_type": "scaffold",
                "score_mode": "fusion",
                "metric": "mrr",
                "estimate": 0.6,
                "ci_lower_95": 0.4,
                "ci_upper_95": 0.8,
            },
        ]
    ).to_csv(tmp_path / "benchmark_mode_comparison_v3.csv", index=False)
    pd.DataFrame(
        [
            {"applicability_domain_flag": "in_domain"},
            {"applicability_domain_flag": "out_of_domain"},
        ]
    ).to_csv(tmp_path / "benchmark_target_scores_by_split_v3.csv", index=False)
    pd.DataFrame(
        columns=["query_id", "target_class", "absolute_rank_shift"]
    ).to_csv(tmp_path / "chemical_evidence_disagreements_v3.csv", index=False)
    pd.DataFrame(
        columns=["n", "mean_calibrated_probability", "observed_fraction_correct"]
    ).to_csv(tmp_path / "scoring_model_reliability_v3.csv", index=False)

    status = generate_v3_figures(tmp_path).set_index("figure")

    assert status.loc["benchmark_mode_comparison", "status"] == "created"
    assert status.loc["applicability_domain", "status"] == "created"
    assert status.loc["rank_disagreement", "status"] == "unavailable_no_evaluable_rows"
    assert status.loc["calibration_reliability", "status"] == "unavailable_no_evaluable_rows"
    assert (tmp_path / "figures" / "benchmark_mode_comparison_v3.png").is_file()
    assert (tmp_path / "figures" / "applicability_domain_v3.png").is_file()
