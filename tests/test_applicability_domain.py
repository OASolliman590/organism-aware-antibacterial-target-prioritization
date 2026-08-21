from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline.applicability_domain import (
    assign_applicability_domain,
    shortlist_with_applicability_domain,
)
from pipeline.config import load_config


def test_ad_flags_use_only_declared_thresholds_and_keep_usrcat_continuous() -> None:
    source = pd.DataFrame(
        [
            {"target_class": "in", "ecfp4_max": 0.5, "usrcat_max": 0.8, "score": 0.1},
            {"target_class": "near", "ecfp4_max": 0.3, "usrcat_max": 0.7, "score": 0.2},
            {"target_class": "out", "ecfp4_max": 0.2, "usrcat_max": 0.6, "score": 0.9},
            {"target_class": "missing", "ecfp4_max": np.nan, "usrcat_max": np.nan, "score": 0.8},
        ]
    )

    annotated = assign_applicability_domain(source, load_config()).set_index(
        "target_class"
    )

    assert annotated.loc["in", "applicability_domain_flag"] == "in_domain"
    assert annotated.loc["near", "applicability_domain_flag"] == "near_domain"
    assert annotated.loc["out", "applicability_domain_flag"] == "out_of_domain"
    assert annotated.loc["missing", "applicability_domain_flag"] == (
        "unassessable_missing_tanimoto"
    )
    assert np.isclose(
        annotated.loc["in", "ad_nearest_reference_usrcat_distance"], 0.2
    )
    assert "no_calibrated_cutoff" in annotated.loc["in", "ad_usrcat_status"]
    assert annotated["score"].equals(source.set_index("target_class")["score"])


def test_out_of_domain_is_demoted_without_changing_score() -> None:
    frame = pd.DataFrame(
        [
            {"query_id": "q", "target_class": "out", "ad_shortlist_eligible": False, "score": 0.9},
            {"query_id": "q", "target_class": "near", "ad_shortlist_eligible": True, "score": 0.2},
        ]
    )
    shortlist = shortlist_with_applicability_domain(
        frame, group_columns=["query_id"], score_column="score", top_n=1
    )

    assert shortlist.iloc[0].target_class == "near"
    assert frame.loc[frame.target_class == "out", "score"].iloc[0] == 0.9
