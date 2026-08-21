"""Applicability-domain annotations without score rescaling or imputation."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

try:
    from pipeline.config import ProjectConfig
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig


def assign_applicability_domain(
    predictions: pd.DataFrame, config: ProjectConfig
) -> pd.DataFrame:
    """Add Tanimoto-threshold AD flags and continuous USRCAT NN distance.

    The specification declares Tanimoto thresholds but no calibrated USRCAT
    threshold. Consequently USRCAT remains an independent continuous AD field and
    is never forced through an invented cutoff.
    """

    required = {"ecfp4_max", "usrcat_max"}
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise ValueError(f"Predictions are missing AD evidence: {missing}")
    annotated = predictions.copy()
    tanimoto = pd.to_numeric(annotated["ecfp4_max"], errors="coerce").where(
        lambda values: np.isfinite(values)
    )
    usrcat = pd.to_numeric(annotated["usrcat_max"], errors="coerce").where(
        lambda values: np.isfinite(values)
    )
    in_threshold = float(config.value("applicability_domain.tanimoto_in"))
    out_threshold = float(config.value("applicability_domain.tanimoto_out"))
    annotated["ad_nearest_reference_tanimoto"] = tanimoto
    annotated["ad_nearest_reference_usrcat_similarity"] = usrcat
    annotated["ad_nearest_reference_usrcat_distance"] = 1.0 - usrcat
    annotated["ad_tanimoto_flag"] = np.select(
        [
            tanimoto.isna(),
            tanimoto >= in_threshold,
            tanimoto < out_threshold,
        ],
        ["unassessable_missing_tanimoto", "in_domain", "out_of_domain"],
        default="near_domain",
    )
    annotated["ad_usrcat_status"] = np.where(
        usrcat.notna(),
        "continuous_similarity_available_no_calibrated_cutoff",
        "unassessable_missing_usrcat",
    )
    annotated["applicability_domain_flag"] = annotated["ad_tanimoto_flag"]
    annotated["applicability_domain_flag_basis"] = (
        "configured ECFP4 Tanimoto thresholds; USRCAT retained continuously"
    )
    annotated["ad_tanimoto_in_threshold"] = in_threshold
    annotated["ad_tanimoto_out_threshold"] = out_threshold
    annotated["ad_shortlist_eligible"] = annotated["ad_tanimoto_flag"].isin(
        ["in_domain", "near_domain"]
    )
    annotated["ad_shortlist_discount_policy"] = (
        "ordering_only_no_score_rescaling"
    )
    return annotated


def shortlist_with_applicability_domain(
    predictions: pd.DataFrame,
    *,
    group_columns: Sequence[str],
    score_column: str,
    top_n: int,
) -> pd.DataFrame:
    """Rank AD-eligible rows first while preserving the original numeric score."""

    required = {*group_columns, "ad_shortlist_eligible", score_column}
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise ValueError(f"AD shortlist is missing fields: {missing}")
    ordered = predictions.sort_values(
        [*group_columns, "ad_shortlist_eligible", score_column],
        ascending=[*[True] * len(group_columns), False, False],
        kind="mergesort",
    )
    return ordered.groupby(list(group_columns), sort=False).head(top_n)
