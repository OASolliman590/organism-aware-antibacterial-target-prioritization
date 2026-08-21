"""Audit cases where adding 3D evidence changes a target's within-query rank."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


IDENTITY_COLUMNS = ["query_id", "target_class"]
RANK_COLUMNS = [
    "chemical_evidence_rank_2d",
    "chemical_evidence_rank_v3",
    "rank_shift_2d_to_v3",
    "absolute_rank_shift",
    "rank_change_direction",
]


def build_disagreement_report(
    scores: pd.DataFrame,
    *,
    component_columns: Sequence[str],
    minimum_absolute_rank_shift: float,
) -> pd.DataFrame:
    """Return every target whose fused rank differs materially from its 2D rank.

    Positive ``rank_shift_2d_to_v3`` means the target was promoted (its rank
    number decreased) after adding 3D/pharmacophore evidence. Exact score ties
    receive their average rank. Rows missing either score cannot establish a
    disagreement and are excluded, while their missingness remains in the source
    score output.
    """

    required = {
        *IDENTITY_COLUMNS,
        "chemical_evidence_score",
        "chemical_evidence_score_v3",
        *component_columns,
    }
    missing = sorted(required - set(scores.columns))
    if missing:
        raise ValueError(f"Score frame is missing disagreement fields: {missing}")
    if minimum_absolute_rank_shift <= 0:
        raise ValueError("minimum_absolute_rank_shift must be positive")
    if scores.duplicated(IDENTITY_COLUMNS).any():
        raise ValueError("Score frame contains duplicate query/target rows")

    compared = scores.copy()
    score_2d = pd.to_numeric(
        compared["chemical_evidence_score"], errors="coerce"
    ).where(lambda values: np.isfinite(values))
    score_v3 = pd.to_numeric(
        compared["chemical_evidence_score_v3"], errors="coerce"
    ).where(lambda values: np.isfinite(values))
    compared["chemical_evidence_rank_2d"] = score_2d.groupby(
        compared["query_id"], sort=False
    ).rank(method="average", ascending=False, na_option="keep")
    compared["chemical_evidence_rank_v3"] = score_v3.groupby(
        compared["query_id"], sort=False
    ).rank(method="average", ascending=False, na_option="keep")
    compared["rank_shift_2d_to_v3"] = (
        compared["chemical_evidence_rank_2d"]
        - compared["chemical_evidence_rank_v3"]
    )
    compared["absolute_rank_shift"] = compared["rank_shift_2d_to_v3"].abs()
    compared["rank_change_direction"] = np.select(
        [
            compared["rank_shift_2d_to_v3"] > 0,
            compared["rank_shift_2d_to_v3"] < 0,
        ],
        ["promoted_by_v3", "demoted_by_v3"],
        default="unchanged",
    )
    compared["material_rank_shift_threshold"] = float(
        minimum_absolute_rank_shift
    )
    compared["disagreement_definition"] = (
        "absolute within-query rank shift: v2 chemical evidence vs v3 rank fusion"
    )
    report = compared[
        compared["absolute_rank_shift"] >= float(minimum_absolute_rank_shift)
    ].copy()

    preferred = [
        *IDENTITY_COLUMNS,
        *(["dataset_scope"] if "dataset_scope" in report else []),
        "chemical_evidence_score",
        "chemical_evidence_score_v3",
        *RANK_COLUMNS,
        *component_columns,
        "fusion_component_count",
        "fusion_missing_components",
        "material_rank_shift_threshold",
        "disagreement_definition",
    ]
    ordered = list(dict.fromkeys(column for column in preferred if column in report))
    report = report[ordered]
    return report.sort_values(
        ["query_id", "absolute_rank_shift", "target_class"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
