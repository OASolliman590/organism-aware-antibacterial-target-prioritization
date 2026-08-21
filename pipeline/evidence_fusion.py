"""Transparent fusion of independent 2D, 3D, and pharmacophore evidence.

The fused value is a deterministic ranking score.  It is not a calibrated
probability.  Missing component measurements remain missing and are exposed by
coverage fields rather than filled with synthetic values.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

try:
    from pipeline.config import ProjectConfig
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig


KEY_COLUMNS = ["query_id", "target_class"]


def _validate_unique_keys(frame: pd.DataFrame, label: str) -> None:
    missing = [column for column in KEY_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing key columns: {missing}")
    duplicate = frame.duplicated(KEY_COLUMNS, keep=False)
    if duplicate.any():
        examples = frame.loc[duplicate, KEY_COLUMNS].head(3).to_dict("records")
        raise ValueError(f"{label} has duplicate query/target rows: {examples}")


def merge_evidence_frames(
    chem2d: pd.DataFrame, chem3d: pd.DataFrame
) -> pd.DataFrame:
    """Outer-join evidence without discarding unmatched rows from either layer."""

    _validate_unique_keys(chem2d, "chem2d")
    _validate_unique_keys(chem3d, "chem3d")
    overlapping = sorted(
        (set(chem2d.columns) & set(chem3d.columns)) - set(KEY_COLUMNS)
    )
    if overlapping:
        raise ValueError(f"Evidence columns occur in both frames: {overlapping}")
    return chem2d.merge(
        chem3d,
        on=KEY_COLUMNS,
        how="outer",
        validate="one_to_one",
        sort=True,
    )


def reciprocal_rank_fusion(
    evidence: pd.DataFrame,
    *,
    components: Sequence[str],
    reciprocal_rank_constant: float,
) -> pd.DataFrame:
    """Add normalized Reciprocal Rank Fusion fields per query.

    Higher component values rank first within each query. Exact ties receive
    their average rank. Each available component contributes ``1 / (k + rank)``.
    The sum is divided by the fixed ideal contribution from every configured
    component, so missing evidence is not silently rewarded. A row with no
    available component remains missing rather than becoming a synthetic zero.
    """

    selected = list(components)
    if not selected or len(selected) != len(set(selected)):
        raise ValueError("fusion components must be a non-empty unique list")
    absent = [column for column in selected if column not in evidence.columns]
    if absent:
        raise ValueError(f"Evidence frame is missing fusion components: {absent}")
    if reciprocal_rank_constant <= 0:
        raise ValueError("reciprocal_rank_constant must be positive")
    _validate_unique_keys(evidence, "evidence")

    fused = evidence.copy()
    contribution_columns: list[str] = []
    for component in selected:
        values = pd.to_numeric(fused[component], errors="coerce")
        finite = values.where(np.isfinite(values))
        rank_column = f"{component}_fusion_rank"
        contribution_column = f"{component}_fusion_contribution"
        fused[rank_column] = finite.groupby(fused["query_id"], sort=False).rank(
            method="average", ascending=False, na_option="keep"
        )
        fused[contribution_column] = 1.0 / (
            float(reciprocal_rank_constant) + fused[rank_column]
        )
        contribution_columns.append(contribution_column)

    available = fused[selected].apply(
        lambda column: pd.to_numeric(column, errors="coerce")
    )
    available = available.apply(lambda column: np.isfinite(column))
    fused["fusion_component_count"] = available.sum(axis=1).astype(int)
    fused["fusion_missing_components"] = available.apply(
        lambda row: ";".join(
            component for component in selected if not bool(row[component])
        ),
        axis=1,
    )
    raw = fused[contribution_columns].sum(axis=1, min_count=1)
    ideal = len(selected) / (float(reciprocal_rank_constant) + 1.0)
    fused["chemical_evidence_score_v3"] = (raw / ideal).where(
        fused["fusion_component_count"] > 0
    )
    fused["chemical_evidence_score_v3_is_probability"] = False
    fused["fusion_method"] = "reciprocal_rank_fusion"
    fused["fusion_reciprocal_rank_constant"] = float(reciprocal_rank_constant)
    fused["fusion_configured_component_count"] = len(selected)
    return fused


def fuse_evidence(
    chem2d: pd.DataFrame, chem3d: pd.DataFrame, config: ProjectConfig
) -> pd.DataFrame:
    """Merge evidence frames and apply the config-selected documented fusion."""

    mode = str(config.value("run.fusion_mode"))
    if mode != "rank_fusion":
        raise ValueError(
            f"Unsupported fusion mode {mode!r}; v3 currently specifies rank_fusion"
        )
    return reciprocal_rank_fusion(
        merge_evidence_frames(chem2d, chem3d),
        components=list(config.value("fusion.components")),
        reciprocal_rank_constant=float(
            config.value("fusion.reciprocal_rank_constant")
        ),
    )
