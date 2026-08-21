from __future__ import annotations

from copy import deepcopy

import numpy as np
import pandas as pd

from pipeline.benchmark_v3 import (
    aggregate_metrics_with_bootstrap,
    query_level_metrics,
)
from pipeline.config import ProjectConfig, load_config


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["benchmark"]["bootstrap_n"] = 200
    return ProjectConfig(path=base.path, root=base.root, data=data)


def _score_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"query_id": "q1", "target_class": "A", "score": 0.9, "is_active": 1},
            {"query_id": "q1", "target_class": "B", "score": 0.1, "is_active": 0},
            {"query_id": "q2", "target_class": "A", "score": 0.1, "is_active": 1},
            {"query_id": "q2", "target_class": "B", "score": 0.9, "is_active": 0},
            {"query_id": "q3", "target_class": "B", "score": 0.5, "is_active": 0},
            {"query_id": "q3", "target_class": "C", "score": 0.4, "is_active": 0},
        ]
    )


def test_query_metrics_cover_auroc_bedroc_ef_mrr_and_coverage() -> None:
    metrics = query_level_metrics(_score_rows(), score_column="score").set_index(
        "query_id"
    )

    assert metrics.loc["q1", "auroc"] == 1.0
    assert metrics.loc["q2", "auroc"] == 0.0
    assert np.isnan(metrics.loc["q3", "auroc"])
    assert metrics.loc["q1", "mrr"] == 1.0
    assert metrics.loc["q2", "mrr"] == 0.5
    assert metrics.loc["q3", "mrr"] == 0.0
    assert metrics["coverage"].mean() == 2 / 3
    assert {"bedroc_alpha_20_0", "bedroc_alpha_80_5", "ef_1pct", "ef_5pct"}.issubset(
        metrics.columns
    )


def test_every_aggregate_has_deterministic_ci_n_and_split_provenance() -> None:
    config = _config()
    query_metrics = query_level_metrics(_score_rows(), score_column="score")
    provenance = pd.DataFrame(
        [
            {
                "status": "available",
                "n_references_input": 10,
                "n_references_after_split": 8,
                "n_removed_close_analogue": 1,
                "n_removed_same_scaffold": 1,
                "n_removed_target_family": 0,
                "n_removed_post_cutoff": 0,
                "n_removed_missing_date": 0,
            }
        ]
    )

    first = aggregate_metrics_with_bootstrap(
        query_metrics,
        split_type="scaffold",
        score_mode="2d",
        split_provenance=provenance,
        config=config,
    )
    second = aggregate_metrics_with_bootstrap(
        query_metrics,
        split_type="scaffold",
        score_mode="2d",
        split_provenance=provenance,
        config=config,
    )

    pd.testing.assert_frame_equal(first, second, check_exact=True)
    assert set(first.metric) == {
        "auroc",
        "bedroc_alpha_20_0",
        "bedroc_alpha_80_5",
        "ef_1pct",
        "ef_5pct",
        "mrr",
        "coverage",
    }
    assert first[["ci_lower_95", "ci_upper_95"]].notna().all().all()
    assert (first.n_queries == 3).all()
    assert (first.bootstrap_unit == "query_id").all()
    assert (first.snapshot_id == config.value("snapshots.snapshot_id")).all()
    assert (first.n_removed_close_analogue == 1).all()
