from __future__ import annotations

from pathlib import Path

import pandas as pd

from pipeline.baseline_external import (
    parse_pidginv4_predictions,
    probe_pidginv4,
    run_or_load_pidginv4,
)
from pipeline.config import load_config


def test_missing_pidgin_runtime_models_and_map_are_pending_not_fabricated() -> None:
    config = load_config()
    status = probe_pidginv4(config)
    queries = pd.DataFrame(
        [{"query_id": "q", "canonical_smiles": "CCOC(=O)N"}]
    )
    scores, run_status = run_or_load_pidginv4(queries, config)

    assert status["status"] == "pending_unavailable_prerequisites"
    assert not status["python2_runtime_available"]
    assert not status["model_pickles_available"]
    assert not status["target_mapping_available"]
    assert scores.empty
    assert run_status["score_is_calibrated_probability"] is False


def test_parser_requires_explicit_versioned_target_mapping(tmp_path: Path) -> None:
    predictions = tmp_path / "predictions.tsv"
    predictions.write_text(
        "target_id\tq1\tq2\nP0ABQ4\t0.8\t0.2\nUNMAPPED\t0.9\t0.9\n",
        encoding="utf-8",
    )
    mapping = tmp_path / "mapping.csv"
    pd.DataFrame(
        [
            {
                "baseline_target_id": "P0ABQ4",
                "target_class": "DHFR",
                "mapping_source": "reviewed mapping fixture",
                "mapping_version": "fixture-v1",
            }
        ]
    ).to_csv(mapping, index=False)

    parsed = parse_pidginv4_predictions(
        predictions,
        mapping,
        code_commit="abc",
        model_release="model-v1",
    )

    assert len(parsed) == 2
    assert set(parsed.target_class) == {"DHFR"}
    assert set(parsed.query_id) == {"q1", "q2"}
    assert not parsed.baseline_score_is_calibrated_probability.any()
    assert "UNMAPPED" not in set(parsed.baseline_target_id)
