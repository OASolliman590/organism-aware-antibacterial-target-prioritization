"""Freeze the tracked, public-only v2 outputs before additive v3 work.

The golden fingerprints deliberately exclude ignored private query structures and
compound-specific results. They are regenerated in memory from the public benchmark
and reference data using the unmodified v2 scoring path.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from pipeline import open_target_discovery_v2 as v2


ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "golden" / "v2_public_outputs.json"


def _sha256_csv(frame: pd.DataFrame) -> str:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _public_v2_outputs() -> dict[str, pd.DataFrame]:
    refs = v2.load_refs()
    quality = pd.read_csv(v2.QUALITY) if v2.QUALITY.exists() else pd.DataFrame()
    ontology = v2.load_ontology()
    compatibility = pd.read_csv(v2.COMPAT) if v2.COMPAT.exists() else pd.DataFrame()

    benchmark = pd.read_csv(v2.BENCH)
    scores = []
    for _, row in benchmark.iterrows():
        molecule = v2.mol(row.canonical_smiles)
        assert molecule is not None, f"Invalid tracked benchmark SMILES for {row.drug}"
        query = {
            "query_id": row.drug,
            "query_name": row.drug,
            "mol": molecule,
            "fp": v2.fp(molecule),
            "maccs": v2.maccs(molecule),
            "source": "eskape_benchmark",
            **row.to_dict(),
        }
        score = v2.score_query(
            query,
            refs,
            quality,
            compatibility,
            ontology,
            exclude_close=True,
        )
        if not score.empty:
            scores.append(score)

    benchmark_scores = pd.concat(scores, ignore_index=True)
    coverage = pd.DataFrame(
        {
            "target_class": sorted(refs),
            "n_reference_ligands": [len(refs[key]) for key in sorted(refs)],
        }
    )
    return {
        "v2_benchmark_open_target_scores.csv": benchmark_scores,
        "v2_open_target_reference_coverage.csv": coverage,
    }


def test_public_v2_outputs_match_pre_v3_golden() -> None:
    expected = json.loads(GOLDEN.read_text(encoding="utf-8"))["outputs"]
    actual = _public_v2_outputs()

    assert set(actual) == set(expected)
    for name, frame in actual.items():
        assert len(frame) == expected[name]["rows"], name
        assert _sha256_csv(frame) == expected[name]["sha256_lf"], name
