from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from rdkit import Chem

from pipeline import open_target_discovery_v2 as v2


FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _load_fixture():
    query_rows = json.loads((FIXTURES / "queries.json").read_text(encoding="utf-8"))
    reference_rows = json.loads(
        (FIXTURES / "references.json").read_text(encoding="utf-8")
    )
    queries = []
    for row in query_rows:
        molecule = v2.mol(row["canonical_smiles"])
        assert molecule is not None
        queries.append(
            {
                "query_id": row["drug"],
                "source": "tiny_public_fixture",
                "mol": molecule,
                "fp": v2.fp(molecule),
                "maccs": v2.maccs(molecule),
                **row,
            }
        )
    references = {}
    for target_class, rows in reference_rows.items():
        references[target_class] = []
        for row in rows:
            molecule = Chem.MolFromSmiles(row["canonical_smiles"])
            assert molecule is not None
            references[target_class].append(
                {
                    **row,
                    "_mol": molecule,
                    "_fp": v2.fp(molecule),
                    "_maccs": v2.maccs(molecule),
                    "_smi": Chem.MolToSmiles(molecule),
                }
            )
    return queries, references


def test_tiny_public_fixture_scores_deterministically_without_quality_imputation() -> None:
    queries, references = _load_fixture()
    missing_quality = pd.DataFrame(columns=["target_class"])
    empty = pd.DataFrame()

    first = [
        v2.score_query(query, references, missing_quality, empty, empty)
        for query in queries
    ]
    second = [
        v2.score_query(query, references, missing_quality, empty, empty)
        for query in queries
    ]

    for left, right in zip(first, second):
        pd.testing.assert_frame_equal(left, right, check_exact=True)
        assert set(left["target_class"]) == {"GyrB", "DHFR", "FabI"}
        assert left["chemical_evidence_score"].between(0, 1).all()
        assert (left["reference_quality_score"] == 0.0).all()
        assert (left["chemical_quality_adjusted_score"] == 0.0).all()
        assert set(left["reference_quality_grade"]) == {"insufficient"}
