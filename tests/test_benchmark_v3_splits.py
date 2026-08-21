from __future__ import annotations

from copy import deepcopy

import pandas as pd
from rdkit import Chem

from pipeline.benchmark_v3 import generate_splits
from pipeline.config import ProjectConfig, load_config


def _config() -> ProjectConfig:
    base = load_config()
    data = deepcopy(base.data)
    data["benchmark"]["analogue_exclusion_threshold"] = 0.85
    data["benchmark"]["time_cutoff"] = "2018-01-01"
    return ProjectConfig(path=base.path, root=base.root, data=data)


def _record(smiles: str, *, year=None):
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    row = {"canonical_smiles": smiles, "_mol": molecule}
    if year is not None:
        row["year"] = year
    return row


def test_all_splits_apply_real_exclusions_and_log_counts() -> None:
    query = _record("CC1=CC=CC=C1", year=2020)
    query.update({"query_id": "q", "target_class": "known_alias"})
    references = {
        "target_a": [
            _record("CC1=CC=CC=C1", year=2010),
            _record("CCC1=CC=CC=C1", year=2019),
        ],
        "target_b": [
            _record("CCOC(=O)N", year=2010),
            _record("CCS(=O)(=O)N"),
        ],
    }
    ontology = pd.DataFrame(
        [
            {
                "target_class": "target_a",
                "target_family": "family_a",
                "benchmark_aliases": "known_alias",
            },
            {
                "target_class": "target_b",
                "target_family": "family_b",
                "benchmark_aliases": "other_alias",
            },
        ]
    )

    results, provenance = generate_splits(
        [query], references, ontology, _config()
    )

    assert {result.split_type for result in results} == {
        "target_family",
        "scaffold",
        "temporal",
    }
    assert len(provenance) == 3
    assert provenance["analogue_leakage_guard_passed"].all()
    assert (
        provenance["max_remaining_query_reference_tanimoto"].dropna() < 0.85
    ).all()
    target = provenance.set_index("split_type")
    assert target.loc["target_family", "n_removed_target_family"] == 2
    assert target.loc["scaffold", "n_removed_same_scaffold"] == 2
    assert target.loc["temporal", "n_removed_post_cutoff"] == 1
    assert target.loc["temporal", "n_removed_missing_date"] == 1


def test_temporal_split_is_pending_without_query_date_not_imputed() -> None:
    query = _record("CCOC(=O)N")
    query.update({"query_id": "q", "target_class": "target_a"})
    references = {"target_a": [_record("CCS(=O)(=O)N", year=2010)]}
    ontology = pd.DataFrame(
        [{"target_class": "target_a", "target_family": "family_a"}]
    )

    results, provenance = generate_splits(
        [query], references, ontology, _config()
    )
    temporal = next(result for result in results if result.split_type == "temporal")
    row = provenance[provenance.split_type == "temporal"].iloc[0]

    assert temporal.references == {}
    assert row.status == "pending_missing_query_date"
    assert pd.isna(row.query_date_value)
