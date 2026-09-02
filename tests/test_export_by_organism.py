from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline.compound_manifest import CompoundManifest
from pipeline.export_by_organism import PREDICTIONS, export_by_organism

MANIFEST = CompoundManifest(
    assignments={"A-1": "Escherichia coli", "A-2": "Staphylococcus aureus"}
)


def _write_predictions(directory: Path) -> pd.DataFrame:
    rows = []
    for compound, base in (("A-1", 0.30), ("A-2", 0.20)):
        for organism in ("Escherichia coli", "Staphylococcus aureus"):
            for index, target in enumerate(("GyrB", "MurC")):
                rows.append(
                    {
                        "query_id": compound,
                        "organism": organism,
                        "target_class": target,
                        "overall_priority_score": base - 0.05 * index,
                        "confidence_class": "Moderate" if index == 0 else "Low",
                        "species_transfer_score": 0.8,
                        "pocket_evidence_score": 0.6,
                        "biological_priority_score": 0.5,
                        "clinical_translation_score": 0.4,
                        "sequence_mapping_status": "mapped",
                        "best_reference_molecule": "CHEMBL1",
                    }
                )
    frame = pd.DataFrame(rows)
    frame.to_csv(directory / PREDICTIONS, index=False)
    return frame


def test_one_directory_per_organism_with_both_views(tmp_path: Path) -> None:
    _write_predictions(tmp_path)

    index = export_by_organism(tmp_path, manifest=MANIFEST)

    assert set(index["organism"]) == {"Escherichia coli", "Staphylococcus aureus"}
    for _, row in index.iterrows():
        assert Path(str(row["targets_predicted_csv"])).is_file()
        assert Path(str(row["target_summary_csv"])).is_file()
        assert row["n_rows"] == 4
        assert row["n_compounds"] == 2
        assert row["n_target_classes"] == 2
    assert (tmp_path / "by_organism" / "index.csv").is_file()


def test_rows_are_partitioned_by_organism_and_ranked(tmp_path: Path) -> None:
    _write_predictions(tmp_path)

    export_by_organism(tmp_path, manifest=MANIFEST)

    exported = pd.read_csv(
        tmp_path / "by_organism" / "escherichia_coli" / "targets_predicted.csv"
    )
    assert set(exported["organism"]) == {"Escherichia coli"}
    scores = exported["overall_priority_score"].tolist()
    assert scores == sorted(scores, reverse=True)


def test_manifest_assignment_is_flagged_per_organism(tmp_path: Path) -> None:
    _write_predictions(tmp_path)

    export_by_organism(tmp_path, manifest=MANIFEST)

    coli = pd.read_csv(
        tmp_path / "by_organism" / "escherichia_coli" / "targets_predicted.csv"
    )
    aureus = pd.read_csv(
        tmp_path / "by_organism" / "staphylococcus_aureus" / "targets_predicted.csv"
    )
    # A-1 was prepared against E. coli only, so it is marked there and nowhere else.
    assert set(coli.loc[coli["manifest_assigned"], "query_id"]) == {"A-1"}
    assert set(aureus.loc[aureus["manifest_assigned"], "query_id"]) == {"A-2"}


def test_export_is_a_projection_not_a_recomputation(tmp_path: Path) -> None:
    source = _write_predictions(tmp_path)

    export_by_organism(tmp_path, manifest=MANIFEST)

    rebuilt = pd.concat(
        [
            pd.read_csv(path)
            for path in sorted(
                (tmp_path / "by_organism").glob("*/targets_predicted.csv")
            )
        ],
        ignore_index=True,
    )
    assert len(rebuilt) == len(source)
    key = ["query_id", "organism", "target_class"]
    merged = source.merge(rebuilt, on=key, suffixes=("_source", "_export"))
    assert len(merged) == len(source)
    pd.testing.assert_series_equal(
        merged["overall_priority_score_source"],
        merged["overall_priority_score_export"],
        check_names=False,
    )


def test_summary_names_the_best_compound_per_target(tmp_path: Path) -> None:
    _write_predictions(tmp_path)

    export_by_organism(tmp_path, manifest=MANIFEST)

    summary = pd.read_csv(
        tmp_path / "by_organism" / "escherichia_coli" / "target_summary.csv"
    ).set_index("target_class")
    # A-1 outscores A-2 on every target in the fixture.
    assert summary.loc["GyrB", "best_compound"] == "A-1"
    assert summary.loc["GyrB", "max_overall_priority_score"] == pytest.approx(0.30)
    assert summary.loc["GyrB", "n_compound_hypotheses"] == 2
    assert bool(summary.loc["GyrB", "best_compound_is_manifest_assigned"]) is True
    # Ranked by the strongest hypothesis, so GyrB precedes MurC.
    ordered = pd.read_csv(
        tmp_path / "by_organism" / "escherichia_coli" / "target_summary.csv"
    )["target_class"].tolist()
    assert ordered == ["GyrB", "MurC"]


def test_missing_predictions_table_is_reported_clearly(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="run the pipeline"):
        export_by_organism(tmp_path, manifest=MANIFEST)


def test_export_without_a_manifest_marks_nothing_assigned(tmp_path: Path) -> None:
    _write_predictions(tmp_path)

    export_by_organism(tmp_path, manifest=None)

    exported = pd.read_csv(
        tmp_path / "by_organism" / "escherichia_coli" / "targets_predicted.csv"
    )
    assert not exported["manifest_assigned"].any()
