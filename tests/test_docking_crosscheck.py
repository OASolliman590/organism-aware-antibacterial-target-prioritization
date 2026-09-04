from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline.docking_crosscheck import (
    PREDICTIONS,
    DockingCrosscheckError,
    DockingSpec,
    build_agreement,
    build_coverage,
    build_crosscheck,
    collapse_to_target_class,
    run_crosscheck,
    undocked_top_targets,
)

ORGANISM = "Klebsiella pneumoniae"
SPEC = DockingSpec(
    organism=ORGANISM,
    target_aliases={"9L5X_FabI": "FabI", "2OV5_KPC": "Beta-lactamase_class_A",
                    "4ZBE_KPC": "Beta-lactamase_class_A"},
    excluded_targets={"BAD_5EIX": "receptor was DNA only"},
    control_ligand_prefixes=("native_",),
)
HEAVY = {"OX-11": 32, "T2Z14": 43, "native_TCL": 17, "native_NXL": 17}


def _predictions() -> pd.DataFrame:
    rows = []
    scores = {"FabI": 0.043, "Beta-lactamase_class_A": 0.0, "MurC": 0.094, "GyrB": 0.080}
    for compound in ("OX-11", "T2Z14"):
        for target, score in scores.items():
            rows.append(
                {
                    "query_id": compound,
                    "organism": ORGANISM,
                    "target_class": target,
                    "overall_priority_score": score,
                    "chemical_evidence_score": score / 2,
                    "confidence_class": "Low" if score else "Insufficient",
                    "species_transfer_score": 0.6,
                    "pocket_evidence_score": 0.6,
                }
            )
    return pd.DataFrame(rows)


def _docking() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"target_folder": "9L5X_FabI", "ligand": "OX-11",
             "docking_affinity_kcal_mol": -9.53},
            {"target_folder": "9L5X_FabI", "ligand": "T2Z14",
             "docking_affinity_kcal_mol": -9.86},
            {"target_folder": "9L5X_FabI", "ligand": "native_TCL",
             "docking_affinity_kcal_mol": -7.71},
            {"target_folder": "2OV5_KPC", "ligand": "OX-11",
             "docking_affinity_kcal_mol": -6.42},
            {"target_folder": "2OV5_KPC", "ligand": "T2Z14",
             "docking_affinity_kcal_mol": -6.04},
            {"target_folder": "4ZBE_KPC", "ligand": "OX-11",
             "docking_affinity_kcal_mol": -6.54},
            {"target_folder": "4ZBE_KPC", "ligand": "T2Z14",
             "docking_affinity_kcal_mol": -5.63},
            {"target_folder": "4ZBE_KPC", "ligand": "native_NXL",
             "docking_affinity_kcal_mol": -5.48},
            {"target_folder": "BAD_5EIX", "ligand": "OX-11",
             "docking_affinity_kcal_mol": -5.31},
            {"target_folder": "UNKNOWN_RECEPTOR", "ligand": "OX-11",
             "docking_affinity_kcal_mol": -8.0},
        ]
    )


def test_coverage_separates_mapped_excluded_and_unmapped() -> None:
    coverage = build_coverage(_docking(), SPEC, _predictions()).set_index(
        "docking_receptor"
    )

    assert coverage.loc["9L5X_FabI", "mapping_status"] == "mapped"
    assert coverage.loc["BAD_5EIX", "mapping_status"] == "excluded_by_config"
    # The whole point: a receptor nobody mapped must surface, not vanish.
    assert coverage.loc["UNKNOWN_RECEPTOR", "mapping_status"] == "unmapped_no_alias"


def test_unmapped_receptors_are_left_out_of_the_comparison() -> None:
    crosscheck = build_crosscheck(_predictions(), _docking(), SPEC, heavy_atoms=HEAVY)

    assert "UNKNOWN_RECEPTOR" not in set(crosscheck.target_folder)
    assert "BAD_5EIX" not in set(crosscheck.target_folder)


def test_native_controls_are_not_treated_as_query_compounds() -> None:
    crosscheck = build_crosscheck(_predictions(), _docking(), SPEC, heavy_atoms=HEAVY)

    assert set(crosscheck.query_id) == {"OX-11", "T2Z14"}


def test_ligand_efficiency_compares_against_the_native_control() -> None:
    crosscheck = build_crosscheck(_predictions(), _docking(), SPEC, heavy_atoms=HEAVY)

    match = crosscheck[
        (crosscheck.query_id == "OX-11") & (crosscheck.target_class == "FabI")
    ]
    assert len(match) == 1
    row = match.iloc[0]
    assert row.ligand_efficiency == pytest.approx(9.53 / 32, rel=1e-6)
    assert row.control_ligand_efficiency == pytest.approx(7.71 / 17, rel=1e-6)
    # Beats the native on raw affinity but loses once size is accounted for:
    # exactly the confound this column exists to expose.
    assert bool(row.beats_control_on_affinity) is True
    assert bool(row.beats_control_on_efficiency) is False


def test_repeated_receptors_for_one_target_class_are_collapsed() -> None:
    crosscheck = build_crosscheck(_predictions(), _docking(), SPEC, heavy_atoms=HEAVY)
    collapsed = collapse_to_target_class(crosscheck)

    kpc = crosscheck[crosscheck.target_class == "Beta-lactamase_class_A"]
    assert len(kpc) == 4, "both KPC structures should survive in the detail table"
    assert len(collapsed) == 4, "two compounds x two target classes"
    # The better (more negative) of the two KPC structures is the one kept.
    best = collapsed[
        (collapsed.query_id == "OX-11")
        & (collapsed.target_class == "Beta-lactamase_class_A")
    ]
    assert best.docking_affinity_kcal_mol.iloc[0] == pytest.approx(-6.54)


def test_agreement_is_computed_without_pseudo_replication() -> None:
    crosscheck = build_crosscheck(_predictions(), _docking(), SPEC, heavy_atoms=HEAVY)

    agreement = build_agreement(crosscheck).set_index("scope")

    # Two target classes per compound after collapsing, four rows pooled: too
    # few to correlate, and it must say so rather than invent a coefficient.
    assert agreement.loc["OX-11", "n_shared_targets"] == 2
    assert agreement.loc["OX-11", "status"] == "unavailable_too_few_comparable_targets"
    assert pd.isna(agreement.loc["OX-11", "spearman_rho"])


def test_undocked_top_targets_names_what_was_never_tested() -> None:
    crosscheck = build_crosscheck(_predictions(), _docking(), SPEC, heavy_atoms=HEAVY)

    undocked = undocked_top_targets(
        _predictions(), crosscheck, SPEC, top_n=4
    ).set_index("target_class")

    assert bool(undocked.loc["MurC", "docked"]) is False
    assert bool(undocked.loc["GyrB", "docked"]) is False
    assert bool(undocked.loc["FabI", "docked"]) is True
    assert undocked.loc["MurC", "pipeline_rank"] == 1


def test_missing_predictions_table_is_reported(tmp_path: Path) -> None:
    _docking().to_csv(tmp_path / "dock.csv", index=False)

    with pytest.raises(DockingCrosscheckError, match="run the pipeline"):
        run_crosscheck(tmp_path, tmp_path / "dock.csv", SPEC)


def test_nothing_mapping_is_an_error_not_a_silent_empty_result(tmp_path: Path) -> None:
    _predictions().to_csv(tmp_path / PREDICTIONS, index=False)
    _docking().to_csv(tmp_path / "dock.csv", index=False)
    blind = DockingSpec(organism=ORGANISM, target_aliases={}, excluded_targets={})

    # The prior campaign's own cross-check wrote an empty correlation file and a
    # "matched_rows: 0" report, which reads as "no disagreement found". Refuse.
    with pytest.raises(DockingCrosscheckError, match="no docked receptor mapped"):
        run_crosscheck(tmp_path, tmp_path / "dock.csv", blind)


def test_run_writes_all_three_tables(tmp_path: Path) -> None:
    _predictions().to_csv(tmp_path / PREDICTIONS, index=False)
    _docking().to_csv(tmp_path / "dock.csv", index=False)

    tables = run_crosscheck(
        tmp_path, tmp_path / "dock.csv", SPEC, heavy_atoms=HEAVY
    )

    assert set(tables) == {"coverage", "crosscheck", "agreement", "undocked_top_targets"}
    for name in ("docking_pipeline_coverage.csv", "docking_pipeline_crosscheck.csv",
                 "docking_pipeline_agreement.csv"):
        assert (tmp_path / name).is_file()


def test_unknown_organism_is_an_error(tmp_path: Path) -> None:
    _predictions().to_csv(tmp_path / PREDICTIONS, index=False)
    _docking().to_csv(tmp_path / "dock.csv", index=False)
    spec = DockingSpec(organism="Nonexistent species",
                       target_aliases=SPEC.target_aliases,
                       excluded_targets=SPEC.excluded_targets)

    with pytest.raises(DockingCrosscheckError, match="no prediction rows"):
        run_crosscheck(tmp_path, tmp_path / "dock.csv", spec)
