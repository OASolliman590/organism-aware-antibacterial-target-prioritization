from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline.figures_suite import (
    STATUS_CREATED,
    STATUS_DEPENDENCY,
    STATUS_NO_ROWS,
    STATUS_SOURCE_MISSING,
    PANELS,
    SuiteContext,
    generate_figure_suite,
)


def _predictions_frame(compounds: tuple[str, ...] = ("A-1", "A-2")) -> pd.DataFrame:
    rows = []
    for compound in compounds:
        for organism in ("Escherichia coli", "Staphylococcus aureus"):
            for index, target in enumerate(("GyrB", "MurC")):
                rows.append(
                    {
                        "query_id": compound,
                        "organism": organism,
                        "target_class": target,
                        "chemical_quality_adjusted_score": 0.3 + 0.1 * index,
                        "species_transfer_score": 0.8,
                        "pocket_evidence_score": 0.6,
                        "biological_priority_score": 0.5,
                        "overall_priority_score": 0.25 + 0.05 * index,
                        "confidence_class": "Moderate" if index == 0 else "Low",
                    }
                )
    return pd.DataFrame(rows)


def _write_full_run(directory: Path) -> None:
    _predictions_frame().to_csv(
        directory / "v2_open_target_predictions_by_organism.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "query_id": "A-1",
                "target_class": "GyrB",
                "ecfp4_max_fusion_contribution": 0.016,
                "maccs_max_fusion_contribution": 0.012,
                "overall_priority_score": 0.3,
            },
            {
                "query_id": "A-2",
                "target_class": "MurC",
                "ecfp4_max_fusion_contribution": 0.011,
                "maccs_max_fusion_contribution": 0.009,
                "overall_priority_score": 0.2,
            },
        ]
    ).to_csv(directory / "open_target_predictions_by_organism_v3.csv", index=False)
    pd.DataFrame(
        [
            {
                "query_id": "A-1",
                "target_class": "GyrB",
                "rank_shift_2d_to_v3": 5.0,
                "absolute_rank_shift": 5.0,
            },
            {
                "query_id": "A-2",
                "target_class": "MurC",
                "rank_shift_2d_to_v3": -4.0,
                "absolute_rank_shift": 4.0,
            },
        ]
    ).to_csv(directory / "chemical_evidence_disagreements_v3.csv", index=False)
    pd.DataFrame(
        [
            {
                "scenario": "s1",
                "scenario_type": "final_biology_weight_perturbation",
                "mean_kendall_tau": 0.98,
                "kendall_tau_ci_lower_95": 0.96,
                "kendall_tau_ci_upper_95": 0.99,
                "mean_rbo": 0.99,
                "rbo_ci_lower_95": 0.97,
                "rbo_ci_upper_95": 1.0,
            }
        ]
    ).to_csv(directory / "final_ranking_sensitivity_v3.csv", index=False)
    pd.DataFrame(
        [
            {
                "query_id": "A-1",
                "target_class": "GyrB",
                "bootstrap_stability_score": 0.42,
                "empirical_decoy_p_value": 0.03,
                "bootstrap_top1_probability": 0.8,
            }
        ]
    ).to_csv(directory / "v2_uncertainty_private.csv", index=False)
    pd.DataFrame(
        [
            {
                "split": "scaffold",
                "mode": "ensemble_score",
                "top1_enrichment_over_random": 11.3,
                "top3_enrichment_over_random": 4.2,
            },
            {
                "split": "scaffold",
                "mode": "prevalence_baseline",
                "top1_enrichment_over_random": 3.4,
                "top3_enrichment_over_random": 1.9,
            },
        ]
    ).to_csv(directory / "benchmark_v2_summary.csv", index=False)


def test_every_panel_renders_from_a_complete_run(tmp_path: Path) -> None:
    _write_full_run(tmp_path)

    status = generate_figure_suite(tmp_path).set_index("figure")

    table_backed = [panel.name for panel in PANELS if panel.sources]
    for name in table_backed:
        assert status.loc[name, "status"] == STATUS_CREATED, name
        assert Path(str(status.loc[name, "output"])).is_file()


def test_missing_tables_are_reported_not_invented(tmp_path: Path) -> None:
    status = generate_figure_suite(tmp_path).set_index("figure")

    for panel in PANELS:
        if not panel.sources:
            continue
        assert status.loc[panel.name, "status"] == STATUS_SOURCE_MISSING
        assert pd.isna(status.loc[panel.name, "output"])
    assert not list((tmp_path / "figures_suite").glob("*.png"))


def test_present_but_unevaluable_tables_are_reported(tmp_path: Path) -> None:
    _write_full_run(tmp_path)
    # A table that exists with the right header but carries no scored rows must
    # not produce an empty panel.
    pd.DataFrame(columns=["query_id", "target_class", "rank_shift_2d_to_v3"]).to_csv(
        tmp_path / "chemical_evidence_disagreements_v3.csv", index=False
    )

    status = generate_figure_suite(tmp_path).set_index("figure")

    assert status.loc["rank_shift_disagreements", "status"] == STATUS_NO_ROWS
    assert pd.isna(status.loc["rank_shift_disagreements", "output"])
    assert not (tmp_path / "figures_suite" / "rank_shift_disagreements.png").exists()


def test_chemical_space_reports_missing_input_without_raising(tmp_path: Path) -> None:
    _write_full_run(tmp_path)

    status = generate_figure_suite(
        tmp_path, context=SuiteContext(private_compounds=tmp_path / "absent.sdf")
    ).set_index("figure")

    assert status.loc["chemical_space", "status"] in {
        STATUS_SOURCE_MISSING,
        STATUS_DEPENDENCY,
    }


def test_panel_names_and_outputs_are_unique() -> None:
    names = [panel.name for panel in PANELS]
    outputs = [panel.output for panel in PANELS]

    assert len(names) == len(set(names))
    assert len(outputs) == len(set(outputs))


@pytest.mark.parametrize("panel", PANELS, ids=lambda panel: panel.name)
def test_every_panel_documents_itself(panel) -> None:
    assert panel.description
    assert panel.output.endswith(".png")
