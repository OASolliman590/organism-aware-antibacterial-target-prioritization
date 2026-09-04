from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline.compound_manifest import CompoundManifest
from pipeline.figures_suite import (
    STATUS_CREATED,
    STATUS_DEPENDENCY,
    STATUS_NO_ROWS,
    STATUS_NOT_APPLICABLE,
    STATUS_NOT_TARGET_SCOPED,
    STATUS_SOURCE_MISSING,
    PANELS,
    SuiteContext,
    generate_figure_suite,
    generate_per_organism_suites,
    organism_slug,
)

MANIFEST = CompoundManifest(
    assignments={"A-1": "Escherichia coli", "A-2": "Staphylococcus aureus"}
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

    status = generate_figure_suite(
        tmp_path, context=SuiteContext(manifest=MANIFEST)
    ).set_index("figure")

    table_backed = [panel.name for panel in PANELS if panel.sources]
    for name in table_backed:
        assert status.loc[name, "status"] == STATUS_CREATED, name
        assert Path(str(status.loc[name, "output"])).is_file()


def test_concordance_needs_a_manifest_and_says_so_without_one(tmp_path: Path) -> None:
    _write_full_run(tmp_path)

    status = generate_figure_suite(tmp_path).set_index("figure")

    # Design intent comes only from the manifest; with none supplied the panel
    # must report that rather than inferring an assignment from the scores.
    assert status.loc["assignment_concordance", "status"] == STATUS_SOURCE_MISSING


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


def test_per_organism_suites_land_in_their_own_directories(tmp_path: Path) -> None:
    _write_full_run(tmp_path)
    organisms = ["Escherichia coli", "Staphylococcus aureus"]

    status = generate_per_organism_suites(
        tmp_path,
        context=SuiteContext(manifest=MANIFEST),
        organisms=organisms,
    )

    assert set(status["organism"]) == set(organisms)
    for organism in organisms:
        directory = tmp_path / "figures_suite" / "by_organism" / organism_slug(organism)
        assert directory.is_dir()
        assert list(directory.glob("*.png"))


def test_cross_organism_panels_are_skipped_per_organism(tmp_path: Path) -> None:
    _write_full_run(tmp_path)
    cross_only = {panel.name for panel in PANELS if panel.cross_organism_only}
    assert cross_only, "expected at least one cross-organism panel"

    status = generate_per_organism_suites(
        tmp_path,
        context=SuiteContext(manifest=MANIFEST),
        organisms=["Escherichia coli"],
    ).set_index("figure")

    for name in cross_only:
        # Filtered to one organism these panels have nothing to compare against,
        # so they must be recorded as skipped rather than drawn from the
        # surviving rows.
        assert status.loc[name, "status"] == STATUS_NOT_APPLICABLE
        assert pd.isna(status.loc[name, "output"])
        assert not (
            tmp_path
            / "figures_suite"
            / "by_organism"
            / organism_slug("Escherichia coli")
            / f"{name}.png"
        ).exists()


def test_per_organism_suite_only_uses_that_organisms_rows(tmp_path: Path) -> None:
    _write_full_run(tmp_path)

    status = generate_per_organism_suites(
        tmp_path,
        context=SuiteContext(manifest=MANIFEST),
        organisms=["Klebsiella pneumoniae"],
    ).set_index("figure")

    # No row in the fixture carries this organism, so organism-aware panels have
    # nothing to draw and must say so.
    assert status.loc["compound_target_priority", "status"] == STATUS_NO_ROWS


def test_per_organism_suite_covers_only_that_organisms_compounds(
    tmp_path: Path,
) -> None:
    _write_full_run(tmp_path)

    generate_per_organism_suites(
        tmp_path,
        context=SuiteContext(manifest=MANIFEST),
        organisms=["Escherichia coli"],
    )

    # A-1 is the only compound prepared against E. coli, so the E. coli suite is
    # about A-1; A-2's rows for that organism are out of scope.
    scoped = SuiteContext(
        manifest=MANIFEST, organism="Escherichia coli", restrict_to_assigned=True
    )
    assert scoped.scoped_compounds() == ["A-1"]


def test_scoped_suites_drop_the_now_redundant_assignment_marker() -> None:
    scoped = SuiteContext(
        manifest=MANIFEST, organism="Escherichia coli", restrict_to_assigned=True
    )

    # Every compound in a scoped suite is assigned, so marking them all says
    # nothing; the footnote carries the scope instead.
    assert scoped.label("A-1") == "A-1"
    assert "Limited to the compounds prepared against this organism" in (
        scoped.assignment_note()
    )
    assert "A-1" in scoped.assignment_note()


def test_unscoped_suites_still_mark_assigned_compounds() -> None:
    context = SuiteContext(manifest=MANIFEST, organism="Escherichia coli")

    assert context.scoped_compounds() is None
    assert context.label("A-1") == "A-1*"
    assert context.label("A-2") == "A-2"
    assert context.labels(["A-1", "A-2"]) == ["A-1*", "A-2"]
    assert "A-1" in context.assignment_note()


def test_scoping_needs_both_a_manifest_and_an_organism() -> None:
    assert SuiteContext(restrict_to_assigned=True).scoped_compounds() is None
    assert (
        SuiteContext(
            restrict_to_assigned=True, organism="Escherichia coli"
        ).scoped_compounds()
        is None
    )
    assert (
        SuiteContext(manifest=MANIFEST, organism="Escherichia coli").scoped_compounds()
        is None
    )


def test_labels_are_unmarked_without_a_manifest_or_organism() -> None:
    assert SuiteContext().label("A-1") == "A-1"
    assert SuiteContext(manifest=MANIFEST).label("A-1") == "A-1"
    assert SuiteContext(manifest=MANIFEST).assignment_note() == ""


def test_organism_slug_is_filesystem_safe() -> None:
    assert organism_slug("Klebsiella pneumoniae") == "klebsiella_pneumoniae"


def test_panel_names_and_outputs_are_unique() -> None:
    names = [panel.name for panel in PANELS]
    outputs = [panel.output for panel in PANELS]

    assert len(names) == len(set(names))
    assert len(outputs) == len(set(outputs))


@pytest.mark.parametrize("panel", PANELS, ids=lambda panel: panel.name)
def test_every_panel_documents_itself(panel) -> None:
    assert panel.description
    assert panel.output.endswith(".png")


def test_target_scoped_suite_covers_only_those_target_classes(tmp_path: Path) -> None:
    _write_full_run(tmp_path)

    status = generate_figure_suite(
        tmp_path,
        context=SuiteContext(
            manifest=MANIFEST,
            organism="Escherichia coli",
            restrict_to_assigned=True,
            restrict_to_targets=("GyrB",),
        ),
        figure_dirname="scoped",
    ).set_index("figure")

    assert status.loc["compound_target_priority", "status"] == STATUS_CREATED

    # A target class the run never scored leaves the organism-aware panels with
    # nothing to draw, rather than silently widening back to every target.
    empty = generate_figure_suite(
        tmp_path,
        context=SuiteContext(
            manifest=MANIFEST,
            organism="Escherichia coli",
            restrict_to_assigned=True,
            restrict_to_targets=("NotATarget",),
        ),
        figure_dirname="scoped_empty",
    ).set_index("figure")
    assert empty.loc["compound_target_priority", "status"] == STATUS_NO_ROWS


def test_run_level_panels_are_skipped_when_scoped_to_targets(tmp_path: Path) -> None:
    _write_full_run(tmp_path)
    run_level = {panel.name for panel in PANELS if panel.run_level}
    assert run_level, "expected at least one run-level panel"

    status = generate_figure_suite(
        tmp_path,
        context=SuiteContext(organism="Escherichia coli", restrict_to_targets=("GyrB",)),
        figure_dirname="scoped2",
    ).set_index("figure")

    for name in run_level:
        # A target filter does not change these, so shipping them in a
        # target-scoped folder would imply they said something about it.
        assert status.loc[name, "status"] == STATUS_NOT_TARGET_SCOPED
        assert pd.isna(status.loc[name, "output"])


def test_run_level_panels_still_render_when_not_target_scoped(tmp_path: Path) -> None:
    _write_full_run(tmp_path)

    status = generate_figure_suite(
        tmp_path, context=SuiteContext(manifest=MANIFEST), figure_dirname="unscoped"
    ).set_index("figure")

    for panel in PANELS:
        if panel.run_level:
            assert status.loc[panel.name, "status"] == STATUS_CREATED


def test_scope_label_names_organism_and_target_count() -> None:
    context = SuiteContext(
        organism="Klebsiella pneumoniae", restrict_to_targets=("FabI", "LpxH")
    )

    assert "Klebsiella pneumoniae" in context.scope_label()
    assert "2 tested targets" in context.scope_label()
