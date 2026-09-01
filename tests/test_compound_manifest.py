from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline.compound_manifest import (
    CompoundManifest,
    ManifestError,
    load_manifest,
    resolve_group,
)

ORGANISMS = [
    "Klebsiella pneumoniae",
    "Bacillus cereus",
    "Escherichia coli",
    "Staphylococcus aureus",
]
ALIASES = {
    "K_Pneumonia": "Klebsiella pneumoniae",
    "Bacillus_cereus": "Bacillus cereus",
    "E_Coli": "Escherichia coli",
    "MRSA": "Staphylococcus aureus",
}


def _write_manifest(path: Path, rows: list[tuple[str, str]]) -> Path:
    pd.DataFrame(rows, columns=["compound_code", "microbe_group"]).to_csv(
        path, index=False
    )
    return path


def test_aliases_resolve_abbreviated_and_misspelled_groups() -> None:
    assert (
        resolve_group("MRSA", aliases=ALIASES, organism_names=ORGANISMS)
        == "Staphylococcus aureus"
    )
    # An already-canonical name needs no alias.
    assert (
        resolve_group(
            "Escherichia coli", aliases=ALIASES, organism_names=ORGANISMS
        )
        == "Escherichia coli"
    )


def test_unknown_group_is_reported_not_guessed() -> None:
    # "Acinitobacter" without an alias entry must not be fuzzy-matched onto a
    # configured organism: a wrong genus would mislabel every figure.
    assert (
        resolve_group("Acinitobacter", aliases=ALIASES, organism_names=ORGANISMS)
        is None
    )


def test_alias_pointing_outside_configured_organisms_is_an_error() -> None:
    with pytest.raises(ManifestError, match="not in organisms.names"):
        resolve_group(
            "X", aliases={"X": "Yersinia pestis"}, organism_names=ORGANISMS
        )


def test_load_manifest_collects_assignments_and_unresolved_groups(
    tmp_path: Path,
) -> None:
    path = _write_manifest(
        tmp_path / "manifest.csv",
        [
            ("BI-1", "Bacillus_cereus"),
            ("BI-6", "Bacillus_cereus"),
            ("OX-11", "K_Pneumonia"),
            ("ZZ-9", "Unlisted_group"),
        ],
    )

    manifest = load_manifest(path, aliases=ALIASES, organism_names=ORGANISMS)

    assert manifest.assignments == {
        "BI-1": "Bacillus cereus",
        "BI-6": "Bacillus cereus",
        "OX-11": "Klebsiella pneumoniae",
    }
    assert manifest.unresolved_groups == ("Unlisted_group",)
    assert "ZZ-9" not in manifest.assignments


def test_manifest_missing_columns_is_an_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    pd.DataFrame([{"compound_code": "BI-1"}]).to_csv(path, index=False)

    with pytest.raises(ManifestError, match="missing required column"):
        load_manifest(path, aliases=ALIASES, organism_names=ORGANISMS)


def test_compounds_for_and_is_assigned() -> None:
    manifest = CompoundManifest(
        assignments={
            "BI-6": "Bacillus cereus",
            "BI-1": "Bacillus cereus",
            "OX-11": "Klebsiella pneumoniae",
        }
    )

    assert manifest.compounds_for("Bacillus cereus") == ["BI-1", "BI-6"]
    assert manifest.organisms() == ["Bacillus cereus", "Klebsiella pneumoniae"]
    assert manifest.is_assigned("BI-1", "Bacillus cereus")
    assert not manifest.is_assigned("BI-1", "Escherichia coli")
    # No organism in scope means nothing is marked as assigned.
    assert not manifest.is_assigned("BI-1", None)


def test_repository_manifest_aliases_cover_the_shipped_config() -> None:
    from pipeline.config import load_config

    config = load_config()
    organisms = list(config.value("organisms.names"))
    for group, organism in dict(config.value("organisms.manifest_aliases")).items():
        assert resolve_group(
            group, aliases={group: organism}, organism_names=organisms
        ) == organism
