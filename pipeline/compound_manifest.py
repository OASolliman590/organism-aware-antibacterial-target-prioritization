"""Read the private compound manifest and resolve its organism vocabulary.

The manifest records which organism each novel-series compound was prepared
against. That assignment is design intent, not evidence: it is used to focus and
annotate figures, never to weight a score. Mapping the manifest's own microbe
labels onto the configured organism names is done through an explicit alias
table in ``config.yaml`` rather than by fuzzy matching, so that an abbreviated
or misspelled group name resolves visibly and reviewably or not at all.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

COMPOUND_COLUMN = "compound_code"
GROUP_COLUMN = "microbe_group"


class ManifestError(ValueError):
    """Raised when a manifest cannot be resolved against the configuration."""


@dataclass(frozen=True)
class CompoundManifest:
    """Compound-to-organism assignments plus the rows that did not resolve."""

    assignments: dict[str, str]
    unresolved_groups: tuple[str, ...] = ()

    def __bool__(self) -> bool:
        return bool(self.assignments)

    def compounds_for(self, organism: str) -> list[str]:
        """Compounds assigned to ``organism``, sorted for stable figure order."""

        return sorted(
            compound
            for compound, assigned in self.assignments.items()
            if assigned == organism
        )

    def organisms(self) -> list[str]:
        return sorted(set(self.assignments.values()))

    def is_assigned(self, compound: str, organism: str | None) -> bool:
        """Whether ``compound`` was prepared against ``organism``."""

        if organism is None:
            return False
        return self.assignments.get(compound) == organism


def resolve_group(
    group: str, *, aliases: dict[str, str], organism_names: list[str]
) -> str | None:
    """Map one manifest microbe group onto a configured organism name."""

    text = str(group).strip()
    if text in organism_names:
        return text
    if text in aliases:
        resolved = aliases[text]
        if resolved not in organism_names:
            raise ManifestError(
                f"organisms.manifest_aliases maps {text!r} to {resolved!r}, "
                "which is not in organisms.names"
            )
        return resolved
    return None


def load_manifest(
    path: Path, *, aliases: dict[str, str], organism_names: list[str]
) -> CompoundManifest:
    """Load ``path`` and resolve its groups; unknown groups are reported, not guessed."""

    frame = pd.read_csv(path)
    missing = {COMPOUND_COLUMN, GROUP_COLUMN} - set(frame.columns)
    if missing:
        raise ManifestError(
            f"manifest {path} is missing required column(s): {sorted(missing)}"
        )

    assignments: dict[str, str] = {}
    unresolved: list[str] = []
    for _, row in frame.iterrows():
        compound = str(row[COMPOUND_COLUMN]).strip()
        group = row[GROUP_COLUMN]
        if not compound or pd.isna(group):
            continue
        resolved = resolve_group(group, aliases=aliases, organism_names=organism_names)
        if resolved is None:
            unresolved.append(str(group).strip())
            continue
        assignments[compound] = resolved
    return CompoundManifest(
        assignments=assignments, unresolved_groups=tuple(sorted(set(unresolved)))
    )


def load_from_config(config) -> CompoundManifest:
    """Load the manifest declared in the run configuration, if it is present.

    A missing manifest is not an error: it is private input that need not exist
    for a public run, and callers degrade to unfocused figures.
    """

    try:
        path = config.path_for("compound_manifest")
    except Exception:
        return CompoundManifest(assignments={})
    if not path.is_file():
        return CompoundManifest(assignments={})
    try:
        aliases = dict(config.value("organisms.manifest_aliases"))
    except Exception:
        aliases = {}
    return load_manifest(
        path,
        aliases=aliases,
        organism_names=list(config.value("organisms.names")),
    )
