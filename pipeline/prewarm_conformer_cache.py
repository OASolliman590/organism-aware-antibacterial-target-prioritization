"""Deterministically prewarm the ETKDGv3 cache with fixed worker scheduling.

Workers only change execution scheduling. Each molecule is still embedded with the
single-threaded, seeded parameters recorded in its conformer-cache manifest.
No identifiers or structures are written to the aggregate run-status output.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

import pandas as pd
from rdkit import Chem

try:
    from pipeline.chem3d_matching import generate_conformer_ensemble
    from pipeline.config import ProjectConfig, load_config
except ModuleNotFoundError:  # direct ``python pipeline/<script>.py`` execution
    from chem3d_matching import generate_conformer_ensemble
    from config import ProjectConfig, load_config


def canonical_smiles(molecule: Chem.Mol) -> str:
    """Return the exact identity used to deduplicate cache-prewarm work."""

    return Chem.MolToSmiles(
        Chem.RemoveHs(Chem.Mol(molecule)), canonical=True, isomericSmiles=True
    )


def prewarm_molecules(
    molecules: Iterable[Chem.Mol],
    config: ProjectConfig,
    *,
    cache_dir: str | Path | None = None,
    workers: int | None = None,
) -> pd.DataFrame:
    """Generate one cache record per unique structure and report honest statuses."""

    unique: dict[str, Chem.Mol] = {}
    invalid_count = 0
    input_count = 0
    for molecule in molecules:
        input_count += 1
        if molecule is None:
            invalid_count += 1
            continue
        unique.setdefault(canonical_smiles(molecule), Chem.Mol(molecule))
    ordered = [unique[key] for key in sorted(unique)]
    worker_count = int(
        config.value("chem3d.prewarm_workers") if workers is None else workers
    )
    if worker_count < 1:
        raise ValueError("workers must be at least one")

    def generate(molecule: Chem.Mol) -> dict[str, object]:
        try:
            ensemble = generate_conformer_ensemble(
                molecule, config, cache_dir=cache_dir
            )
            return {
                "status": ensemble.status,
                "cache_hit": bool(ensemble.cache_hit),
                "n_conformers": int(ensemble.n_conformers),
            }
        except Exception as exc:  # preserve the gap; never synthesize an ensemble
            return {
                "status": f"exception:{type(exc).__name__}",
                "cache_hit": False,
                "n_conformers": 0,
            }

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        outcomes = list(executor.map(generate, ordered))
    rows = [
        {
            "prewarm_workers": worker_count,
            "n_input_records": input_count,
            "n_unique_valid_structures": len(ordered),
            "n_invalid_structures": invalid_count,
            "n_cache_hits": sum(bool(row["cache_hit"]) for row in outcomes),
            "n_cache_misses": sum(not bool(row["cache_hit"]) for row in outcomes),
            "n_with_conformers": sum(int(row["n_conformers"]) > 0 for row in outcomes),
            "n_without_conformers": sum(int(row["n_conformers"]) == 0 for row in outcomes),
            "status_counts": ";".join(
                f"{status}:{sum(row['status'] == status for row in outcomes)}"
                for status in sorted({str(row["status"]) for row in outcomes})
            ),
            "scheduling_only": True,
            "per_molecule_num_threads": int(config.value("chem3d.num_threads")),
        }
    ]
    return pd.DataFrame(rows)


def load_run_molecules(config: ProjectConfig) -> list[Chem.Mol]:
    """Load exactly the reference and query structures used by the configured run."""

    try:
        from pipeline import open_target_discovery_v2 as discovery
    except ModuleNotFoundError:
        import open_target_discovery_v2 as discovery

    molecules = [
        record["_mol"]
        for records in discovery.load_refs().values()
        for record in records
        if record.get("_mol") is not None
    ]
    molecules.extend(
        query["mol"]
        for query in discovery.load_queries(config.path_for("private_compounds"))
    )
    benchmark_queries = discovery.load_queries(config.path_for("benchmark_structures"))
    if benchmark_queries:
        molecules.extend(query["mol"] for query in benchmark_queries)
    elif config.path_for("benchmark").is_file():
        benchmark = pd.read_csv(config.path_for("benchmark"))
        for smiles in benchmark.get("canonical_smiles", pd.Series(dtype=str)):
            molecule = discovery.mol(smiles)
            if molecule is not None:
                molecules.append(molecule)
    return molecules


def main() -> None:
    config = load_config()
    status = prewarm_molecules(load_run_molecules(config), config)
    output = config.path_for("results") / "conformer_cache_prewarm_status_v3.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    status.to_csv(output, index=False)
    print(status.to_string(index=False))


if __name__ == "__main__":
    main()
