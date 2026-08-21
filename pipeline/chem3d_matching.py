"""Deterministic 3D conformer generation and matching evidence.

This module never substitutes synthetic scores for failed embeddings or unavailable
force-field parameters. Such cases carry an explicit status and remain missing.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

try:
    from pipeline.config import ProjectConfig
except ModuleNotFoundError:  # direct module execution/import compatibility
    from config import ProjectConfig


@dataclass(frozen=True)
class ConformerEnsemble:
    molecule: Chem.Mol | None
    relative_energies_kcal: tuple[float | None, ...]
    optimization_statuses: tuple[int | None, ...]
    cache_key: str
    cache_path: Path | None
    cache_hit: bool
    status: str
    n_embedded: int

    @property
    def n_conformers(self) -> int:
        return 0 if self.molecule is None else self.molecule.GetNumConformers()


def _canonical_smiles(molecule: Chem.Mol) -> str:
    no_hydrogens = Chem.RemoveHs(Chem.Mol(molecule))
    return Chem.MolToSmiles(no_hydrogens, canonical=True, isomericSmiles=True)


def _conformer_parameters(config: ProjectConfig) -> dict[str, Any]:
    chem3d = config.value("chem3d")
    return {
        "method": "ETKDGv3",
        "force_field": "MMFF94",
        "seed": int(config.value("seeds.conformer")),
        "n_confs": int(chem3d["n_confs"]),
        "prune_rms": float(chem3d["prune_rms"]),
        "energy_window_kcal": float(chem3d["energy_window_kcal"]),
        "max_iterations": int(chem3d["max_iterations"]),
        "num_threads": int(chem3d["num_threads"]),
        "rdkit_version": rdBase.rdkitVersion,
    }


def conformer_cache_key(molecule: Chem.Mol, config: ProjectConfig) -> str:
    payload = {
        "canonical_smiles": _canonical_smiles(molecule),
        "parameters": _conformer_parameters(config),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(canonical).hexdigest()


def _cache_root(config: ProjectConfig, cache_dir: str | Path | None) -> Path:
    selected = Path(cache_dir) if cache_dir is not None else Path(
        str(config.value("chem3d.cache_dir"))
    )
    return selected.resolve() if selected.is_absolute() else (config.root / selected).resolve()


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def _load_cached(
    binary_path: Path, manifest_path: Path, cache_key: str
) -> ConformerEnsemble | None:
    if not binary_path.is_file() or not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("cache_key") != cache_key:
            return None
        molecule = Chem.Mol(binary_path.read_bytes())
        energies = tuple(manifest["relative_energies_kcal"])
        statuses = tuple(manifest["optimization_statuses"])
        if molecule.GetNumConformers() != len(energies) or len(energies) != len(statuses):
            return None
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return None
    return ConformerEnsemble(
        molecule=molecule,
        relative_energies_kcal=energies,
        optimization_statuses=statuses,
        cache_key=cache_key,
        cache_path=binary_path,
        cache_hit=True,
        status=str(manifest["status"]),
        n_embedded=int(manifest["n_embedded"]),
    )


def _filter_by_energy(
    molecule: Chem.Mol,
    optimization: list[tuple[int, float]],
    energy_window_kcal: float,
) -> tuple[Chem.Mol, tuple[float, ...], tuple[int, ...]]:
    finite = [float(energy) for _, energy in optimization if math.isfinite(energy)]
    minimum = min(finite)
    keep_indices = [
        index
        for index, (_, energy) in enumerate(optimization)
        if math.isfinite(energy) and float(energy) <= minimum + energy_window_kcal
    ]
    filtered = Chem.Mol(molecule)
    filtered.RemoveAllConformers()
    relative: list[float] = []
    statuses: list[int] = []
    for index in keep_indices:
        filtered.AddConformer(molecule.GetConformer(index), assignId=True)
        status, energy = optimization[index]
        relative.append(float(energy) - minimum)
        statuses.append(int(status))
    return filtered, tuple(relative), tuple(statuses)


def generate_conformer_ensemble(
    molecule: Chem.Mol,
    config: ProjectConfig,
    *,
    cache_dir: str | Path | None = None,
) -> ConformerEnsemble:
    """Generate a seeded ETKDGv3/MMFF94 ensemble or load its exact binary cache."""

    if molecule is None:
        raise ValueError("molecule must not be None")
    parameters = _conformer_parameters(config)
    cache_key = conformer_cache_key(molecule, config)
    root = _cache_root(config, cache_dir)
    binary_path = root / f"{cache_key}.rdkit.bin"
    manifest_path = root / f"{cache_key}.json"
    cached = _load_cached(binary_path, manifest_path, cache_key)
    if cached is not None:
        return cached

    embedded = Chem.AddHs(Chem.Mol(molecule), addCoords=True)
    embed_parameters = AllChem.ETKDGv3()
    embed_parameters.randomSeed = parameters["seed"]
    embed_parameters.pruneRmsThresh = parameters["prune_rms"]
    embed_parameters.numThreads = parameters["num_threads"]
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(
            embedded,
            numConfs=parameters["n_confs"],
            params=embed_parameters,
        )
    )
    n_embedded = len(conformer_ids)
    if not conformer_ids:
        return ConformerEnsemble(
            molecule=None,
            relative_energies_kcal=(),
            optimization_statuses=(),
            cache_key=cache_key,
            cache_path=None,
            cache_hit=False,
            status="embedding_failed",
            n_embedded=0,
        )

    if not AllChem.MMFFHasAllMoleculeParams(embedded):
        ensemble = embedded
        relative_energies: tuple[float | None, ...] = tuple(
            None for _ in conformer_ids
        )
        optimization_statuses: tuple[int | None, ...] = tuple(
            None for _ in conformer_ids
        )
        status = "mmff94_parameters_unavailable_unfiltered"
    else:
        optimization = list(
            AllChem.MMFFOptimizeMoleculeConfs(
                embedded,
                numThreads=parameters["num_threads"],
                maxIters=parameters["max_iterations"],
                mmffVariant="MMFF94",
            )
        )
        if not optimization or not any(math.isfinite(x[1]) for x in optimization):
            return ConformerEnsemble(
                molecule=None,
                relative_energies_kcal=(),
                optimization_statuses=(),
                cache_key=cache_key,
                cache_path=None,
                cache_hit=False,
                status="mmff94_energy_failed",
                n_embedded=n_embedded,
            )
        ensemble, relative_energies, optimization_statuses = _filter_by_energy(
            embedded, optimization, parameters["energy_window_kcal"]
        )
        status = "ok"

    manifest = {
        "schema_version": 1,
        "cache_key": cache_key,
        "canonical_smiles": _canonical_smiles(molecule),
        "parameters": parameters,
        "status": status,
        "n_embedded": n_embedded,
        "n_kept": ensemble.GetNumConformers(),
        "relative_energies_kcal": list(relative_energies),
        "optimization_statuses": list(optimization_statuses),
    }
    _atomic_write(binary_path, ensemble.ToBinary())
    _atomic_write(
        manifest_path,
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    round_tripped = _load_cached(binary_path, manifest_path, cache_key)
    if round_tripped is None:
        raise RuntimeError(f"Conformer cache verification failed: {binary_path}")
    return ConformerEnsemble(
        molecule=round_tripped.molecule,
        relative_energies_kcal=round_tripped.relative_energies_kcal,
        optimization_statuses=round_tripped.optimization_statuses,
        cache_key=cache_key,
        cache_path=binary_path,
        cache_hit=False,
        status=round_tripped.status,
        n_embedded=n_embedded,
    )
