"""Fetch and curate ChEMBL reference ligands for v2.1 target subtypes.

The script uses a versioned alias configuration plus ChEMBL target search. It keeps
raw target/activity provenance, filters molecular activities by assay confidence and
valid numeric potency, standardizes structures conservatively with RDKit, and writes
one reference-ligand JSON per target subtype. It never reads private compounds.

Network use is bounded and reproducible: target IDs discovered by search are cached,
manual IDs remain in the configuration, and a local cache can be reused with
CHEMBL_V21_OFFLINE=1 after the first successful retrieval.
"""
from __future__ import annotations
from pathlib import Path
import json, math, os, re, time
from typing import Any
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(os.environ.get("PROJECT_ROOT", Path(__file__).resolve().parents[1]))
CONFIG = ROOT / "data" / "chembl_target_aliases_v21.json"
ONTO = ROOT / "data" / "target_subtype_ontology_v21.csv"
OUT = ROOT / "data" / "reference_ligands"
CACHE = ROOT / "data" / "reference_quality" / "chembl_v21_cache"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)
BASE = "https://www.ebi.ac.uk/chembl/api/data"
OFFLINE = os.environ.get("CHEMBL_V21_OFFLINE", "0") == "1"
DISCOVER = os.environ.get("CHEMBL_V21_DISCOVER", "1") == "1"
SUBTYPE_FILTER = {x.strip() for x in os.environ.get("CHEMBL_V21_SUBTYPES", "").split(",") if x.strip()}
MIN_PCHEMBL = float(os.environ.get("CHEMBL_V21_MIN_PCHEMBL", "5.0"))
MAX_PER_SUBTYPE = int(os.environ.get("CHEMBL_V21_MAX_PER_SUBTYPE", "4000"))
HTTP_TIMEOUT = int(os.environ.get("CHEMBL_V21_HTTP_TIMEOUT", "30"))

session = requests.Session()
retry = Retry(total=4, backoff_factor=0.7, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
session.mount("https://", HTTPAdapter(max_retries=retry))


def get_json(url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    if OFFLINE:
        raise RuntimeError("offline mode: network retrieval disabled")
    r = session.get(url, params=params, timeout=(10, HTTP_TIMEOUT))
    r.raise_for_status()
    time.sleep(0.12)
    return r.json()


def cached_json(path: Path, url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text())
    data = get_json(url, params)
    path.write_text(json.dumps(data, indent=2))
    return data


def pchembl(activity: dict[str, Any]) -> float | None:
    raw = activity.get("pchembl_value")
    try:
        if raw not in (None, ""):
            return float(raw)
    except (TypeError, ValueError):
        pass
    try:
        value = float(activity.get("standard_value"))
        if value <= 0:
            return None
        unit = str(activity.get("standard_units", "")).strip().lower().replace("μ", "µ")
        molar = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "µm": 1e-6, "nm": 1e-9}.get(unit)
        if molar is None:
            return None
        return -math.log10(value * molar)
    except (TypeError, ValueError):
        return None


def standardize(smiles: str) -> tuple[str, str] | None:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol)
        # Remove salts/solvents but preserve stereochemistry, charge, beta-lactams,
        # covalent warheads, and tautomer-sensitive ribosome pharmacophores.
        mol = rdMolStandardize.FragmentParent(mol)
        Chem.SanitizeMol(mol)
        can = Chem.MolToSmiles(mol, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(mol)
        return can, inchikey
    except Exception:
        return None


def target_search(term: str) -> list[dict[str, Any]]:
    safe = re.sub(r"[^A-Za-z0-9]+", "_", term.lower()).strip("_")
    path = CACHE / f"search_{safe}.json"
    try:
        data = cached_json(path, f"{BASE}/target/search.json", {"q": term, "limit": 1000})
    except Exception as exc:
        print("WARN target search", term, repr(exc))
        return []
    return data.get("targets", data.get("target", [])) or []


def target_metadata(tid: str) -> dict[str, Any]:
    path = CACHE / f"target_{tid}.json"
    try:
        return cached_json(path, f"{BASE}/target/{tid}.json")
    except Exception as exc:
        print("WARN target metadata", tid, repr(exc))
        return {"target_chembl_id": tid}


def discover_ids(cfg: dict[str, Any]) -> list[str]:
    ids = set(cfg.get("manual_target_ids", []))
    if not DISCOVER:
        return sorted(ids)
    for term in cfg.get("search_terms", []):
        for hit in target_search(term):
            tid = hit.get("target_chembl_id")
            name = f"{hit.get('pref_name', '')} {hit.get('organism', '')}".lower()
            # Do not import mammalian/viral searches into antibacterial subtype sets.
            if tid and not any(x in name for x in ["human", "homo sapiens", "mouse", "virus", "fungus", "yeast"]):
                ids.add(tid)
    return sorted(ids)


def subtype_accepts(subtype: str, meta: dict[str, Any]) -> bool:
    name = str(meta.get("pref_name", "")).lower()
    if subtype.startswith("PBP"):
        return "penicillin-binding protein" in name or "pbp" in name
    if subtype == "Beta-lactamase_class_A":
        return any(k in name for k in ["ctx-m", "tem", "shv", "kpc", "pse", "ges", "class a"])
    if subtype == "Beta-lactamase_class_B":
        return "metallo" in name or any(k in name for k in ["vim", "ndm", "imp", "ind", "gob", "l1"])
    if subtype == "Beta-lactamase_class_C":
        return any(k in name for k in ["class c", "ampc", "adc", "acc"])
    if subtype == "Beta-lactamase_class_D":
        return "class d" in name or "oxa" in name
    if subtype.startswith("30S"):
        return any(k in name for k in ["30s", "16s", "ribosomal", "ribosome", "ribosomal subunit"])
    if subtype.startswith("50S"):
        return any(k in name for k in ["50s", "23s", "ribosomal", "ribosome", "ribosomal subunit"])
    return True


def assay_quality(a: dict[str, Any]) -> tuple[str, bool]:
    conf = a.get("confidence_score")
    try: conf = int(conf)
    except (TypeError, ValueError): conf = 0
    at = str(a.get("assay_type", "")).upper()
    standard = str(a.get("standard_type", "")).upper()
    valid = str(a.get("data_validity_comment", "") or "").lower()
    numeric = a.get("standard_value") not in (None, "") or a.get("pchembl_value") not in (None, "")
    direct = at in {"B", "F"} or standard in {"KD", "KI", "IC50", "EC50", "XC50", "POTENCY"}
    if not numeric or valid not in {"", "none", "nan"}:
        return "exclude", False
    if conf >= 8 and direct: return "A", True
    if conf >= 6 and direct: return "B", True
    if conf >= 4: return "C", True
    return "D", False


def fetch_activities(tid: str) -> list[dict[str, Any]]:
    path = CACHE / f"activities_{tid}.json"
    if path.exists():
        return json.loads(path.read_text())
    all_rows: list[dict[str, Any]] = []
    offset = 0
    while True:
        data = get_json(f"{BASE}/activity.json", {"target_chembl_id": tid, "limit": 1000, "offset": offset})
        rows = data.get("activities", []) or []
        if not rows: break
        all_rows.extend(rows)
        meta = data.get("page_meta", {})
        total = int(meta.get("total_count", offset + len(rows)))
        offset += len(rows)
        if offset >= total or len(rows) < 1000: break
    path.write_text(json.dumps(all_rows))
    return all_rows


def curate(subtype: str, cfg: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ids = discover_ids(cfg)
    (CACHE / f"ids_{subtype}.json").write_text(json.dumps(ids, indent=2))
    rows: list[dict[str, Any]] = []
    for tid in ids:
        meta = target_metadata(tid)
        if not subtype_accepts(subtype, meta):
            continue
        try: activities = fetch_activities(tid)
        except Exception as exc:
            print("WARN activities", subtype, tid, repr(exc)); continue
        for a in activities:
            smi = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            pc = pchembl(a)
            grade, use = assay_quality(a)
            if not smi or not mid or pc is None or pc < MIN_PCHEMBL:
                continue
            norm = standardize(smi)
            if norm is None:
                continue
            can, inchikey = norm
            rows.append({
                "molecule_chembl_id": mid, "canonical_smiles": smi,
                "canonical_smiles_standardized": can, "standardized_inchikey": inchikey,
                "pchembl_value": round(pc, 3), "standard_type": a.get("standard_type"),
                "standard_value": a.get("standard_value"), "standard_units": a.get("standard_units"),
                "assay_type": a.get("assay_type"), "confidence_score": a.get("confidence_score"),
                "data_validity_comment": a.get("data_validity_comment"),
                "target_chembl_id": tid, "target_pref_name": meta.get("pref_name", ""),
                "target_organism": meta.get("organism", ""), "target_type": meta.get("target_type", ""),
                "document_chembl_id": a.get("document_chembl_id"),
                "assay_chembl_id": a.get("assay_chembl_id"), "reference_source": "ChEMBL",
                "source_release": os.environ.get("CHEMBL_RELEASE", "runtime_api"),
                "quality_grade": grade, "quality_included": use,
            })
    # Keep each molecule once per subtype, retaining the best potency but preserving
    # assay-count metadata needed to audit weak or contradictory classes.
    grouped: dict[str, list[dict[str, Any]]] = {}
    for r in rows: grouped.setdefault(r["standardized_inchikey"], []).append(r)
    out = []
    for key, group in grouped.items():
        group.sort(key=lambda r: r["pchembl_value"], reverse=True)
        best = dict(group[0])
        best["n_assays_for_standardized_structure"] = len(group)
        best["assay_quality_grades"] = ";".join(sorted({str(x["quality_grade"]) for x in group}))
        best["n_target_chembl_ids"] = len({x["target_chembl_id"] for x in group})
        out.append(best)
    out.sort(key=lambda r: r["pchembl_value"], reverse=True)
    return out[:MAX_PER_SUBTYPE], rows


def main() -> None:
    cfg = json.loads(CONFIG.read_text())
    manifest = []
    for subtype, spec in cfg.items():
        if SUBTYPE_FILTER and subtype not in SUBTYPE_FILTER:
            continue
        try:
            curated, raw = curate(subtype, spec)
        except Exception as exc:
            print("FAILED", subtype, repr(exc)); curated, raw = [], []
        (OUT / f"ref_ligands_{subtype}.json").write_text(json.dumps(curated, indent=2))
        manifest.append({"target_subtype": subtype, "n_curated_ligands": len(curated), "n_raw_qualifying_records": len(raw), "manual_target_ids": spec.get("manual_target_ids", []), "offline": OFFLINE, "min_pchembl": MIN_PCHEMBL})
        print(subtype, "curated", len(curated), "raw", len(raw))
    (CACHE / "manifest_v21.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__": main()
