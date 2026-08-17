"""Target scoring engine: ligand-based similarity + pharmacophore/SAR evidence.

Reads:
  - user compounds (from SDF directory)
  - ChEMBL reference ligand sets per target class (data/ref_ligands_*.json)

Computes for each (compound, target_class):
  - max ECFP4 Tanimoto to any reference ligand
  - mean of top-5 ECFP4 Tanimoto
  - MACCS key Tanimoto (best)
  - top reference ligand identity + its pChEMBL (potency context)
  - pharmacophore flags per class

Writes scores.csv and supporting JSONs into results/.
"""
import json, glob, os, math
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, rdMolDescriptors
import numpy as np
import pandas as pd

from pathlib import Path
PROJECT_ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
WORK = str(PROJECT_ROOT)
DATA = str(PROJECT_ROOT / 'data')
RES = str(PROJECT_ROOT / 'results')
os.makedirs(RES, exist_ok=True)
REF_DATA = Path(DATA) / 'reference_ligands'
if not REF_DATA.exists():
    REF_DATA = Path(DATA)

SDF_DIR = str(PROJECT_ROOT / 'data' / 'compounds')

# ---------------------------------------------------------------- compound set ---
def load_compounds():
    rows = []
    for f in [f"{SDF_DIR}/compounds_normalized.sdf"]:
        base = os.path.basename(f)
        for m in Chem.SDMolSupplier(f, sanitize=True, removeHs=True):
            if m is None:
                continue
            name = m.GetProp('_Name') if m.HasProp('_Name') else os.path.splitext(base)[0]
            Chem.SanitizeMol(m)
            smiles = Chem.MolToSmiles(m)
            rows.append({
                "compound": name, "canonical_smiles": smiles, "source_file": base,
                "mol": m,
                "mw": Descriptors.MolWt(m),
                "logp": Descriptors.MolLogP(m),
                "tpsa": rdMolDescriptors.CalcTPSA(m),
                "hbd": rdMolDescriptors.CalcNumHBD(m),
                "hba": rdMolDescriptors.CalcNumHBA(m),
                "rotbonds": rdMolDescriptors.CalcNumRotatableBonds(m),
                "heavy_atoms": m.GetNumHeavyAtoms(),
            })
    df = pd.DataFrame(rows)
    df = df.sort_values("compound").reset_index(drop=True)
    return df


def get_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def tanimoto(fp1, fp2):
    return AllChem.DataStructs.TanimotoSimilarity(fp1, fp2)


# -------------------------------------------------- pharmacophore / SAR flags ---
# SMILES patterns (SMARTS) describing class-defining features.
PHARM = {
    "DHPS": {
        "sulfonamide_aniline": "[SX3](=O)(=O)[NX3;H1,H2][a]",
        "sulfonamide_heteroaryl": "[SX3](=O)(=O)[NX3;H1,H2]c1cscn1",
        "paba_mimic": "[NX3;H1,H2][cH0][cH0]1[c,n,s]",
    },
    "DHFR": {
        "pterin_like_pyrimidine_N": "[nX2]1c[nX2]c(N)n1",  # diaminopyrimidine core of trimethoprim
        "aroyl_hydrazone": "C(=O)NN=Cc",
    },
    "FabI": {
        "phenol_like_oh_aromatic": "[OX2H][a]",
        "thioether_carbonyl": "[SX2]C(=O)",
    },
    "FtsZ": {
        "heteroaryl_scaffold": "c1n[n,c]nc1",
    },
    "Gyr": {
        "gyrase_pharmacophore": "c1n[n,c]c(N)c1",
        "xanthine_dione": "O=C1[n]c(=O)n",
    },
    "MurA": {
        "electrophile": "[C;$(C=[O,S])]",
    },
    "LpxC": {
        "hydroxamate": "C(=O)N[OX2H]",
    },
    "LpxH": {
        "pyridine_amine": "c1ccncc1[NX3]",
    },
    "PBP2a": {
        "beta_lactam": "C1C(=O)N1",
        "carboxylate": "[CX3](=O)[OX1H0-,OX2H1]",
    },
}


def flag_molecule(mol, flags):
    out = {}
    for k, smarts in flags.items():
        try:
            pat = Chem.MolFromSmarts(smarts)
            out[k] = int(mol.HasSubstructMatch(pat))
        except Exception:
            out[k] = 0
    return out


def class_pharm_flags(cls_name):
    return PHARM.get(cls_name, {})


# ---------------------------------------------------------------- scoring ----
def score_compounds(comp_df, out_prefix="scores"):
    results = []
    top_refs = {}  # per class: list of dicts for figure support

    # precompute compound fingerprints
    comp_fps = {row["compound"]: get_fp(row["mol"]) for _, row in comp_df.iterrows()}
    comp_maccs = {row["compound"]: MACCSkeys.GenMACCSKeys(row["mol"]) for _, row in comp_df.iterrows()}
    comp_pharm = {row["compound"]: flag_molecule(row["mol"],
                       {"xanthine": "O=C1[nX3]c(=O)[nX3]",
                        "benzimidazole": "c1c[nH]c2ccccc12",
                        "benzothiazole": "c1c2ccccc2nc1",
                        "sulfonamide": "[SX3](=O)(=O)[NX3]",
                        "triazolopyrimidine": "[n,c]1[n,c]n[n,c]1",
                        "oxadiazole": "c1nnoc1",
                        "hydrazone": "C(=O)NN=C",
                        "aryl_halide": "[a][F,Cl,Br]"})
                          for _, row in comp_df.iterrows()}

    classes = sorted(set(f.split("_")[-1].replace(".json", "")
                         for f in glob.glob(str(REF_DATA / 'ref_ligands_*.json'))
                         if "summary" not in f))

    for cls in classes:
        refs = json.load(open(REF_DATA / f"ref_ligands_{cls}.json"))
        if not refs:
            print(f"[skip] {cls}: no reference ligands")
            continue
        print(f"[score] {cls}: {len(refs)} reference ligands")
        ref_fps = []
        ref_maccs = []
        for r in refs:
            mol = Chem.MolFromSmiles(r["canonical_smiles"])
            ref_fps.append(get_fp(mol))
            ref_maccs.append(MACCSkeys.GenMACCSKeys(mol))
        # keep as plain lists (numpy object arrays iterate incorrectly)

        for _, row in comp_df.iterrows():
            c = row["compound"]
            fp = comp_fps[c]
            maccs = comp_maccs[c]
            # ECFP4 similarities
            sim = np.array([tanimoto(fp, r) for r in ref_fps])
            order = np.argsort(-sim)
            best_idx = order[0]
            top5_mean = sim[order[:5]].mean()
            maccs_sim = np.array([tanimoto(maccs, r) for r in ref_maccs])
            maccs_idx = int(np.argmax(maccs_sim)) if len(maccs_sim) else 0
            maccs_best = float(maccs_sim[maccs_idx]) if len(maccs_sim) else 0.0
            best = refs[best_idx]
            results.append({
                "compound": c, "target_class": cls,
                "ecfp4_max_tanimoto": float(sim[best_idx]),
                "ecfp4_top5_mean": float(top5_mean),
                "maccs_max_tanimoto": maccs_best,
                "best_ref_molecule": best["molecule_chembl_id"],
                "best_ref_pchembl": best["pchembl_value"],
                "best_ref_ic50_nm": float(best["standard_value"]),
                "best_ref_target": best["target_chembl_id"],
                "best_ref_organism": best["organism"],
                "n_references": len(refs),
                "ecfp4_2nd_tanimoto": float(sim[order[1]]) if len(order) > 1 else 0.0,
                "n_refs_above_0.4": int((sim >= 0.4).sum()),
                "n_refs_above_0.5": int((sim >= 0.5).sum()),
            })
        # Store potency-representative references for chemical-space figures.
        top_refs[cls] = sorted(refs, key=lambda x: -float(x.get('pchembl_value', 0)))[:10]

    scores = pd.DataFrame(results)
    scores.to_csv(f"{RES}/{out_prefix}.csv", index=False)
    json.dump(top_refs, open(f"{RES}/top_refs.json", "w"))
    json.dump(comp_pharm, open(f"{RES}/compound_pharm_flags.json", "w"), indent=1)
    comp_df.drop(columns=["mol"]).to_csv(f"{RES}/compound_properties.csv", index=False)
    print(f"Saved {RES}/{out_prefix}.csv with {len(scores)} rows")
    return scores


if __name__ == "__main__":
    comp_df = load_compounds()
    print(f"Loaded {len(comp_df)} compounds")
    scores = score_compounds(comp_df)
