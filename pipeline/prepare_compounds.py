"""Prepare the complete compound set.

The uploaded README.md contains the BI-1 MOL record, while the other 11 compounds
are separate SDF files. This script extracts BI-1, validates all structures with RDKit,
and writes a normalized 12-compound SDF plus a compound manifest.
"""
from pathlib import Path
import os
import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

PROJECT_ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
UPLOAD = Path(os.environ.get('INPUT_DIR', PROJECT_ROOT / 'inputs'))
OUT = PROJECT_ROOT / 'data' / 'compounds'
OUT.mkdir(parents=True, exist_ok=True)


def read_first_mol(path):
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=True)
    mols = [m for m in supplier if m is not None]
    return mols[0] if mols else None


def extract_readme_record(name):
    text = (UPLOAD / 'README.md').read_text(errors='replace')
    start = text.find(name + '\n')
    if start < 0:
        raise ValueError(f'{name} not found in README.md')
    end = text.find('$$$$', start)
    if end < 0:
        end = text.find('M  END', start)
        if end < 0:
            raise ValueError(f'end of {name} record not found')
        end += len('M  END')
    block = text[start:end]
    mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=True, strictParsing=False)
    if mol is None:
        raise ValueError(f'RDKit could not parse {name} MOL block')
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    mol.SetProp('_Name', name)
    return mol


names = ['BI-1', 'BI-6', 'OX-11', 'T2Z14', 'T2Z5', 'T2Z6', 'T2Z9',
         'X1V11', 'X1V19', 'X1V20', 'X1V26', 'X1V9']
writer = Chem.SDWriter(str(OUT / 'compounds_normalized.sdf'))
manifest = []
for name in names:
    # README.md is the authoritative name-to-structure source. Some uploaded SDF
    # filenames do not match their internal _Name fields, so do not use filenames
    # for assignment.
    path = UPLOAD / f'{name}.sdf'
    mol = extract_readme_record(name)
    if mol is None:
        raise ValueError(f'Could not parse {name}')
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    writer.write(mol)
    manifest.append({
        'compound': name,
        'source': 'README.md embedded MOL record',
        'canonical_smiles': smiles,
        'mw': round(Descriptors.MolWt(mol), 3),
        'logp': round(Descriptors.MolLogP(mol), 3),
        'tpsa': round(rdMolDescriptors.CalcTPSA(mol), 3),
        'hbd': rdMolDescriptors.CalcNumHBD(mol),
        'hba': rdMolDescriptors.CalcNumHBA(mol),
        'heavy_atoms': mol.GetNumHeavyAtoms(),
        'rings': rdMolDescriptors.CalcNumRings(mol),
    })
writer.close()
with (OUT / 'compound_manifest.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=manifest[0].keys())
    writer.writeheader()
    writer.writerows(manifest)
print(f'Prepared {len(manifest)} compounds: {", ".join(x["compound"] for x in manifest)}')
print(f'Wrote {OUT / "compounds_normalized.sdf"}')
print(f'Wrote {OUT / "compound_manifest.csv"}')
for row in manifest:
    print(row['compound'], row['mw'], row['tpsa'], row['canonical_smiles'])

if len(manifest) != 12:
    raise SystemExit('Expected 12 compounds')
