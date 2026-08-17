"""Build v2 reference-ligand quality and decoy evidence tables.

Positive ChEMBL records are kept separate from negative evidence. Cross-target
molecules are labelled as *decoys* and are never treated as experimentally inactive;
this prevents a common but serious target-prediction error.
"""
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT=Path(__file__).resolve().parents[1]
REF=ROOT/'data'/'reference_ligands'
ONTO=ROOT/'data'/'target_ontology_v2.csv'
OUT=ROOT/'data'/'reference_quality'
OUT.mkdir(parents=True,exist_ok=True)

ontology=pd.read_csv(ONTO)
all_rows=[]
for path in sorted(REF.glob('ref_ligands_*.json')):
    if path.stem.endswith('summary'): continue
    target=path.stem.replace('ref_ligands_','')
    try: records=json.loads(path.read_text())
    except Exception as e:
        print('WARN',path,e); continue
    for rec in records:
        smi=rec.get('canonical_smiles','')
        mol=Chem.MolFromSmiles(smi) if smi else None
        if not mol: continue
        can=Chem.MolToSmiles(mol)
        scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        pchembl=rec.get('pchembl_value','')
        try: pchembl=float(pchembl)
        except: pchembl=np.nan
        all_rows.append({**rec,'reference_class':target,'canonical_smiles_rdkit':can,'scaffold_smiles':scaffold,
                         'mw':Descriptors.MolWt(mol),'logp':Descriptors.MolLogP(mol),'hba':Lipinski.NumHAcceptors(mol),
                         'hbd':Lipinski.NumHDonors(mol),'valid_structure':True,'pchembl_numeric':pchembl})
refs=pd.DataFrame(all_rows)
refs.to_csv(OUT/'reference_ligands_v2.csv',index=False)

rows=[]
for _,ann in ontology.iterrows():
    target=ann.target_class
    g=refs[refs.reference_class==target].copy()
    n=len(g); n_scaf=g.scaffold_smiles.nunique() if n else 0
    p=g.pchembl_numeric.dropna()
    assay=g.standard_type.fillna('unknown').value_counts().to_dict() if n else {}
    org=g.organism.fillna('unknown').value_counts().head(8).to_dict() if n else {}
    target_ids=g.target_chembl_id.fillna('unknown').nunique() if n else 0
    min_ref=float(ann.min_reference_ligands)
    if n==0: grade='insufficient'
    elif n < min_ref: grade='low'
    elif n_scaf < max(3,min(10,n//5)): grade='moderate_redundancy'
    else: grade='usable'
    rows.append({
        'target_class':target,'n_valid_ligands':n,'n_unique_scaffolds':n_scaf,
        'scaffold_fraction':n_scaf/max(1,n),'n_chembl_targets':target_ids,
        'n_with_pchembl':len(p),'median_pchembl':float(p.median()) if len(p) else np.nan,
        'assay_types':json.dumps(assay,sort_keys=True),'organism_counts':json.dumps(org,sort_keys=True),
        'minimum_reference_ligands':min_ref,'quality_grade':grade,
        'clinical_status':ann.clinical_status,'target_role':ann.target_role,
        'sequence_mapping_required':ann.sequence_mapping_required,'pocket_mapping_required':ann.pocket_mapping_required,
        'evidence_warning':ann.notes,
    })
quality=pd.DataFrame(rows)
quality.to_csv(OUT/'target_reference_quality_v2.csv',index=False)

# Decoys are sampled from public ligands assigned to other classes, matched only by
# coarse molecular-weight bins. They are labelled decoys, not inactive measurements.
rng=np.random.default_rng(20260817)
decoy_rows=[]
for target in sorted(refs.reference_class.unique()):
    pos=refs[refs.reference_class==target]
    pool=refs[refs.reference_class!=target].copy()
    if len(pos)==0 or len(pool)==0: continue
    pos_bins=pd.cut(pos.mw,bins=[0,250,350,450,600,800,10000],include_lowest=True)
    for bin_value,group in pos.groupby(pos_bins,observed=True):
        candidates=pool[pd.cut(pool.mw,bins=[0,250,350,450,600,800,10000],include_lowest=True)==bin_value]
        if len(candidates)==0: continue
        n=min(len(group),len(candidates),25)
        pick=candidates.sample(n=n,random_state=20260817+len(decoy_rows))
        for _,r in pick.iterrows():
            decoy_rows.append({'query_target_class':target,'decoy_reference_class':r.reference_class,
                               'molecule_chembl_id':r.get('molecule_chembl_id',''),'canonical_smiles':r.canonical_smiles_rdkit,
                               'decoy_label':'cross_target_decoy_not_inactive','mw':r.mw,'logp':r.logp,
                               'scaffold_smiles':r.scaffold_smiles})
pd.DataFrame(decoy_rows).drop_duplicates(['query_target_class','molecule_chembl_id']).to_csv(OUT/'cross_target_decoys_v2.csv',index=False)
print('reference rows',len(refs),'quality rows',len(quality),'decoy rows',len(decoy_rows))
print(quality[['target_class','n_valid_ligands','n_unique_scaffolds','quality_grade']].to_string(index=False))
