"""Calibrate v2 uncertainty with reference-set bootstrap stability and empirical decoy p-values.

These are uncertainty/stability measures, not probabilities of target binding. Bootstrap
sampling quantifies dependence on the public reference set; decoy p-values quantify
whether the target class is more similar than unrelated public target classes.
"""
from pathlib import Path
import os, json
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1]))
REF=ROOT/'data'/'reference_ligands'; RES=ROOT/'results'
N_BOOT=100; SEED=20260817

def mol(s):
    m=Chem.MolFromSmiles(str(s)) if s and str(s)!='nan' else None
    return Chem.RemoveHs(m) if m else None
def fp(m): return AllChem.GetMorganGenerator(radius=2,fpSize=2048).GetFingerprint(m)
def maccs(m): return MACCSkeys.GenMACCSKeys(m)
def sim(a,b): return float(DataStructs.TanimotoSimilarity(a,b))
def chem_score(ecfp_max,top5,maccs_max):
    return float(np.clip(.50*np.clip((ecfp_max-.10)/.55,0,1)+.25*np.clip((top5-.08)/.45,0,1)+.15*np.clip((maccs_max-.10)/.70,0,1),0,1))

def load_refs():
    refs={}
    for p in sorted(REF.glob('ref_ligands_*.json')):
        if p.stem.endswith('summary'): continue
        cls=p.stem.replace('ref_ligands_',''); cls='GyrB' if cls=='Gyr' else cls
        seen=set()
        for r in json.loads(p.read_text()):
            m=mol(r.get('canonical_smiles'))
            if m is None: continue
            smi=Chem.MolToSmiles(m)
            if smi in seen: continue
            seen.add(smi); refs.setdefault(cls,[]).append((fp(m),maccs(m)))
    return refs

def load_private():
    out=[]; path=ROOT/'data'/'compounds'/'compounds_normalized.sdf'
    if not path.exists(): return out
    for m in Chem.SDMolSupplier(str(path),removeHs=True):
        if m is None: continue
        name=m.GetProp('_Name') if m.HasProp('_Name') else 'private_compound'
        out.append((name,fp(m),maccs(m)))
    return out

def class_score(qfp,qmac,items,rng):
    if not items: return 0.0
    idx=rng.integers(0,len(items),size=len(items))
    ef=np.array([sim(qfp,items[i][0]) for i in idx]); mc=np.array([sim(qmac,items[i][1]) for i in idx])
    order=np.argsort(-ef); top=ef[order[:min(5,len(ef))]].mean()
    return chem_score(float(ef.max()),float(top),float(mc.max()))

def main():
    refs=load_refs(); queries=load_private(); rng=np.random.default_rng(SEED); rows=[]
    for qid,qfp,qmac in queries:
        boot_scores={c:[] for c in refs}
        for _ in range(N_BOOT):
            for c,items in refs.items(): boot_scores[c].append(class_score(qfp,qmac,items,rng))
        classes=list(refs)
        for c,items in refs.items():
            arr=np.array(boot_scores[c]); ranks=[]
            for b in range(N_BOOT): ranks.append(1+sum(boot_scores[x][b]>boot_scores[c][b] for x in classes if x!=c))
            # observed max similarity and empirical target-class-vs-decoy-class p-value
            target_max=max([sim(qfp,x[0]) for x in items],default=0.0)
            decoy_class_max=[max([sim(qfp,x[0]) for x in other],default=0.0) for x,other in refs.items() if x!=c]
            p=(1+sum(v>=target_max for v in decoy_class_max))/(1+len(decoy_class_max)) if decoy_class_max else 1.0
            rank_std=float(np.std(ranks)); top3=float(np.mean(np.array(ranks)<=3)); top1=float(np.mean(np.array(ranks)<=1))
            stability=float(.5*top3+.5*(1/(1+rank_std)))
            rows.append({'query_id':qid,'target_class':c,'bootstrap_n':N_BOOT,'bootstrap_chemical_mean':float(arr.mean()),'bootstrap_chemical_std':float(arr.std()),
                         'bootstrap_chemical_q025':float(np.quantile(arr,.025)),'bootstrap_chemical_q975':float(np.quantile(arr,.975)),
                         'bootstrap_rank_mean':float(np.mean(ranks)),'bootstrap_rank_std':rank_std,'bootstrap_top1_probability':top1,
                         'bootstrap_top3_probability':top3,'bootstrap_stability_score':stability,'empirical_decoy_p_value':p,
                         'uncertainty_note':'bootstrap reference-set stability; empirical decoy p is target-class specificity, not binding probability'})
    df=pd.DataFrame(rows); df.to_csv(RES/'v2_uncertainty_private.csv',index=False)
    # Merge into private ranked output and downgrade only when instability is strong.
    ranked_path=RES/'v2_open_target_predictions_by_organism.csv'
    if ranked_path.exists():
        ranked=pd.read_csv(ranked_path); ranked=ranked.merge(df,on=['query_id','target_class'],how='left')
        def recal(row):
            old=row.get('confidence_class','Insufficient')
            st=row.get('bootstrap_stability_score',0); p=row.get('empirical_decoy_p_value',1)
            if pd.isna(st): return old
            if old=='High' and (st<.55 or p>.35): return 'Moderate'
            if old=='Moderate' and (st<.35 or p>.60): return 'Low'
            return old
        ranked['confidence_class_calibrated']=ranked.apply(recal,axis=1)
        ranked['calibration_note']='Confidence downgraded when bootstrap target ranking was unstable or empirical decoy specificity was weak; not a binding probability.'
        ranked.to_csv(RES/'v2_open_target_predictions_by_organism_calibrated.csv',index=False)
        ranked.sort_values(['organism','query_id','overall_priority_score'],ascending=[True,True,False]).groupby(['organism','query_id']).head(10).to_csv(RES/'v2_open_target_shortlist_by_organism_calibrated.csv',index=False)
    print('Wrote',len(df),'uncertainty rows')
    print(df.groupby('target_class')[['bootstrap_stability_score','empirical_decoy_p_value']].mean().sort_values('bootstrap_stability_score',ascending=False).head(15).to_string())

if __name__=='__main__': main()
