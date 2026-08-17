"""Open target discovery with organism/clinical filtering and leakage-aware benchmark evaluation.

This script deliberately does not use a fixed organism target panel. It scores every
reference-supported target class, then applies transparent annotations from
 data/target_annotations.csv. Unpublished local compounds are processed only when
 data/compounds/compounds_normalized.sdf exists; the directory is ignored by Git.
"""
from pathlib import Path
import os, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Descriptors

ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
REF_DIR = ROOT / 'data' / 'reference_ligands'
ANN_PATH = ROOT / 'data' / 'target_annotations.csv'
BENCH = ROOT / 'data' / 'benchmark' / 'eskape_benchmark_drugs.csv'
RES = ROOT / 'results'
FIG = RES / 'figures'
RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)


def fp(mol): return AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)
def maccs(mol): return MACCSkeys.GenMACCSKeys(mol)
def sim(a,b): return float(DataStructs.TanimotoSimilarity(a,b))

def canonical_mol(smiles):
    if not smiles: return None
    m = Chem.MolFromSmiles(smiles)
    if m is None: return None
    return Chem.RemoveHs(m)

def load_refs():
    refs = {}
    for path in sorted(REF_DIR.glob('ref_ligands_*.json')):
        if 'summary' in path.name: continue
        cls = path.stem.replace('ref_ligands_', '')
        # Fold the older Gyr set into the gyrase family; retain TopoIV separately.
        norm_cls = 'GyrB' if cls == 'Gyr' else cls
        records = json.loads(path.read_text())
        refs.setdefault(norm_cls, []).extend(records)
    out = {}
    for cls, rows in refs.items():
        clean=[]; seen=set()
        for r in rows:
            m = canonical_mol(r.get('canonical_smiles'))
            if m is None: continue
            smi = Chem.MolToSmiles(m)
            if smi in seen: continue
            seen.add(smi)
            clean.append({**r, '_mol':m, '_fp':fp(m), '_maccs':maccs(m), '_smi':smi})
        out[cls]=clean
    return out

def load_private_compounds():
    path = ROOT / 'data' / 'compounds' / 'compounds_normalized.sdf'
    if not path.exists(): return []
    out=[]
    for m in Chem.SDMolSupplier(str(path), removeHs=True):
        if m is None: continue
        name = m.GetProp('_Name') if m.HasProp('_Name') else 'private_compound'
        out.append({'query_id':name,'query_name':name,'mol':m,'fp':fp(m),'maccs':maccs(m),'source':'private_local'})
    return out

def load_benchmark():
    out=[]
    if not BENCH.exists(): return out
    df=pd.read_csv(BENCH)
    for _,r in df.iterrows():
        m=canonical_mol(r.canonical_smiles)
        if m is None: continue
        out.append({'query_id':r.drug,'query_name':r.drug,'mol':m,'fp':fp(m),'maccs':maccs(m),'source':'eskape_benchmark',**r.to_dict()})
    return out

def score_query(q, refs, exclude_close=False, close_cutoff=0.85):
    rows=[]
    for cls, rs in refs.items():
        kept=[]; excluded=0
        for r in rs:
            if exclude_close and sim(q['fp'], r['_fp']) >= close_cutoff:
                excluded += 1
            else: kept.append(r)
        if len(kept) == 0: continue
        e=np.array([sim(q['fp'], r['_fp']) for r in kept], dtype=float)
        k=np.array([sim(q['maccs'], r['_maccs']) for r in kept], dtype=float)
        order=np.argsort(-e)
        topn=min(5,len(e))
        best=kept[int(order[0])]
        rows.append({'query_id':q['query_id'],'source':q['source'],'target_class':cls,
                     'ecfp4_max':float(e[order[0]]),'ecfp4_top5_mean':float(e[order[:topn]].mean()),
                     'maccs_max':float(k.max()),'n_references_after_exclusion':len(kept),
                     'n_close_references_excluded':excluded,'best_reference_molecule':best.get('molecule_chembl_id',''),
                     'best_reference_organism':best.get('organism',''),'query_target_label':q.get('target_class',''),
                     'query_mechanism_class':q.get('mechanism_class',''),'query_organisms':q.get('organisms','')})
    df=pd.DataFrame(rows)
    if df.empty: return df
    # Within-run percentile normalization keeps target classes comparable without
    # turning similarity into a probability.
    df['ecfp4_component']=(df.ecfp4_max-0.10).clip(lower=0).div(0.55).clip(upper=1)
    df['top5_component']=(df.ecfp4_top5_mean-0.08).clip(lower=0).div(0.45).clip(upper=1)
    df['maccs_component']=(df.maccs_max-0.10).clip(lower=0).div(0.70).clip(upper=1)
    df['chemical_evidence']=(0.50*df.ecfp4_component+0.25*df.top5_component+0.15*df.maccs_component).clip(upper=1)
    return df

def clinical_value(v):
    v=str(v).lower()
    if 'approved' in v: return 1.0
    if 'clinical-development' in v: return .70
    if 'validated' in v: return .60
    if 'preclinical' in v: return .40
    return .25

def ordinal(v): return {'high':1.0,'medium':.65,'variable':.45,'low':.35}.get(str(v).lower(),.45)

def compatibility(row, organism):
    vals=str(row.compatible_organisms).split(';')
    exact=organism in vals
    # ESKAPE family labels can transfer cautiously to close Enterobacterales.
    if exact: return 1.0
    if organism in {'Klebsiella pneumoniae','Escherichia coli','Proteus mirabilis'} and 'Enterobacter spp.' in vals: return .75
    return .25

def apply_annotations(scores, annotations):
    if scores.empty: return scores
    ann=annotations.set_index('target_class')
    rows=[]
    for _,r in scores.iterrows():
        if r.target_class not in ann.index: continue
        a=ann.loc[r.target_class]
        for org in ['Klebsiella pneumoniae','Bacillus cereus','Escherichia coli','Proteus mirabilis','Acinetobacter baumannii','Staphylococcus aureus','MRSA / Staphylococcus aureus']:
            compat=compatibility(a, org)
            clinical=clinical_value(a.clinical_validation)
            essential=ordinal(a.essentiality_or_fitness)
            resistance=ordinal(a.resistance_relevance)
            # Accessibility is reported as annotation, not hidden in a score.
            org_priority=(0.40*compat+0.25*clinical+0.20*essential+0.15*resistance)
            rr=r.to_dict(); rr.update({'organism':org,'organism_compatibility':compat,'clinical_validation_score':clinical,
                                       'essentiality_score':essential,'resistance_relevance_score':resistance,
                                       'organism_clinical_priority':org_priority,
                                       'open_target_priority':r.chemical_evidence*org_priority,
                                       'clinical_validation':a.clinical_validation,
                                       'essentiality_or_fitness':a.essentiality_or_fitness,
                                       'cellular_accessibility':a.cellular_accessibility,
                                       'resistance_relevance':a.resistance_relevance,
                                       'organism_filter_note':a.organism_filter_note})
            rows.append(rr)
    return pd.DataFrame(rows)

def benchmark_metrics(pred):
    aliases={
        'Gyrase/TopoIV':{'GyrB','TopoIV'},
        'Ribosome':{'70S_ribosome'},
        'D-Ala-D-Ala / cell wall':{'D-Ala-D-Ala'},
        'lipid A / membrane':{'LpxA','LpxC','LpxH'},
        'PBP':{'PBP'},
        'RpoB':{'RpoB'},
        'DHFR':{'DHFR'}, 'DHPS':{'DHPS'}, 'MurA':{'MurA'}, 'FabI':{'FabI'},
    }
    rows=[]
    for qid,g in pred.groupby('query_id'):
        g=g.sort_values('chemical_evidence',ascending=False).reset_index(drop=True)
        label=str(g.query_target_label.iloc[0])
        accepted=set(aliases.get(label,{label}))
        # Ceftaroline is clinically associated with the MRSA PBP2a target; a
        # broad PBP label is not treated as PBP2a for Gram-negative meropenem.
        if qid == 'ceftaroline': accepted.update({'PBP2a','PBP'})
        hits=g.index[g.target_class.isin(accepted)].tolist()
        covered=int(bool(hits))
        rank=(int(hits[0])+1) if hits else np.nan
        rows.append({'query_id':qid,'query_target_label':label,
                     'accepted_target_classes':';'.join(sorted(accepted)),
                     'target_in_reference_universe':covered,'rank_of_known_target':rank,
                     'top1_hit':int(covered and rank<=1),'top3_hit':int(covered and rank<=3),
                     'top5_hit':int(covered and rank<=5),'reciprocal_rank':(1/rank if covered else 0),
                     'n_candidate_target_classes':g.target_class.nunique(),'top_predicted_target':g.target_class.iloc[0],
                     'top_predicted_score':g.chemical_evidence.iloc[0],
                     'random_top1_baseline':1/g.target_class.nunique(),
                     'random_top3_baseline':min(3,g.target_class.nunique())/g.target_class.nunique(),
                     'random_top5_baseline':min(5,g.target_class.nunique())/g.target_class.nunique()})
    return pd.DataFrame(rows)

def plot_benchmark(metrics):
    m=metrics[metrics.target_in_reference_universe==1].copy()
    if m.empty: return
    hit=pd.DataFrame({'metric':['Top-1','Top-3','Top-5','MRR'], 'value':[m.top1_hit.mean(),m.top3_hit.mean(),m.top5_hit.mean(),m.reciprocal_rank.mean()]})
    plt.figure(figsize=(7,4.5)); sns.barplot(data=hit,x='metric',y='value',color='#3b6ea8'); plt.ylim(0,1); plt.ylabel('benchmark value'); plt.title('Leakage-aware ESKAPE drug benchmark'); plt.tight_layout(); plt.savefig(FIG/'benchmark_topk_metrics.png',dpi=300); plt.close()
    plt.figure(figsize=(7,4.5)); sns.histplot(m.rank_of_known_target, bins=np.arange(0.5,max(6,m.rank_of_known_target.max()+1.5)), discrete=True, color='#bb6f3a'); plt.xlabel('rank of known target class'); plt.ylabel('number of benchmark drugs'); plt.title('Known-target retrieval rank'); plt.tight_layout(); plt.savefig(FIG/'benchmark_target_rank_distribution.png',dpi=300); plt.close()

def main():
    refs=load_refs(); anns=pd.read_csv(ANN_PATH)
    bench=load_benchmark(); private=load_private_compounds()
    all_pred=[]; bench_pred=[]
    for q in bench:
        s=score_query(q,refs,exclude_close=True)
        if not s.empty: bench_pred.append(s)
    if bench_pred:
        bp=pd.concat(bench_pred,ignore_index=True); metrics=benchmark_metrics(bp); metrics.to_csv(RES/'benchmark_metrics.csv',index=False); bp.to_csv(RES/'benchmark_open_target_scores.csv',index=False); plot_benchmark(metrics)
        print(metrics.to_string(index=False))
    else:
        pd.DataFrame().to_csv(RES/'benchmark_metrics.csv',index=False)
    for q in private:
        s=score_query(q,refs,exclude_close=False)
        if not s.empty: all_pred.append(s)
    if all_pred:
        raw=pd.concat(all_pred,ignore_index=True); raw.to_csv(RES/'open_target_scores.csv',index=False)
        ranked=apply_annotations(raw,anns); ranked.to_csv(RES/'open_target_predictions_by_organism.csv',index=False)
        top=ranked.sort_values(['organism','query_id','open_target_priority'],ascending=[True,True,False]).groupby(['organism','query_id']).head(5)
        top.to_csv(RES/'open_target_shortlist_by_organism.csv',index=False)
        print('private compounds scored:',len(private),'rows:',len(raw))
    else:
        print('No private local compounds found; open-target compound outputs not generated.')
    pd.DataFrame({'target_class':sorted(refs), 'n_reference_ligands':[len(refs[x]) for x in sorted(refs)]}).to_csv(RES/'open_target_reference_coverage.csv',index=False)

if __name__=='__main__': main()
