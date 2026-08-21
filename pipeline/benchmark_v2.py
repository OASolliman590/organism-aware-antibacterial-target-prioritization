"""V2 benchmark for open antibacterial target discovery.

Evaluates close-analogue exclusion, exact Bemis-Murcko scaffold exclusion, ECFP4-only,
MACCS-only, ensemble, random and prevalence baselines. Target-family and temporal
splits are explicit diagnostics; if the public records lack the required metadata,
the result is reported as unavailable rather than fabricated.
"""
import json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles

try:
    from pipeline.config import load_config
except ModuleNotFoundError:  # direct ``python pipeline/<script>.py`` execution
    from config import load_config


CONFIG = load_config()
ROOT = CONFIG.root
REF = CONFIG.path_for("reference_ligands")
BENCH = CONFIG.path_for("benchmark")
SUBTYPE_ONTO = CONFIG.path_for("target_subtype_ontology")
RES = CONFIG.path_for("results"); FIG=RES/'figures'; RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
CHEM2D = CONFIG.value("chem2d")
V2_BENCHMARK = CONFIG.value("v2_benchmark")

ALIASES={
 'Gyrase/TopoIV':{'GyrB','TopoIV'},'Ribosome':{'70S_ribosome','30S_ribosome','50S_ribosome'},
 'D-Ala-D-Ala / cell wall':{'D-Ala-D-Ala'},'lipid A / membrane':{'LpxA','LpxC','LpxH'},
 'PBP':{'PBP','PBP2a'},'RpoB':{'RpoB'},'DHFR':{'DHFR'},'DHPS':{'DHPS'},'MurA':{'MurA'},'FabI':{'FabI'},
 'LeuRS':{'LeuRS'},'MurC':{'MurC'},'FtsZ':{'FtsZ'},'LpxC':{'LpxC'},'membrane':{'Membrane'},
 'Beta-lactamase':{'Beta-lactamase','Beta-lactamase_class_A','Beta-lactamase_class_B','Beta-lactamase_class_C','Beta-lactamase_class_D'},
}

def mol(s):
    m=Chem.MolFromSmiles(str(s)) if s and str(s)!='nan' else None
    return Chem.RemoveHs(m) if m else None
def fp(m):
    return AllChem.GetMorganGenerator(
        radius=int(CHEM2D["fingerprint_radius"]),
        fpSize=int(CHEM2D["fingerprint_bits"]),
    ).GetFingerprint(m)
def maccs(m): return MACCSkeys.GenMACCSKeys(m)
def sim(a,b): return float(DataStructs.TanimotoSimilarity(a,b))
def load_refs():
    out={}
    for p in sorted(REF.glob('ref_ligands_*.json')):
        if p.stem.endswith('summary'): continue
        cls=p.stem.replace('ref_ligands_',''); cls='GyrB' if cls=='Gyr' else cls
        for r in json.loads(p.read_text()):
            m=mol(r.get('canonical_smiles'))
            if m: out.setdefault(cls,[]).append({**r,'_mol':m,'_fp':fp(m),'_maccs':maccs(m),'_scaffold':MurckoScaffoldSmiles(mol=m)})
    # remove duplicate structures within target class
    for cls,rows in out.items():
        seen=set(); clean=[]
        for r in rows:
            s=Chem.MolToSmiles(r['_mol'])
            if s not in seen: seen.add(s); clean.append(r)
        out[cls]=clean
    return out

def load_bench():
    df=pd.read_csv(BENCH); out=[]
    for _,r in df.iterrows():
        m=mol(r.canonical_smiles)
        if m: out.append({**r.to_dict(),'_mol':m,'_fp':fp(m),'_maccs':maccs(m),'_scaffold':MurckoScaffoldSmiles(mol=m)})
    return out

def score(q,refs,split):
    close_cutoff = float(CHEM2D["close_analogue_cutoff"])
    ensemble_weights = V2_BENCHMARK["ensemble_weights"]
    rows=[]
    for cls,rs in refs.items():
        keep=[]; excluded=0
        for r in rs:
            close=sim(q['_fp'],r['_fp'])>=close_cutoff
            same_scaf=(q['_scaffold'] and r['_scaffold'] and q['_scaffold']==r['_scaffold'])
            remove=close if split=='close_analogue' else same_scaf if split=='scaffold' else False
            if remove: excluded+=1
            else: keep.append(r)
        if not keep: continue
        ef=np.array([sim(q['_fp'],r['_fp']) for r in keep]); mc=np.array([sim(q['_maccs'],r['_maccs']) for r in keep])
        rows.append({'query_id':q['drug'],'target_class':cls,'ecfp4_score':float(ef.max()),'maccs_score':float(mc.max()),
                     'ensemble_score':float(float(ensemble_weights["ecfp4"])*ef.max()+float(ensemble_weights["maccs"])*mc.max()),'n_refs_after_split':len(keep),'n_excluded':excluded})
    return pd.DataFrame(rows)

def accepted(label):
    label=str(label)
    acc=set(ALIASES.get(label,{label}))
    if SUBTYPE_ONTO.exists():
        ont=pd.read_csv(SUBTYPE_ONTO)
        rows=ont[ont.parent_target_class.isin(acc)] if 'parent_target_class' in ont.columns else pd.DataFrame()
        acc.update(rows.target_class.tolist())
    return acc
def evaluate(pred,mode):
    rows=[]
    for qid,g in pred.groupby('query_id'):
        g=g.sort_values(mode,ascending=False).reset_index(drop=True)
        label=str(g.query_target_label.iloc[0]); acc=accepted(label)
        hits=g.index[g.target_class.isin(acc)].tolist(); covered=int(bool(hits)); rank=int(hits[0]+1) if hits else np.nan
        n=len(g); random1=1/n if n else np.nan; random3=min(3,n)/n if n else np.nan; random5=min(5,n)/n if n else np.nan
        rows.append({'query_id':qid,'query_target_label':label,'mode':mode,'target_in_reference_universe':covered,'rank_of_known_target':rank,
                     'top1_hit':int(covered and rank<=1),'top3_hit':int(covered and rank<=3),'top5_hit':int(covered and rank<=5),
                     'reciprocal_rank':1/rank if covered else 0,'n_candidate_target_classes':n,'random_top1_baseline':random1,
                     'random_top3_baseline':random3,'random_top5_baseline':random5,'top_predicted_target':g.target_class.iloc[0],
                     'top_predicted_score':g[mode].iloc[0]})
    return pd.DataFrame(rows)

def prevalence_order(refs): return sorted(refs,key=lambda x:sum(len(refs[x]) for _ in [0]),reverse=True)
def evaluate_prevalence(bench,refs,split):
    order=prevalence_order(refs); rows=[]
    for q in bench:
        acc=accepted(q.get('target_class','')); hits=[i+1 for i,c in enumerate(order) if c in acc]; rank=hits[0] if hits else np.nan; n=len(order); covered=int(bool(hits))
        rows.append({'query_id':q['drug'],'query_target_label':q.get('target_class',''),'mode':'prevalence_baseline','target_in_reference_universe':covered,'rank_of_known_target':rank,
                     'top1_hit':int(covered and rank<=1),'top3_hit':int(covered and rank<=3),'top5_hit':int(covered and rank<=5),'reciprocal_rank':1/rank if covered else 0,
                     'n_candidate_target_classes':n,'random_top1_baseline':1/n,'random_top3_baseline':min(3,n)/n,'random_top5_baseline':min(5,n)/n,'split':split})
    return pd.DataFrame(rows)

def summary(metrics):
    rows=[]
    for (split,mode),g in metrics.groupby(['split','mode']):
        covered=g[g.target_in_reference_universe==1]
        rows.append({'split':split,'mode':mode,'n_queries':len(g),'n_covered':int(g.target_in_reference_universe.sum()),'coverage':g.target_in_reference_universe.mean(),
                     'top1_recall_covered':covered.top1_hit.mean() if len(covered) else np.nan,'top3_recall_covered':covered.top3_hit.mean() if len(covered) else np.nan,
                     'top5_recall_covered':covered.top5_hit.mean() if len(covered) else np.nan,'mrr_all':g.reciprocal_rank.mean(),
                     'random_top1':g.random_top1_baseline.mean(),'random_top3':g.random_top3_baseline.mean(),'random_top5':g.random_top5_baseline.mean(),
                     'top1_enrichment_over_random':(covered.top1_hit.mean()/g.random_top1_baseline.mean()) if len(covered) else np.nan,
                     'top3_enrichment_over_random':(covered.top3_hit.mean()/g.random_top3_baseline.mean()) if len(covered) else np.nan})
    return pd.DataFrame(rows)

def main():
    refs=load_refs(); bench=load_bench(); all_rows=[]; query_rows=[]
    for split in V2_BENCHMARK["splits"]:
        for q in bench:
            p=score(q,refs,split)
            if p.empty: continue
            p['query_target_label']=q.get('target_class',''); p['split']=split; query_rows.append(p)
            for mode in ['ecfp4_score','maccs_score','ensemble_score']:
                m=evaluate(pd.DataFrame([{**r,'query_target_label':q.get('target_class','')} for r in p.to_dict('records')]),mode)
                m['split']=split; all_rows.append(m)
        prev=evaluate_prevalence(bench,refs,split); all_rows.append(prev)
    metrics=pd.concat(all_rows,ignore_index=True); metrics.to_csv(RES/'benchmark_v2_query_metrics.csv',index=False)
    qdf=pd.concat(query_rows,ignore_index=True); qdf.to_csv(RES/'benchmark_v2_target_scores_by_split.csv',index=False)
    sm=summary(metrics); sm.to_csv(RES/'benchmark_v2_summary.csv',index=False)
    # Metadata readiness checks are explicit; no temporal result is fabricated if years are absent.
    has_year=any('year' in r for cls in refs.values() for r in cls)
    status=pd.DataFrame([{'evaluation':'target_family_holdout','status':'diagnostic_only','reason':'Holding out the only reference universe for a target family makes direct retrieval impossible; use family-related transfer sets in a future benchmark.'},
                         {'evaluation':'temporal_split','status':'available' if has_year else 'not_available','reason':'Reference activity records lack standardized assay-year metadata in the cached public files.' if not has_year else 'Activity year metadata present; implement chronological cutoff.'},
                         {'evaluation':'species_holdout','status':'partial','reason':'The curated benchmark contains organism scope labels, but public reference records are not uniformly species-resolved.'}])
    status.to_csv(RES/'benchmark_v2_split_status.csv',index=False)
    plt.figure(figsize=(9,5)); plot=sm[sm['mode'].isin(['ecfp4_score','maccs_score','ensemble_score','prevalence_baseline'])].copy(); plot['metric_label']=plot['mode'].map({'ecfp4_score':'ECFP4','maccs_score':'MACCS','ensemble_score':'Ensemble','prevalence_baseline':'Prevalence'})
    sns.barplot(data=plot,x='split',y='top3_recall_covered',hue='metric_label'); plt.ylim(0,1.05); plt.ylabel('Top-3 recall among covered queries'); plt.title('V2 benchmark: scaffold and close-analogue leakage controls'); plt.tight_layout(); plt.savefig(FIG/'benchmark_v2_split_comparison.png',dpi=300); plt.close()
    plt.figure(figsize=(9,5)); plot=sm[sm['mode'].isin(['ensemble_score','prevalence_baseline'])].copy(); plot['metric_label']=plot['mode'].map({'ensemble_score':'Ensemble','prevalence_baseline':'Prevalence'}); plot=pd.melt(plot,id_vars=['split','metric_label'],value_vars=['top1_enrichment_over_random','top3_enrichment_over_random'],var_name='metric',value_name='enrichment'); sns.barplot(data=plot,x='metric',y='enrichment',hue='metric_label'); plt.axhline(1,color='black',lw=1); plt.ylabel('Enrichment over random'); plt.title('V2 benchmark enrichment over random retrieval'); plt.tight_layout(); plt.savefig(FIG/'benchmark_v2_enrichment.png',dpi=300); plt.close()
    print(sm.to_string(index=False)); print('\nSplit status:\n',status.to_string(index=False))

if __name__=='__main__': main()
