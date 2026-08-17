"""V2 open target discovery for protected local compounds and public benchmarks.

The output keeps chemical evidence, reference quality, species transfer, biological
priority, anti-target risk, and overall prioritization as separate auditable fields.
Cross-target molecules are decoys for specificity analysis only; they are not called
experimentally inactive.
"""
from pathlib import Path
import os, json, math
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1]))
REF=ROOT/'data'/'reference_ligands'
ONTO=ROOT/'data'/'target_ontology_v2.csv'
QUALITY=ROOT/'data'/'reference_quality'/'target_reference_quality_v2.csv'
COMPAT=ROOT/'data'/'species_targets'/'species_target_compatibility.csv'
CARD_SUM=ROOT/'data'/'resistance_v2'/'card_resistance_family_summary_v2.csv'
CARD_SNP=ROOT/'data'/'resistance_v2'/'card_snp_family_summary_v2.csv'
CARD_SNP_ORG=ROOT/'data'/'resistance_v2'/'card_snp_organism_family_summary_v2.csv'
STRUCT_SUM=ROOT/'data'/'structures_v2'/'rcsb_structure_summary_v2.csv'
ANN_OLD=ROOT/'data'/'target_annotations.csv'
BENCH=ROOT/'data'/'benchmark'/'eskape_benchmark_drugs.csv'
RES=ROOT/'results'; RES.mkdir(exist_ok=True)

ORGANISMS=['Klebsiella pneumoniae','Bacillus cereus','Escherichia coli','Proteus mirabilis','Acinetobacter baumannii','Staphylococcus aureus']
ORGANISM_ALIASES={'MRSA / Staphylococcus aureus':'Staphylococcus aureus'}
GRAM_NEG={'Klebsiella pneumoniae','Escherichia coli','Proteus mirabilis','Acinetobacter baumannii'}
GRAM_POS={'Bacillus cereus','Staphylococcus aureus'}


def fp(m): return AllChem.GetMorganGenerator(radius=2,fpSize=2048).GetFingerprint(m)
def maccs(m): return MACCSkeys.GenMACCSKeys(m)
def sim(a,b): return float(DataStructs.TanimotoSimilarity(a,b))
def mol(s):
    x=Chem.MolFromSmiles(str(s)) if s and str(s)!='nan' else None
    return Chem.RemoveHs(x) if x else None


def load_refs():
    refs={}
    for p in sorted(REF.glob('ref_ligands_*.json')):
        if p.stem.endswith('summary'): continue
        cls=p.stem.replace('ref_ligands_',''); cls='GyrB' if cls=='Gyr' else cls
        for r in json.loads(p.read_text()):
            m=mol(r.get('canonical_smiles'))
            if m is None: continue
            smi=Chem.MolToSmiles(m)
            refs.setdefault(cls,[]).append({**r,'_mol':m,'_fp':fp(m),'_maccs':maccs(m),'_smi':smi})
    out={}
    for cls,rows in refs.items():
        seen=set(); out[cls]=[]
        for r in rows:
            if r['_smi'] in seen: continue
            seen.add(r['_smi']); out[cls].append(r)
    return out


def load_queries(path):
    out=[]
    if not path.exists(): return out
    for m in Chem.SDMolSupplier(str(path),removeHs=True):
        if m is None: continue
        name=m.GetProp('_Name') if m.HasProp('_Name') else 'private_compound'
        out.append({'query_id':name,'query_name':name,'mol':m,'fp':fp(m),'maccs':maccs(m),'source':'private_local'})
    return out


def quality_score(grade):
    return {'usable':1.0,'moderate_redundancy':0.85,'low':0.60,'insufficient':0.0}.get(str(grade),0.0)

def clinical_value(x):
    x=str(x).lower()
    if 'approved' in x: return 1.0
    if 'clinical' in x: return .75
    if 'validated' in x: return .65
    if 'preclinical' in x: return .40
    return .20

def ordinal(x): return {'high':1.0,'medium':.65,'variable':.45,'low':.30}.get(str(x).lower(),.45)

def access_value(x): return {'cytosolic':.70,'periplasmic':.90,'outer_membrane':.55,'membrane':.45,'extracellular_or_periplasmic':.85}.get(str(x).lower(),.50)

def scope_score(scope,org,target):
    s=str(scope).lower()
    if target=='PBP2a' and org!='Staphylococcus aureus': return 0.05
    if 'mrsa' in s and org=='Staphylococcus aureus': return 1.0
    if 'gram-negative' in s and org in GRAM_NEG: return 1.0
    if 'gram-positive' in s and org in GRAM_POS: return 1.0
    if 'broad' in s or 'bacteria' in s: return .85
    if 'staphylococci' in s and org=='Staphylococcus aureus': return 1.0
    if org in s: return 1.0
    return .25

def anti_target_annotation(target):
    # Annotation-only risk prior; no claim of measured human off-target activity.
    high={'DHFR':('human DHFR homolog; inspect selectivity','high'),'LeuRS':('human LARS homolog; inspect selectivity','high'),'70S_ribosome':('mitochondrial translation selectivity','medium'),'30S_ribosome':('mitochondrial translation selectivity','medium'),'50S_ribosome':('mitochondrial translation selectivity','medium')}
    if target in high: note,r=high[target]; return (0.75 if r=='high' else 0.50,note,'annotation_only')
    return (0.10,'no direct human orthologue risk assigned in ontology; safety remains untested','annotation_only')

def validation_plan(target,role):
    if target in {'GyrB','TopoIV','FtsZ','DHFR','FabI','FabH','MurA','MurC','MurE','LpxA','LpxC','LpxH','RpoB','LeuRS','PBP2a'}:
        return 'purified-protein inhibition or binding; species-orthologue assay; resistant-mutant or complementation test'
    if 'ribosome' in target or target=='D-Ala-D-Ala':
        return 'target-complex binding or biochemical translation assay; cell-based target-dependence experiment'
    return 'phenotypic assay followed by mechanism-specific orthogonal validation'

def score_query(q,refs,quality,compat,ontology,exclude_close=False,cutoff=.85):
    rows=[]
    for cls,rs in refs.items():
        kept=[]; excluded=0
        for r in rs:
            if exclude_close and sim(q['fp'],r['_fp'])>=cutoff: excluded+=1
            else: kept.append(r)
        if not kept: continue
        e=np.array([sim(q['fp'],r['_fp']) for r in kept]); k=np.array([sim(q['maccs'],r['_maccs']) for r in kept])
        order=np.argsort(-e); topn=min(5,len(e)); best=kept[int(order[0])]
        # Cross-target molecules are decoys for specificity only.
        other=[]
        for other_cls,other_rs in refs.items():
            if other_cls==cls: continue
            other.extend(other_rs)
        decoy_max=max([sim(q['fp'],r['_fp']) for r in other],default=0.0)
        margin=float(e[order[0]]-decoy_max)
        qrow=quality[quality.target_class==cls]
        if len(qrow): qr=qrow.iloc[0]; qscore=quality_score(qr.quality_grade); nref=int(qr.n_valid_ligands); nscaf=int(qr.n_unique_scaffolds); grade=qr.quality_grade
        else: qscore=0.0; nref=len(kept); nscaf=0; grade='insufficient'
        ec=(float(e[order[0]])-.10)/.55; top=(float(e[order[:topn]].mean())-.08)/.45; mc=(float(k.max())-.10)/.70
        chem=float(np.clip(.50*np.clip(ec,0,1)+.25*np.clip(top,0,1)+.15*np.clip(mc,0,1),0,1))
        specificity=float(np.clip(.5+.5*margin/.25,0,1))
        chem_quality=float(chem*qscore*specificity)
        rows.append({'query_id':q['query_id'],'source':q['source'],'target_class':cls,
                     'ecfp4_max':float(e[order[0]]),'ecfp4_top5_mean':float(e[order[:topn]].mean()),'maccs_max':float(k.max()),
                     'cross_target_decoy_max':decoy_max,'target_specificity_margin':margin,'target_specificity_score':specificity,
                     'chemical_evidence_score':chem,'reference_quality_score':qscore,'chemical_quality_adjusted_score':chem_quality,
                     'reference_quality_grade':grade,'n_references_after_exclusion':len(kept),'n_close_references_excluded':excluded,
                     'n_unique_scaffolds':nscaf,'best_reference_molecule':best.get('molecule_chembl_id',''),'best_reference_organism':best.get('organism',''),
                     'query_target_label':q.get('target_class',''),'query_mechanism_class':q.get('mechanism_class',''),'query_organisms':q.get('organisms','')})
    return pd.DataFrame(rows)

def apply_biology(raw,ontology,compat,card_summary=None,snp_summary=None,snp_org=None,struct_summary=None):
    if raw.empty: return raw
    ann=ontology.set_index('target_class'); rows=[]
    for _,r in raw.iterrows():
        if r.target_class not in ann.index: continue
        a=ann.loc[r.target_class]
        for org in ORGANISMS:
            c=compat[(compat.organism==org)&(compat.target_class==r.target_class)] if not compat.empty else pd.DataFrame()
            transfer=float(c.iloc[0].species_transfer_score) if len(c) and pd.notna(c.iloc[0].species_transfer_score) else 0.0
            mapping_status=c.iloc[0].sequence_status if len(c) else 'no_mapping_record'
            scope=scope_score(a.organism_scope,org,r.target_class)
            clinical=clinical_value(a.clinical_status); essential=ordinal(a.essentiality_level); access=access_value(a.cellular_localization); resistance=ordinal(a.resistance_relevance)
            resistance_family={'GyrB':'GyrB/TopoIV resistance','TopoIV':'GyrB/TopoIV resistance','RpoB':'RpoB resistance','PBP2a':'PBP2a','DHFR':'DHFR resistance','DHPS':'DHFR/DHPS resistance','D-Ala-D-Ala':'D-Ala-D-Ala resistance','70S_ribosome':'Ribosome resistance','30S_ribosome':'Ribosome resistance','50S_ribosome':'Ribosome resistance','LpxA':'Lipid-A/envelope resistance','LpxC':'Lipid-A/envelope resistance','LpxH':'Lipid-A/envelope resistance','MurA':'MurA-pathway resistance','Beta-lactamase':'Beta-lactamase'}.get(r.target_class,'')
            card_models=0; snp_rows=0; org_snp_rows=0
            if card_summary is not None and resistance_family:
                z=card_summary[card_summary.target_resistance_family==resistance_family]
                card_models=int(z.n_models.iloc[0]) if len(z) else 0
            if snp_summary is not None and resistance_family:
                z=snp_summary[snp_summary.resistance_family==resistance_family]
                snp_rows=int(z.n_snp_rows.iloc[0]) if len(z) else 0
            if snp_org is not None and resistance_family:
                z=snp_org[(snp_org.resistance_family==resistance_family)&(snp_org.organism.str.contains(org.split()[0],case=False,na=False))]
                org_snp_rows=int(z.n_snp_rows.sum()) if len(z) else 0
            card_context=float(np.clip(.50*(1 if card_models else 0)+.50*(1 if snp_rows else 0),0,1))
            if struct_summary is not None:
                z=struct_summary[struct_summary.target_class==r.target_class]
                struct_candidates=int(z.n_search_candidates.iloc[0]) if len(z) else 0; co_crystal=int(z.n_with_co_crystal_ligand.iloc[0]) if len(z) else 0
            else: struct_candidates=0; co_crystal=0
            pocket=float(1.0 if co_crystal else .60 if struct_candidates else 0.0)
            biological=.25*scope+.20*clinical+.20*essential+.15*access+.10*resistance+.10*card_context
            anti,anti_note,anti_status=anti_target_annotation(r.target_class)
            overall=float(r.chemical_quality_adjusted_score*(.65+.35*transfer)*(.75+.25*pocket)*(.50+.50*biological)*(1-.20*anti))
            reasons=[]
            if r.reference_quality_grade in {'low','insufficient'}: reasons.append('reference coverage limited')
            if mapping_status!='mapped' or transfer<.5: reasons.append('species sequence mapping unresolved or weak')
            if r.target_specificity_score<.5: reasons.append('similarity is not target-specific versus cross-target decoys')
            if anti>=.5: reasons.append('human-homologue or mitochondrial selectivity risk is annotation-only')
            if pocket==0: reasons.append('no RCSB co-crystal/pocket evidence in bounded public catalog')
            if overall>=.50 and r.chemical_quality_adjusted_score>=.45 and transfer>=.70: conf='High'
            elif overall>=.25 and r.chemical_quality_adjusted_score>=.25: conf='Moderate'
            elif r.chemical_evidence_score>=.20: conf='Low'
            else: conf='Insufficient'
            rr=r.to_dict(); rr.update({'organism':org,'species_transfer_score':transfer,'sequence_mapping_status':mapping_status,
                'organism_scope_score':scope,'clinical_priority_score':clinical,'essentiality_score':essential,'cellular_access_score':access,
                'resistance_relevance_score':resistance,'card_resistance_context_score':card_context,'card_model_count':card_models,'card_snp_row_count':snp_rows,'organism_specific_snp_row_count':org_snp_rows,
                'rcsb_structure_candidate_count':struct_candidates,'rcsb_co_crystal_ligand_count':co_crystal,'pocket_evidence_score':pocket,'biological_priority_score':biological,'anti_target_risk_score':anti,
                'anti_target_evidence_status':anti_status,'anti_target_note':anti_note,'overall_priority_score':overall,
                'confidence_class':conf,'uncertainty_reasons':'; '.join(reasons) if reasons else 'none',
                'recommended_validation':validation_plan(r.target_class,a.target_role),'clinical_status':a.clinical_status,
                'target_role':a.target_role,'organism_scope':a.organism_scope,'cellular_localization':a.cellular_localization,
                'resistance_relevance':a.resistance_relevance})
            rows.append(rr)
    return pd.DataFrame(rows)

def main():
    refs=load_refs(); quality=pd.read_csv(QUALITY) if QUALITY.exists() else pd.DataFrame()
    ontology=pd.read_csv(ONTO); compat=pd.read_csv(COMPAT) if COMPAT.exists() else pd.DataFrame()
    card_summary=pd.read_csv(CARD_SUM) if CARD_SUM.exists() else pd.DataFrame()
    snp_summary=pd.read_csv(CARD_SNP) if CARD_SNP.exists() else pd.DataFrame()
    snp_org=pd.read_csv(CARD_SNP_ORG) if CARD_SNP_ORG.exists() else pd.DataFrame(columns=['organism','resistance_family','n_snp_rows'])
    struct_summary=pd.read_csv(STRUCT_SUM) if STRUCT_SUM.exists() else pd.DataFrame()
    private=load_queries(ROOT/'data'/'compounds'/'compounds_normalized.sdf')
    bench=load_queries(ROOT/'data'/'benchmark'/'benchmark_structures.sdf')
    if not bench and BENCH.exists():
        bdf=pd.read_csv(BENCH)
        for _,b in bdf.iterrows():
            m=mol(b.canonical_smiles)
            if m: bench.append({'query_id':b.drug,'query_name':b.drug,'mol':m,'fp':fp(m),'maccs':maccs(m),'source':'eskape_benchmark',**b.to_dict()})
    private_scores=[]
    for q in private:
        s=score_query(q,refs,quality,compat,ontology,exclude_close=False)
        if not s.empty: private_scores.append(s)
    if private_scores:
        raw=pd.concat(private_scores,ignore_index=True); raw.to_csv(RES/'v2_open_target_scores_private.csv',index=False)
        ranked=apply_biology(raw,ontology,compat,card_summary,snp_summary,snp_org,struct_summary); ranked.to_csv(RES/'v2_open_target_predictions_by_organism.csv',index=False)
        ranked.sort_values(['organism','query_id','overall_priority_score'],ascending=[True,True,False]).groupby(['organism','query_id']).head(10).to_csv(RES/'v2_open_target_shortlist_by_organism.csv',index=False)
        ranked[ranked.confidence_class.isin(['High','Moderate'])].sort_values('overall_priority_score',ascending=False).to_csv(RES/'v2_validation_priority_candidates.csv',index=False)
        print('private queries',len(private),'raw rows',len(raw),'ranked rows',len(ranked))
    else: print('No protected private structures available; public-only run completed.')
    bench_scores=[]
    for q in bench:
        s=score_query(q,refs,quality,compat,ontology,exclude_close=True)
        if not s.empty: bench_scores.append(s)
    if bench_scores:
        bp=pd.concat(bench_scores,ignore_index=True); bp.to_csv(RES/'v2_benchmark_open_target_scores.csv',index=False)
        # Benchmarks are summarized in the dedicated v2 benchmark script.
        print('benchmark queries',len(bench),'rows',len(bp))
    pd.DataFrame({'target_class':sorted(refs),'n_reference_ligands':[len(refs[k]) for k in sorted(refs)]}).to_csv(RES/'v2_open_target_reference_coverage.csv',index=False)

if __name__=='__main__': main()
