from pathlib import Path
import os, pandas as pd

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1])); RES=ROOT/'results'

def assay_plan(target, role, localization, clinical):
    role=str(role).lower(); loc=str(localization).lower(); target=str(target)
    direct='purified bacterial target inhibition or binding assay; species-orthologue comparison'
    cellular='MIC and time-kill; resistant-mutant selection; complementation or target-rescue experiment'
    if 'complex' in role or 'ribosome' in target.lower(): direct='purified target-complex binding or biochemical translation assay; species-specific ribosome comparison'
    if 'resistance' in role or 'beta-lactamase' in target.lower(): direct='enzyme inhibition assay with clinically relevant resistance variants; inhibitor-rescue panel'
    if 'membrane' in loc or 'envelope' in loc: cellular='MIC with permeability/efflux controls; membrane or envelope phenotype assay; resistant-mutant selection'
    if target in {'PBP','PBP2a'}: direct='PBP transpeptidase acylation/inhibition assay; species-appropriate PBP panel'
    if target in {'D-Ala-D-Ala'}: direct='ligand-binding or cell-wall precursor incorporation assay; precursor competition'
    safety='human orthologue/mitochondrial selectivity screen where homologous or mitochondrial risk is annotated'
    return direct,cellular,safety

def main():
    p=RES/'v2_open_target_predictions_by_organism.csv'
    df=pd.read_csv(p)
    # Top three hypotheses per private compound-organism, preserving uncertainty fields.
    top=df.sort_values(['query_id','organism','overall_priority_score'],ascending=[True,True,False]).groupby(['query_id','organism']).head(3).copy()
    top['priority_rank']=top.groupby(['query_id','organism'])['overall_priority_score'].rank(method='first',ascending=False).astype(int)
    plans=top.apply(lambda r: assay_plan(r.target_class,r.target_role,r.cellular_localization,r.clinical_status),axis=1,result_type='expand')
    plans.columns=['direct_target_assay','cellular_causality_assay','safety_followup']
    top=pd.concat([top.reset_index(drop=True),plans.reset_index(drop=True)],axis=1)
    keep=['query_id','organism','priority_rank','target_class','overall_priority_score','confidence_class','chemical_quality_adjusted_score','target_specificity_score','species_transfer_score','pocket_evidence_score','card_resistance_context_score','organism_specific_snp_row_count','anti_target_risk_score','uncertainty_reasons','recommended_validation','direct_target_assay','cellular_causality_assay','safety_followup']
    top[keep].to_csv(RES/'v2_compound_organism_validation_plan.csv',index=False)
    # Target-class summary for experimental ordering.
    g=df.groupby('target_class').agg(n_hypotheses=('overall_priority_score','size'),mean_priority=('overall_priority_score','mean'),max_priority=('overall_priority_score','max'),mean_chemical=('chemical_quality_adjusted_score','mean'),mean_transfer=('species_transfer_score','mean'),mean_pocket=('pocket_evidence_score','mean'),max_card_snp_rows=('card_snp_row_count','max'),max_anti_target_risk=('anti_target_risk_score','max')).reset_index()
    meta=df.sort_values('overall_priority_score',ascending=False).drop_duplicates('target_class')[['target_class','target_role','clinical_status','cellular_localization','recommended_validation']]
    g=g.merge(meta,on='target_class',how='left')
    g['recommended_sequence']='1) purified target/orthologue assay; 2) MIC and time-kill; 3) resistant-mutant or complementation study; 4) permeability/efflux controls where relevant'
    g.sort_values(['mean_priority','max_priority'],ascending=False).to_csv(RES/'v2_target_class_validation_plan.csv',index=False)
    print('compound-organism plans',len(top),'target-class plans',len(g))

if __name__=='__main__': main()
