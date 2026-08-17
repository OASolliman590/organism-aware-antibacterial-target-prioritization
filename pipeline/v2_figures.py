from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1])); RES=ROOT/'results'; FIG=RES/'figures_v2'; FIG.mkdir(parents=True,exist_ok=True)
sns.set_theme(style='whitegrid',context='talk')


def savefig(name):
    plt.tight_layout(); plt.savefig(FIG/name,dpi=220,bbox_inches='tight'); plt.close()

def main():
    p=RES/'v2_open_target_predictions_by_organism.csv'
    if not p.exists(): raise FileNotFoundError(p)
    df=pd.read_csv(p)
    # Top target evidence heatmap by organism, averaged over private compounds.
    m=df.groupby(['organism','target_class'])['overall_priority_score'].mean().unstack(fill_value=0)
    plt.figure(figsize=(18,7)); sns.heatmap(m,cmap='mako',annot=False,linewidths=.25,cbar_kws={'label':'Mean overall priority'}); plt.xlabel('Open target class'); plt.ylabel('Bacterial organism'); plt.xticks(rotation=60,ha='right'); savefig('v2_overall_priority_heatmap.png')
    # ECFP/MACCS/biological/overall evidence decomposition for the strongest target per organism.
    top=df.sort_values(['organism','overall_priority_score'],ascending=[True,False]).groupby('organism').head(1).copy()
    cols=['chemical_quality_adjusted_score','species_transfer_score','pocket_evidence_score','biological_priority_score','overall_priority_score']
    z=top.set_index('organism')[cols].rename(columns={'chemical_quality_adjusted_score':'Chemical','species_transfer_score':'Species transfer','pocket_evidence_score':'Pocket','biological_priority_score':'Biology','overall_priority_score':'Overall'})
    z.plot(kind='bar',figsize=(16,7),color=['#3b82f6','#14b8a6','#f59e0b','#8b5cf6','#ef4444']); plt.ylim(0,1); plt.ylabel('Score'); plt.xlabel('Organism'); plt.xticks(rotation=35,ha='right'); plt.legend(ncol=3,loc='upper center',bbox_to_anchor=(.5,1.20)); savefig('v2_top_evidence_decomposition.png')
    # Confidence distribution per organism.
    c=pd.crosstab(df.organism,df.confidence_class).reindex(columns=['High','Moderate','Low','Insufficient'],fill_value=0)
    c.plot(kind='bar',stacked=True,figsize=(16,7),color=['#16a34a','#eab308','#f97316','#9ca3af']); plt.ylabel('Number of compound-target hypotheses'); plt.xlabel('Organism'); plt.xticks(rotation=35,ha='right'); plt.legend(title='Confidence'); savefig('v2_confidence_by_organism.png')
    # Sequence transfer heatmap, deduplicated target classes.
    if (ROOT/'data/species_targets/species_target_compatibility.csv').exists():
        s=pd.read_csv(ROOT/'data/species_targets/species_target_compatibility.csv'); s['species_transfer_score']=pd.to_numeric(s.species_transfer_score,errors='coerce').fillna(0)
        sm=s.pivot_table(index='organism',columns='target_class',values='species_transfer_score',aggfunc='max',fill_value=0)
        plt.figure(figsize=(18,7)); sns.heatmap(sm,cmap='viridis',vmin=0,vmax=1,cbar_kws={'label':'Species transfer score'}); plt.xlabel('Target class'); plt.ylabel('Organism'); plt.xticks(rotation=60,ha='right'); savefig('v2_species_transfer_heatmap.png')
    # Resistance/pocket context.
    ctx=df.groupby('target_class').agg(card_models=('card_model_count','max'),card_snp_rows=('card_snp_row_count','max'),co_crystal=('rcsb_co_crystal_ligand_count','max'),structure_candidates=('rcsb_structure_candidate_count','max')).sort_values('card_snp_rows',ascending=False).head(20)
    ax=ctx[['card_snp_rows','card_models','co_crystal']].plot(kind='bar',figsize=(18,7),logy=True,color=['#dc2626','#2563eb','#f59e0b']); ax.set_ylabel('Count (log scale)'); ax.set_xlabel('Target class'); plt.xticks(rotation=60,ha='right'); plt.legend(['CARD SNP rows','CARD models','RCSB co-crystal entries']); savefig('v2_resistance_structure_context.png')
    # Uncertainty: chemical quality vs calibrated overall score, colored by confidence.
    plt.figure(figsize=(11,8)); sns.scatterplot(data=df.sample(min(3000,len(df)),random_state=7),x='chemical_quality_adjusted_score',y='overall_priority_score',hue='confidence_class',style='organism',alpha=.65); plt.xlabel('Quality-adjusted chemical evidence'); plt.ylabel('Overall priority'); plt.title('Chemical evidence versus organism-aware priority'); savefig('v2_uncertainty_priority_scatter.png')
    # v2.1 parallel chemical-versus-clinical-translational diagnostic.
    if {'chemical_hypothesis_score','clinical_translation_score'}.issubset(df.columns):
        plt.figure(figsize=(12,8)); sns.scatterplot(data=df.sample(min(3000,len(df)),random_state=7),x='chemical_hypothesis_score',y='clinical_translation_score',hue='confidence_class',style='organism',alpha=.65)
        plt.xlim(0,1); plt.ylim(0,1); plt.xlabel('Chemical hypothesis score'); plt.ylabel('Clinical-translational score'); plt.title('v2.1: chemical compatibility versus clinical translation'); savefig('v21_chemical_vs_clinical_translation.png')
        hm=df.groupby(['organism','target_class'])['clinical_translation_score'].mean().unstack(fill_value=0)
        plt.figure(figsize=(18,7)); sns.heatmap(hm,cmap='crest',vmin=0,vmax=1,cbar_kws={'label':'Mean clinical-translational score'}); plt.xlabel('Target class / subtype'); plt.ylabel('Bacterial organism'); plt.xticks(rotation=60,ha='right'); savefig('v21_clinical_translation_heatmap.png')
    # Benchmark split/baseline figure if available.
    bp=RES/'benchmark_v2_summary.csv'
    if bp.exists():
        b=pd.read_csv(bp)
        if {'split','mode','top1_recall_covered'}.issubset(b.columns):
            plot=b[b['mode'].isin(['ecfp4_score','maccs_score','ensemble_score','prevalence_baseline'])].copy()
            plot['method']=plot['mode'].map({'ecfp4_score':'ECFP4','maccs_score':'MACCS','ensemble_score':'Ensemble','prevalence_baseline':'Prevalence'})
            plt.figure(figsize=(14,7)); sns.barplot(data=plot,x='split',y='top1_recall_covered',hue='method'); plt.ylim(0,1); plt.ylabel('Top-1 retrieval among covered queries'); plt.xlabel('Evaluation split'); plt.xticks(rotation=25,ha='right'); plt.legend(bbox_to_anchor=(1.02,1),loc='upper left'); savefig('v2_benchmark_baselines.png')
    print('created',len(list(FIG.glob('*.png'))),'figures in',FIG)

if __name__=='__main__': main()
