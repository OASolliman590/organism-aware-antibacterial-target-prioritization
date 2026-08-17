"""Summarize organism-panel rankings and prior-assignment agreement."""
from pathlib import Path
import os
import pandas as pd
import numpy as np

WORK = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
RES = WORK / 'results'

panel = pd.read_csv(RES / 'organism_panel_rankings.csv')
short = pd.read_csv(RES / 'organism_target_shortlist.csv')
ranked = pd.read_csv(RES / 'ranked_target_evidence.csv')
props = pd.read_csv(RES / 'compound_properties.csv')

# One row per organism/compound/target with top two plus all panel rank details.
rows = []
for (org, comp), q in panel.groupby(['organism', 'compound']):
    q = q.sort_values('panel_rank')
    top = q.iloc[0]
    second = q.iloc[1]
    rows.append({
        'organism': org,
        'compound': comp,
        'top_target': top.target_class,
        'top_score': top.organism_adjusted_score,
        'top_ecfp4': top.ecfp4_max_tanimoto,
        'top_maccs': top.maccs_max_tanimoto,
        'top_sar_bonus': top.sar_bonus,
        'second_target': second.target_class,
        'second_score': second.organism_adjusted_score,
        'score_margin': top.organism_adjusted_score - second.organism_adjusted_score,
        'top_ref_ligand': top.best_ref_molecule,
        'top_ref_pchembl': top.best_ref_pchembl,
        'top_ref_organism': top.best_ref_organism,
        'top_n_references': top.n_references,
        'tpsa': top.tpsa,
        'gram_negative_tpsa_caution': int(org in {'Klebsiella pneumoniae','Escherichia coli','Proteus mirabilis','Acinetobacter baumannii'} and top.tpsa > 130),
    })
summary = pd.DataFrame(rows)
summary.to_csv(RES / 'organism_compound_top_predictions.csv', index=False)

# Target-level mean scores within each organism panel.
target_means = panel.groupby(['organism','target_class'], as_index=False).agg(
    mean_score=('organism_adjusted_score','mean'), median_score=('organism_adjusted_score','median'),
    max_score=('organism_adjusted_score','max'), n_compounds=('compound','nunique'),
    mean_ecfp4=('ecfp4_max_tanimoto','mean'), mean_maccs=('maccs_max_tanimoto','mean'),
    mean_sar_bonus=('sar_bonus','mean'), n_references=('n_references','first'))
target_means = target_means.sort_values(['organism','mean_score'], ascending=[True,False])
target_means.to_csv(RES / 'organism_target_mean_scores.csv', index=False)

# Compound series summary.
series = {'BI-1':'BI-series','BI-6':'BI-series','OX-11':'OX/T2Z-series','T2Z14':'OX/T2Z-series','T2Z5':'OX/T2Z-series','T2Z6':'OX/T2Z-series','T2Z9':'OX/T2Z-series','X1V11':'X1V-series','X1V19':'X1V-series','X1V20':'X1V-series','X1V26':'X1V-series','X1V9':'X1V-series'}
props['series'] = props.compound.map(series)
props.to_csv(RES / 'compound_properties_with_series.csv', index=False)

print('=== ORGANISM MEAN TARGET SCORES ===')
for org, q in target_means.groupby('organism'):
    print('\n' + org)
    print(q[['target_class','mean_score','max_score','mean_ecfp4','mean_maccs','mean_sar_bonus','n_references']].head(5).to_string(index=False, float_format=lambda x: f'{x:.3f}'))
print('\n=== SELECTED TOP PREDICTIONS ===')
print(summary.sort_values(['organism','compound']).to_string(index=False, float_format=lambda x: f'{x:.3f}'))
print('\nRows:', len(summary), 'ranked:', len(ranked), 'props:', len(props))
