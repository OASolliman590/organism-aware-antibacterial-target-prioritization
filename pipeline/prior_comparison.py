"""Compare prior docking/MD target assignments to scored target classes and panels."""
from pathlib import Path
import os
import pandas as pd

WORK = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
RES = WORK / 'results'
CSV = Path(os.environ.get('PRIOR_CSV', WORK / 'inputs' / 'final_target_list_per_bacterium.csv'))

prior = pd.read_csv(CSV)
scored = {'DHFR','DHPS','FabH','FabI','FtsZ','Gyr','GyrB','LpxA','LpxC','LpxH','MurA','MurC','PBP2a','TopoIV'}
original_panel = {
    'Klebsiella pneumoniae': {'DHFR','FabI','TopoIV','MurA','LpxH'},
    'Bacillus cereus': {'FabI','DHFR','FtsZ','GyrB','MurA'},
    'Escherichia coli': {'MurA','DHFR','FabI','GyrB','LpxC'},
    'Proteus mirabilis': {'MurA','DHFR','FabI','GyrB','LpxC'},
    'Acinetobacter baumannii': {'PBP2a','FabI','GyrB','LpxA','MurC'},
    'MRSA / Staphylococcus aureus': {'PBP2a','DHFR','FabI','FtsZ','GyrB'},
}

def map_class(name):
    n = name.lower()
    if 'beta-lactamase' in n or 'blaz' in n or 'ampc' in n or 'oxa-23' in n or 'kpc-2' in n:
        return 'Beta-lactamase (unscored)'
    if 'fts' in n: return 'FtsZ'
    if 'mura' in n: return 'MurA'
    if 'gyrase' in n: return 'GyrB'
    if 'fabh' in n or '3-oxoacyl' in n: return 'FabH'
    if 'leurs' in n: return 'LeuRS (unscored)'
    return 'unmapped / unscored'

def expand(row):
    cls = map_class(row.final_target)
    in_universe = cls in scored
    in_panel = cls in original_panel.get(row.source_bacterium, set())
    return pd.Series({
        'source_bacterium': row.source_bacterium,
        'compound_scope': row.compound_scope,
        'prior_target': row.final_target,
        'recommended_pdb': row.recommended_pdb,
        'match_status': row.match_status,
        'priority': row.priority,
        'mapped_scored_class': cls,
        'in_scored_ligand_universe': int(in_universe),
        'in_original_organism_panel': int(in_panel),
        'validation_interpretation': ('direct scored overlap' if in_universe and in_panel else 'scored family but outside original panel' if in_universe else 'not scored; retain as orthogonal docking/MD hypothesis'),
    })
comparison = prior.apply(expand, axis=1)
comparison.to_csv(RES / 'prior_target_comparison.csv', index=False)
print(comparison.to_string(index=False))
print('\nSummary by interpretation:')
print(comparison.groupby('validation_interpretation').size().to_string())
