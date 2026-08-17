from pathlib import Path
import os, pandas as pd
ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1])); RES=ROOT/'results'

def main():
    p=pd.read_csv(RES/'v2_open_target_predictions_by_organism.csv')
    s=pd.read_csv(ROOT/'data/species_targets/species_target_compatibility.csv')
    print('predictions',len(p),'queries',p.query_id.nunique(),'organisms',p.organism.nunique(),'targets',p.target_class.nunique())
    print('\nconfidence\n',p.groupby('confidence_class').size().to_string())
    print('\nmean priority by target\n',p.groupby('target_class').overall_priority_score.mean().sort_values(ascending=False).head(20).to_string())
    print('\nmean priority by organism/top3\n',p.sort_values('overall_priority_score',ascending=False).groupby('organism').head(3)[['organism','query_id','target_class','overall_priority_score','confidence_class','species_transfer_score','pocket_evidence_score']].to_string(index=False))
    print('\ncompatibility mapping status\n',s.groupby(['organism','sequence_status']).size().unstack(fill_value=0).to_string())
    print('\nbenchmark\n',pd.read_csv(RES/'benchmark_v2_summary.csv').to_string(index=False))
    print('\nvalidation plans',len(pd.read_csv(RES/'v2_compound_organism_validation_plan.csv')),len(pd.read_csv(RES/'v2_target_class_validation_plan.csv')))
    print('\ncard family summary\n',pd.read_csv(ROOT/'data/resistance_v2/card_snp_family_summary_v2.csv').sort_values('n_snp_rows',ascending=False).head(12).to_string(index=False))

if __name__=='__main__': main()
