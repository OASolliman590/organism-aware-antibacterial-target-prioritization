from pathlib import Path
import os
import pandas as pd

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1]))
RES=ROOT/'results'
m=pd.read_csv(RES/'benchmark_metrics.csv')
covered=m[m.target_in_reference_universe==1].copy()
rows=[]
for label,g in m.groupby('query_target_label'):
    c=g[g.target_in_reference_universe==1]
    rows.append({'query_target_label':label,'n_queries':len(g),'n_covered':len(c),
                 'coverage_fraction':len(c)/len(g) if len(g) else 0,
                 'top1_recall_covered':c.top1_hit.mean() if len(c) else None,
                 'top3_recall_covered':c.top3_hit.mean() if len(c) else None,
                 'top5_recall_covered':c.top5_hit.mean() if len(c) else None,
                 'mrr_covered':c.reciprocal_rank.mean() if len(c) else None})
summary=pd.DataFrame(rows)
overall={'query_target_label':'ALL','n_queries':len(m),'n_covered':len(covered),
         'coverage_fraction':len(covered)/len(m),
         'top1_recall_covered':covered.top1_hit.mean(),
         'top3_recall_covered':covered.top3_hit.mean(),
         'top5_recall_covered':covered.top5_hit.mean(),
         'mrr_covered':covered.reciprocal_rank.mean()}
summary=pd.concat([summary,pd.DataFrame([overall])],ignore_index=True)
summary.to_csv(RES/'benchmark_summary.csv',index=False)
print(summary.to_string(index=False))
