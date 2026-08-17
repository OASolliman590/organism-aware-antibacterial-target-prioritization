"""Compute transparent species-specific sequence compatibility metrics.

This module avoids claiming that a target is conserved merely because the gene name
matches. It performs a simple global protein alignment within each target class,
uses the longest mapped sequence as the reference, and reports aligned identity,
coverage, and a conservative transfer score. The metrics are for target mapping and
not binding prediction.
"""
from pathlib import Path
import csv, math
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'data'/'species_targets'/'species_target_proteins.csv'
OUT=ROOT/'data'/'species_targets'


def nw_identity(a,b,match=1,mismatch=-1,gap=-1):
    a=str(a); b=str(b)
    if not a or not b: return 0.0,0.0,0
    n,m=len(a),len(b)
    dp=[[0]*(m+1) for _ in range(n+1)]
    tr=[[None]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=i*gap; tr[i][0]='U'
    for j in range(1,m+1): dp[0][j]=j*gap; tr[0][j]='L'
    for i in range(1,n+1):
        for j in range(1,m+1):
            vals=[dp[i-1][j-1]+(match if a[i-1]==b[j-1] else mismatch),dp[i-1][j]+gap,dp[i][j-1]+gap]
            k=max(range(3),key=lambda x:vals[x]); dp[i][j]=vals[k]; tr[i][j]=('D','U','L')[k]
    i,j=n,m; aa=[]; bb=[]
    while i or j:
        t=tr[i][j]
        if t=='D': aa.append(a[i-1]); bb.append(b[j-1]); i-=1; j-=1
        elif t=='U': aa.append(a[i-1]); bb.append('-'); i-=1
        else: aa.append('-'); bb.append(b[j-1]); j-=1
    aa=''.join(reversed(aa)); bb=''.join(reversed(bb))
    aligned=sum(x!='-' and y!='-' for x,y in zip(aa,bb))
    ident=sum(x==y and x!='-' for x,y in zip(aa,bb))
    coverage=aligned/max(1,min(n,m))
    return ident/max(1,aligned),coverage,len(aa)


df=pd.read_csv(IN)
mapped=df[df.status=='mapped'].copy()
rows=[]
for target, grp in mapped.groupby('target_class'):
    grp=grp.sort_values(['reviewed','length','accession'],ascending=[True,False,True])
    ref=grp.iloc[0]
    for _,r in df[df.target_class==target].iterrows():
        if r.status!='mapped':
            rows.append({**r.to_dict(),'reference_accession':ref.accession,'reference_organism':ref.organism,'pairwise_identity':'','alignment_coverage':'','species_transfer_score':0.0,'sequence_status':'unmapped'})
            continue
        ident,cov,alen=nw_identity(r.sequence,ref.sequence)
        # Transfer score penalizes short/partial mappings and is not a binding probability.
        score=max(0.0,min(1.0,ident*math.sqrt(max(0.0,min(1.0,cov)))))
        rows.append({**r.to_dict(),'reference_accession':ref.accession,'reference_organism':ref.organism,'pairwise_identity':ident,'alignment_coverage':cov,'species_transfer_score':score,'sequence_status':'mapped'})

out=OUT/'species_target_compatibility.csv'
pd.DataFrame(rows).to_csv(out,index=False)
print(f'Wrote {len(rows)} compatibility rows to {out}')
print(pd.DataFrame(rows).groupby(['target_class','sequence_status']).size().to_string())
