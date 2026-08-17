from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from umap import UMAP

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1]))
RES=ROOT/'results'; FIG=RES/'figures'; FIG.mkdir(parents=True,exist_ok=True)
PRIVATE=ROOT/'data'/'compounds'/'compounds_normalized.sdf'

sns.set_theme(style='whitegrid',context='notebook')

def load_private():
    rows=[]
    if not PRIVATE.exists(): return rows
    for m in Chem.SDMolSupplier(str(PRIVATE),removeHs=True):
        if m is None: continue
        name=m.GetProp('_Name') if m.HasProp('_Name') else 'compound'
        rows.append((name,m))
    return rows

def fp(m): return AllChem.GetMorganGenerator(radius=2,fpSize=2048).GetFingerprint(m)
def mf(m): return MACCSkeys.GenMACCSKeys(m)

def heatmap(path,col,title,outfile,cmap='viridis'):
    d=pd.read_csv(path)
    if d.empty: return
    p=d.pivot_table(index='query_id',columns='target_class',values=col,aggfunc='max')
    p=p.fillna(0)
    plt.figure(figsize=(max(10,0.45*p.shape[1]),max(5,0.38*p.shape[0])))
    sns.heatmap(p,cmap=cmap,vmin=0,vmax=1,linewidths=.15,linecolor='white',cbar_kws={'label':col.replace('_',' ')})
    plt.title(title); plt.xlabel('open target class'); plt.ylabel('private query compound'); plt.xticks(rotation=45,ha='right'); plt.tight_layout(); plt.savefig(FIG/outfile,dpi=300); plt.close()

# 1. Private compound x open target similarity heatmaps.
heatmap(RES/'open_target_scores.csv','ecfp4_max','Open target discovery: ECFP4 nearest-neighbour similarity','open_target_ecfp4_heatmap.png','mako')
heatmap(RES/'open_target_scores.csv','maccs_max','Open target discovery: MACCS maximum similarity','open_target_maccs_heatmap.png','rocket')

# 2. Organism-specific target filtering: chemical evidence vs clinical priority.
p=RES/'open_target_predictions_by_organism.csv'
if p.exists():
    d=pd.read_csv(p)
    # Use one private compound at a time to avoid hiding query-specific evidence.
    top=d.sort_values('open_target_priority',ascending=False).groupby(['organism','query_id']).head(5)
    organisms=sorted(top.organism.unique())
    fig,axes=plt.subplots(2,3,figsize=(16,9),sharex=True,sharey=True)
    for ax,org in zip(axes.flat,organisms):
        g=top[top.organism==org]
        sns.scatterplot(data=g,x='chemical_evidence',y='organism_clinical_priority',hue='target_class',style='query_id',s=85,ax=ax,legend=False)
        ax.set_title(org); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xlabel('chemical target evidence'); ax.set_ylabel('organism/clinical priority')
    handles,labels=axes.flat[0].get_legend_handles_labels()
    if handles: fig.legend(handles,labels,loc='center left',bbox_to_anchor=(1.00,.5),title='target class',fontsize=8)
    fig.suptitle('Open target discovery followed by organism-specific clinical filtering',y=1.02)
    fig.tight_layout(); fig.savefig(FIG/'open_target_organism_clinical_filter.png',dpi=300,bbox_inches='tight'); plt.close(fig)

# 3. Private compounds and benchmark drugs in ECFP4 UMAP chemical space.
items=[]
for name,m in load_private(): items.append({'id':name,'source':'private compounds','mol':m})
bench=ROOT/'data'/'benchmark'/'eskape_benchmark_drugs.csv'
if bench.exists():
    b=pd.read_csv(bench)
    for _,r in b.iterrows():
        m=Chem.MolFromSmiles(r.canonical_smiles)
        if m is not None: items.append({'id':r.drug,'source':'benchmark drugs','mol':m})
if len(items)>=5:
    X=[]
    for x in items:
        bit=np.zeros(2048,dtype=np.float32); DataStructs.ConvertToNumpyArray(fp(x['mol']),bit); X.append(bit)
    emb=UMAP(n_neighbors=min(10,len(items)-1),min_dist=.25,metric='jaccard',random_state=17).fit_transform(np.asarray(X))
    df=pd.DataFrame({'id':[x['id'] for x in items],'source':[x['source'] for x in items],'x':emb[:,0],'y':emb[:,1]})
    plt.figure(figsize=(11,8)); sns.scatterplot(data=df,x='x',y='y',hue='source',style='source',s=75,palette={'private compounds':'#b64b4b','benchmark drugs':'#356fa8'})
    for _,r in df[df.source=='private compounds'].iterrows(): plt.text(r.x+.03,r.y+.03,r.id,fontsize=8)
    plt.title('ECFP4 chemical space: private compounds versus public benchmark drugs'); plt.xlabel('UMAP-1'); plt.ylabel('UMAP-2'); plt.tight_layout(); plt.savefig(FIG/'open_target_ecfp4_chemical_space.png',dpi=300); plt.close()

print('Generated revised open-target figures in',FIG)
