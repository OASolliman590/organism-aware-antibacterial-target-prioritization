"""Fetch public ChEMBL reference ligands for additional clinical target families.

The target IDs were selected from ChEMBL target search results and are recorded in
 data/additional_target_hits.json. The resulting files are public reference data;
no unpublished compound structures are used.
"""
from pathlib import Path
import json, math, time
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data'/'reference_ligands'; OUT.mkdir(parents=True,exist_ok=True)
TARGETS={
 'RpoB':['CHEMBL1852','CHEMBL2006'],
 'LeuRS':['CHEMBL4295566','CHEMBL4295595','CHEMBL4662935','CHEMBL4105844'],
 'PBP':['CHEMBL3512','CHEMBL1813','CHEMBL4269','CHEMBL2363020','CHEMBL3309036'],
 'D-Ala-D-Ala':['CHEMBL1956','CHEMBL4523951','CHEMBL2030'],
 '70S_ribosome':['CHEMBL2363965','CHEMBL2364022','CHEMBL2363135','CHEMBL2363853','CHEMBL2364096'],
 'Beta-lactamase':['CHEMBL2026','CHEMBL4114','CHEMBL1293246','CHEMBL3499','CHEMBL1667700','CHEMBL3562180','CHEMBL1667679'],
}

def pchembl(a):
    if a.get('pchembl_value') not in (None,''):
        try:return float(a['pchembl_value'])
        except:pass
    st=a.get('standard_type','')
    try:
        val=float(a.get('standard_value'))
        unit=str(a.get('standard_units','')).lower()
        if val<=0:return None
        if unit=='nm': val=val/1000
        if unit not in {'nm','um','µm','mm'} and unit!='': return None
        return -math.log10(val*1e-6) if unit in {'um','µm'} else -math.log10(val*1e-9) if unit=='nm' else -math.log10(val*1e-3)
    except:return None

def fetch_target(tid):
    out=[]; offset=0
    while True:
        url='https://www.ebi.ac.uk/chembl/api/data/activity.json'
        r=requests.get(url,params={'target_chembl_id':tid,'limit':1000,'offset':offset},timeout=60)
        r.raise_for_status(); d=r.json(); acts=d.get('activities',[])
        if not acts: break
        out.extend(acts)
        meta=d.get('page_meta',{})
        if len(acts)<1000 or offset+len(acts)>=int(meta.get('total_count',offset+len(acts))): break
        offset += len(acts); time.sleep(.2)
    return out

for cls,ids in TARGETS.items():
    records=[]; seen=set()
    for tid in ids:
        try: acts=fetch_target(tid)
        except Exception as e:
            print('FAILED',cls,tid,repr(e)); continue
        for a in acts:
            if a.get('canonical_smiles') in (None,'') or a.get('molecule_chembl_id') in (None,''): continue
            pc=pchembl(a)
            if pc is None or pc<4.5: continue
            smi=a.get('canonical_smiles')
            key=(a.get('molecule_chembl_id'),tid)
            if key in seen: continue
            seen.add(key)
            records.append({'molecule_chembl_id':a.get('molecule_chembl_id'),'canonical_smiles':smi,
                            'pchembl_value':round(pc,3),'standard_type':a.get('standard_type'),
                            'standard_value':a.get('standard_value'),'target_chembl_id':tid,
                            'organism':a.get('target_organism') or a.get('organism') or ''})
    # Deduplicate same molecule/SMILES, preserve best pChEMBL and one target context.
    best={}
    for r in records:
        key=(r['molecule_chembl_id'],r['canonical_smiles'])
        if key not in best or r['pchembl_value']>best[key]['pchembl_value']: best[key]=r
    out=sorted(best.values(),key=lambda x:-x['pchembl_value'])[:2000]
    (OUT/f'ref_ligands_{cls}.json').write_text(json.dumps(out,indent=2))
    print(cls,'records',len(out))
