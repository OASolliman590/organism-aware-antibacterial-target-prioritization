"""Build a conservative RCSB structure/pocket availability catalog.

Text search hits are marked as candidates, not as validated species-specific structures.
A co-crystallized non-solvent ligand is recorded as active-site evidence for future user
controlled docking. This module does not download or prepare PDB files and does not dock.
"""
from pathlib import Path
import time, requests, pandas as pd

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data'/'structures_v2'; OUT.mkdir(parents=True,exist_ok=True)
TERMS={
'DHFR':'dihydrofolate reductase bacterial','DHPS':'dihydropteroate synthase bacterial','FabI':'enoyl ACP reductase FabI bacterial',
'FabH':'beta-ketoacyl ACP synthase III FabH','FtsZ':'FtsZ bacterial','GyrB':'DNA gyrase B bacterial ATPase',
'TopoIV':'topoisomerase IV ParE ParC bacterial','RpoB':'RNA polymerase beta subunit RpoB bacterial','LeuRS':'leucyl tRNA synthetase bacterial',
'MurA':'UDP-N-acetylglucosamine enolpyruvyl transferase MurA','MurC':'UDP-N-acetylmuramoyl tripeptide synthetase MurC',
'MurE':'UDP-N-acetylmuramoyl tripeptide synthetase MurE','LpxA':'UDP-3-O-acyl-glucosamine N-acyltransferase LpxA',
'LpxC':'UDP-3-O-acyl-N-acetylglucosamine deacetylase LpxC','LpxH':'UDP-2,3-diacylglucosamine hydrolase LpxH',
'PBP2a':'penicillin binding protein 2a PBP2a','PBP':'penicillin binding protein bacterial','70S_ribosome':'bacterial 70S ribosome antibiotic','D-Ala-D-Ala':'D-Ala-D-Ala ligase bacterial',
'Beta-lactamase':'bacterial beta lactamase','30S_ribosome':'bacterial 30S ribosome antibiotic','50S_ribosome':'bacterial 50S ribosome antibiotic'}
SOLVENT={'HOH','DOD','WAT','SOL','SO4','PO4','CL','NA','K','MG','CA','GOL','EDO','PEG','DMS','ACT','FMT','MES','TRS','MPD'}

def search(term):
    payload={'query':{'type':'terminal','service':'full_text','parameters':{'value':term}},'return_type':'entry','request_options':{'paginate':{'start':0,'rows':20},'results_content_type':['experimental']}}
    r=requests.post('https://search.rcsb.org/rcsbsearch/v2/query',json=payload,timeout=10)
    r.raise_for_status(); return [x.get('identifier') for x in r.json().get('result_set',[]) if x.get('identifier')]

def entry(pdb):
    r=requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}',timeout=8); r.raise_for_status(); return r.json()

def main():
    rows=[]
    out_path=OUT/'rcsb_structure_candidates_v2.csv'
    for cls,term in TERMS.items():
        ids=[]; error=''
        try: ids=search(term)
        except Exception as e: error=repr(e)
        for pdb in ids[:2]:
            row={'target_class':cls,'search_term':term,'pdb_id':pdb,'source_url':f'https://www.rcsb.org/structure/{pdb}','search_status':'text_search_candidate','search_error':error}
            try:
                d=entry(pdb); info=d.get('rcsb_entry_info',{}); lig_raw=info.get('nonpolymer_bound_components') or []
                lig=[]
                for x in lig_raw:
                    if isinstance(x,dict): lig.append(str(x.get('comp_id') or x.get('name') or x.get('id') or ''))
                    else: lig.append(str(x))
                lig=[x for x in lig if x.upper() not in SOLVENT]
                exptl=d.get('exptl') or []
                methods=[]
                for x in exptl if isinstance(exptl,list) else [exptl]:
                    if isinstance(x,dict) and x.get('method'): methods.append(str(x.get('method')))
                row.update({'experimental_method':';'.join(methods),
                            'resolution_A':info.get('resolution_combined',[None])[0] if info.get('resolution_combined') else None,
                            'nonpolymer_ligands':';'.join(lig),'co_crystal_ligand_present':bool(lig),
                            'active_site_evidence':'co-crystallized non-solvent ligand reported by RCSB entry' if lig else 'no non-solvent ligand reported in entry metadata',
                            'structure_metadata_status':'retrieved'})
            except Exception as e: row.update({'structure_metadata_status':'not_retrieved','structure_error':repr(e),'co_crystal_ligand_present':False})
            rows.append(row); time.sleep(.05)
        # checkpoint after every target class so a remote timeout never loses all progress
        pd.DataFrame(rows).to_csv(out_path,index=False)
    df=pd.DataFrame(rows)
    summ=df.groupby('target_class').agg(n_search_candidates=('pdb_id','nunique'),n_with_co_crystal_ligand=('co_crystal_ligand_present','sum')).reset_index() if len(df) else pd.DataFrame(columns=['target_class','n_search_candidates','n_with_co_crystal_ligand'])
    summ.to_csv(OUT/'rcsb_structure_summary_v2.csv',index=False)
    print(summ.to_string(index=False))

if __name__=='__main__': main()
