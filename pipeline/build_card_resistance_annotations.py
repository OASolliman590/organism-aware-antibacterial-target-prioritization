"""Parse the current public CARD release into target-family resistance annotations.

The mapping is intentionally conservative and mechanism-level. CARD is a resistance
resource, not a direct assay for the user's compounds; outputs are therefore used as
resistance-risk/context annotations, never as binding evidence.
"""
from pathlib import Path
import json, re
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
CARD=ROOT/'data'/'card'/'data'/'card.json'
OUT=ROOT/'data'/'resistance_v2'; OUT.mkdir(parents=True,exist_ok=True)


def text_of(m):
    cats=m.get('ARO_category',{})
    bits=[m.get('ARO_name',''),m.get('ARO_description',''),m.get('model_name',''),m.get('model_description','')]
    for c in cats.values() if isinstance(cats,dict) else []:
        bits += [c.get('category_aro_name',''),c.get('category_aro_description',''),c.get('category_aro_class_name','')]
    return ' '.join(str(x) for x in bits).lower()

def family(text):
    # Order matters: specific target-alteration and major enzyme families first.
    if re.search(r'mec[-_]?a|pennicillin.binding protein 2a|pbp2a',text): return 'PBP2a'
    if re.search(r'beta[- ]?lactamase|bla[a-z]|carbapenemase|extended[- ]spectrum beta',text): return 'Beta-lactamase'
    if re.search(r'gyr[a-b]|par[c-e]|quinolone resistance|topoisomerase',text): return 'GyrB/TopoIV resistance'
    if re.search(r'fol[a-p]|dihydrofolate reductase|trimethoprim resistance',text): return 'DHFR resistance'
    if re.search(r'folp|dihydropteroate synthase|sulfonamide resistance',text): return 'DHPS resistance'
    if re.search(r'rpo[bcd]|rifampicin resistance|rifamycin resistance',text): return 'RpoB resistance'
    if re.search(r'van[a-z]|d-ala-d-ala|glycopeptide resistance',text): return 'D-Ala-D-Ala resistance'
    if re.search(r'erm\(|23s rRNA|50s ribosom|30s ribosom|ribosomal protein|tetracycline resistance',text): return 'Ribosome resistance'
    if re.search(r'porin|omp[a-z]|outer membrane permeability',text): return 'Porin resistance'
    if re.search(r'efflux|acr[a-z]|mex[a-z]|ade[a-z]|oqx[a-z]|nor[a-z]',text): return 'Efflux resistance'
    if re.search(r'lpx|arn|pmr|mcr|lipid a|colistin resistance',text): return 'Lipid-A/envelope resistance'
    if re.search(r'fos[a-z]|fosfomycin resistance',text): return 'MurA-pathway resistance'
    return 'Other resistance determinant'


def main():
    if not CARD.exists(): raise FileNotFoundError(CARD)
    data=json.loads(CARD.read_text())
    rows=[]
    for key,m in data.items():
        if not isinstance(m,dict) or 'model_id' not in m: continue
        txt=text_of(m); cats=m.get('ARO_category',{})
        drug_classes=[]; mechanisms=[]
        for c in cats.values() if isinstance(cats,dict) else []:
            n=str(c.get('category_aro_name',''))
            k=str(c.get('category_aro_class_name',''))
            if k=='Drug Class': drug_classes.append(n)
            if k=='Resistance Mechanism': mechanisms.append(n)
        seqs=m.get('model_sequences',{}).get('sequence',{})
        taxa=[]
        for s in seqs.values() if isinstance(seqs,dict) else []:
            tax=s.get('NCBI_taxonomy',{}).get('NCBI_taxonomy_name','') if isinstance(s,dict) else ''
            if tax: taxa.append(tax)
        rows.append({'model_id':m.get('model_id'),'model_name':m.get('model_name'),'model_type':m.get('model_type'),
                     'aro_accession':m.get('ARO_accession'),'aro_name':m.get('ARO_name'),'target_resistance_family':family(txt),
                     'drug_classes':';'.join(sorted(set(drug_classes))),'resistance_mechanisms':';'.join(sorted(set(mechanisms))),
                     'reference_taxa':';'.join(sorted(set(taxa))),'n_model_sequences':len(seqs),
                     'has_snp_data':bool(m.get('model_param',{}).get('snp_data') or m.get('model_sequences',{}).get('snp_data') or m.get('snp_data')),
                     'source':'CARD latest data archive','source_url':'https://card.mcmaster.ca/latest/data','card_version':data.get('_version','')})
    df=pd.DataFrame(rows); df.to_csv(OUT/'card_model_annotations_v2.csv',index=False)
    summary=df.groupby('target_resistance_family').agg(n_models=('model_id','count'),n_aro_terms=('aro_accession','nunique'),n_with_snp_data=('has_snp_data','sum')).reset_index()
    summary.to_csv(OUT/'card_resistance_family_summary_v2.csv',index=False)
    print(summary.sort_values('n_models',ascending=False).to_string(index=False))

if __name__=='__main__': main()
