from pathlib import Path
import re
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SNP=ROOT/'data'/'card'/'data'/'snps.txt'
OUT=ROOT/'data'/'resistance_v2'; OUT.mkdir(parents=True,exist_ok=True)

def family(text):
    t=text.lower()
    if re.search(r'mec[-_]?a|pbp2a',t): return 'PBP2a'
    if re.search(r'gyr[a-b]|par[c-e]|topoisomerase|fluoroquinolone',t): return 'GyrB/TopoIV resistance'
    if re.search(r'rpo[bcd]|rifamp',t): return 'RpoB resistance'
    if re.search(r'fol[a-p]|trimethoprim|sulfonamide',t): return 'DHFR/DHPS resistance'
    if re.search(r'rps|rrs|23s|ribosom|tetracycline|macrolide',t): return 'Ribosome resistance'
    if re.search(r'van[a-z]|glycopeptide|d-ala',t): return 'D-Ala-D-Ala resistance'
    if re.search(r'porin|omp[a-z]|efflux|acr|mex|ade|oqx',t): return 'Porin/Efflux resistance'
    if re.search(r'lpx|mcr|lipid a|colistin',t): return 'Lipid-A/envelope resistance'
    if re.search(r'fos[a-z]|fosfomycin',t): return 'MurA-pathway resistance'
    return 'Other resistance determinant'

def main():
    rows=[]
    for line in SNP.read_text(errors='ignore').splitlines():
        if not line or line.startswith('#'): continue
        parts=line.rstrip('\n').split('\t')
        if len(parts)<8: continue
        model_id,description,model_type,variant_type,variant,model_name,confidence,pmid=parts[:8]
        text=description+' '+model_name
        organism='unresolved'
        for candidate in ['Acinetobacter baumannii','Bacillus subtilis','Bacillus cereus','Escherichia coli','Klebsiella pneumoniae','Proteus mirabilis','Staphylococcus aureus','Staphylococcus epidermidis','Pseudomonas aeruginosa']:
            if candidate.lower() in text.lower(): organism=candidate; break
        rows.append({'model_id':model_id,'description':description,'model_type':model_type,'variant_type':variant_type,'variant':variant,
                     'model_name':model_name,'confidence':confidence,'pmid':pmid,'resistance_family':family(text),'organism':organism,
                     'source':'CARD current snps.txt','source_url':'https://card.mcmaster.ca/latest/data'})
    df=pd.DataFrame(rows); df.to_csv(OUT/'card_snp_annotations_v2.csv',index=False)
    summary=df.groupby('resistance_family').agg(n_snp_rows=('variant','count'),n_models=('model_id','nunique'),n_organism_descriptions=('description','nunique')).reset_index()
    summary.to_csv(OUT/'card_snp_family_summary_v2.csv',index=False)
    by_org=df.groupby(['organism','resistance_family']).agg(n_snp_rows=('variant','count'),n_models=('model_id','nunique')).reset_index()
    by_org.to_csv(OUT/'card_snp_organism_family_summary_v2.csv',index=False)
    print(summary.sort_values('n_snp_rows',ascending=False).to_string(index=False))

if __name__=='__main__': main()
