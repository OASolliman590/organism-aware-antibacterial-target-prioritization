"""Fetch public PubChem structures for the ESKAPE benchmark drugs."""
from pathlib import Path
import csv, json, time, argparse
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'benchmark'
OUT.mkdir(parents=True, exist_ok=True)

# Target labels are mechanism-level and intentionally conservative.
# Reference URLs are retained for every row for auditability.
rows = [
    ('ciprofloxacin','Gyrase/TopoIV','gyrase/topoisomerase inhibitor','Klebsiella pneumoniae;Acinetobacter baumannii;Escherichia coli;Proteus mirabilis;MRSA / Staphylococcus aureus'),
    ('levofloxacin','Gyrase/TopoIV','gyrase/topoisomerase inhibitor','Klebsiella pneumoniae;Acinetobacter baumannii;Escherichia coli;Proteus mirabilis;MRSA / Staphylococcus aureus'),
    ('moxifloxacin','Gyrase/TopoIV','gyrase/topoisomerase inhibitor','Klebsiella pneumoniae;Escherichia coli;Proteus mirabilis;MRSA / Staphylococcus aureus'),
    ('trimethoprim','DHFR','folate-pathway inhibitor','Klebsiella pneumoniae;Escherichia coli;Proteus mirabilis;MRSA / Staphylococcus aureus'),
    ('sulfamethoxazole','DHPS','folate-pathway inhibitor','Klebsiella pneumoniae;Escherichia coli;Proteus mirabilis;MRSA / Staphylococcus aureus'),
    ('fosfomycin','MurA','peptidoglycan precursor inhibitor','Escherichia coli;Klebsiella pneumoniae;Proteus mirabilis;MRSA / Staphylococcus aureus'),
    ('triclosan','FabI','enoyl-ACP reductase inhibitor','Escherichia coli;Klebsiella pneumoniae;Proteus mirabilis;Acinetobacter baumannii'),
    ('ceftaroline','PBP','PBP/transpeptidase inhibitor','MRSA / Staphylococcus aureus'),
    ('rifampicin','RpoB','RNA polymerase inhibitor','MRSA / Staphylococcus aureus;Klebsiella pneumoniae;Escherichia coli'),
    ('linezolid','Ribosome','protein-synthesis inhibitor','MRSA / Staphylococcus aureus;Enterococcus faecium'),
    ('vancomycin','D-Ala-D-Ala / cell wall','cell-wall precursor binder','MRSA / Staphylococcus aureus;Enterococcus faecium'),
    ('daptomycin','membrane','membrane-active antibiotic','MRSA / Staphylococcus aureus;Enterococcus faecium'),
    ('doxycycline','Ribosome','protein-synthesis inhibitor','Klebsiella pneumoniae;Escherichia coli;Proteus mirabilis;Acinetobacter baumannii;MRSA / Staphylococcus aureus'),
    ('tigecycline','Ribosome','protein-synthesis inhibitor','Klebsiella pneumoniae;Acinetobacter baumannii;MRSA / Staphylococcus aureus'),
    ('meropenem','PBP','beta-lactam cell-wall inhibitor','Klebsiella pneumoniae;Acinetobacter baumannii;Escherichia coli;Proteus mirabilis'),
    ('colistin','lipid A / membrane','membrane-active antibiotic','Klebsiella pneumoniae;Acinetobacter baumannii;Escherichia coli;Proteus mirabilis'),
]

parser=argparse.ArgumentParser()
parser.add_argument('--refresh',action='store_true',help='refresh the public PubChem structure cache')
args=parser.parse_args()
if (OUT/'eskape_benchmark_drugs.csv').exists() and not args.refresh:
    print(f'Using cached benchmark structures in {OUT}')
    raise SystemExit(0)

out=[]
for name, target, mechanism, organisms in rows:
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON'
    r = requests.get(url, timeout=45)
    r.raise_for_status()
    props = r.json()['PropertyTable']['Properties'][0]
    out.append({
        'drug': name,
        'target_class': target,
        'mechanism_class': mechanism,
        'organisms': organisms,
        'canonical_smiles': props.get('ConnectivitySMILES') or props.get('CanonicalSMILES'),
        'isomeric_smiles': props.get('SMILES') or props.get('IsomericSMILES'),
        'inchikey': props.get('InChIKey',''),
        'pubchem_url': url,
        'target_label_source': 'Mechanism-level curated label; see benchmark README and primary label sources',
    })
    time.sleep(0.25)

with open(OUT/'eskape_benchmark_drugs.csv','w',newline='') as f:
    w=csv.DictWriter(f, fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
with open(OUT/'eskape_benchmark_sources.json','w') as f:
    json.dump({'description':'Public PubChem structures with conservative mechanism-level target labels for leakage-aware benchmarking.','rows':out},f,indent=2)
print(f'Wrote {len(out)} benchmark drugs to {OUT}')
