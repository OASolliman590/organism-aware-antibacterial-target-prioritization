"""Refresh public PubChem benchmark structures into a new dated snapshot."""
import csv, io, json, time
import requests

try:
    from pipeline.config import load_config
    from pipeline.provenance import utc_now
    from pipeline.snapshots import require_refresh_output, write_text_exclusive
except ModuleNotFoundError:  # direct ``python pipeline/<script>.py`` execution
    from config import load_config
    from provenance import utc_now
    from snapshots import require_refresh_output, write_text_exclusive


CONFIG = load_config()
CSV_OUT = CONFIG.path_for("benchmark")
JSON_OUT = CSV_OUT.with_name("eskape_benchmark_sources.json")
TIMEOUT = int(CONFIG.value("refresh.pubchem.request_timeout_seconds"))

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

def main():
    # All destinations are validated before the first network request.
    require_refresh_output(CONFIG, CSV_OUT)
    require_refresh_output(CONFIG, JSON_OUT)
    queried_at = utc_now()
    out=[]
    for name, target, mechanism, organisms in rows:
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON'
        response = requests.get(url, timeout=TIMEOUT)
        response.raise_for_status()
        props = response.json()['PropertyTable']['Properties'][0]
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

    csv_buffer = io.StringIO(newline="")
    writer=csv.DictWriter(csv_buffer, fieldnames=out[0].keys())
    writer.writeheader(); writer.writerows(out)
    write_text_exclusive(CONFIG, CSV_OUT, csv_buffer.getvalue())
    write_text_exclusive(
        CONFIG,
        JSON_OUT,
        json.dumps(
            {
                'description':'Public PubChem structures with conservative mechanism-level target labels for leakage-aware benchmarking.',
                'queried_at_utc': queried_at,
                'source_version': None,
                'rows':out,
            },
            indent=2,
        ) + "\n",
    )
    print(f'Wrote {len(out)} benchmark drugs to {CSV_OUT.parent}')


if __name__ == "__main__":
    main()
