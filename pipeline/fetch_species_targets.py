"""Fetch species-specific bacterial target proteins from UniProt REST in batches.

The output is an auditable mapping layer, not a claim of target engagement. The
mapper uses canonical gene symbols first and protein-name fallback terms second,
because bacterial resistance alleles and fragmented annotations often omit the
canonical gene symbol.
"""
import csv, io, json, time
import requests
from requests.exceptions import RequestException

try:
    from pipeline.config import load_config
    from pipeline.provenance import utc_now
    from pipeline.snapshots import require_refresh_output, refresh_snapshot_root, write_text_exclusive
except ModuleNotFoundError:  # direct ``python pipeline/<script>.py`` execution
    from config import load_config
    from provenance import utc_now
    from snapshots import require_refresh_output, refresh_snapshot_root, write_text_exclusive


CONFIG = load_config()
CSV_OUT = CONFIG.path_for("species_proteins")
FASTA_OUT = CONFIG.path_for("species_fasta")
METADATA_OUT = CONFIG.path_for("species_metadata")
UNIPROT_CONFIG = CONFIG.value("refresh.uniprot")

ORGANISMS = {
    "Bacillus cereus": 1396,
    "Klebsiella pneumoniae": 573,
    "Escherichia coli": 562,
    "Proteus mirabilis": 584,
    "Acinetobacter baumannii": 470,
    "Staphylococcus aureus": 1280,
}

TARGET_GENES = {
    "DHFR": "folA", "FabI": "fabI", "FabH": "fabH", "FtsZ": "ftsZ",
    "GyrB": "gyrB", "TopoIV": "parE", "RpoB": "rpoB", "LeuRS": "leuS",
    "MurA": "murA", "MurC": "murC", "MurE": "murE", "LpxA": "lpxA",
    "LpxC": "lpxC", "LpxH": "lpxH", "PBP2a": "mecA",
    "PBP1A": "mrcA", "PBP1B": "mrcB", "PBP2": "pbpA", "PBP2B": "pbp2b",
    "PBP2X": "pbp2x", "PBP3": "ftsI", "PBP4": "dacB",
}
GENE_ALIASES = {
    ("Staphylococcus aureus", "FabI"): ["fabI", "fabL"],
    ("Bacillus cereus", "FabI"): ["fabL", "fabI"],
}
PROTEIN_TERMS = {
    "DHFR": ["dihydrofolate reductase"],
    "FabI": ["enoyl-ACP reductase", "enoyl-[acyl-carrier-protein] reductase"],
    "FabH": ["3-oxoacyl-ACP synthase III", "beta-ketoacyl-ACP synthase III"],
    "FtsZ": ["cell division protein FtsZ"],
    "GyrB": ["DNA gyrase subunit B"],
    "TopoIV": ["DNA topoisomerase IV subunit B"],
    "RpoB": ["DNA-directed RNA polymerase subunit beta", "RNA polymerase beta subunit"],
    "LeuRS": ["leucine--tRNA ligase", "leucyl-tRNA synthetase"],
    "MurA": ["UDP-N-acetylglucosamine enolpyruvyl transferase", "UDP-N-acetylglucosamine 1-carboxyvinyltransferase"],
    "MurC": ["UDP-N-acetylmuramate--L-alanine ligase"],
    "MurE": ["UDP-N-acetylmuramoyl-L-alanine--D-glutamate ligase"],
    "LpxA": ["UDP-N-acetylglucosamine acyltransferase"],
    "LpxC": ["UDP-3-O-[3-hydroxymyristoyl]glucosamine N-acyltransferase"],
    "LpxH": ["UDP-2,3-diacylglucosamine hydrolase"],
    "PBP2a": ["penicillin-binding protein 2a", "PBP2a"],
    "PBP1A": ["penicillin-binding protein 1A", "PBP1A"],
    "PBP1B": ["penicillin-binding protein 1B", "PBP1B"],
    "PBP2": ["penicillin-binding protein 2", "PBP2"],
    "PBP2B": ["penicillin-binding protein 2B", "PBP2B"],
    "PBP2X": ["penicillin-binding protein 2X", "PBP2X"],
    "PBP3": ["penicillin-binding protein 3", "PBP3"],
    "PBP4": ["penicillin-binding protein 4", "PBP4"],
}
FIELDS = "accession,id,protein_name,organism_name,organism_id,gene_names,length,sequence,reviewed,annotation_score"


def blank_row(organism, taxid, target, gene, status="not_found"):
    return {"organism": organism, "taxid": taxid, "target_class": target,
            "gene_query": gene, "accession": "", "protein_id": "", "protein_name": "",
            "organism_name": "", "organism_id": "", "gene_names": "", "length": "",
            "sequence": "", "reviewed": "", "annotation_score": "", "mapping_method": "",
            "source_url": "", "status": status}


def fetch_batch(organism, taxid):
    genes = sorted(set(TARGET_GENES.values()) | {a for k, aliases in GENE_ALIASES.items() if k[0] == organism for a in aliases})
    protein_terms = sorted(set(term for terms in PROTEIN_TERMS.values() for term in terms))
    pieces = [f"gene:{g}" for g in genes] + [f'protein_name:"{term}"' for term in protein_terms]
    query = f"organism_id:{taxid} AND (" + " OR ".join(pieces) + ")"
    url = "https://rest.uniprot.org/uniprotkb/search"
    try:
        r = requests.get(url, params={"query": query, "format": "tsv", "fields": FIELDS, "size": 500}, timeout=(10, int(UNIPROT_CONFIG["request_timeout_seconds"])))
        r.raise_for_status()
    except RequestException as e:
        print(f"WARN {organism}: {e}")
        return [], "", "retrieval_failed"
    lines = r.text.splitlines()
    if len(lines) <= 1:
        return [], r.url, "empty_response"
    header = lines[0].split("\t")
    records = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]
    return records, r.url, "retrieved"


def choose_record(records, target, genes):
    exact=[]
    for rec in records:
        tokens = set((rec.get("Gene Names", "") or "").replace(";", " ").split())
        if any(g in tokens for g in genes):
            exact.append(("gene", rec))
    if exact:
        candidates=[r for _,r in exact]
        candidates.sort(key=lambda x: (x.get("Reviewed", "") != "reviewed", -int(x.get("Length", "0") or 0), x.get("Entry", "")))
        return candidates[0], "gene"
    terms=[t.casefold() for t in PROTEIN_TERMS.get(target, [])]
    named=[r for r in records if any(t in (r.get("Protein names", "") or "").casefold() for t in terms)]
    if named:
        named.sort(key=lambda x: (x.get("Reviewed", "") != "reviewed", -int(x.get("Length", "0") or 0), x.get("Entry", "")))
        return named[0], "protein_name"
    return None, ""


def main():
    refresh_snapshot_root(CONFIG)
    if str(UNIPROT_CONFIG["source_release"]).lower() == "unrecorded":
        raise RuntimeError(
            "refresh.uniprot.source_release must name the UniProt release before refresh"
        )
    for path in (CSV_OUT, FASTA_OUT, METADATA_OUT):
        require_refresh_output(CONFIG, path)

    queried_at = utc_now()
    rows=[]
    retrieval_statuses = {}
    for organism, taxid in ORGANISMS.items():
        records, source_url, retrieval_status = fetch_batch(organism, taxid)
        retrieval_statuses[organism] = retrieval_status
        print(f"{organism}: {len(records)} candidate records")
        for target, primary_gene in TARGET_GENES.items():
            genes = GENE_ALIASES.get((organism, target), [primary_gene])
            rec, method = choose_record(records, target, genes)
            if rec is None:
                missing_status = (
                    "not_found" if retrieval_status == "retrieved" else retrieval_status
                )
                rows.append(blank_row(organism, taxid, target, "/".join(genes), missing_status))
            else:
                rows.append({"organism": organism, "taxid": taxid, "target_class": target,
                             "gene_query": "/".join(genes), "accession": rec.get("Entry", ""),
                             "protein_id": rec.get("Entry Name", ""), "protein_name": rec.get("Protein names", ""),
                             "organism_name": rec.get("Organism", ""), "organism_id": rec.get("Organism (ID)", ""),
                             "gene_names": rec.get("Gene Names", ""), "length": rec.get("Length", ""),
                             "sequence": rec.get("Sequence", ""), "reviewed": rec.get("Reviewed", ""),
                             "annotation_score": rec.get("Annotation", ""), "mapping_method": method,
                             "source_url": source_url, "status": "mapped"})
        time.sleep(0.25)

    csv_buffer = io.StringIO(newline="")
    writer = csv.DictWriter(csv_buffer, fieldnames=rows[0].keys())
    writer.writeheader(); writer.writerows(rows)
    write_text_exclusive(CONFIG, CSV_OUT, csv_buffer.getvalue())
    fasta = "".join(
        f">{row['organism'].replace(' ', '_')}|{row['target_class']}|{row['accession']}\n{row['sequence']}\n"
        for row in rows
        if row["sequence"]
    )
    write_text_exclusive(CONFIG, FASTA_OUT, fasta)
    write_text_exclusive(
        CONFIG,
        METADATA_OUT,
        json.dumps(
            {
                "organisms": ORGANISMS,
                "target_genes": TARGET_GENES,
                "protein_name_terms": PROTEIN_TERMS,
                "n_rows": len(rows),
                "n_mapped": sum(r['status']=='mapped' for r in rows),
                "source": "UniProt REST API",
                "source_release": UNIPROT_CONFIG["source_release"],
                "queried_at_utc": queried_at,
                "retrieval_status_by_organism": retrieval_statuses,
            },
            indent=2,
        ) + "\n",
    )
    print(f"Wrote {len(rows)} species-target rows; mapped {sum(r['status']=='mapped' for r in rows)}")


if __name__ == "__main__":
    main()
