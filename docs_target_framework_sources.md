# Target framework sources and design notes

## Sources consulted

1. World Health Organization. **WHO bacterial priority pathogens list, 2024**. The list covers 24 pathogens across 15 families and is intended to guide antibacterial R&D and public-health prioritization. https://www.who.int/publications/i/item/9789240093461
2. De Oliveira DMP et al. **Antimicrobial Resistance in ESKAPE Pathogens**. Clinical Microbiology Reviews. 2020;33(3):e00181-19. doi:10.1128/CMR.00181-19. https://pmc.ncbi.nlm.nih.gov/articles/PMC7227449/
3. ChEMBL database and web services. Public bioactivity and target annotations used for ligand-reference sets. https://www.ebi.ac.uk/chembl/
4. RDKit documentation. Molecular standardization, fingerprints, descriptors, and substructure operations. https://www.rdkit.org/

## Design decisions

The revised workflow separates **chemical target evidence** from **organism/clinical priority**. The first layer is open over all target classes for which reference-ligand evidence is available; it does not select from a user-defined panel. The second layer annotates each chemically plausible target with clinical validation, essentiality/fitness relevance, accessibility, resistance relevance, and organism compatibility.

Target classes with clinically validated antibacterial chemistry include DNA gyrase/topoisomerase (fluoroquinolones), DHFR/DHPS (trimethoprim/sulfonamides), MurA (fosfomycin), RpoB (rifamycins), ribosomal 30S/50S targets, PBPs/transpeptidases (beta-lactams), D-Ala-D-Ala/cell-wall precursors (glycopeptides), LeuRS (mupirocin, mainly staphylococcal topical use), membrane/lipid targets (daptomycin or polymyxins), and beta-lactamases as clinically important resistance targets. FabI, FtsZ, LpxA, LpxC, LpxH, MurC, MurE, and related enzymes are retained as chemically valid discovery targets but are annotated as preclinical or limited clinical-translation classes rather than falsely labelled clinically validated.

The organism filter must distinguish direct growth targets from resistance targets. A beta-lactamase can be a high clinical-priority target for a resistant isolate without being essential for bacterial viability. PBP2a is a high-priority target for MRSA but should not receive the same organism priority for Acinetobacter baumannii. Gram-negative outer-membrane or envelope targets require accessibility and permeability annotations rather than an assumption that biochemical binding implies cellular activity.

The ESKAPE benchmark uses known antibacterial drugs with conservative mechanism-level labels. Evaluation must be leakage-aware: the query drug and close analogues are excluded from reference evidence for that fold. Performance is reported only for target classes with sufficient positive examples; one-example target labels are reported as coverage limitations, not as failed predictions.


## v2 public sequence and resistance resources

5. UniProt Proteomes/API documentation: https://www.uniprot.org/api-documentation/proteomes. UniProt provides searchable proteome records and programmatic access suitable for retrieving bacterial reference proteomes and protein features; the v2 mapping layer will record the organism identifier, protein accession, sequence, and provenance.
6. Comprehensive Antibiotic Resistance Database (CARD): https://card.mcmaster.ca/ and primary citation https://academic.oup.com/nar/article/51/D1/D690/6764414. CARD organizes curated resistance determinants and associated antibiotics through the Antibiotic Resistance Ontology and provides reference sequences, mutation records, and AMR detection models. The v2 workflow will use CARD-derived annotations as resistance evidence, not as direct proof of compound binding.
7. RCSB Protein Data Bank antimicrobial-resistance resources: https://pdb101.rcsb.org/browse/antimicrobial-resistance. PDB structures will be used only to annotate target-pocket availability and structural precedent; docking and MD execution remain separate downstream analyses.

## V2 resistance and structure sources

CARD automated downloads were retrieved from https://card.mcmaster.ca/latest/data and https://card.mcmaster.ca/latest/ontology on 2026-08-17. The current CARD page reports release 4.0.2 (2026-08-12) and provides ARO ontology, model metadata, protein sequences, and `snps.txt`. The local parser uses `card.json` for model-level resistance-family context and `snps.txt` for explicit resistance variants; these are contextual resistance annotations, not direct evidence that a project compound binds a target.

RCSB structure candidates were queried through `https://search.rcsb.org/rcsbsearch/v2/query` and entry metadata through `https://data.rcsb.org/rest/v1/core/entry/{pdb_id}`. The structure catalog stores text-search candidates, metadata retrieval status, non-solvent bound components, and a co-crystallized-ligand active-site flag. Text-search candidates are not treated as validated species-specific structures; user-controlled docking remains a separate downstream activity.
