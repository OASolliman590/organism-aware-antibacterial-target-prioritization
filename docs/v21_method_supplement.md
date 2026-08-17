# v2.1 ChEMBL reference-universe expansion and dual target prioritization

## Scope

Version 2.1 extends the v2 open-target-discovery workflow with a mechanism- and binding-site-aware reference layer for clinically important antibacterial targets that were previously absent or too broad in the public ChEMBL reference universe. The expansion is deliberately conservative: a target subtype is added to the molecular reference set only when ChEMBL target metadata, assay metadata, and a valid molecular structure support the assignment.

The expanded ontology separates PBPs, Ambler beta-lactamase classes, and selected 30S/50S ribosomal binding-site classes. Parent-family labels remain available for reporting and benchmark acceptance, but chemical similarity is calculated at the narrowest supported subtype. Beta-lactamases are annotated as resistance-modifying hypotheses rather than assumed direct bactericidal targets.

## ChEMBL retrieval and provenance

`pipeline/fetch_chembl_reference_subtypes_v21.py` consumes `data/chembl_target_aliases_v21.json`. The configuration stores manual ChEMBL target-ID seeds and search terms for each subtype. The fetcher records target IDs, target names, organisms, target types, assay IDs, document IDs, activity types, values, units, pChEMBL values, confidence scores, validity comments, source release, and standardized structure identifiers.

The retrieval is cacheable under `data/reference_quality/chembl_v21_cache/`. `CHEMBL_V21_DISCOVER=0` restricts a run to the manually verified IDs; `CHEMBL_V21_SUBTYPES` limits the target subtypes; `CHEMBL_V21_OFFLINE=1` prevents network access after cache creation; and `CHEMBL_V21_MAX_PER_SUBTYPE` bounds output size. The public pipeline must preserve an empty subtype as missing evidence rather than filling it with simulated or cross-family records.

## Activity-quality grading

The fetcher retains numeric potency records with pChEMBL ≥ 5.0 by default. Each record receives an assay-quality grade:

| Grade | Interpretation | Molecular similarity use |
|---|---|---|
| A | Direct single-protein or direct RNA/complex assignment, high target confidence, interpretable biochemical or binding assay | Full weight |
| B | Direct complex or homologous assignment with moderate-to-high confidence | Reduced weight |
| C | Protein-family assignment with clear family identity but uncertain subtype | Parent-family evidence only where explicitly allowed |
| D | Phenotypic or low-confidence assignment without defensible molecular resolution | Excluded from subtype similarity, retained for audit |

Replicate activities are grouped by standardized InChIKey. The strongest retained record is accompanied by the number of assays, the number of target IDs, and the observed quality-grade set. The original ChEMBL structure and the standardized structure are both retained for provenance.

## Structure standardization

RDKit removes explicit hydrogens and selects a parent fragment while preserving stereochemistry, formal charge, beta-lactam rings, covalent warheads, and tautomer-sensitive ribosomal pharmacophores. Structures that fail sanitization or InChIKey generation are excluded from molecular similarity and recorded as failed curation. Distinct stereoisomers are not collapsed. The canonical isomeric SMILES and standardized InChIKey are stored alongside the original ChEMBL representation.

## Species transfer

The UniProt mapper now includes PBP orthologue queries where canonical gene symbols are available, including `mrcA`, `mrcB`, `pbpA`, `pbp2b`, `pbp2x`, `ftsI`, and `dacB`. If a subtype has no organism-specific mapping but its parent class has a mapping, the scorer uses the parent mapping and records `organism_transfer_source=parent_class`. This is a transfer prior, not evidence that the compound binds the organism protein.

## Parallel scoring outputs

The v2.1 scorer retains the original overall organism-aware priority but adds two explicit fields:

`chemical_hypothesis_score` is the reference-quality-adjusted chemical evidence for the target subtype, including ECFP4/MACCS similarity and cross-target specificity.

`clinical_translation_score` is a separate prior combining clinical precedent, organism scope, essentiality, cellular accessibility, pocket/co-crystal context, and resistance context. For beta-lactamase subtypes, CARD resistance context contributes positively because the target is explicitly a resistance-modifying hypothesis. For direct antibacterial targets, resistance burden is retained as a translational caution.

The outputs therefore distinguish chemical compatibility from clinical development relevance. The pipeline does not interpret either score as a binding probability, MIC prediction, target-engagement measurement, or evidence that an unpublished compound is an antibiotic drug.

## Benchmark leakage control

The ESKAPE benchmark continues to exclude exact and close analogues at ECFP4 Tanimoto similarity ≥0.85 and evaluates an exact Bemis–Murcko scaffold split separately. Subtype references are evaluated under the same exclusion rules as parent-class references. Parent-level benchmark labels accept the corresponding v2.1 subtype labels, but subtype retrieval is reported separately when sufficient subtype reference coverage exists.

## Current release limitation

During the v2.1 implementation run, the public ChEMBL REST service returned timeouts and HTTP 500 responses for target/activity requests. The new ontology, retrieval code, cache controls, curation rules, scoring integration, benchmark acceptance mapping, figures, and documentation were implemented and tested. No missing PBP, beta-lactamase, or ribosomal subtype records were fabricated. Consequently, the current numerical benchmark remains the v2 reference benchmark until the ChEMBL refresh completes successfully.

## Recommended refresh acceptance criteria

A successful refresh should produce at least ten quality A/B ligands and at least five unique scaffolds for a subtype before it is treated as a stable subtype-specific ranking source. The benchmark should then be rerun, reporting coverage, retrieval, MRR, enrichment, parent-family retrieval, and site/mechanism-subtype retrieval separately. Any subtype that remains below the threshold should remain visible as low-confidence or parent-family-only evidence.

## References

[1]: https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services ChEMBL Data Web Services documentation.

[2]: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions ChEMBL activity and assay confidence guidance.

[3]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5131749/ Silver, Appropriate Targets for Antibacterial Drugs.
