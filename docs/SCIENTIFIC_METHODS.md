# Scientific Methods Reporting

## Organism-aware antibacterial target prioritization from small-molecule chemical evidence

### Purpose and reporting scope

This document provides a manuscript-ready description of the computational methodology implemented in this repository. It is intended for use in a Methods section, computational supplement, preregistration, thesis chapter, or reproducibility appendix.

The framework performs **open target discovery** from small molecules by comparing each query compound with a broad public reference-ligand universe and subsequently adding organism-specific biological, resistance, structural, and translational annotations. The method is designed to generate **ranked, experimentally testable target hypotheses**. It does **not** establish direct binding, target engagement, antibacterial efficacy, or clinical activity. Similarity-derived scores and heuristic prioritization scores are therefore reported as ranking or decision-support quantities rather than probabilities unless a field is explicitly produced by a held-out calibration procedure.

The current configured organism panel comprises:

- *Klebsiella pneumoniae*
- *Bacillus cereus*
- *Escherichia coli*
- *Proteus mirabilis*
- *Acinetobacter baumannii*
- *Staphylococcus aureus*

This should be described as a six-species antibacterial panel rather than as the canonical ESKAPE set.

---

## 1. Study design

The computational workflow separates evidence into two major stages.

First, a **chemical target-discovery stage** scores each query compound against all target classes for which public reference-ligand evidence is available. No user-selected docking panel is used to restrict the initial search. Chemical evidence is derived from two-dimensional molecular similarity, three-dimensional conformer similarity, alignment-based molecular shape and feature overlap, pharmacophore similarity, target specificity relative to cross-target controls, and the quality of the public reference set.

Second, an **organism-aware prioritization stage** annotates each chemically supported target with species-transfer evidence, organism scope, essentiality, cellular accessibility, clinical precedent, resistance context, structure/pocket precedent, and annotation-only anti-target risk. The chemical and biological layers are retained as separate output fields before an overall prioritization score is calculated.

Missing measurements are retained as missing or explicitly unavailable wherever the v3 workflow supports that state. They are not replaced by synthetic values, averages, or fabricated evidence.

---

## 2. Reproducibility, configuration, and provenance

The pipeline is executed from the repository root using:

```bash
python run_pipeline.py --config config.yaml
```

`config.yaml` is the declarative source for scoring parameters, thresholds, paths, analysis modes, and random seeds. The primary global seed is `20260817`. Separate seeds are used for conformer generation (`20240601`), bootstrap procedures (`20260817`), model fitting (`20240601`), figure sampling (`7`), and UMAP (`17`).

The default configuration uses a frozen public-data snapshot (`v2-public-baseline-2ed4684`) and sets `run.refresh_external_data: false`. The runner verifies the snapshot before analysis. Live API refresh is therefore not part of the default analytical execution.

For each run, the workflow writes `results/run_manifest.json` and an immutable copy under `results/run_manifests/`. The manifest records the configuration SHA-256 hash, snapshot metadata, random seeds, Git commit and branch when available, tracked working-tree status, Python version and executable, installed package versions, command-line invocation, and start/completion timestamps.

The pinned public snapshot records the following provenance state:

- ChEMBL reference ligands: source release and acquisition date were not retained in the legacy v2 files and remain unknown.
- PubChem benchmark: query date was not retained.
- UniProt species-target data: release and query date were not retained.
- CARD-derived annotations: CARD release 4.0.2 is recorded.
- RCSB PDB structure catalog: query date was not retained.
- Repository-curated target ontologies are versioned by Git.

These missing provenance fields must be reported as unknown rather than inferred from Git history.

The locked computational environment includes Python 3.12 and repository-pinned packages including RDKit 2026.3.5, NumPy 2.5.2, pandas 3.0.5, SciPy 1.18.0, scikit-learn 1.9.0, matplotlib 3.11.1, and UMAP-learn 0.5.12.

---

## 3. Target universe and biological ontology

The search space is defined by repository-curated target ontologies rather than by a compound-specific target list. The base ontology includes direct viability targets, protein complexes, cell-wall precursor targets, resistance determinants, phenotypic mechanisms, permeability/resistance mechanisms, and a safety anti-target class. A v2.1 subtype ontology adds mechanism- or site-resolved subclasses including PBP, beta-lactamase, and ribosomal subtypes where public reference evidence supports that resolution.

Examples of represented target families include folate metabolism, fatty-acid synthesis, DNA gyrase/topoisomerase, RNA polymerase, aminoacyl-tRNA synthetases, peptidoglycan biosynthesis, lipid-A biosynthesis, PBPs, beta-lactamases, ribosomal targets, membrane-associated mechanisms, efflux, and porin-related resistance context.

Target annotations include target role, mechanism granularity, clinical status, essentiality, cellular localization, resistance relevance, organism scope, benchmark aliases, minimum reference-ligand requirements, and requirements for sequence or pocket mapping.

Target classes lacking usable reference ligands are retained in the ontology and are not silently removed from the conceptual target universe. In the v3 output path, unavailable chemical components are represented as unavailable evidence rather than as a positive signal.

---

## 4. Public reference-ligand data and reference quality

Public target-associated ligands are stored under `data/reference_ligands/` and are derived from ChEMBL records. During loading, structures are parsed with RDKit, explicit hydrogens are removed for comparison, canonical SMILES are generated, and duplicate structures within a target class are removed on canonical SMILES.

Reference quality is summarized per target class using:

- number of valid ligands;
- number and fraction of unique Bemis-Murcko scaffolds;
- number of represented ChEMBL targets;
- number of records with pChEMBL values;
- median pChEMBL where available;
- assay-type distribution; and
- source-organism distribution.

Reference quality is classified using target-specific minimum ligand counts and scaffold diversity. In the current configuration, targets below their minimum reference count are classified as `low`, targets with no valid ligands as `insufficient`, targets meeting count requirements but with limited scaffold diversity as `moderate_redundancy`, and sufficiently diverse sets as `usable`.

Reference quality is preserved as an independent evidence field and contributes multiplicatively to the quality-adjusted chemical score.

---

## 5. Two-dimensional molecular representations and chemical evidence

Each sanitized molecule is represented using:

1. a Morgan/ECFP4 fingerprint generated with radius 2 and 2,048 bits; and
2. the standard RDKit MACCS structural keys.

Tanimoto similarity is used for both fingerprint types.

For each query-target class pair, the workflow calculates:

- maximum ECFP4 similarity to the target reference set;
- mean ECFP4 similarity among the five most similar target references (`top_k = 5`); and
- maximum MACCS similarity.

The legacy v2 chemical evidence score is calculated after parameterized linear normalization and clipping of each component to [0,1]. The configured normalized components are:

```text
ECFP4_max_normalized = clip((ECFP4_max - 0.10) / 0.55, 0, 1)
ECFP4_top5_normalized = clip((ECFP4_top5_mean - 0.08) / 0.45, 0, 1)
MACCS_max_normalized = clip((MACCS_max - 0.10) / 0.70, 0, 1)
```

The configured weighted score is:

```text
chemical_evidence_score_v2
    = 0.50 * ECFP4_max_normalized
    + 0.25 * ECFP4_top5_normalized
    + 0.15 * MACCS_max_normalized
```

The implementation clips the final value to [0,1]. The listed weights sum to 0.90; the code does not renormalize them to 1.0, so this should be reported exactly as implemented.

---

## 6. Target-specificity control

To estimate whether a chemical similarity signal is specific to a target class rather than broadly shared across the reference universe, the maximum ECFP4 similarity to ligands from all other target classes is calculated as a cross-target control.

The target-specificity margin is:

```text
target_specificity_margin
    = best_target_ECFP4_similarity
    - best_cross_target_ECFP4_similarity
```

The configured specificity transform is:

```text
target_specificity_score
    = clip(0.50 + 0.50 * margin / 0.25, 0, 1)
```

Cross-target reference ligands are used only as **specificity controls**. They are not interpreted as experimentally inactive compounds against the query target.

The v2 quality-adjusted chemical score is:

```text
chemical_quality_adjusted_score
    = chemical_evidence_score_v2
    * reference_quality_score
    * target_specificity_score
```

---

## 7. Deterministic three-dimensional conformer generation

V3 adds deterministic conformer-based chemical evidence. For each sanitized molecule, RDKit ETKDGv3 is used with a fixed conformer seed. The default configuration requests 30 conformers, uses an RMS pruning threshold of 0.5 Å, and uses one RDKit thread per molecule to preserve deterministic behavior.

Hydrogens are added before embedding. When MMFF94 parameters are available, conformers are optimized with MMFF94 for up to 500 iterations and conformers with finite energies within 10 kcal/mol of the lowest-energy conformer are retained. If MMFF94 parameters are unavailable, the ensemble is retained without energy filtering and is explicitly marked accordingly. Failed embedding or failed finite MMFF energy calculation is stored as a failure state and does not generate synthetic 3D similarity.

Conformers are cached as RDKit binary molecules. The cache key is SHA-256 over the canonical isomeric SMILES and the complete conformer-generation parameter set, including the RDKit version. This prevents a cached conformer ensemble generated with different parameters or software versions from being reused silently.

The public configuration uses eight prewarming workers and eight independent scoring workers. Worker scheduling does not alter molecular seeds, component values, or row order.

---

## 8. Three-dimensional shape and pharmacophore evidence

### 8.1 USRCAT

For each query-reference conformer pair, USRCAT is used as an alignment-free three-dimensional shape-plus-feature measure. Reference-level scores are retained, and target-level summaries include the maximum score and a configured top-k aggregate.

### 8.2 O3A alignment and shape evidence

References within each target class are first ranked by USRCAT. Only the top 25 references enter the more expensive alignment-based stage.

RDKit O3A is applied using MMFF94 atom typing. After alignment, shape similarity is calculated as:

```text
shape_similarity = 1 - ShapeTanimotoDist
```

The O3A color/feature score is retained separately. Shape and color values are bounded to [0,1]. O3A failures are recorded as failures or missing evidence rather than converted to numeric similarity.

### 8.3 Two-dimensional pharmacophore similarity

An alignment-free pharmacophore field is calculated with the RDKit Gobbi/Poppe feature-pair fingerprint and Tanimoto similarity. This is stored as `pharmacophore_2d_gobbi_sim_max`.

### 8.4 Three-dimensional pharmacophore overlap

For successfully O3A-aligned conformer pairs, RDKit `ChemicalFeatures` is used with the pinned `BaseFeatures.fdef` feature definition. Feature overlap is scored with RDKit `FeatMaps` using:

- same-family feature matching;
- Gaussian profile;
- 2.5 Å cutoff radius;
- width 1.0;
- unit feature weights; and
- ignored feature directions.

For each directed comparison, the best same-family feature match is selected independently for each probe feature and normalized by the number of probe features. The final pair score is the mean of query-to-reference and reference-to-query directed coverage:

```text
pharmacophore_3d_similarity
    = mean(directed_coverage_query_to_reference,
           directed_coverage_reference_to_query)
```

A valid pair with no spatial overlap may receive zero. In contrast, absence of detectable features, feature-extraction failure, or failure to obtain a valid aligned comparison remains missing evidence.

---

## 9. V3 evidence fusion

V3 combines seven configured chemical components using normalized Reciprocal Rank Fusion (RRF):

1. `ecfp4_max`
2. `maccs_max`
3. `usrcat_max`
4. `o3a_shape_tanimoto_max`
5. `o3a_color_max`
6. `pharmacophore_2d_gobbi_sim_max`
7. `pharmacophore_3d_sim_max`

For each query, target classes are ranked separately within each component. Higher values rank first, and exact ties receive their average rank.

With RRF constant `k = 60`, the raw fused score is:

```text
RRF_raw(q,t) = sum_i [ I_i / (k + rank_i(q,t)) ]
```

where `I_i` indicates whether component `i` is available.

The normalized score is:

```text
chemical_evidence_score_v3
    = RRF_raw / (N_configured_components / (k + 1))
```

The fixed denominator uses all seven configured components. Consequently, missing measurements reduce achievable fused support and are not rewarded by denominator shrinkage. A row with no valid chemical components remains missing.

The v3 quality-adjusted chemical score is:

```text
chemical_quality_adjusted_score_v3
    = chemical_evidence_score_v3
    * reference_quality_score
    * target_specificity_score
```

Neither the v3 fused score nor its quality-adjusted form is interpreted as a binding probability.

---

## 10. Species-specific sequence transfer

Species-target transfer is assessed using mapped bacterial protein sequences. Within each target class, the mapped sequence selected as the reference is the longest sequence after deterministic sorting of reviewed status, sequence length, and accession.

A simple global Needleman-Wunsch alignment is performed with match = +1, mismatch = -1, and gap = -1. Pairwise identity is calculated over aligned non-gap residues, and coverage is calculated as the number of aligned residue pairs divided by the shorter input sequence length.

The species-transfer score is:

```text
species_transfer_score
    = identity * sqrt(coverage)
```

with the result bounded to [0,1]. Unmapped targets receive a transfer score of 0 and an explicit unmapped status. This score is used as a mapping/transfer heuristic and is not a binding probability.

---

## 11. Biological and organism-aware annotation

For every chemically scored target, the pipeline generates an organism-specific row for each configured species.

The following evidence fields are calculated or attached independently:

- organism-scope compatibility;
- clinical-priority score;
- essentiality score;
- cellular-accessibility score;
- resistance-relevance score;
- CARD resistance-context score;
- structure/pocket evidence;
- species-transfer score; and
- annotation-only anti-target risk.

### 11.1 Organism scope

The ontology specifies broad bacterial, Gram-positive, Gram-negative, MRSA-specific, or other target scope. Exact/group-compatible scope scores receive 1.0, broad bacterial scope 0.85, default unmatched scope 0.25, and PBP2a outside *S. aureus* receives 0.05 in the current configuration.

### 11.2 Clinical, essentiality, and accessibility mappings

Clinical evidence is mapped to configured ordinal values: approved = 1.0, clinical = 0.75, validated = 0.65, preclinical = 0.40, otherwise 0.20.

General ordinal biological annotations are mapped as high = 1.0, medium = 0.65, variable = 0.45, low = 0.30, default = 0.45.

Cellular accessibility is mapped as cytosolic = 0.70, periplasmic = 0.90, outer membrane = 0.55, membrane = 0.45, extracellular/periplasmic = 0.85, otherwise 0.50.

### 11.3 Resistance context

CARD-derived resistance evidence is used as contextual biological information only. Presence of a matching CARD model contributes 0.5 and presence of corresponding SNP rows contributes 0.5, bounded to [0,1]. CARD counts are not interpreted as prevalence and do not provide evidence that a project compound binds the target.

### 11.4 Structure and pocket precedent

RCSB PDB metadata are used to identify bounded structure candidates and whether co-crystallized non-solvent ligands are present. The configured pocket evidence is 1.0 when a co-crystallized ligand is present, 0.60 when a structure candidate is present without such a ligand, and 0.0 when no structure evidence is available.

This layer provides structural precedent only. The repository does not equate RCSB structure availability with compound-specific docking support.

### 11.5 Annotation-only anti-target risk

Selected target classes receive an early selectivity warning based on obvious human-homologue or mitochondrial considerations. These annotations are explicitly labelled `annotation_only` and are not toxicity predictions or experimentally measured off-target activities.

---

## 12. Biological priority and overall organism-aware score

The configured biological-priority score is:

```text
biological_priority_score
    = 0.25 * organism_scope
    + 0.20 * clinical_priority
    + 0.20 * essentiality
    + 0.15 * cellular_accessibility
    + 0.10 * resistance_relevance
    + 0.10 * CARD_context
```

The overall organism-aware priority score is then calculated as:

```text
overall_priority_score
    = chemical_quality
    * (0.65 + 0.35 * species_transfer)
    * (0.75 + 0.25 * pocket_evidence)
    * (0.50 + 0.50 * biological_priority)
    * (1 - 0.20 * anti_target_risk)
```

where `chemical_quality` is the configured quality-adjusted chemical score for the selected analysis path.

This formulation is a transparent heuristic prioritization function. It is not fitted to antibacterial outcomes and is not a probability of target engagement or efficacy.

A separate clinical-translational score is also generated from clinical status, organism scope, essentiality, accessibility, pocket evidence, and resistance/CARD context. This score is retained separately from the chemical hypothesis score.

---

## 13. Confidence classification and uncertainty reasons

The workflow assigns heuristic confidence classes after ranking.

High confidence requires:

- overall priority >= 0.50;
- chemical-quality score >= 0.45; and
- species-transfer score >= 0.70.

Moderate confidence requires:

- overall priority >= 0.25; and
- chemical-quality score >= 0.25.

Rows with weaker prioritization but chemical evidence >= 0.20 may be labelled Low; otherwise they are Insufficient.

Uncertainty reasons are recorded explicitly when reference coverage is limited, sequence mapping is unresolved or weak, specificity against cross-target controls is poor, anti-target risk is elevated, or no RCSB structural precedent is available.

These confidence labels are heuristic categories and must not be described as statistically calibrated confidence probabilities.

---

## 14. Applicability-domain assessment

Applicability domain is estimated primarily from the nearest-reference ECFP4 Tanimoto similarity. The configured rules are:

```text
in domain:    Tanimoto >= 0.40
near domain:  0.25 <= Tanimoto < 0.40
out of domain:Tanimoto < 0.25
unassessable: nearest-reference Tanimoto missing
```

Nearest-reference USRCAT similarity and distance are retained as independent continuous 3D domain measurements, but no calibrated USRCAT threshold is imposed.

Out-of-domain and unassessable rows are ordered after eligible rows in organism-specific shortlists while their underlying numeric evidence scores are left unchanged.

---

## 15. Leakage-controlled benchmarking

Public antibacterial benchmark compounds with curated mechanism-level target labels are evaluated under leakage-aware splits.

An ECFP4 analogue guard removes reference molecules with Tanimoto similarity >= 0.85 to the query in every supported v3 split. After filtering, the implementation verifies that no retained reference violates this threshold.

The v3 split types are:

### 15.1 Target-family holdout

References belonging to the query's labelled target family are withheld using the pinned ontology. If the query label cannot be mapped to a target family, the result is marked pending rather than inferred.

### 15.2 Scaffold split

References sharing an exact non-empty Bemis-Murcko scaffold with the query are removed in addition to the close-analogue guard.

### 15.3 Temporal split

The configured cutoff is 2018-01-01. Only pre-cutoff references are eligible and only dated post-cutoff queries may form the test set. Records with missing dates are excluded. The pinned legacy data lack sufficient date provenance for a complete temporal benchmark, so unavailable temporal results must remain pending.

---

## 16. Benchmark metrics and statistical analysis

Benchmarking reports:

- AUROC;
- BEDROC with alpha = 20.0 and 80.5;
- enrichment factor at 1% and 5%;
- mean reciprocal rank (MRR); and
- coverage.

Exact score ties receive average ranks. The enrichment-factor implementation handles a tied decision boundary fractionally.

Metrics are calculated per query and then aggregated. Nonparametric bootstrap confidence intervals resample query identifiers rather than individual target rows, thereby preserving within-query dependence. The default benchmark bootstrap uses 1,000 resamples with the configured bootstrap seed. Aggregate results report point estimates, 95% confidence intervals, total and evaluable query counts, seed, snapshot identifier, and split-removal provenance.

V3 compares three chemical modes:

1. legacy 2D-only evidence;
2. 3D-only evidence comprising USRCAT, O3A shape, O3A color, and aligned 3D pharmacophore; and
3. seven-component fused evidence.

The 2D Gobbi pharmacophore field is excluded from the declared 3D-only mode.

Undefined metrics remain unavailable rather than being replaced by zeros.

---

## 17. Sensitivity analysis

Because the current pinned benchmark does not provide the data partitions required for a publishable trained and calibrated combiner, the active justification path uses sensitivity analysis of the transparent heuristic scores.

For chemical RRF, each implicit equal component weight is multiplied by 0.50, 0.75, 1.25, and 1.50; individual components and grouped evidence layers are also removed in turn. Ranking stability is summarized with Kendall tau and rank-biased overlap (RBO).

For the final organism-aware ranking, configured coefficients in the specificity transform, reference-quality mapping, transfer, pocket, biological, anti-target, and biological-subweight terms are perturbed or ablated. Confidence intervals are derived by resampling complete organism-query ranked lists.

Sensitivity outputs remain non-probabilistic.

---

## 18. External baseline and calibration status

PIDGINv4 is the configured external baseline. The code is pinned to commit `df0f6068a8aa16e2278e3779a1ad5e6d552731dc`, with model release DOI `10.6084/m9.figshare.19108382.v1`. The adapter requires the original Python 2.7 runtime, model files, applicability-domain resources, and an explicit mapping from PIDGIN targets to the repository ontology.

These prerequisites are not present in the pinned default analytical environment. Consequently, PIDGIN head-to-head outputs are correctly reported as pending rather than simulated.

The repository also contains a regularized logistic-regression plus Platt-calibration path, but publishable fitting requires non-overlapping training, calibration, and test query sets. Rows missing selected features are excluded rather than imputed. Because the current pinned benchmark does not contain the required valid partitions, no calibrated probability should be reported from the default snapshot.

---

## 19. Output interpretation

The primary scientific output is a ranked set of organism-target hypotheses with a transparent decomposition of the supporting evidence. The recommended interpretation is:

- chemical similarity supports **target plausibility**, not binding proof;
- sequence transfer supports **orthologue compatibility**, not equal ligand affinity;
- essentiality and accessibility support **biological relevance**, not whole-cell activity;
- CARD annotations provide **resistance context**, not target engagement;
- RCSB structures provide **structural precedent**, not docking evidence; and
- the overall score supports **prioritization**, not probability of success.

The strongest next-step validation for a high-ranking direct target is an orthogonal sequence of biochemical or binding assays, species-orthologue testing, whole-cell susceptibility with permeability/efflux controls, resistant-mutant or complementation experiments, and selectivity testing. Docking and molecular dynamics should be treated as downstream structural analyses after a target survives chemical, organism, and assay-level review.

---

## 20. Limitations

The principal limitations are the heterogeneity of public bioactivity assays, unequal reference-ligand coverage among target classes, potential scaffold and chemical-space bias, incomplete provenance for several legacy public datasets, strain dependence of species mapping, limited interpretability of curated ordinal biological scores, and the inability of chemical similarity alone to capture permeability, efflux, metabolism, expression, or compound exposure.

The v3 3D stage is dependent on successful conformer generation and MMFF/O3A compatibility. Large, highly flexible, macrocyclic, or peptide-like molecules may therefore have partial or unavailable 3D evidence. O3A is evaluated only for the USRCAT shortlist, which should be acknowledged when interpreting 3D sensitivity analyses.

Temporal benchmarking, property-matched decoys, external PIDGIN comparison, and a publishable calibrated combiner are currently incomplete for the frozen legacy snapshot. These are explicit data/runtime limitations and should not be filled by simulated values.

---

## 21. Data and software resources

The workflow uses or references the following public resources:

- ChEMBL for public target-associated bioactivity records: https://www.ebi.ac.uk/chembl/
- UniProt for species-specific protein mapping: https://www.uniprot.org/
- Comprehensive Antibiotic Resistance Database (CARD) for resistance annotations: https://card.mcmaster.ca/
- RCSB Protein Data Bank for structure and co-crystal metadata: https://www.rcsb.org/
- PubChem for benchmark-structure retrieval: https://pubchem.ncbi.nlm.nih.gov/
- RDKit for molecular standardization, fingerprints, conformer generation, shape alignment, and pharmacophore calculations: https://www.rdkit.org/
- WHO bacterial priority pathogens list (2024): https://www.who.int/publications/i/item/9789240093461
- De Oliveira DMP et al. Antimicrobial Resistance in ESKAPE Pathogens. *Clin Microbiol Rev.* 2020;33(3):e00181-19. doi:10.1128/CMR.00181-19.

---

## 22. Suggested concise manuscript Methods paragraph

A reproducible organism-aware target-prioritization workflow was used to generate antibacterial target hypotheses from small molecules without restricting the initial search to a predefined docking panel. Query compounds were standardized with RDKit and compared with ChEMBL-derived reference ligands using Morgan/ECFP4 (radius 2, 2,048 bits), MACCS, deterministic ETKDGv3/MMFF94 conformer ensembles, USRCAT, O3A shape/color overlap, and both Gobbi/Poppe 2D and O3A-aligned 3D pharmacophore similarity. Seven chemical evidence components were integrated by normalized Reciprocal Rank Fusion (k = 60), while target-specificity, reference quality, and applicability-domain information were retained separately. Chemically supported targets were subsequently annotated for species-specific sequence transfer, organism scope, essentiality, cellular accessibility, clinical precedent, CARD resistance context, RCSB structure/pocket precedent, and annotation-only anti-target risk across six configured bacterial species. The final organism-aware priority score was a transparent heuristic combination of chemical quality, species transfer, pocket evidence, biological priority, and anti-target adjustment and was not interpreted as a probability. Benchmarking used an ECFP4 close-analogue exclusion threshold of 0.85 and target-family, Bemis-Murcko scaffold, and temporal split logic, with AUROC, BEDROC, EF, MRR, coverage, and query-level bootstrap confidence intervals. All parameters, random seeds, software versions, snapshot identifiers, configuration hashes, and run metadata were recorded for reproducibility.