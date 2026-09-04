# Organism-Aware Antibacterial Open Target Discovery (v3)

[![CI](https://github.com/OASolliman590/organism-aware-antibacterial-target-prioritization/actions/workflows/ci.yml/badge.svg)](https://github.com/OASolliman590/organism-aware-antibacterial-target-prioritization/actions/workflows/ci.yml)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3120/)
[![RDKit](https://img.shields.io/badge/cheminformatics-RDKit-3838ff.svg)](https://www.rdkit.org/)

This repository contains a reproducible cheminformatics framework for **open target discovery** from small molecules, followed by organism-specific biological and clinical annotation. It is designed for antibacterial research involving the ESKAPE pathogens and related organisms.

> **Data-protection policy:** unpublished compound structures, compound names, docking/MD outputs, and compound-specific predictions are intentionally excluded from this GitHub repository. Authorized local inputs belong in `inputs/`, which is ignored by Git. The public repository contains code, public benchmark data, public reference data, and documentation only.

## Scientific objective

The pipeline does not force a compound into a preselected target panel. It first searches a broad public target universe using chemical evidence, then annotates each chemically plausible target with species-specific protein mapping, sequence transfer, clinical and essentiality context, cellular accessibility, resistance biology, structure-pocket precedent, anti-target risk, and uncertainty. The output is a ranked set of **testable target hypotheses**, not proof of binding or target engagement.

The workflow separates the following evidence layers:

| Layer | Examples | Scientific interpretation |
|---|---|---|
| Chemical | ECFP4/Morgan, MACCS, nearest-neighbour consistency | Ligand-space support |
| 3D shape | ETKDGv3, USRCAT, O3A shape/color | Conformer-dependent shape/feature support |
| Pharmacophore | O3A-aligned BaseFeatures 3D overlap; Gobbi/Poppe 2D feature-pair similarity | Separate spatial and alignment-free scaffold-independent feature support |
| Specificity | Cross-target decoys, target-specificity margin | Whether the signal is target-specific |
| Reference quality | Reference count, scaffold diversity, evidence grade | Stability of public support |
| Species transfer | UniProt accession, mapping status, sequence similarity/coverage | Transferability to the organism protein |
| Biology | Essentiality/fitness, clinical status, accessibility | Biological prioritization, not binding proof |
| Resistance | CARD models, CARD SNP rows, organism-linked mutation counts | Resistance context and target-variation risk |
| Structure | RCSB candidates, co-crystallized ligand flag | Future docking/pocket precedent; no docking performed |
| Safety | Human-homologue/mitochondrial annotations | Early selectivity warning |
| Uncertainty | Bootstrap stability, decoy calibration, reference coverage | Reliability of each hypothesis |

## Repository layout

| Path | Purpose |
|---|---|
| `pipeline/` | Public workflow modules for preparation, open target scoring, benchmarking, calibration, sequence mapping, resistance parsing, structure cataloguing, figures, and reporting |
| `data/reference_ligands/` | Public ChEMBL-derived reference-ligand records used by the target classes |
| `data/benchmark/` | Public ESKAPE-active antibacterial benchmark metadata and structures |
| `data/target_ontology_v2.csv` | Auditable parent target ontology separating direct targets, complexes, pathways, resistance determinants, and phenotypic mechanisms |
| `data/target_subtype_ontology_v21.csv` | v2.1 mechanism/site subtypes for PBPs, beta-lactamases, 30S, and 50S ribosomal targets |
| `data/chembl_target_aliases_v21.json` | Versioned ChEMBL target aliases and manually verified target-ID seeds for v2.1 retrieval |
| `data/species_targets/` | Public species-target mapping and sequence-transfer metadata when present in the repository release |
| `data/resistance_v2/` | CARD-derived public resistance-family and SNP summaries when present in the repository release |
| `data/structures_v2/` | Bounded RCSB structure/pocket metadata; no docking inputs or private structures |
| `inputs/` | Local-only research inputs; ignored by Git and represented publicly only by `.gitkeep` |
| `data/compounds/` | Local-only normalized compound structures; ignored by Git |
| `results/` | Local-only private prediction tables, figures, and reports; ignored by Git |
| `run_pipeline.py` | Local end-to-end runner |
| `config.yaml` | Single declarative source for paths, thresholds, modes, and seeds |
| `requirements.lock` | Exact transitive Python environment |

The full v3 protocol, formulas, failure semantics, and current data gaps are in
[`docs/v3_method_supplement.md`](docs/v3_method_supplement.md). V3 preserves the
v2 engine and adds deterministic 3D/pharmacophore evidence, query-local rank
fusion, target/scaffold/temporal split machinery, bootstrap uncertainty,
sensitivity analysis, and applicability-domain annotations. Scores are not
probabilities unless a held-out calibrated field explicitly says so.

## Local execution with authorized unpublished inputs

Copy authorized input structures into the local `inputs/` directory. Do not commit them. Install dependencies and run:

```bash
python -m pip install -r requirements.lock
python run_pipeline.py --config config.yaml
```

The private workflow writes normalized structures and compound-specific outputs to ignored directories. Public-only benchmark and annotation modules can be run without private inputs when their public data are present.

The default config is analysis-only: it verifies the immutable public snapshot in
`data/snapshots/SNAPSHOT_VERSIONS.json`, never calls a live external API, and writes a
provenance manifest for the run. Live refresh requires a separate config with
`run.refresh_external_data: true`, a date-prefixed snapshot ID, explicit ChEMBL/CARD/
UniProt releases, and every refreshed output path located under
`data/snapshots/<snapshot-id>/`. Refresh writers use create-only files and refuse to
overwrite an existing snapshot. A raw refresh is marked `analysis_ready: false` until
all derived layers have been rebuilt and independently pinned.

## V2 target-discovery workflow

The method has six conceptual stages:

1. **Structure standardization.** Molecules are sanitized, canonicalized, and represented using ECFP4/Morgan and MACCS fingerprints. Explicit hydrogens are removed for fingerprint comparison.
2. **Open chemical discovery.** Each query is compared against all target classes with public reference support rather than a user-defined target panel.
3. **Specificity and quality correction.** Cross-target decoys, target-specificity margins, reference count, scaffold diversity, and evidence grade are retained as separate fields.
4. **Species and biological annotation.** Organism-specific UniProt mapping, sequence transfer, clinical status, essentiality/fitness, localization/accessibility, and resistance relevance are attached after chemical scoring.
5. **Resistance and structure context.** CARD model/SNP context and bounded RCSB structure/co-crystal metadata are added without treating them as direct binding evidence.
6. **Uncertainty and validation planning.** Bootstrap stability, decoy calibration, mapping status, and evidence coverage generate confidence classes and assay recommendations.

A simplified priority summary is:

```text
chemical_quality_adjusted_evidence
    = chemical_similarity × target_specificity × reference_quality

organism_aware_priority
    = chemical_quality_adjusted_evidence
    × species_transfer
    × pocket_evidence
    × biological_priority
    × anti_target_adjustment
```

These are decision-support scores, not calibrated binding probabilities. A prediction with unresolved sequence mapping or poor target specificity is explicitly downgraded or classified as insufficient.

## V2.1 reference-universe expansion

v2.1 adds a reproducible target-subtype layer to the **same broad target universe for all six organisms**. PBPs are only one additive family among the existing targets and newly represented beta-lactamase and ribosomal mechanism/site classes. PBPs are separated into PBP1A, PBP1B, PBP2, PBP2B, PBP2X, PBP3, and PBP4. Beta-lactamases are separated by Ambler class, and ribosomal records are represented by 30S aminoglycoside/tetracycline sites and 50S macrolide, oxazolidinone, pleuromutilin, and chloramphenicol sites when ChEMBL target and assay metadata support that resolution. No organism is assigned a PBP-only panel, and no target family receives a hard-coded ranking bonus.

The new `pipeline/fetch_chembl_reference_subtypes_v21.py` performs bounded ChEMBL target discovery, activity retrieval, pChEMBL filtering, confidence-aware assay grading, conservative RDKit structure standardization, InChIKey deduplication, and provenance retention. The retrieval is cacheable and supports `CHEMBL_V21_DISCOVER=0`, `CHEMBL_V21_SUBTYPES=...`, `CHEMBL_V21_OFFLINE=1`, and `CHEMBL_V21_MAX_PER_SUBTYPE=...` for controlled reruns. A failed or temporarily unavailable ChEMBL endpoint must not be replaced by simulated data; the manifest and empty subtype outputs preserve the missingness for later refresh.

v2.1 reports two linked rankings across the complete available target universe. The **chemical hypothesis score** answers which target or binding-site subtype is chemically compatible with the compound. The **clinical-translational score** combines clinical precedent, organism scope, essentiality, accessibility, structure precedent, and resistance context. Species-specific mapping and organism scope modify each target independently; they do not replace the broad universe with the original user docking panel. These scores are deliberately retained separately rather than collapsed into a single claim of antibiotic activity.

Because the ChEMBL REST service was temporarily unavailable during this implementation run, the v2.1 code and ontology were integrated and tested, but no new PBP, beta-lactamase, or ribosomal subtype activity records were fabricated. The same limitation applies equally to every newly added subtype; it is not a reason to prioritize PBPs over other targets. The current numerical benchmark therefore remains the v2 public-reference benchmark until the cached subtype retrieval completes successfully. This limitation is recorded in the local fetch logs and is an explicit reproducibility safeguard.

## Benchmarking

The benchmark uses public antibacterial drugs with mechanism-level target labels. Query molecules and close analogues are excluded from the reference set at ECFP4 similarity ≥0.85, and exact Bemis–Murcko scaffold exclusion is evaluated separately. Results report coverage, top-1/top-3/top-5 retrieval, MRR, enrichment over random, prevalence baselines, split status, and per-query ranks.

The benchmark currently supports close-analogue and scaffold diagnostics. The scaffold result is reported transparently; if the available benchmark contains no additional same-scaffold exclusions, identical summary values are not interpreted as proof of scaffold-level generalization. Target-family holdout, temporal, and species-holdout results are reported as unavailable or partial when the public metadata do not support a valid evaluation.

## Resistance and structure annotations

CARD-derived annotations describe resistance determinants and mutations associated with target families. They are not treated as evidence that a private compound binds a target. RCSB PDB metadata identify structure candidates and whether non-solvent co-crystallized ligands are present. Co-crystallized ligands are used to identify a defensible future docking site; PDBQT preparation, docking, and molecular dynamics remain downstream user-controlled steps.

## Visualization suite

`pipeline/figures_suite.py` renders the compound-facing view of a completed run into
`results/figures_suite/`, using the shared visual language defined in
`pipeline/figure_style.py`. It runs as part of `run_pipeline.py`, or on its own against an
existing `results/` directory:

```bash
python pipeline/figures_suite.py
```

Like every other pipeline module it reads its configuration from `config.yaml`, or from the
path in the `OATP_CONFIG` environment variable when one is set.

| Panel | Question it answers | Source table |
| --- | --- | --- |
| `compound_target_priority` | Which organism–target pairs does each private compound favour? | `v2_open_target_predictions_by_organism.csv` |
| `organism_target_atlas` | Which target classes carry the study overall? | `v2_open_target_predictions_by_organism.csv` |
| `evidence_decomposition` | Which evidence layer drives each compound's best hypothesis? | `v2_open_target_predictions_by_organism.csv` |
| `confidence_profile` | How is confidence distributed by compound and organism? | `v2_open_target_predictions_by_organism.csv` |
| `fusion_contribution` | What does each fusion component contribute to the v3 score? | `open_target_predictions_by_organism_v3.csv` |
| `rank_shift_disagreements` | Where does 3D fusion disagree with 2D-only ranking? | `chemical_evidence_disagreements_v3.csv` |
| `ranking_stability` | How far does the ranking move under weight perturbation? | `final_ranking_sensitivity_v3.csv` |
| `uncertainty_landscape` | Which hypotheses are both stable and decoy-separated? | `v2_uncertainty_private.csv` |
| `benchmark_enrichment` | How much better than random is retrieval? | `benchmark_v2_summary.csv` |
| `chemical_space` | Do the private compounds sit near known benchmark drugs? | private SDF plus `data/benchmark/eskape_benchmark_drugs.csv` |
| `assignment_concordance` | Does each compound score best against the organism it was made for? | manifest plus `v2_open_target_predictions_by_organism.csv` |

### Per-organism suites

When a private compound manifest is present at `paths.compound_manifest` (default
`inputs/compound_manifest.csv`, gitignored), the suite additionally renders one focused set
per configured organism under `results/figures_suite/by_organism/<organism>/`.

The manifest needs a `compound_code` and a `microbe_group` column. Its microbe labels are
mapped onto `organisms.names` through the explicit `organisms.manifest_aliases` table in
`config.yaml` — never by fuzzy matching, because the manifest uses abbreviations and at
least one misspelled genus, and a silently wrong resolution would mislabel every
per-organism figure. A group with no alias entry is left unassigned and reported.

Each per-organism suite covers **the compounds the manifest prepared against that organism**,
scored against it: the `klebsiella_pneumoniae` folder is about the Klebsiella compounds, and
a footnote on each figure names the set. Because the target axis and the colour scale are
recomputed from only those compounds, a target that leads for them is not buried under a
stronger compound designed for a different organism.

The cross-series comparison is not lost — it lives in the overall suite, whose per-compound
panels already show every compound against every organism. Setting `restrict_to_assigned`
to false on the `SuiteContext`, or supplying no manifest, falls back to the whole series
scored against one organism, with assigned compounds marked `*` in the axis labels.

Panels that exist to compare organisms against each other (`organism_target_atlas`,
`assignment_concordance`) are skipped in the per-organism suites and recorded as
`not_applicable_for_single_organism`, since filtering to one organism would leave them
drawing a misleading figure from the surviving rows.

The manifest records design intent, not evidence: it is used to focus and annotate figures,
and never to weight a score.

Every panel is drawn only from a completed run table. A panel with no source file, no
evaluable rows, or an unavailable optional dependency writes no image and instead records
the reason in `results/figure_suite_status.csv` as `unavailable_source_missing`,
`unavailable_no_evaluable_rows`, or `unavailable_dependency_unavailable`. A missing figure
must always be explainable from that table; empty or placeholder panels are never emitted.
`chemical_space` needs RDKit and UMAP, whose compiled extensions are unavailable on some
hosts, so it degrades to a dependency status rather than failing the run.

## Tying a docking campaign to the ranking

`pipeline/docking_crosscheck.py` joins an external docking/MD campaign to the pipeline
ranking, so a run can be reported against the targets that were actually tested:

```bash
python pipeline/docking_crosscheck.py
```

It is not part of `run_pipeline.py`, because it needs a private campaign summary that a
public run will not have. Point `paths.docking_summary` at that table and declare the
receptor mapping in the `docking` block of `config.yaml`.

It writes three tables to `results/` and one figure:

| Output | Answers |
| --- | --- |
| `docking_pipeline_coverage.csv` | Which docked receptors map onto a target class, which are excluded and why, and which map to nothing |
| `docking_pipeline_crosscheck.csv` | Per compound and target: docked affinity, ligand efficiency, the native co-crystal control, and the pipeline priority, rank and confidence |
| `docking_pipeline_agreement.csv` | Spearman rank agreement between affinity and priority, pooled and per compound |

Three deliberate choices make the comparison honest rather than flattering:

- **The receptor mapping is declared, never inferred.** An unmapped receptor is reported as
  `unmapped_no_alias`, and a campaign where *nothing* maps raises instead of writing an
  empty correlation. An empty agreement file reads as "no disagreement found", which is the
  opposite of what it means.
- **Repeated structures for one target class are collapsed** before correlating. Two crystal
  forms of the same enzyme share one pipeline priority, so keeping both would pseudo-replicate
  that value and inflate both the sample size and the coefficient.
- **Coverage is reported next to agreement.** A correlation over five shared targets says
  little without knowing which of the pipeline's strongest hypotheses were never docked, so
  the undocked top targets are listed alongside it.

Docking affinity and pipeline priority are different quantities — a pose in one structure
against target-class plausibility — so only their *rank* agreement is compared, never their
values. Because docking scores grow with molecular size, affinity is also reported per heavy
atom against the native co-crystal ligand, which is the control that says whether an affinity
is good for a molecule of that size.

## Limitations

Reference-ligand assays are heterogeneous, target-family coverage is uneven, and whole-cell antibacterial activity depends on permeability, efflux, expression, metabolism, and target accessibility. UniProt mapping is also strain-dependent; unresolved mapping is uncertainty, not evidence of target absence. CARD record counts reflect curation and database density rather than resistance prevalence in a user’s isolates. Anti-target fields are early alerts rather than toxicity predictions.

The strongest next validation is a staged experiment: purified target or target-complex assay, species-orthologue comparison, MIC/time-kill with permeability and efflux controls, resistant-mutant selection, and complementation or target rescue. Docking and MD should be performed only for targets that survive chemical, species, and assay-level review.

## Public sources

- [WHO bacterial priority pathogens list, 2024](https://www.who.int/publications/i/item/9789240093461)
- [ESKAPE antimicrobial-resistance review](https://pmc.ncbi.nlm.nih.gov/articles/PMC7227449/)
- [ChEMBL](https://www.ebi.ac.uk/chembl/)
- [RDKit](https://www.rdkit.org/)
- [UniProt API documentation](https://www.uniprot.org/api-documentation/proteomes)
- [CARD downloads and ontology](https://card.mcmaster.ca/download)
- [RCSB PDB Search API](https://search.rcsb.org/)
- [RCSB PDB Data API](https://data.rcsb.org/)

## Reproducibility

The main v2.1 modules are `pipeline/open_target_discovery_v2.py`, `pipeline/benchmark_v2.py`, `pipeline/fetch_chembl_reference_subtypes_v21.py`, `pipeline/calibrate_uncertainty_v2.py`, `pipeline/fetch_species_targets.py`, `pipeline/sequence_compatibility.py`, `pipeline/build_reference_quality.py`, `pipeline/build_card_resistance_annotations.py`, `pipeline/parse_card_snps_v2.py`, `pipeline/fetch_structure_catalog_v2.py`, `pipeline/v2_figures.py`, including `v21_chemical_vs_clinical_translation.png` and `v21_clinical_translation_heatmap.png`, `pipeline/build_validation_plan_v2.py`, `pipeline/figures_suite.py` with
`pipeline/figure_style.py`, and `pipeline/summarize_v2.py`.

All unpublished inputs and compound-specific outputs remain local until the compounds are published and the author explicitly approves a new repository release. The frozen
legacy v2 snapshot records CARD 4.0.2; the original ChEMBL release, UniProt release/query
date, PubChem query date, and RCSB query date were not retained by v2 and are therefore
recorded as explicit provenance gaps rather than inferred.
