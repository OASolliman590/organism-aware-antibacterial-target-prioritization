# Open organism-aware antibacterial target prediction

This repository contains a reproducible cheminformatics framework for **open target discovery** from small-molecule structures, followed by organism-specific biological and clinical filtering. It is designed for antibacterial research involving ESKAPE organisms and related pathogens.

> **Data-protection policy:** unpublished compound structures, compound names, docking/MD outputs, and compound-specific predictions are intentionally excluded from this GitHub repository. Place authorized local inputs in `inputs/`; that directory is ignored by Git.

## Scientific objective

The pipeline does not force a compound into a preselected target panel. It first generates chemically plausible target hypotheses from a broad reference universe, then applies organism-specific filters based on target family, bacterial essentiality, cellular localization, clinical precedent, resistance biology, and target availability. The output is a ranked hypothesis set for docking, molecular dynamics, biochemical assays, and cellular validation—not proof of target engagement.

The framework uses RDKit ECFP4/Morgan and MACCS fingerprints, nearest-neighbour reference-ligand evidence, target-family aggregation, structural/SAR flags, chemical-space visualization, and leakage-aware validation on known antibacterial drugs. Target scores report their evidence components separately so that target ranking is auditable rather than presented as a calibrated probability.

## Repository layout

| Path | Purpose |
|---|---|
| `pipeline/` | Compound preparation, reference-ligand scoring, open target discovery, organism filtering, benchmarking, figures, and report generation |
| `data/reference_ligands/` | Public reference-ligand records used by the current target-family modules |
| `inputs/` | Local-only research inputs; ignored by Git and represented only by `.gitkeep` |
| `data/compounds/` | Local-only normalized compound structures; ignored by Git |
| `results/` | Local-only compound-specific outputs; ignored by Git |
| `run_pipeline.py` | Local end-to-end runner |
| `requirements.txt` | Python dependencies |

## Local execution with authorized unpublished inputs

Copy authorized input files into `inputs/`, including the structure files and any private target metadata. The input directory is not committed or uploaded by the standard Git workflow.

```bash
python -m pip install -r requirements.txt
python run_pipeline.py
```

Set `PROJECT_ROOT` and `INPUT_DIR` when running from another location. The pipeline writes compound-specific structures, tables, figures, and reports to ignored local directories.

## Open target-discovery design

The revised workflow is intentionally two-stage. **Stage A: chemically valid target discovery.** Compounds are compared against a broad target-family reference universe, including clinically relevant antibacterial target families such as beta-lactamases, DNA gyrase/topoisomerase IV, FabI/FabH, Mur enzymes, FtsZ, DHFR/DHPS, Lpx enzymes, aminoacyl-tRNA synthetases, and other families when reference coverage is adequate. Candidates are retained only when similarity, nearest-neighbour consistency, fingerprint agreement, and structural plausibility meet transparent criteria.

**Stage B: organism-specific clinical filtering.** Each chemically plausible target is annotated for organism relevance, essentiality or fitness impact, subcellular accessibility, clinical validation, resistance relevance, and evidence quality. This stage ranks target hypotheses for a specific organism without pretending that the organism itself determines ligand binding.

The general score is:

```text
chemical_target_evidence =
    weighted ECFP4 nearest-neighbour evidence
  + MACCS agreement
  + reference-consistency evidence
  + target-family SAR plausibility

organism_clinical_priority =
    chemical_target_evidence
  × organism relevance
  × essentiality/fitness evidence
  × accessibility/clinical evidence
```

The factors are reported separately. They are not calibrated probabilities, and a high clinical priority cannot rescue chemically implausible ligand evidence.

## Benchmarking

The repository includes an evaluation framework for known antibacterial drugs with curated target labels. Benchmarking is leakage-aware: a drug or near-duplicate analogue must not be allowed to contribute its own target evidence to the reference set used to predict it. Performance is reported with top-1/top-3/top-5 recall, mean reciprocal rank, mean percentile rank, enrichment over random target ranking, and per-target-family confusion or retrieval plots.

Benchmark drugs should be sourced from public, citable resources such as ChEMBL, DrugBank records where available, FDA labels, and primary literature. Their target labels must distinguish direct molecular targets from resistance mechanisms, phenotypic mechanisms, and broad target-family annotations.

## Reproducibility and limitations

Reference-ligand data are heterogeneous and target-family coverage is uneven. Fingerprint similarity is a ligand-space hypothesis generator, not a binding assay. Cross-organism target records, whole-complex measurements, resistance determinants, efflux, permeability, expression, and target-site sequence variation require separate annotation or experimental validation. Sparse target classes should be marked low confidence rather than hidden.

All unpublished inputs and compound-specific outputs remain local until the compounds are publicly released and the author explicitly approves a new repository release.

## References

1. [RDKit documentation](https://www.rdkit.org/)
2. [ChEMBL database and API](https://www.ebi.ac.uk/chembl/)
3. [UMAP documentation](https://umap-learn.readthedocs.io/)
4. [WHO bacterial priority pathogens](https://www.who.int/publications/i/item/9789240093461)
