from pathlib import Path
import os
import pandas as pd

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1]))
RES=ROOT/'results'
OUT=RES/'revised_open_target_report.md'

def md(df, cols=None, n=10):
    if cols: df=df[cols]
    return df.head(n).to_markdown(index=False)

metrics=pd.read_csv(RES/'benchmark_metrics.csv')
summary=pd.read_csv(RES/'benchmark_summary.csv')
covered=metrics[metrics.target_in_reference_universe==1]
allrow=summary[summary.query_target_label=='ALL'].iloc[0]
short=pd.read_csv(RES/'open_target_shortlist_by_organism.csv')
agg=(short.groupby(['organism','target_class'],as_index=False)
     .agg(mean_open_target_priority=('open_target_priority','mean'),
          max_chemical_evidence=('chemical_evidence','max'),
          n_compound_appearances=('query_id','nunique'),
          mean_clinical_priority=('organism_clinical_priority','mean')))
agg=agg.sort_values(['organism','mean_open_target_priority'],ascending=[True,False])

org_table=agg.groupby('organism',group_keys=False).head(5)
bench_cols=['query_id','query_target_label','accepted_target_classes','target_in_reference_universe','rank_of_known_target','top1_hit','top3_hit','top5_hit','reciprocal_rank','top_predicted_target']

text=f'''# Revised open organism-aware antibacterial target-discovery report

## Executive summary

This revision changes the workflow from **ranking within a fixed, user-selected target panel** to a two-stage procedure. First, the chemistry is compared against every target family represented in the public reference-ligand universe. Second, chemically plausible targets are filtered and prioritized for each organism using explicit annotations for clinical validation, essentiality or fitness, cellular accessibility, and resistance relevance. The organism filter therefore modifies biological priority; it does not manufacture chemical evidence.

The unpublished compound structures and compound-specific results are retained locally and are excluded from GitHub. The public repository contains the code, public reference data, benchmark structures/labels, target annotations, and documentation only. This separation should remain in place until the compounds are published.

A public benchmark of **16 known antibacterial drugs** produced **14/16 target-label coverage** in the current reference universe. Among covered labels, top-1 retrieval was **{allrow.top1_recall_covered:.1%}**, top-3 retrieval **{allrow.top3_recall_covered:.1%}**, top-5 retrieval **{allrow.top5_recall_covered:.1%}**, and mean reciprocal rank **{allrow.mrr_covered:.3f}**. These values are encouraging for a hypothesis generator but are not a clinical validation claim: the benchmark is small, target-class labels are mechanism-level, and two drug mechanisms were outside current target-reference coverage.

## 1. Data protection and reproducibility

The following categories are intentionally absent from GitHub: unpublished SDF/MOL structures, compound names associated with the private project, compound-specific fingerprints and descriptors, private target lists, docking/MD coordinates and scores, and private target-prediction tables or figures. The local `inputs/`, `data/compounds/`, and `results/` directories are ignored by Git. The pipeline can still process authorized local inputs when they are present.

Public reference-ligand records were retrieved from ChEMBL and benchmark structures from PubChem. The target annotation table is versioned separately from compound data and records clinical validation and organism-fit notes. The WHO bacterial priority-pathogen list and the ESKAPE AMR review provide the public clinical context for prioritization [1,2].

## 2. Revised methodology

### 2.1 Stage A: open chemical target discovery

Each query molecule is standardized with RDKit and represented using Morgan/ECFP4 and MACCS fingerprints. For every target class with reference ligands, the pipeline computes the maximum ECFP4 Tanimoto similarity, the mean of the five best ECFP4 similarities, and the maximum MACCS similarity. Reference molecules with ECFP4 similarity at least 0.85 are excluded for benchmark queries to prevent exact or close-analogue leakage.

The chemical evidence score is a transparent composite of ECFP4 nearest-neighbour evidence, top-five consistency, and MACCS agreement. It is a ranking score, not a probability and not a predicted binding affinity. Target classes with sparse or absent reference chemistry are reported as low-coverage hypotheses rather than silently interpreted as negative predictions.

### 2.2 Stage B: organism-specific clinical filtering

The chemically supported target candidates are annotated by target family. Approved or clinically validated families include gyrase/topoisomerase, DHFR/DHPS, MurA, RpoB, ribosome, PBPs, glycopeptide cell-wall precursor recognition, and selected resistance targets. FabI, FtsZ, MurC/MurE, LpxA/LpxC/LpxH, and related families remain valid discovery targets but are labelled preclinical or limited-translation where appropriate.

For each organism, the clinical-priority annotation combines organism compatibility, clinical validation, essentiality/fitness evidence, and resistance relevance. Accessibility and organism-specific notes are retained as fields rather than hidden in an overconfident scalar score. Beta-lactamases and PBP2a are treated as resistance or pathogen-specific targets; they are not assigned universal viability-target status.

### 2.3 Relationship to prior docking and molecular dynamics

The prior docking/MD work for *K. pneumoniae*, *A. baumannii*, and MRSA is retained as an **orthogonal structural-validation layer**. It should be compared against the new open-target output after chemical discovery, not used to define the target universe or to train the fingerprint score. Agreement supports prioritization; disagreement is scientifically useful and should trigger sequence, pocket, permeability, and assay review rather than automatic rejection.

## 3. Public ESKAPE benchmark

The benchmark contains 16 real antibacterial drugs with PubChem structures and conservative mechanism-level labels. It includes fluoroquinolone gyrase/topoisomerase chemistry, trimethoprim/DHFR, sulfamethoxazole/DHPS, fosfomycin/MurA, triclosan/FabI, ceftaroline/PBP, rifampicin/RpoB, ribosome-active drugs, vancomycin/D-Ala-D-Ala, daptomycin/membrane action, meropenem/PBP, and colistin/lipid-A or membrane action. When a drug mechanism was too broad for a single protein, the label was kept at the appropriate target-family level.

| Metric | Value |
|---|---:|
| Benchmark drugs | 16 |
| Target labels represented by the reference universe | {int(allrow.n_covered)}/{int(allrow.n_queries)} ({allrow.coverage_fraction:.1%}) |
| Top-1 recall among covered labels | {allrow.top1_recall_covered:.1%} |
| Top-3 recall among covered labels | {allrow.top3_recall_covered:.1%} |
| Top-5 recall among covered labels | {allrow.top5_recall_covered:.1%} |
| Mean reciprocal rank among covered labels | {allrow.mrr_covered:.3f} |
| Close-analogue leakage cutoff | ECFP4 Tanimoto ≥ 0.85 |

The detailed query-level results are shown below.

{md(metrics,bench_cols,20)}

The benchmark also includes two important negative controls. Broad PBP mechanisms were only partially covered because the current ChEMBL PBP reference set is sparse, and general membrane action is not a single protein target. These are **coverage limitations**, not evidence that the open method disproves those mechanisms. Fosfomycin/MurA and vancomycin/D-Ala-D-Ala were represented but ranked below the top prediction in this ligand-reference snapshot, demonstrating that the method still requires larger, better-curated, target-specific reference sets.

## 4. Private-compound open-target results

The private run includes all 12 normalized compounds locally. The full compound-by-target scores, MACCS and ECFP4 heatmaps, chemical-space map, and organism-clinical filtering plot remain in the ignored local `results/` directory. The main shortlist is produced by aggregating the top five open-target hypotheses per compound and organism; it should be interpreted as a starting point for docking or biochemical testing, not as a definitive target assignment.

### Aggregate shortlist by organism

{md(org_table,['organism','target_class','mean_open_target_priority','max_chemical_evidence','n_compound_appearances','mean_clinical_priority'],60)}

The shortlist is broader than the original five-target panels. It includes clinically validated families such as gyrase/topoisomerase, DHFR/DHPS, MurA, RpoB, and ribosome where ligand evidence exists, while retaining preclinical but biologically credible families such as FtsZ, MurC, LpxC, LpxH, and FabI. PBP2a remains organism-specific to MRSA, and Gram-negative lipid-A targets are filtered through permeability and target-distribution notes.

## 5. Recommended decision rule for docking and MD

Docking and MD should be performed only after the open chemical score and organism annotations have been inspected separately. A practical selection rule is to choose candidates with at least moderate nearest-neighbour evidence, either MACCS agreement or a coherent ECFP4 neighbour cluster, and a biologically plausible organism annotation. At least one clinically validated target and one mechanistically novel/preclinical target should be retained when their chemical evidence is comparable, because this prevents the pipeline from collapsing all hypotheses into familiar drug targets.

For *K. pneumoniae*, *A. baumannii*, and MRSA, the prior docking/MD target list should now be reclassified into three categories: direct open-target overlap, target-family overlap, and orthogonal structural hypothesis. Only the first two categories should be described as convergent evidence; the third remains valuable but requires explicit structural justification.

## 6. Limitations

Fingerprint similarity is a ligand-space method and cannot establish binding, selectivity, target engagement, or cellular activity. ChEMBL activity records differ in assay format, organism, construct, endpoint, and data quality. Target classes with few references can be unstable, and the current benchmark is not large enough for statistically precise generalization estimates. Organism-level clinical filtering is an evidence annotation, not a bacterial susceptibility model. It does not model permeability, porins, efflux, expression, target mutation, resistance-gene context, protein abundance, or pharmacokinetics.

The next scientific improvement should be a curated, target-family-specific reference set with sequence-aware bacterial target mapping and explicit exclusion of cytotoxic/human off-target references. The next validation layer should combine ligand-based prediction with docking/MD, target-specific biochemical assays, and organism-level MIC or growth-rescue experiments.

## 7. Public references

1. [WHO bacterial priority pathogens list, 2024](https://www.who.int/publications/i/item/9789240093461)
2. [De Oliveira et al., Antimicrobial Resistance in ESKAPE Pathogens, Clinical Microbiology Reviews (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7227449/)
3. [ChEMBL database](https://www.ebi.ac.uk/chembl/)
4. [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial)
5. [RDKit documentation](https://www.rdkit.org/)

## Local output files

`results/benchmark_metrics.csv`, `results/benchmark_summary.csv`, `results/open_target_scores.csv`, `results/open_target_predictions_by_organism.csv`, `results/open_target_shortlist_by_organism.csv`, `results/figures/open_target_ecfp4_heatmap.png`, `results/figures/open_target_maccs_heatmap.png`, `results/figures/open_target_ecfp4_chemical_space.png`, `results/figures/open_target_organism_clinical_filter.png`, `results/figures/benchmark_topk_metrics.png`, and `results/figures/benchmark_target_rank_distribution.png` are generated locally and are intentionally ignored by Git until publication approval.
'''
OUT.write_text(text)
print(OUT)
