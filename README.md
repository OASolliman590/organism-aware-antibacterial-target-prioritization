# Organism-aware ligand-based antibacterial target prioritization

This repository implements an auditable cheminformatics workflow for prioritizing antibacterial target hypotheses for 12 user-provided compounds across six bacterial organisms: *Bacillus cereus*, *Klebsiella pneumoniae*, *Escherichia coli*, *Proteus mirabilis*, *Acinetobacter baumannii*, and MRSA/*Staphylococcus aureus*.

> **Interpretation:** the output is a ligand-based target-prioritization analysis, not proof of target engagement. It is designed to improve target-panel selection before docking, molecular dynamics, biochemical assays, or cellular validation.

## Scientific design

The workflow combines four evidence layers. First, RDKit ECFP4 fingerprints are used to calculate the maximum Tanimoto similarity and the mean similarity of the five closest active reference ligands for each candidate target class. Second, MACCS-key similarity is calculated as an interpretable complementary fingerprint rather than treating one fingerprint as definitive. Third, class-relevant structural flags are recorded, including sulfonamides, beta-lactams, hydroxamates, diaminopyrimidine-like motifs, xanthine-like motifs, hydrazones, benzimidazoles, benzothiazoles, and aryl halides. Fourth, the raw ligand evidence is adjusted by a small, explicit organism/target prior. These priors are deliberately modest and are not probabilities.

The target reference sets were assembled from ChEMBL activity records during the analysis. The repository preserves the downloaded JSON reference sets used by the scoring run. Because the public data are uneven across target classes, the number of references and the best-reference assay context are reported for every score. LpxA, PBP2a, and some Mur-family classes have sparse or heterogeneous public data; they must therefore be interpreted more cautiously.

The principal evidence score is:

```text
raw_ligand_evidence = 0.50 * normalized_max_ECFP4
                    + 0.25 * normalized_top5_mean_ECFP4
                    + 0.15 * normalized_max_MACCS
                    + SAR_bonus

organism_adjusted_score = raw_ligand_evidence * organism_target_prior
```

The normalized fingerprint components are transparent affine rescalings for ranking, not calibrated probabilities. A target is labelled high, moderate, or low evidence using the raw ligand score, while organism-adjusted rankings are used only within the requested organism panel.

## Current target panels

| Organism | Panel used for ranking |
|---|---|
| *K. pneumoniae* | DHFR, FabI, TopoIV, MurA, LpxH |
| *B. cereus* | FabI, DHFR, FtsZ, GyrB, MurA |
| *E. coli* | MurA, DHFR, FabI, GyrB, LpxC |
| *P. mirabilis* | MurA, DHFR, FabI, GyrB, LpxC |
| *A. baumannii* | PBP2a, FabI, GyrB, LpxA, MurC |
| MRSA/*S. aureus* | PBP2a, DHFR, FabI, FtsZ, GyrB |

The current panel is kept separate from the prior docking/MD CSV. This distinction is important because several prior hypotheses—KPC-2, OXA-23, BlaZ, AmpC, and LeuRS—are not represented by the current ChEMBL ligand-reference universe. They are retained as orthogonal hypotheses rather than silently treated as negative results.

## Input files

Place the following files in `inputs/` before running the workflow:

* `README.md` containing the authoritative named MOL records;
* `final_target_list_per_bacterium.csv` containing prior docking/MD target assignments;
* the available SDF files and/or `all_compounds_combined.sdf`.

The preparation script intentionally uses the named records in `README.md` as the authoritative name-to-structure mapping because the uploaded SDF filenames did not consistently match their internal `_Name` fields. It writes a normalized 12-compound SDF and manifest.

## Execution

```bash
python -m pip install -r requirements.txt
python pipeline/prepare_compounds.py
python pipeline/target_scoring.py
python pipeline/rank_and_figures.py
python pipeline/summarize_results.py
python pipeline/prior_comparison.py
```

The same workflow can be run from a clone with:

```bash
export PROJECT_ROOT="$PWD"
export INPUT_DIR="$PWD/inputs"
python pipeline/prepare_compounds.py
python pipeline/target_scoring.py
python pipeline/rank_and_figures.py
python pipeline/summarize_results.py
python pipeline/prior_comparison.py
```

## Main outputs

| Output | Description |
|---|---|
| `results/organism_target_shortlist.csv` | Top two target hypotheses per compound within each organism panel, with justification text |
| `results/organism_panel_rankings.csv` | Full organism-panel rankings, raw evidence, organism prior, and adjusted score |
| `results/ranked_target_evidence.csv` | Full compound-by-target-class ligand evidence table |
| `results/organism_target_mean_scores.csv` | Mean, median, maximum, and component scores by organism and target class |
| `results/prior_target_comparison.csv` | Reconciliation of prior docking/MD targets with the scored target universe and original panels |
| `results/compound_target_ecfp4_heatmap.png` | Maximum ECFP4 similarity heatmap |
| `results/compound_target_maccs_heatmap.png` | Maximum MACCS-key similarity heatmap |
| `results/compound_maccs_similarity_heatmap.png` | Pairwise MACCS similarity of the user compounds |
| `results/tmap_like_ecfp4_reference_map.png` | UMAP/Jaccard map of user compounds and representative reference ligands |
| `results/organism_panel_target_rankings.png` | Organism-adjusted top-two target rankings |
| `results/prior_docking_md_overlap.png` | Prior assignment overlap classification |
| `results/physchem_permeability_space.png` | MW/logP/TPSA space with a Gram-negative permeability caution flag |

## Interpretation of the current analysis

The ligand-only component is strongest for the benzimidazole-thioacetamide BI series against the GyrB reference set, and for several T2Z/OX-11 compounds against Mur-ligase or GyrB-related reference chemistry. The X1V sulfonamide-benzothiazole series is chemically coherent in the map and is not automatically assigned to DHPS because the current organism panels do not include DHPS; DHPS is included as an orthogonal rescue target in the full target-class score table because its sulfonamide/PABA rationale is mechanistically relevant.

The current panel ranking should not be interpreted as a universal claim that GyrB is the target for every organism. The repeated GyrB ranking reflects the composition of the available reference sets and the fact that all organisms are evaluated against the same compound set. The organism prior changes rank ordering modestly, but it cannot correct a target universe that omits beta-lactamases or LeuRS. The prior CSV therefore remains important: beta-lactamase and LeuRS hypotheses should be tested by dedicated protein-level docking, sequence-appropriate structures, and biochemical assays rather than discarded from the basis of this ligand-only screen.

## Validation against prior docking/MD work

The prior-target comparison currently identifies four direct scored overlaps: *B. cereus*–FtsZ, *E. coli*–MurA, *P. mirabilis*–MurA, and MRSA–DNA gyrase mapped to the GyrB class. Three additional prior hypotheses belong to scored families but fall outside the original organism panels: *K. pneumoniae*–FtsZ, *E. coli*–FtsZ, and *A. baumannii*–FabH. Eight prior hypotheses are not scored because they are beta-lactamase or LeuRS hypotheses. Those eight should be treated as orthogonal validation targets, not failed predictions.

## Limitations and next steps

This workflow does not model protein sequence divergence, target expression, cellular uptake, efflux, permeability, resistance determinants, or polypharmacology. ChEMBL assays are heterogeneous, and target-family records can contain cross-organism or whole-complex measurements. The public ligand sets are also imbalanced: some target classes have hundreds of references, while LpxA has only one retained reference and PBP2a has a small set. The map is TMAP-like rather than an exact TMAP implementation: it uses UMAP with Jaccard distance on ECFP4 bit vectors when UMAP is available.

The recommended next experimental/computational stage is to retain two or three hypotheses per compound: one high-scoring ligand-based target, one biologically important prior docking/MD target, and one orthogonal mechanistic control. For the beta-lactamase hypotheses, use organism-matched proteins and consistent catalytic-state preparation. For the GyrB/MurC/FabI hypotheses, use sequence-matched proteins and repeat docking with protonation, cofactors, and conserved-water treatment standardized across targets. Finally, test target engagement with purified-enzyme inhibition or thermal-shift assays before interpreting cellular antibacterial activity as target-specific.

## References

1. RDKit: cheminformatics and machine-learning software. <https://www.rdkit.org/>
2. ChEMBL database and web services. <https://www.ebi.ac.uk/chembl/>
3. Morgan/ECFP fingerprints are used here for similarity ranking; MACCS keys provide an orthogonal structural-key representation.
4. UMAP is used only for visualization of fingerprint chemical space. It is not used as a predictive model.
