# v2 evidence schema and quality gates

## Objective

The v2 pipeline must distinguish what is supported by ligand chemistry from what is biologically desirable and from what is experimentally established. No composite score may be interpreted as a probability of target engagement.

## Evidence layers

| Layer | Fields | Interpretation |
|---|---|---|
| Chemical | ECFP4 maximum, ECFP4 top-five mean, MACCS maximum, scaffold distance, pharmacophore/substructure flags | Similarity to reference ligands of one target family |
| Reference quality | Number of unique references, assay types, target confidence, organism/construct, pChEMBL distribution | Stability and transferability of chemical evidence |
| Negative evidence | Inactive target-tested molecules, unrelated-target activity, decoy similarity, reference-set coverage | Whether evidence is target-selective rather than generic chemical similarity |
| Species mapping | Organism protein accession, sequence identity, active-site conservation, target complex, orthologue status | Whether the reference target maps to the specific organism |
| Pocket compatibility | PDB/AlphaFold structure, binding-site residues, pocket similarity, ligandable-site status | Whether a plausible binding site exists in the organism-specific protein |
| Biology | Essentiality/fitness, growth phenotype, expression, pathway redundancy | Biological desirability of target engagement |
| Cellular access | Localization, outer-membrane exposure, porin/efflux risk, compound polarity and ionization | Probability that a compound can reach the target in cells |
| Resistance | Resistance-gene prevalence, target mutation frequency, bypass pathways | Risk of rapid loss of activity or organism-specific failure |
| Safety | Human orthologue similarity, human off-target similarity, cytotoxicity alerts | Selectivity and translational risk |
| Validation | Biochemical assay, genetic rescue, resistant mutant, MIC/time-kill, target engagement | Experimental confirmation status |

## Target-quality gates

A target family is eligible for a **moderate-confidence chemical prediction** only when it has at least 30 unique reference ligands, at least one direct biochemical or binding assay layer, and a target annotation that maps to a defined bacterial protein or complex. Classes with fewer than 30 references can still appear as exploratory hypotheses but must carry a low-confidence label and a reference-sparsity penalty.

Phenotypic mechanisms such as general membrane disruption, efflux, and porin loss are not treated as protein targets unless a mechanism-specific reference set is available. Resistance targets such as beta-lactamases and PBP2a are retained but are labelled separately from direct viability targets.

## Composite outputs

The pipeline produces four independent outputs:

1. `chemical_evidence_score`: ligand-space evidence only.
2. `species_transfer_score`: sequence and pocket compatibility for the organism-specific target.
3. `biological_priority_score`: essentiality, clinical validation, accessibility, and resistance annotations.
4. `overall_priority_score`: a transparent decision-support product, never a binding probability.

Every result also receives a confidence class:

- **High:** strong chemical evidence, adequate reference support, species mapping available, and no major accessibility contradiction.
- **Moderate:** chemically supported but one of reference size, species mapping, or accessibility evidence is incomplete.
- **Low:** sparse reference set, weak similarity, unresolved orthologue mapping, or strong cellular-access concern.
- **Insufficient:** no appropriate reference universe or mechanism-specific evidence.

## Benchmark requirements

Benchmark evaluation must exclude the query and close analogues, report scaffold and target-family splits, retain negative controls, and compare against random, prevalence, ECFP4-only, and MACCS-only baselines. Results must be reported with confidence intervals when the sample size permits them.
