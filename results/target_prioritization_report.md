# Organism-aware target prioritization for 12 antibacterial compounds

**Analysis type:** RDKit ligand-based target prioritization with organism-panel constraints, MACCS/ECFP4 consensus, structural-feature flags, and prior docking/MD reconciliation.

> **Scope statement.** These results prioritize testable hypotheses. They do not demonstrate target engagement, biochemical inhibition, or organism-specific antibacterial mechanism.

## Executive summary

The uploaded structure set contains **12 validated compounds**. The authoritative compound-to-name mapping was reconstructed from the named records embedded in `README.md` because the filenames of the supplied SDF files did not consistently match their internal names. The resulting normalized set is BI-1, BI-6, OX-11, T2Z14, T2Z5, T2Z6, T2Z9, X1V11, X1V19, X1V20, X1V26, and X1V9.

The ligand-based analysis is strongest for the BI series against GyrB-like reference chemistry. The X1V compounds form a coherent sulfonamide-benzothiazole chemical series, but the available organism panels do not contain DHPS; DHPS was therefore retained as an orthogonal target class in the full score matrix rather than silently omitted. The T2Z/OX-11 compounds occupy a separate heteroaromatic/xanthine-like chemical-space region and show mixed Mur-ligase, GyrB, FabI, or TopoIV-related similarity depending on the organism panel.

The main scientific conclusion is conservative: **GyrB, MurC/MurA, FabI, TopoIV, LpxC, FtsZ, and the prior beta-lactamase hypotheses should be treated as a ranked hypothesis set rather than a single universal target.** The repeated GyrB ranking is partly a consequence of shared ligand-reference chemistry and must not be misreported as proof that all 12 compounds inhibit GyrB.

## 1. Data preparation and input integrity

RDKit successfully parsed and sanitized 12 unique named structures. Explicit hydrogens were removed before fingerprint calculation. The compound manifest is available at `results/compound_properties_with_series.csv`.

| Series | Compounds | Chemical interpretation |
|---|---|---|
| BI-series | BI-1, BI-6 | benzimidazole-thioacetamide |
| OX/T2Z-series | OX-11, T2Z14, T2Z5, T2Z6, T2Z9 | xanthine/heteroaromatic hydrazone |
| X1V-series | X1V11, X1V19, X1V20, X1V26, X1V9 | sulfonamide-benzothiazole hydrazone |

The supplied `final_target_list_per_bacterium.csv` was retained unchanged and compared to the ligand-scored target universe. Targets outside that universe were explicitly labelled as orthogonal prior hypotheses rather than treated as negative predictions.

## 2. Methodology

For each candidate target class, active ChEMBL reference ligands were collected from target-specific or target-family records. The data assembly retained the canonical SMILES, ChEMBL molecule identifier, assay type, standard value, derived or reported pChEMBL value, target identifier, and reference organism. Sparse classes were preserved but flagged by their reference count.

Each user compound was compared with every reference ligand using 2048-bit ECFP4/Morgan fingerprints and MACCS keys. The reported ECFP4 metrics are maximum similarity and mean similarity of the five closest reference ligands. The MACCS metric is the maximum key similarity. Class-relevant substructure flags provide a separately inspectable SAR layer.

The raw evidence score is a weighted combination of normalized ECFP4 maximum similarity, top-five ECFP4 mean similarity, MACCS similarity, and a small class-specific SAR bonus. A modest organism/target prior was then applied only inside each requested target panel. For *A. baumannii*, PBP2a was down-weighted because MecA/PBP2a is a MRSA resistance determinant and is biologically mismatched to the requested *A. baumannii* panel; this is a panel-quality correction, not a claim that all *A. baumannii* PBPs are irrelevant.

References: [RDKit](https://www.rdkit.org/), [ChEMBL](https://www.ebi.ac.uk/chembl/), and [UMAP](https://umap-learn.readthedocs.io/).

## 3. Organism-panel target scores

The table below reports mean organism-adjusted score across all 12 compounds. The companion CSV contains medians, maxima, ECFP4 and MACCS components, SAR bonuses, and reference counts.

| Organism | Rank | Target class | Mean adjusted score | Max score | Mean ECFP4 | Mean MACCS | Reference ligands |
|---|---:|---|---:|---:|---:|---:|---:|
| Klebsiella pneumoniae | 1 | TopoIV | 0.312 | 0.437 | 0.220 | 0.632 | 190 |
| Klebsiella pneumoniae | 2 | FabI | 0.269 | 0.306 | 0.206 | 0.599 | 331 |
| Klebsiella pneumoniae | 3 | LpxH | 0.256 | 0.312 | 0.174 | 0.688 | 252 |
| Klebsiella pneumoniae | 4 | MurA | 0.246 | 0.301 | 0.197 | 0.545 | 94 |
| Klebsiella pneumoniae | 5 | DHFR | 0.208 | 0.297 | 0.169 | 0.582 | 232 |
| Bacillus cereus | 1 | GyrB | 0.394 | 0.511 | 0.266 | 0.732 | 253 |
| Bacillus cereus | 2 | FabI | 0.269 | 0.306 | 0.206 | 0.599 | 331 |
| Bacillus cereus | 3 | MurA | 0.246 | 0.301 | 0.197 | 0.545 | 94 |
| Bacillus cereus | 4 | FtsZ | 0.234 | 0.345 | 0.185 | 0.515 | 67 |
| Bacillus cereus | 5 | DHFR | 0.208 | 0.297 | 0.169 | 0.582 | 232 |
| Escherichia coli | 1 | GyrB | 0.394 | 0.511 | 0.266 | 0.732 | 253 |
| Escherichia coli | 2 | LpxC | 0.274 | 0.324 | 0.207 | 0.560 | 600 |
| Escherichia coli | 3 | FabI | 0.269 | 0.306 | 0.206 | 0.599 | 331 |
| Escherichia coli | 4 | MurA | 0.246 | 0.301 | 0.197 | 0.545 | 94 |
| Escherichia coli | 5 | DHFR | 0.208 | 0.297 | 0.169 | 0.582 | 232 |
| Proteus mirabilis | 1 | GyrB | 0.394 | 0.511 | 0.266 | 0.732 | 253 |
| Proteus mirabilis | 2 | LpxC | 0.274 | 0.324 | 0.207 | 0.560 | 600 |
| Proteus mirabilis | 3 | FabI | 0.269 | 0.306 | 0.206 | 0.599 | 331 |
| Proteus mirabilis | 4 | MurA | 0.246 | 0.301 | 0.197 | 0.545 | 94 |
| Proteus mirabilis | 5 | DHFR | 0.208 | 0.297 | 0.169 | 0.582 | 232 |
| Acinetobacter baumannii | 1 | GyrB | 0.394 | 0.511 | 0.266 | 0.732 | 253 |
| Acinetobacter baumannii | 2 | MurC | 0.350 | 0.439 | 0.241 | 0.671 | 167 |
| Acinetobacter baumannii | 3 | FabI | 0.269 | 0.306 | 0.206 | 0.599 | 331 |
| Acinetobacter baumannii | 4 | LpxA | 0.146 | 0.248 | 0.148 | 0.515 | 1 |
| Acinetobacter baumannii | 5 | PBP2a | 0.026 | 0.038 | 0.126 | 0.517 | 17 |
| MRSA / Staphylococcus aureus | 1 | GyrB | 0.401 | 0.520 | 0.266 | 0.732 | 253 |
| MRSA / Staphylococcus aureus | 2 | FabI | 0.269 | 0.306 | 0.206 | 0.599 | 331 |
| MRSA / Staphylococcus aureus | 3 | FtsZ | 0.229 | 0.339 | 0.185 | 0.515 | 67 |
| MRSA / Staphylococcus aureus | 4 | DHFR | 0.208 | 0.297 | 0.169 | 0.582 | 232 |
| MRSA / Staphylococcus aureus | 5 | PBP2a | 0.146 | 0.207 | 0.126 | 0.517 | 17 |

### Interpretation by organism

**Klebsiella pneumoniae.** TopoIV is the highest mean-scoring class within the original panel, followed by FabI, MurA, and LpxH. The prior KPC-2 beta-lactamase hypotheses are not scored by the current ligand-reference universe and should remain primary orthogonal hypotheses for the OX-11/T2Z14 subset. FtsZ is a scored family but lies outside the original panel.

**Bacillus cereus.** GyrB is the highest mean-scoring class, followed by FabI, MurA, and FtsZ. FtsZ is the strongest direct overlap with the prior target CSV, while the metallo-beta-lactamase hypothesis remains outside the ligand-scored universe.

**Escherichia coli.** GyrB is the highest mean-scoring class, with LpxC and FabI next. MurA is a lower-ranked ligand hypothesis but a direct prior overlap for T2Z6 and should remain a control target because the prior structure contains the established MurA/fosfomycin context.

**Proteus mirabilis.** GyrB is highest, followed by LpxC, FabI, and MurA. The prior X1V9 beta-lactamase class-C hypothesis is not scored; MurA is a direct scored overlap and should remain the better-controlled cell-wall target.

**Acinetobacter baumannii.** GyrB and MurC are the highest ligand-scored classes. This supports adding MurC/GyrB to the next computational panel, but does not invalidate the prior OXA-23 hypothesis for X1V20 or LeuRS/FabH hypotheses for X1V11/X1V26. PBP2a is strongly down-ranked as an organism-mismatched panel entry.

**MRSA / Staphylococcus aureus.** GyrB is the strongest ligand-scored class, followed by FabI, FtsZ, and DHFR. The prior DNA-gyrase hypothesis maps directly to the GyrB family. The prior BlaZ hypothesis is not scored and should be retained as an orthogonal resistance-mechanism target; PBP2a is biologically important but shows weak ligand-reference evidence for this compound set.

## 4. Shortlisted hypotheses by compound and organism

The full top-two table is available in `results/organism_target_shortlist.csv`. A target is listed as a computational shortlist item when it ranks first or second within the organism-specific panel. The table below shows the highest-ranked target for each compound in each organism.

| Organism | Compound | Top target | Adjusted score | ECFP4 max | MACCS max | Second target | Score margin |
|---|---|---|---:|---:|---:|---|---:|
| Klebsiella pneumoniae | BI-1 | TopoIV | 0.318 | 0.229 | 0.625 | FabI | 0.057 |
| Klebsiella pneumoniae | BI-6 | TopoIV | 0.358 | 0.259 | 0.639 | LpxH | 0.080 |
| Klebsiella pneumoniae | OX-11 | FabI | 0.279 | 0.198 | 0.701 | TopoIV | 0.012 |
| Klebsiella pneumoniae | T2Z14 | FabI | 0.267 | 0.174 | 0.705 | LpxH | 0.017 |
| Klebsiella pneumoniae | T2Z5 | FabI | 0.259 | 0.176 | 0.658 | TopoIV | 0.001 |
| Klebsiella pneumoniae | T2Z6 | TopoIV | 0.259 | 0.162 | 0.689 | MurA | 0.001 |
| Klebsiella pneumoniae | T2Z9 | FabI | 0.301 | 0.200 | 0.705 | MurA | 0.000 |
| Klebsiella pneumoniae | X1V11 | TopoIV | 0.345 | 0.250 | 0.598 | FabI | 0.039 |
| Klebsiella pneumoniae | X1V19 | TopoIV | 0.308 | 0.234 | 0.544 | DHFR | 0.062 |
| Klebsiella pneumoniae | X1V20 | TopoIV | 0.325 | 0.244 | 0.573 | LpxH | 0.058 |
| Klebsiella pneumoniae | X1V26 | TopoIV | 0.437 | 0.301 | 0.689 | LpxH | 0.124 |
| Klebsiella pneumoniae | X1V9 | TopoIV | 0.368 | 0.278 | 0.544 | DHFR | 0.102 |
| Bacillus cereus | BI-1 | GyrB | 0.511 | 0.375 | 0.631 | FtsZ | 0.165 |
| Bacillus cereus | BI-6 | GyrB | 0.436 | 0.316 | 0.644 | FtsZ | 0.117 |
| Bacillus cereus | OX-11 | GyrB | 0.302 | 0.208 | 0.705 | FabI | 0.023 |
| Bacillus cereus | T2Z14 | GyrB | 0.310 | 0.200 | 0.702 | FabI | 0.043 |
| Bacillus cereus | T2Z5 | GyrB | 0.354 | 0.228 | 0.736 | FabI | 0.096 |
| Bacillus cereus | T2Z6 | GyrB | 0.353 | 0.228 | 0.736 | MurA | 0.095 |
| Bacillus cereus | T2Z9 | GyrB | 0.308 | 0.200 | 0.702 | FabI | 0.007 |
| Bacillus cereus | X1V11 | GyrB | 0.425 | 0.282 | 0.742 | FabI | 0.119 |
| Bacillus cereus | X1V19 | GyrB | 0.423 | 0.282 | 0.814 | DHFR | 0.178 |
| Bacillus cereus | X1V20 | GyrB | 0.423 | 0.282 | 0.831 | FabI | 0.165 |
| Bacillus cereus | X1V26 | GyrB | 0.433 | 0.300 | 0.690 | FabI | 0.138 |
| Bacillus cereus | X1V9 | GyrB | 0.445 | 0.294 | 0.845 | DHFR | 0.178 |
| Escherichia coli | BI-1 | GyrB | 0.511 | 0.375 | 0.631 | LpxC | 0.248 |
| Escherichia coli | BI-6 | GyrB | 0.436 | 0.316 | 0.644 | LpxC | 0.137 |
| Escherichia coli | OX-11 | GyrB | 0.302 | 0.208 | 0.705 | LpxC | 0.016 |
| Escherichia coli | T2Z14 | GyrB | 0.310 | 0.200 | 0.702 | LpxC | 0.011 |
| Escherichia coli | T2Z5 | GyrB | 0.354 | 0.228 | 0.736 | LpxC | 0.030 |
| Escherichia coli | T2Z6 | GyrB | 0.353 | 0.228 | 0.736 | LpxC | 0.065 |
| Escherichia coli | T2Z9 | GyrB | 0.308 | 0.200 | 0.702 | FabI | 0.007 |
| Escherichia coli | X1V11 | GyrB | 0.425 | 0.282 | 0.742 | FabI | 0.119 |
| Escherichia coli | X1V19 | GyrB | 0.423 | 0.282 | 0.814 | LpxC | 0.160 |
| Escherichia coli | X1V20 | GyrB | 0.423 | 0.282 | 0.831 | FabI | 0.165 |
| Escherichia coli | X1V26 | GyrB | 0.433 | 0.300 | 0.690 | FabI | 0.138 |
| Escherichia coli | X1V9 | GyrB | 0.445 | 0.294 | 0.845 | DHFR | 0.178 |
| Proteus mirabilis | BI-1 | GyrB | 0.511 | 0.375 | 0.631 | LpxC | 0.248 |
| Proteus mirabilis | BI-6 | GyrB | 0.436 | 0.316 | 0.644 | LpxC | 0.137 |
| Proteus mirabilis | OX-11 | GyrB | 0.302 | 0.208 | 0.705 | LpxC | 0.016 |
| Proteus mirabilis | T2Z14 | GyrB | 0.310 | 0.200 | 0.702 | LpxC | 0.011 |
| Proteus mirabilis | T2Z5 | GyrB | 0.354 | 0.228 | 0.736 | LpxC | 0.030 |
| Proteus mirabilis | T2Z6 | GyrB | 0.353 | 0.228 | 0.736 | LpxC | 0.065 |
| Proteus mirabilis | T2Z9 | GyrB | 0.308 | 0.200 | 0.702 | FabI | 0.007 |
| Proteus mirabilis | X1V11 | GyrB | 0.425 | 0.282 | 0.742 | FabI | 0.119 |
| Proteus mirabilis | X1V19 | GyrB | 0.423 | 0.282 | 0.814 | LpxC | 0.160 |
| Proteus mirabilis | X1V20 | GyrB | 0.423 | 0.282 | 0.831 | FabI | 0.165 |
| Proteus mirabilis | X1V26 | GyrB | 0.433 | 0.300 | 0.690 | FabI | 0.138 |
| Proteus mirabilis | X1V9 | GyrB | 0.445 | 0.294 | 0.845 | DHFR | 0.178 |
| Acinetobacter baumannii | BI-1 | GyrB | 0.511 | 0.375 | 0.631 | MurC | 0.190 |
| Acinetobacter baumannii | BI-6 | GyrB | 0.436 | 0.316 | 0.644 | MurC | 0.156 |
| Acinetobacter baumannii | OX-11 | MurC | 0.439 | 0.290 | 0.753 | GyrB | 0.138 |
| Acinetobacter baumannii | T2Z14 | MurC | 0.384 | 0.241 | 0.744 | GyrB | 0.074 |
| Acinetobacter baumannii | T2Z5 | MurC | 0.370 | 0.235 | 0.707 | GyrB | 0.015 |
| Acinetobacter baumannii | T2Z6 | MurC | 0.411 | 0.260 | 0.728 | GyrB | 0.058 |
| Acinetobacter baumannii | T2Z9 | MurC | 0.420 | 0.263 | 0.765 | GyrB | 0.112 |
| Acinetobacter baumannii | X1V11 | GyrB | 0.425 | 0.282 | 0.742 | FabI | 0.119 |
| Acinetobacter baumannii | X1V19 | GyrB | 0.423 | 0.282 | 0.814 | MurC | 0.058 |
| Acinetobacter baumannii | X1V20 | GyrB | 0.423 | 0.282 | 0.831 | MurC | 0.140 |
| Acinetobacter baumannii | X1V26 | GyrB | 0.433 | 0.300 | 0.690 | MurC | 0.083 |
| Acinetobacter baumannii | X1V9 | GyrB | 0.445 | 0.294 | 0.845 | MurC | 0.155 |
| MRSA / Staphylococcus aureus | BI-1 | GyrB | 0.520 | 0.375 | 0.631 | FtsZ | 0.182 |
| MRSA / Staphylococcus aureus | BI-6 | GyrB | 0.445 | 0.316 | 0.644 | FtsZ | 0.132 |
| MRSA / Staphylococcus aureus | OX-11 | GyrB | 0.307 | 0.208 | 0.705 | FabI | 0.029 |
| MRSA / Staphylococcus aureus | T2Z14 | GyrB | 0.316 | 0.200 | 0.702 | FabI | 0.049 |
| MRSA / Staphylococcus aureus | T2Z5 | GyrB | 0.361 | 0.228 | 0.736 | FabI | 0.102 |
| MRSA / Staphylococcus aureus | T2Z6 | GyrB | 0.360 | 0.228 | 0.736 | FabI | 0.104 |
| MRSA / Staphylococcus aureus | T2Z9 | GyrB | 0.314 | 0.200 | 0.702 | FabI | 0.013 |
| MRSA / Staphylococcus aureus | X1V11 | GyrB | 0.433 | 0.282 | 0.742 | FabI | 0.128 |
| MRSA / Staphylococcus aureus | X1V19 | GyrB | 0.432 | 0.282 | 0.814 | DHFR | 0.186 |
| MRSA / Staphylococcus aureus | X1V20 | GyrB | 0.432 | 0.282 | 0.831 | FabI | 0.173 |
| MRSA / Staphylococcus aureus | X1V26 | GyrB | 0.442 | 0.300 | 0.690 | FabI | 0.147 |
| MRSA / Staphylococcus aureus | X1V9 | GyrB | 0.453 | 0.294 | 0.845 | DHFR | 0.187 |

## 5. Validation against prior docking and molecular dynamics targets

The user-provided prior target list contains 16 rows. The current reconciliation identifies four direct scored overlaps, three scored-family hypotheses outside the original organism panels, and eight hypotheses that are not scored because the current ligand-reference universe does not include beta-lactamases or LeuRS.

| Validation category | Count | Meaning |
|---|---:|---|
| Direct scored overlap | 4 | Prior target maps to a class scored inside that organism panel |
| Scored family outside original panel | 3 | Prior target family is represented, but not in the original requested panel |
| Orthogonal unscored hypothesis | 8 | Dedicated docking/MD hypothesis retained; no ligand-based negative claim |

Direct overlaps are *B. cereus*–FtsZ, *E. coli*–MurA, *P. mirabilis*–MurA, and MRSA–DNA gyrase mapped to GyrB. The three outside-panel scored families are *K. pneumoniae*–FtsZ, *E. coli*–FtsZ, and *A. baumannii*–FabH. The unscored set includes KPC-2, KPC-2–avibactam, AmpC, class-C beta-lactamase, OXA-23, BlaZ, and LeuRS hypotheses.

This validation pattern is scientifically useful: it shows where the ligand-only model agrees with existing structural hypotheses, where the old target panel was incomplete, and where an orthogonal protein-level hypothesis cannot be assessed by the present ligand database.

## 6. Figures and what they show

| Figure | Scientific purpose |
|---|---|
| `results/figures/compound_structure_grid.png` | Visual audit of the 12 input structures and naming. |
| `results/figures/compound_maccs_similarity_heatmap.png` | Pairwise MACCS-key similarity; shows series coherence and chemical redundancy. |
| `results/figures/compound_target_ecfp4_heatmap.png` | Maximum ECFP4 similarity to active reference ligands by target class. |
| `results/figures/compound_target_maccs_heatmap.png` | Orthogonal MACCS-key target-class similarity heatmap. |
| `results/figures/tmap_like_ecfp4_reference_map.png` | UMAP/Jaccard chemical-space map of user compounds and representative reference ligands; labelled TMAP-like, not exact TMAP. |
| `results/figures/organism_panel_target_rankings.png` | Top two organism-panel hypotheses per compound with an explicit moderate-evidence guide line. |
| `results/figures/prior_docking_md_overlap.png` | Counts of prior target assignments that overlap scored families. |
| `results/figures/physchem_permeability_space.png` | MW/logP/TPSA space and a screening-level Gram-negative permeability caution. |

## 7. Limitations and recommended next experiments

The analysis is limited by the quality and coverage of public ligand-reference data. A high fingerprint score is not a calibrated binding probability. Cross-organism ChEMBL records can blur species specificity, whole-complex measurements can blur subunit specificity, and reference-set size is highly unequal. LpxA has only one retained reference ligand and PBP2a has a small set, so rankings for these classes are low-confidence.

The current workflow does not model protein sequence identity, binding-site conservation, active-site protonation, cofactors, membrane permeability, efflux, target expression, or cellular resistance. The Gram-negative TPSA flag is only a screening caution and is not a permeability predictor.

For the next round, retain three hypotheses per compound: the top ligand-based target in the organism panel, the strongest prior docking/MD target even when unscored, and one mechanistic control. Use sequence-matched proteins for GyrB, MurA/MurC, FabI, FtsZ, and LpxC. For KPC-2, OXA-23, AmpC, BlaZ, and class-C beta-lactamases, maintain organism-matched catalytic-state docking and validate with beta-lactamase inhibition assays. For X1V compounds, add DHPS and beta-lactamase family scoring when a curated reference ligand set is available. For MRSA PBP2a, add a target-specific reference set and verify whether the X1V series contains a chemically plausible covalent or noncovalent PBP2a motif before investing in further MD.

## 8. Reproducibility

Run `python run_pipeline.py` after placing the input files in `inputs/`. The generated tables, figures, and report are written to `results/`. The code is intentionally modular: preparation, scoring, organism ranking, figure generation, result summarization, and prior-target comparison are separate scripts.

### References

1. [RDKit documentation](https://www.rdkit.org/).
2. [ChEMBL database and API](https://www.ebi.ac.uk/chembl/).
3. [UMAP documentation](https://umap-learn.readthedocs.io/).
4. Rogers D, Brown RD. Extended-connectivity fingerprints and molecular similarity; the ECFP/Morgan representation is used here as a ligand-space descriptor.
