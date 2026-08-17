# V2 Method Supplement

## 1. Molecular representation

Input molecules are read with RDKit, sanitized, converted to canonical molecular graphs, and stored as canonical SMILES. ECFP4/Morgan fingerprints use radius 2 and 2048 bits. MACCS keys are calculated using RDKit’s standard MACCS implementation. Explicit hydrogens are removed before fingerprinting. Invalid structures are logged and excluded rather than silently repaired.

## 2. Reference evidence

Reference ligands are grouped by `target_class`. Duplicate canonical structures are collapsed within a target class. The quality layer records reference count, unique Murcko-scaffold count, assay/evidence grade, and the number of records remaining after leakage exclusion. Sparse target classes are quality-penalized and labelled as limited coverage.

For a query molecule q and target class t, the primary similarity terms are:

```text
E_ecfp4(q,t) = max_i Tanimoto(ECFP4(q), ECFP4(r_i,t))
E_maccs(q,t) = max_i Tanimoto(MACCS(q), MACCS(r_i,t))
E_mean5(q,t) = mean of the five largest ECFP4 similarities
```

Cross-target decoys are drawn from reference molecules associated with other target classes. The specificity margin is:

```text
margin(q,t) = positive_target_score(q,t) - max_decoy_score(q,t)
```

The target-specificity score maps this margin to a bounded 0–1 value. The quality-adjusted chemical score combines similarity, target specificity, reference quality, and reference coverage. All component values remain in the output table.

## 3. Species transfer

Each target class is mapped to an organism-specific protein when possible. The mapping table records organism, taxonomy identifier, target class, protein accession, sequence, mapping method, and mapping status. Sequence transfer is calculated from the mapped protein relative to the reference target record or within-class reference sequence set. Identity and alignment coverage are retained separately; the transfer score is not treated as an affinity estimate.

An unresolved mapping receives a transfer score of zero and an explicit uncertainty reason. It is not interpreted as target absence.

## 4. Biological and safety layers

Clinical priority, essentiality/fitness, cellular accessibility, and resistance relevance are ontology annotations. They represent desirability and context, not direct molecular evidence. Resistance-target roles are distinguished from direct viability targets. Anti-target fields are qualitative early-alert annotations; they do not constitute an in vitro or in vivo safety result.

## 5. CARD and structural evidence

CARD model counts and SNP counts are used as resistance context. Organism-specific SNP counts are counted only when an organism description can be extracted from the public record. RCSB records are used to annotate candidate structures, experimental method, and the presence of non-solvent co-crystallized components. A co-crystallized ligand is interpreted as active-site precedent for later user-controlled docking. No docking, PDBQT preparation, or molecular dynamics is executed by this v2 workflow.

## 6. Benchmark leakage controls

The public benchmark excludes a reference ligand when ECFP4 similarity to the query is at least 0.85 in the close-analogue split. In the scaffold split, a reference ligand is excluded when the Bemis–Murcko scaffold string matches exactly. Performance is calculated only for queries whose target label appears in the reference-supported candidate universe. Uncovered labels are reported as coverage limitations.

The evaluation reports top-1, top-3, top-5 retrieval, reciprocal rank, random baselines, prevalence baseline, and enrichment over random. Because the current benchmark is small, target-family holdout, temporal, and species-holdout results are reported as unavailable or partial when their prerequisites are absent.

## 7. Calibration and confidence

Bootstrap resampling of target reference ligands estimates score stability. Decoy calibration estimates how frequently an observed target score exceeds cross-target scores. Confidence classes combine chemical quality, specificity, reference coverage, mapping status, species transfer, and bootstrap/decoy stability:

| Class | Meaning |
|---|---|
| Moderate | Initial direct testing is justified, but the prediction is not target confirmation |
| Low | A plausible lower-priority hypothesis with important uncertainty |
| Insufficient | New reference, sequence, structural, or experimental evidence is required before prioritization |

No confidence class is interpreted as a calibrated probability of binding.

## 8. Validation planning

The recommended validation sequence is: purified target or target-complex assay; species-orthologue comparison; MIC/time-kill with permeability and efflux controls; resistant-mutant selection and sequencing; complementation or target rescue; and human-orthologue/mitochondrial selectivity where relevant. Different ontology classes receive different first assays: PBP2a requires transpeptidase acylation/inhibition, beta-lactamase requires enzyme inhibition against resistance variants, ribosome targets require translation or complex-binding assays, and envelope targets require permeability and whole-cell mechanism tests.

## 9. Data protection

Private compound structures, names, and compound-specific results are stored only in ignored local paths. Public commits include code, public benchmark data, public reference data, and documentation. A clean-room clone without authorized local inputs must not be expected to reproduce private predictions.
