# Scientific Reporting and Claim-Control Checklist

Use this checklist when converting pipeline outputs into a manuscript, thesis, preprint, conference abstract, or supplementary methods file.

## 1. Study description

- [ ] Describe the workflow as **open chemical target discovery followed by organism-aware prioritization**.
- [ ] State that the initial chemical search is not restricted to a compound-specific docking panel.
- [ ] Describe the configured organisms explicitly rather than calling the panel the canonical ESKAPE set.
- [ ] State that the output is a ranked set of experimentally testable target hypotheses.

## 2. Claims that are supported

The following phrasing is scientifically consistent with the implementation:

- chemical similarity is **consistent with / supports compatibility with** a target reference-ligand space;
- a target is **prioritized** for experimental evaluation;
- sequence evidence supports **species-transfer plausibility**;
- CARD annotations provide **resistance context**;
- RCSB annotations provide **structural or co-crystal precedent**;
- an applicability-domain flag indicates whether the compound lies near the available reference chemical space;
- a benchmark metric measures retrieval performance under the specified leakage controls.

## 3. Claims that must not be made from this pipeline alone

Do not report any of the following without independent experimental or appropriately calibrated evidence:

- "the compound binds target X";
- "target X is the mechanism of action";
- "the compound inhibits target X";
- "the compound is active against organism Y";
- "the score is the probability of target engagement";
- "the score predicts MIC";
- "CARD proves resistance to the compound";
- "an RCSB structure proves the proposed binding pocket";
- "a sequence-transfer score proves equal binding affinity between orthologues";
- "High confidence" as a statistical confidence level.

The default confidence classes and final priority scores are heuristic categories.

## 4. Mandatory chemical-method reporting

Report all of the following:

- [ ] RDKit version.
- [ ] Morgan/ECFP4 radius = 2.
- [ ] Morgan fingerprint length = 2,048 bits.
- [ ] MACCS keys used as an independent 2D representation.
- [ ] Tanimoto similarity.
- [ ] v2 top-k reference aggregation (`top_k = 5`).
- [ ] Exact v2 component normalization and weights if v2 results are reported.
- [ ] Cross-target specificity margin and transform if quality-adjusted scores are reported.

If reporting v3, also report:

- [ ] ETKDGv3 conformer generation.
- [ ] 30 requested conformers.
- [ ] RMS pruning threshold = 0.5 Å.
- [ ] MMFF94 optimization.
- [ ] maximum optimization iterations = 500.
- [ ] retained energy window = 10 kcal/mol.
- [ ] one RDKit thread per molecule for deterministic embedding.
- [ ] USRCAT.
- [ ] O3A shortlist = top 25 USRCAT references per target class.
- [ ] O3A shape similarity as `1 - ShapeTanimotoDist`.
- [ ] O3A color evidence.
- [ ] Gobbi/Poppe 2D pharmacophore similarity.
- [ ] BaseFeatures/O3A-aligned 3D pharmacophore overlap.
- [ ] 3D pharmacophore radius = 2.5 Å and width = 1.0.
- [ ] seven-component RRF fusion.
- [ ] RRF constant `k = 60`.
- [ ] missing component measurements are not imputed.

## 5. Reference-ligand reporting

- [ ] Identify ChEMBL as the source of public target-associated reference ligands.
- [ ] Report that within-target duplicate structures are removed after RDKit canonicalization.
- [ ] Report reference count and scaffold diversity where target-level conclusions depend on them.
- [ ] Treat cross-target ligands as specificity controls rather than experimentally inactive compounds.
- [ ] Do not infer the ChEMBL source release for the frozen legacy snapshot; it is unrecorded.

## 6. Species-transfer reporting

If species-aware results are used, report:

- [ ] UniProt-derived mapped protein sequences.
- [ ] deterministic global Needleman-Wunsch alignment.
- [ ] match = +1, mismatch = -1, gap = -1.
- [ ] identity calculated over aligned non-gap residue pairs.
- [ ] coverage relative to the shorter sequence.
- [ ] `species_transfer_score = identity * sqrt(coverage)`.
- [ ] unmapped targets are explicitly marked unmapped.

Do not describe this score as a binding-affinity predictor.

## 7. Biological-priority reporting

If the final organism-aware ranking is used, report the component weights rather than only saying "biological context was integrated".

Current biological-priority function:

```text
0.25 organism_scope
+ 0.20 clinical_priority
+ 0.20 essentiality
+ 0.15 cellular_accessibility
+ 0.10 resistance_relevance
+ 0.10 CARD_context
```

Current final prioritization function:

```text
chemical_quality
* (0.65 + 0.35 * species_transfer)
* (0.75 + 0.25 * pocket_evidence)
* (0.50 + 0.50 * biological_priority)
* (1 - 0.20 * anti_target_risk)
```

- [ ] Explicitly call this a heuristic decision-support function.
- [ ] Do not call it trained, calibrated, or probabilistic.

## 8. Resistance reporting

- [ ] Identify CARD release 4.0.2 for the frozen snapshot.
- [ ] Describe CARD model and SNP evidence as contextual resistance annotations.
- [ ] Do not interpret CARD record counts as resistance prevalence.
- [ ] Do not interpret CARD annotations as evidence of compound-target binding.

## 9. Structure reporting

- [ ] Identify RCSB PDB as the structure source.
- [ ] State whether a structure candidate and/or co-crystallized ligand was found.
- [ ] Distinguish structure precedent from compound-specific docking evidence.
- [ ] Do not state that docking or MD was performed by this repository unless a separate downstream analysis was actually run and reported.

## 10. Applicability-domain reporting

If prediction tables are reported, include the reference-domain status:

```text
in domain:     ECFP4 Tanimoto >= 0.40
near domain:   0.25 <= ECFP4 Tanimoto < 0.40
out of domain: ECFP4 Tanimoto < 0.25
unassessable:  nearest-reference Tanimoto missing
```

- [ ] Report nearest-reference USRCAT as a continuous additional 3D domain measure when available.
- [ ] Do not invent a USRCAT cutoff; none is calibrated in the current method.

## 11. Benchmark reporting

For v3 benchmark claims, report:

- [ ] close-analogue exclusion at ECFP4 Tanimoto >= 0.85 for every split;
- [ ] target-family holdout rules;
- [ ] exact non-empty Bemis-Murcko scaffold exclusion;
- [ ] temporal cutoff = 2018-01-01;
- [ ] temporal results are pending when date metadata are unavailable;
- [ ] AUROC;
- [ ] BEDROC alpha = 20 and 80.5;
- [ ] enrichment factor at 1% and 5%;
- [ ] MRR;
- [ ] coverage;
- [ ] tie-aware ranking;
- [ ] query-level bootstrap, `n = 1000`;
- [ ] 95% confidence intervals;
- [ ] number of total and evaluable queries.

A split with insufficient metadata must be reported as unavailable/pending, not as zero performance.

## 12. Model-calibration reporting

The repository contains a regularized logistic-regression/Platt-calibration path, but the default frozen dataset does not provide the required non-overlapping train/calibration/test roles.

- [ ] Do not report calibrated probabilities from the default snapshot.
- [ ] If a future calibrated model is fitted, report the exact training, calibration, and test query partitions.
- [ ] Confirm zero query-ID overlap among those partitions.
- [ ] Report feature-complete row exclusions; selected-feature missing values are not imputed.
- [ ] Report AUROC, Brier score, and reliability analysis on held-out data.

## 13. Reproducibility reporting

Every publication-quality run should archive or provide:

- [ ] Git commit SHA.
- [ ] configuration file.
- [ ] configuration SHA-256.
- [ ] data snapshot identifier.
- [ ] snapshot manifest and hashes.
- [ ] random seeds.
- [ ] run manifest.
- [ ] Python version.
- [ ] package versions / `requirements.lock`.
- [ ] benchmark split provenance.
- [ ] applicability-domain results.
- [ ] explicit missing/unavailable statuses.

For unpublished compound structures, preserve the repository data-protection design: private structures and compound-specific outputs should remain in ignored local paths and should not be committed to the public repository unless disclosure is authorized.

## 14. Known provenance gaps that should be disclosed

For the frozen `v2-public-baseline-2ed4684` snapshot:

- ChEMBL source release: unrecorded.
- ChEMBL acquisition date: unrecorded.
- PubChem benchmark query date: unrecorded.
- UniProt release/query date: unrecorded.
- RCSB query date: unrecorded.
- CARD release: 4.0.2 recorded.

Do not infer missing dates or releases from repository commit history.

## 15. Minimum Results-table columns recommended for publication

For each reported compound-organism-target hypothesis, retain at minimum:

- query/compound identifier;
- organism;
- target class and parent target class;
- chemical evidence score;
- reference quality grade;
- target-specificity score;
- chemical quality-adjusted score;
- nearest-reference ECFP4 similarity;
- applicability-domain flag;
- species-transfer score and mapping status;
- biological-priority score;
- pocket evidence;
- anti-target risk annotation;
- overall priority score;
- confidence class;
- uncertainty reasons; and
- recommended validation.

For v3, also retain the seven component scores, component coverage, RRF ranks/contributions, and fused-score probability flag.

## 16. Recommended validation hierarchy

A computational target hypothesis should ideally be advanced through:

1. purified-target binding or inhibition, or a target-complex biochemical assay;
2. species-orthologue comparison;
3. whole-cell susceptibility with permeability and efflux controls;
4. resistant-mutant selection and target sequencing;
5. complementation, target rescue, or genetic dependence where feasible;
6. selectivity testing against relevant human homologues or mitochondrial liabilities; and
7. downstream docking/MD only when used as structural support rather than as proof of mechanism.

The manuscript should distinguish clearly between completed validation and future recommended experiments.