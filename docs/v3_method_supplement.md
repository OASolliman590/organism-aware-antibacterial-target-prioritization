# V3 Method Supplement

## Scope and interpretation

V3 is additive to the frozen v2 target-discovery path. The golden public v2
outputs are tested before every release. V3 adds conformer-based chemical
evidence, rank fusion, leakage-controlled evaluation, sensitivity analysis, and
applicability-domain annotations. It does not turn a similarity score into a
binding probability. Only an output explicitly named `calibrated_probability`,
produced by a held-out Platt calibration, may be interpreted as a probability.

Every evidence layer remains a separate column. Missing values mean that the
measurement was unavailable; they are not replaced with zeros, averages, or
synthetic observations.

## Reproducibility and provenance

The default command is:

```bash
python run_pipeline.py --config config.yaml
```

The command verifies `data/snapshots/SNAPSHOT_VERSIONS.json` before analysis and
writes `results/run_manifest.json`. The manifest includes the immutable snapshot
identifier and hashes, full config hash, Git commit, branch/dirty state, random
seeds, package versions, and start/completion timestamps. Conformer, bootstrap,
model, and figure seeds are declared separately in `config.yaml`. Live APIs are
refresh-only and cannot overwrite the pinned snapshot.

Before scoring, the runner prewarms the conformer cache with the fixed worker count
in `chem3d.prewarm_workers`. This is scheduling only: every molecule still uses the
single RDKit thread, seed, ETKDGv3/MMFF94 parameters, and cache key recorded in its
own manifest. Unique structures are processed in canonical-SMILES order, and the
aggregate prewarm status reports cache hits, failures, and coverage without writing
identifiers or structures.

The legacy snapshot records CARD 4.0.2. Its original ChEMBL release, UniProt
release/query date, PubChem query date, and RCSB query date were not retained by
v2; these remain explicit `null` provenance gaps and are never inferred from Git
history.

## Molecular representations and conformers

The 2D layer is unchanged: ECFP4 is Morgan radius 2 with 2,048 bits, and MACCS
uses RDKit's standard keys. For v3, each sanitized molecule receives a seeded
ETKDGv3 ensemble. The default configuration requests 30 conformers, prunes at
0.5 Å RMS, optimizes with MMFF94 for at most 500 iterations, and retains
conformers within 10 kcal/mol of the lowest finite MMFF94 energy. Embedding uses
one thread for determinism.

Conformers are cached as RDKit binary molecules. The cache key is SHA-256 over
canonical isomeric SMILES, every conformer parameter, the conformer seed, and the
RDKit version. A failed embedding or unavailable MMFF parameterization yields an
explicit status and missing 3D score; it never yields a fabricated similarity.

## 3D shape and pharmacophore evidence

For each query/reference conformer pair, USRCAT provides an alignment-free
shape-plus-feature score. V3 retains every reference-level score and reports the
maximum and configured top-k mean for each query/target class.

References are ranked by USRCAT within target class. Only the configured top 25
references enter the alignment-based layer. RDKit O3A uses MMFF94 atom typing by
default. After alignment, shape similarity is `1 - ShapeTanimotoDist`; the color
field is RDKit's normalized shape-feature color score. Both remain bounded in
`[0,1]`, with overlay attempts and failures reported separately.

The complementary pharmacophore field uses the RDKit Gobbi/Poppe feature-pair
fingerprint and Tanimoto similarity. This signal is alignment-free; it is kept
separate from O3A color rather than presented as the same measurement.

## Evidence fusion

The default chemical v3 score is normalized Reciprocal Rank Fusion (RRF), applied
within each query over target classes. The configured components are ECFP4,
MACCS, USRCAT, O3A shape, O3A color, and Gobbi pharmacophore similarity. If
`r_i(q,t)` is the rank for component `i`, `k=60`, and `I_i` indicates that the
component exists, then:

```text
RRF_raw(q,t) = Σ_i I_i / (k + r_i(q,t))
chemical_evidence_score_v3 = RRF_raw / (number_configured_components / (k + 1))
```

Exact ties receive their average rank. The fixed denominator means missing
evidence is not rewarded. A row with no component remains missing. Component
ranks and contributions are emitted alongside the fused score. The fused score
is a ranking score, not a probability.

The 2D-vs-v3 disagreement table compares the legacy 2D chemical rank with the v3
rank. The default materiality threshold is an absolute within-query shift of
three rank positions. Every qualifying promotion and demotion is listed with all
source scores. If no case qualifies, the schema-valid report is empty; no
scaffold-hop example is manufactured.

## Leakage-controlled benchmark

The v3 benchmark applies an ECFP4 analogue guard at Tanimoto `>=0.85` in every
split and asserts that no retained reference violates the guard.

- **Target-family:** references belonging to the labelled target family are
  withheld using the pinned ontology. Missing family mappings are pending rows,
  not guessed mappings. A similarity-to-known-reference method may consequently
  have zero positive-target coverage; that failure is retained.
- **Scaffold:** exact nonempty Bemis–Murcko scaffold matches are removed in
  addition to close analogues. The provenance table reports both removal counts.
- **Temporal:** only pre-cutoff reference measurements may train, and only dated
  post-cutoff queries may test. Records with missing dates are excluded. The
  pinned public benchmark lacks the needed dates, so temporal evaluation is
  currently pending rather than simulated.

Metrics are AUROC, BEDROC at alpha 20 and 80.5, EF at 1% and 5%, MRR, and
coverage. Metrics are computed per query and averaged. Exact score ties receive
average ranks; EF handles a tied boundary fractionally. A seeded nonparametric
bootstrap resamples query identifiers, preserving the dependent target rows
inside each query. Every aggregate carries its 95% interval, total/evaluable
query counts, seed, snapshot, and split-removal provenance. Undefined metrics
carry `NaN` bounds and an unavailable status.

The single mode-comparison table contains 2D-only, 3D-plus-pharmacophore-only,
and fused retrieval for all three splits. `performance_vs_2d` is assigned
mechanically as improved, equal, worse, or unavailable. A worse 3D or fused
result is reported unchanged.

## Decoys and external baseline

Cross-target reference ligands remain specificity controls only. They are not
called experimentally inactive. The property-matched decoy loader requires a
versioned record identifier, source dataset/release, matching method, target
class, and valid SMILES. The pinned snapshot has no such artifact. The official
DUD-E arbitrary-ligand workflow requires interactive submission, so T2.3 remains
pending and enrichment is not claimed against substituted negatives.

PIDGINv4 is the selected external-baseline adapter. Code is pinned to commit
`df0f6068a8aa16e2278e3779a1ad5e6d552731dc` and models to DOI
`10.6084/m9.figshare.19108382.v1`. The adapter checks the official Python 2.7
runtime, code commit, model pickle/AD directories, metadata, and a versioned
PIDGIN-target-to-project-ontology map. These prerequisites are absent locally,
so the head-to-head file is empty with a pending status. PIDGINv4 scores are
explicitly marked uncalibrated, consistent with its documentation.

## Combiner calibration and sensitivity

The learned path is standardized regularized logistic regression followed by a
Platt model fit on a disjoint calibration-query set. Train, calibration, and test
query identifiers must not overlap. Rows missing any selected feature are
excluded and counted, never imputed. Held-out AUROC and Brier score use query
bootstrap intervals; the reliability table uses fixed probability bins.

The pinned benchmark does not contain valid temporal results, property-matched
decoys, or declared train/calibration/test roles. Therefore no publishable model
is fit and no probability is emitted. The active justification path is mandatory
sensitivity analysis of the transparent rank fusion.

Sensitivity multiplies each implicit equal RRF weight by 0.50, 0.75, 1.25, and
1.50; removes every component individually; and removes the 2D, shape, and
pharmacophore layers in turn. Mean Kendall tau and finite extrapolated RBO quantify
rank stability. The bootstrap resamples reference ligands within each
query/target class and re-aggregates component maxima. O3A bootstrap variability
is conditional on the originally evaluated USRCAT shortlist, which is retained
as a limitation rather than hidden.

## Applicability domain

Every v3 prediction retains nearest-reference ECFP4 Tanimoto and USRCAT
similarity (`1 - similarity` is also reported as distance). The declared
Tanimoto thresholds are:

| Flag | Rule |
|---|---|
| in domain | Tanimoto >= 0.40 |
| near domain | 0.25 <= Tanimoto < 0.40 |
| out of domain | Tanimoto < 0.25 |
| unassessable | nearest-reference Tanimoto missing |

No calibrated USRCAT cutoff is specified in the research plan, so USRCAT remains
a continuous independent AD measurement. Out-of-domain and unassessable rows are
ordered after eligible rows in organism shortlists, while their numeric evidence
and priority scores remain unchanged.

## Figure and table generation

`pipeline/v3_figures.py` reads only run outputs. It can render the split/mode
metric comparison with confidence intervals, AD coverage, material rank shifts,
and a reliability curve when those tables contain data. Missing or empty inputs
are listed in `results/v3_figure_status.csv`; the script does not draw placeholder
scientific panels. Sensitivity plotting is generated by
`pipeline/sensitivity_analysis.py` from the same reference-bootstrap table.

## Limitations and negative-result policy

The reference universe is assay-heterogeneous and strongly imbalanced across
target classes. Large macrocycles and peptides can fail ETKDG/MMFF or lie outside
the small-molecule reference domain. O3A is limited to the USRCAT shortlist.
Target-family holdout can be structurally unidentifiable for a nearest-reference
classifier. Temporal evaluation, property-matched decoys, the external baseline,
and a scientifically supportable calibrated combiner remain pending data/runtime
gaps in the legacy snapshot.

Neither 3D improvement nor a material disagreement is assumed. The generated
comparison and disagreement tables are authoritative even when they show no
gain, worse retrieval, or no qualifying case. Target hypotheses still require
orthogonal biochemical, genetic, organism-specific, permeability/efflux, and
selectivity experiments.
