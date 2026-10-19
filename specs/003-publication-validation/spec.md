# Feature Spec 003 — Testing & validation program for publication

**Status:** Draft for implementation
**Owner:** (author) → Opus 5.5 spec; execution by coding agents
**Depends on:** v2/v3 engine, `memory/constitution.md`, `docs/v3_method_supplement.md`,
`SCIENTIFIC_VALIDITY_REVIEW.md`, spec 001 (benchmark machinery already built)

## 1. Why this spec exists

The engine, benchmark machinery, calibration path, sensitivity analysis, and reproducibility
spine already exist (spec 001). What does **not** yet exist is a *completed, defensible
validation program* whose results a methods paper can stand on. Today the decisive
evaluations are marked **pending** in the method supplement — temporal split, property-matched
decoys, external baseline, and a calibrated combiner — and the reportable retrieval numbers
come from the easiest split (close-analogue). This spec defines, **in advance of running
anything**, the experiments, datasets, metrics, statistics, and *pre-committed pass/fail
thresholds* that determine whether and how the method is published. It is a pre-registration,
not a results document.

This spec adds **no new scoring science**. It validates the method that exists. If a
validation fails, the paper reports the failure (Constitution III); it does not trigger a
redesign hidden inside the validation.

## 2. The claims the paper is allowed to make

The paper makes **only** claims that a pre-committed experiment in §4 supports. Each claim is
bound to the experiment that can refute it.

| ID | Claim | Supported by | Refuted if |
|----|-------|--------------|------------|
| C1 | The method retrieves the correct target class above chance for compounds **inside** the reference chemical space. | E1 (analogue-controlled retrieval) | AUROC 95% CI includes 0.5, or EF5% CI includes 1.0 |
| C2 | 3D/pharmacophore fusion adds information **over** 2D-only retrieval. | E2 (mode comparison, same split) | Fused AUROC/BEDROC CI not above 2D-only on the pre-chosen primary metric |
| C3 | The method **generalises beyond close analogues** (scaffold hop). | E3 (scaffold-holdout retrieval) | Scaffold-split primary metric CI includes the no-skill value |
| C4 | Ranking conclusions are **stable** to reasonable weight perturbation. | E6 (sensitivity, already built) | Mean Kendall-τ < 0.7 under ±25% perturbation on the primary shortlist |
| C5 | Where a probability is emitted, it is **calibrated**. | E5 (held-out Platt calibration) | Reliability slope outside [0.8, 1.2] or Brier not better than prevalence baseline |
| C6 | The method is **competitive with / complementary to** an established external baseline. | E4 (PIDGINv4 head-to-head) | Pre-registered descriptive comparison only; see §4.4 |

**Claims the paper will NOT make** (hard scope limits, Constitution VII): any claim of
antibacterial activity, target engagement, binding, or clinical relevance; any claim of
target-family generalisation if E7 shows the method cannot (it likely cannot — see §4.7);
any single-number headline metric without its split and CI.

## 3. Primary vs secondary endpoints (decided now, not after seeing results)

- **Primary retrieval metric:** BEDROC(α=20) per query, averaged, with 1000-sample query
  bootstrap 95% CI. Rationale: early recognition is the decision-relevant regime and α=20
  (not 80.5) resolves more than the top slot at this universe size (review §3.3).
- **Primary split for the headline generalisation claim (C3):** scaffold holdout.
- **Secondary metrics:** AUROC, BEDROC(α=80.5), EF1%, EF5%, MRR, coverage — all reported with
  CIs and with their effective `top_n` annotated (review §3.3).
- **Primary calibration endpoint (C5):** reliability-curve slope + Brier vs prevalence.
- Any metric whose per-query evaluation is structurally degenerate for a query (e.g. EF1%
  with `n_candidates` too small to resolve 1%) is reported as `unavailable` for that query
  with the reason, never silently collapsed into the mean.

## 4. The experiments (pre-committed)

### E1 — Analogue-controlled retrieval (supports C1)
Run the existing benchmark with the ≥0.85 ECFP4 analogue guard asserted at runtime. Report
the full metric panel per split. **Pass:** primary metric CI strictly above no-skill on the
analogue-controlled set. This is *interpolation within known chemistry* and the paper must
label it as such.

### E2 — Mode comparison 2D vs 3D-only vs fused (supports C2)
Use `compare_score_modes`. The 3D-only component list is fixed in
`benchmark.three_dimensional_components` (no 2D Gobbi leakage). **Pass for C2:** fused primary
metric CI lies above 2D-only on the *same* split; `performance_vs_2d` is assigned mechanically
and a worse/equal result is reported unchanged (Constitution III). A failure here is a
publishable negative result, not a blocker.

### E3 — Scaffold-holdout generalisation (supports C3, primary)
Exact Bemis–Murcko scaffold removal **plus** the analogue guard; both removal counts reported.
**Pass:** scaffold-split primary metric CI above no-skill. **Guard against the review §2
finding:** the spec requires asserting that the scaffold split actually removed additional
molecules beyond the analogue guard (`n_scaffold_removed_beyond_analogue > 0` for the reported
classes); if it did not, the result is reported as "scaffold generalisation untested on this
benchmark," never as a pass.

### E4 — External baseline head-to-head (supports C6)
PIDGINv4 (`external_baseline` config: pinned commit, figshare model DOI, Python 2.7 runtime,
AD percentile 90). Because PIDGIN scores are uncalibrated, the comparison is **descriptive and
rank-based** on the shared target-mappable subset, pre-registered as such; no significance
claim of superiority is made. **Deliverable:** a head-to-head table over the mappable
intersection + a provenance record of every query/target dropped for lack of mapping. If the
PIDGIN runtime cannot be stood up, E4 is reported `unavailable` with the exact blocker and C6
is dropped from the paper — not substituted with a weaker baseline silently.

### E5 — Calibrated combiner (supports C5)
Standardised L2-logistic + Platt on **disjoint** train/calibration/test query sets
(`scoring_model.py` already enforces disjointness and refuses to emit a probability without
real split roles). **Prerequisite:** E-data must supply valid train/calibration/test role
labels and ≥`model.minimum_class_count` per class per role. **Pass:** reliability slope in
[0.8, 1.2] and held-out Brier below the prevalence-only baseline, both with query-bootstrap
CIs. If the labelled set is too sparse (the current state), E5 is reported `pending_sparse`
via `assess_training_readiness` and C5 is dropped; **no probability is published** in that
case.

### E6 — Sensitivity / ablation (supports C4)
Already implemented (`sensitivity_analysis.py`): per-component RRF weight perturbation
(×0.5/0.75/1.25/1.5), leave-one-component-out, layer removal; Kendall-τ + RBO with bootstrap
CIs; and the final-ranking organism-query grouped perturbation of every priority subweight
(now including the new `transfer.missing_factor`). **Pass:** primary-shortlist mean Kendall-τ
≥ 0.7 under ±25% perturbation. **New requirement from spec 002:** the sweep must include
`transfer.missing_factor` and report ranking stability as it moves across its plausible range
(0.40–0.65).

### E7 — Target-family holdout (honest ceiling, supports the scope limit, not a claim)
Withhold all references of the labelled family. A nearest-reference method is expected to have
near-zero positive coverage. **This experiment exists to quantify and disclose the ceiling,
not to pass.** The paper reports its result verbatim as the method's structural limit
(Constitution III; review §2).

### E8 — Temporal holdout (supports C3 as corroboration if data allow)
Pre-cutoff references train, dated post-cutoff queries test (`time_cutoff: 2018-01-01`).
Records lacking dates are excluded, never dated by inference. The pinned public benchmark
currently lacks the dates, so E8 is reported `pending_temporal_metadata` unless a dated
snapshot is produced. If available, it corroborates C3; if not, its absence is stated.

## 5. Datasets and provenance (Constitution IV)

- Every experiment runs against a **pinned snapshot** (`data/snapshots/SNAPSHOT_VERSIONS.json`)
  and emits a `run_manifest.json` (snapshot id + hashes, git commit, config hash, seeds,
  package versions, timestamps).
- Any new data required (property-matched decoys for E-decoy, dated references for E8, the
  PIDGIN target-ontology map for E4) is acquired via a **refresh-only** config that writes a
  new dated snapshot and never overwrites the pinned one; a raw refresh stays
  `analysis_ready: false` until all derived layers are rebuilt and pinned.
- No fabricated negative, decoy, date, or score is ever written (Constitution I). Missing
  prerequisites produce an `unavailable`/`pending` status row, not a filled value.

## 6. Statistics (decided now)

- All point metrics carry a 1000-sample (`benchmark.bootstrap_n`) **query-level** bootstrap
  95% CI; the bootstrap resamples query identifiers and keeps each query's dependent target
  rows together (and, for the decoy task, resamples target-class clusters). Seeds from
  `config.yaml`.
- Multiplicity: the paper evaluates a small, **pre-registered** set of primary endpoints
  (C1–C6). Secondary metrics are descriptive. No post-hoc metric or split is promoted to
  primary after seeing results.
- A claim passes only if its pre-committed CI criterion in §2/§4 is met. "Trend toward" is not
  a pass.

## 7. Required paper artifacts (what "done" produces)

1. `results/` tables for E1–E8 with CIs, split-removal provenance, and per-query detail.
2. A **figure set**: split×mode metric comparison with CIs; AD coverage; 2D-vs-3D material
   rank-shift (disagreement) panel; calibration reliability curve (only if E5 passes);
   sensitivity τ/RBO panels including the `missing_factor` sweep.
3. A **claims-to-evidence table** (C1–C6 → experiment → pass/fail → figure/table) generated
   from the results, not hand-written.
4. A **reproducibility package manifest**: pinned snapshot id, commit, config hash,
   environment lock, and the exact commands to regenerate every table and figure.
5. A **limitations section** sourced from E7 (family-holdout ceiling), E4/E5/E8 pending
   statuses, assay heterogeneity (review §3.5), and the analogue-split interpretation caveat
   (review §2).

## 8. Out of scope

- No change to the scoring formulae, fusion, or priority weights (except that E6 must sweep
  the already-added `missing_factor`).
- No new modelling method beyond the specified logistic+Platt combiner.
- No wet-lab validation; the paper is a methods/retrospective-validation paper only.

## 9. Pre-registration integrity

This spec is committed to version control **before** the validation runs. The passing
criteria in §2/§4/§6 are fixed at that commit. Any later change to a threshold or endpoint is
a tracked amendment with a rationale, not a silent edit — mirroring the constitution's
amendment rule.
