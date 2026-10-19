# Plan 003 — Publication validation program (HOW)

Maps each experiment in `spec.md` §4 to the modules that run it, the data it needs, and the
output it must emit. Most machinery already exists (spec 001); this plan is mostly
orchestration, data acquisition, pre-committed thresholds, and artifact generation — not new
scoring code.

## Module map (what already exists vs. what is new)

| Experiment | Existing module(s) | New work |
|---|---|---|
| E1 analogue retrieval | `benchmark_v3.py` (`query_level_metrics`, `aggregate_metrics_with_bootstrap`), `summarize_benchmark.py` | thresholds + artifact row |
| E2 mode comparison | `benchmark_v3.py:compare_score_modes`, `add_3d_only_score` | primary-metric gate |
| E3 scaffold holdout | `benchmark_v3.py:make_split`/`generate_splits` (scaffold) | assert `n_scaffold_removed_beyond_analogue > 0` guard |
| E4 external baseline | `baseline_external.py` (PIDGINv4 adapter) | stand up runtime OR record `unavailable` + blocker |
| E5 calibrated combiner | `scoring_model.py` (`fit_calibrated_combiner`, `predict_calibrated`, `reliability_curve`, `assess_training_readiness`) | supply disjoint split roles OR record `pending_sparse` |
| E6 sensitivity | `sensitivity_analysis.py` | add `transfer.missing_factor` to the swept subweights; range 0.40–0.65 |
| E7 family holdout | `benchmark_v3.py` target-family split | report coverage/ceiling verbatim |
| E8 temporal | `benchmark_v3.py` temporal split; `_record_date` | produce dated snapshot OR record `pending_temporal_metadata` |

## Data acquisition (refresh-only, new dated snapshots)

1. **Property-matched decoys** (for the decoy retrieval task feeding E1/E2 robustness and E5
   labels): the loader in `benchmark_decoys.py` requires versioned record id, source
   dataset/release, matching method, matched active id, target class, valid SMILES. Acquire a
   DUD-E / property-matched set via a refresh config; it must pass the loader's validation
   (rejects blank provenance, dupes, off-ontology targets, exact active overlap). If the
   official DUD-E arbitrary-ligand workflow cannot be automated, record T2.3-style pending.
2. **Dated references** (for E8): only if a dated ChEMBL refresh snapshot is produced; else E8
   stays pending.
3. **PIDGIN target-ontology map** (for E4): `data/benchmark/pidginv4_target_map.csv` mapping
   PIDGIN targets → project ontology; versioned.

All three go through the existing refresh path (`run.refresh_external_data: true`, dated
snapshot id, create-only writers). None overwrites the pinned snapshot.

## Thresholds file

Add a single declarative `specs/003-publication-validation/acceptance.yaml` holding the
pre-committed pass/fail numbers from spec §2/§3/§4 (primary metric = BEDROC α20; scaffold
primary split; τ ≥ 0.7 at ±25%; reliability slope [0.8,1.2]; Brier < prevalence). The
claims-to-evidence generator reads this file so the pass/fail verdicts are mechanical and
auditable, never hand-entered.

## Claims-to-evidence generator (new, small)

`pipeline/validation_report.py` (new): reads the E1–E8 result tables + `acceptance.yaml`,
emits `results/claims_to_evidence.csv` (claim id → experiment → observed metric + CI →
pass/fail/unavailable → figure/table ref) and a `results/validation_summary.md`. It performs
**no** statistics itself — it only compares already-computed CIs against the pinned
thresholds. It must mark a claim `unavailable` (not `fail`, not `pass`) when its experiment is
pending, and never invent a verdict for a missing table.

## Figures

Extend `v3_figures.py` only where a panel is missing; it already renders split/mode metric
comparison with CIs, AD coverage, material rank shifts, and reliability. Add the
`missing_factor` sweep panel to the sensitivity output. Every panel reads run tables only and
writes a status row (`unavailable_source_missing` / `_no_evaluable_rows` /
`_dependency_unavailable`) rather than a placeholder.

## Determinism / constitution check

- Seeds for every bootstrap/model from `config.yaml` (IV). No fabricated data; pending ≠ zero
  (I). Only Platt output named a probability, gated on real splits (II). Splits + CIs +
  external baseline per III. AD retained (V). Evidence layers stay separate (VI). No activity
  claim (VII). The pre-registration (spec committed before runs) enforces honest reporting.

## Non-goals

No scoring-formula change; no new metric invented; no promotion of a secondary metric to
primary after results are seen.
