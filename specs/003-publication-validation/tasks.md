# Tasks — Spec 003 (ordered, agent-executable)

Rules: obey `memory/constitution.md`. This spec is a **pre-registration** — commit Phase 0
(thresholds) BEFORE running any experiment, and never edit a threshold after seeing a result
except as a tracked amendment (spec §9). Never fabricate data; a missing prerequisite is an
`unavailable`/`pending` status row, not a filled value. Do not change scoring/fusion/priority
formulae (the only allowed touch is adding `transfer.missing_factor` to the E6 sweep).

## Phase 0 — Pre-registration (do first, commit before any run)
- [ ] **T0.1** Write `specs/003-publication-validation/acceptance.yaml` with the pre-committed
      thresholds from spec §2/§3/§4/§6: primary metric `bedroc_alpha_20_0`; primary
      generalisation split `scaffold`; no-skill values per metric; τ ≥ 0.70 at ±25%;
      reliability slope ∈ [0.8, 1.2]; Brier < prevalence. **AC:** parses; every C1–C6 claim has
      an entry naming its experiment, metric, and pass rule.
- [ ] **T0.2** Commit Phase 0 on its own so the pre-registration predates results. **AC:** a
      commit contains `acceptance.yaml` and this spec, with no experiment outputs.

## Phase 1 — Run the experiments that need no new data
- [ ] **T1.1 (E1)** Run the analogue-controlled benchmark; emit the full metric panel per
      split with 1000× query bootstrap CIs and split-removal provenance. **AC:** tables carry
      CIs, seed, snapshot id; EF/BEDROC cells annotate effective `top_n` (review §3.3).
- [ ] **T1.2 (E2)** Run `compare_score_modes` (2D / 3D-only / fused); the 3D-only list stays
      `benchmark.three_dimensional_components`. **AC:** `performance_vs_2d` assigned
      mechanically; a worse/equal result reported unchanged.
- [ ] **T1.3 (E3)** Run the scaffold-holdout split with the analogue guard; assert
      `n_scaffold_removed_beyond_analogue > 0` for reported classes. **AC:** if the assertion
      fails, the output says "scaffold generalisation untested on this benchmark" and C3 is
      NOT marked pass.
- [ ] **T1.4 (E7)** Run the target-family holdout; report coverage and the primary metric
      verbatim. **AC:** zero/near-zero positive coverage is reported as the method's ceiling,
      not hidden or retried.
- [ ] **T1.5 (E6)** Extend `sensitivity_analysis.py` so the final-ranking sweep includes
      `v2_scoring.overall_factors.transfer.missing_factor` across 0.40–0.65, alongside the
      existing subweights. **AC:** the sensitivity table/figure gains a `missing_factor` row;
      full test suite stays green.

## Phase 2 — Experiments that need acquired data (refresh-only, new dated snapshots)
- [ ] **T2.1 (E8 data)** Attempt a dated-reference refresh snapshot for the temporal split. If
      dates are unavailable, STOP and record `pending_temporal_metadata`. **AC:** either a new
      dated snapshot (create-only, pinned, `analysis_ready` gated) or an explicit pending row;
      no inferred dates.
- [ ] **T2.2 (E8 run)** If T2.1 produced dates, run the temporal split and report; else skip
      with the pending status. **AC:** records-missing-dates excluded, never dated by guess.
- [ ] **T2.3 (E-decoys)** Acquire a property-matched decoy artifact that passes
      `benchmark_decoys.py` validation; run the decoy retrieval task (2D/3D/fusion) with
      target-class-cluster bootstrap. If the official workflow can't be automated, record the
      T2.3 pending status. **AC:** no fabricated negatives; validation rejects bad provenance.
- [ ] **T2.4 (E4)** Stand up the PIDGINv4 adapter (pinned commit, figshare models, Python 2.7,
      AD percentile 90) and produce the descriptive, rank-based head-to-head over the
      target-mappable intersection, plus a dropped-pair provenance table. If the runtime can't
      be built, record `unavailable` + the exact blocker and drop C6. **AC:** no superiority
      significance claim; PIDGIN scores marked uncalibrated.
- [ ] **T2.5 (E5)** Only if T2.2/T2.3 yield valid disjoint train/calibration/test role labels
      with ≥`model.minimum_class_count` per class/role: fit logistic+Platt, emit held-out
      AUROC/Brier with CIs and the reliability curve. Else `assess_training_readiness` reports
      `pending_sparse` and NO probability is published. **AC:** disjointness asserted; only the
      Platt output named a probability.

## Phase 3 — Artifacts
- [ ] **T3.1** Add `pipeline/validation_report.py` generating
      `results/claims_to_evidence.csv` + `results/validation_summary.md` from the E1–E8 tables
      and `acceptance.yaml`. **AC:** verdicts are mechanical (pass/fail/unavailable); a pending
      experiment yields `unavailable`, never a fabricated verdict; unit-tested on a tiny
      fixture.
- [ ] **T3.2** Ensure the figure set (split×mode w/ CIs, AD coverage, 2D-vs-3D rank shift,
      reliability curve iff E5 passed, sensitivity incl. `missing_factor`) renders from run
      tables only, with a status row for any unavailable panel. **AC:** no placeholder
      scientific panels.
- [ ] **T3.3** Write the reproducibility-package manifest (snapshot id, commit, config hash,
      env lock, exact regenerate commands) and a limitations section sourced from E7 + pending
      statuses + review §2/§3.5. **AC:** every table/figure in the paper has a one-command
      regeneration path.

## Gates (run throughout, pinned venv `/tmp/claude-0/venv`)
- `python -m pytest tests/test_v2_golden.py -q` stays green (no scoring change).
- `python -m pytest -q --basetemp=.pytest_tmp` stays green after each phase (add tests for
  T1.5 and T3.1).
- Each experiment's output is traceable to a snapshot id + commit + seed, or carries an
  explicit `unavailable`/`pending` status with its reason.

## Reporting discipline
- A claim is `pass` only if its pre-committed CI criterion in `acceptance.yaml` is met.
- A failed or pending experiment is reported as such in `claims_to_evidence.csv` and the
  paper's limitations; it is never silently dropped and never triggers a hidden threshold
  change. Any threshold change is a tracked amendment to this spec with a rationale.
