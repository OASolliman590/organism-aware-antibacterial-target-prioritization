# Tasks — Spec 001 (ordered, agent-executable)

Rules: obey `memory/constitution.md`. Each task lists acceptance criteria. Do not start a phase until the prior
phase's tasks pass. Never fabricate data (constitution I). Report negative results honestly (constitution III).

## Phase 0 — Reproducibility spine (do first)
- [x] **T0.1** Add `pipeline/config.py` + `config.yaml` (keys per plan §4). All thresholds/seeds/paths read from
      config; remove reliance on ad-hoc env vars for scientific params. **AC:** engine runs from a single config; no
      scientific constant is hard-coded in two places.
- [x] **T0.2** Add `pipeline/provenance.py` emitting `results/run_manifest.json` (data snapshot versions, git commit,
      config hash, seeds, timestamps, package versions). **AC:** manifest written on every run.
- [x] **T0.3** Create `data/snapshots/SNAPSHOT_VERSIONS.json` and repoint loaders to pinned snapshots; live fetch is
      refresh-only and writes a new dated snapshot, never overwrites. **AC:** a reported run cites a snapshot id.
- [x] **T0.4** Add `tests/` with pytest: (a) tiny fixture reference+query set; (b) **golden-file test freezing current
      v2 outputs** so v3 cannot silently change v2 (FR-10); (c) CI workflow running the fast subset. **AC:** `pytest`
      green; golden test fails if v2 numbers change.

## Phase 1 — 3D matching layer (the core new science)
- [x] **T1.1** `pipeline/chem3d_matching.py`: ETKDGv3 conformer generation (seeded, N/energy/RMS pruning per config),
      disk cache keyed by canonical SMILES+params. **AC:** deterministic conformers; cache hit on rerun.
- [x] **T1.2** USRCAT descriptor + score per (query, reference); aggregate `usrcat_max`, `usrcat_top5_mean` per
      (query, target_class). **AC:** fields present for all query×class pairs; runtime acceptable via caching.
- [x] **T1.3** O3A alignment on USRCAT-shortlisted references → `o3a_shape_tanimoto_max`, `o3a_color_max`. **AC:**
      shortlist size from config; values in [0,1]; documented.
- [x] **T1.4** Pharmacophore similarity (RDKit features + Gobbi_Pharm2D) → `pharmacophore_sim_max`. **AC:** field present.
- [ ] **T1.5** `pipeline/evidence_fusion.py`: rank/score fusion → `chemical_evidence_score_v3`, retaining all
      component columns (constitution VI). **AC:** components auditable; fusion deterministic.
- [ ] **T1.6** Wire v3 chemical evidence into `open_target_discovery_v2.py` behind `combiner`/`fusion_mode` config,
      leaving the v2 path intact. **AC:** golden v2 test still green; new `*_v3.csv` outputs produced.
- [ ] **T1.7** Emit a **2D-vs-3D disagreement report** (cases where 3D materially changes a target's rank, e.g. a
      scaffold hop). **AC:** report file lists ≥ the cases found, with per-component scores.

## Phase 2 — Leakage-controlled benchmark
- [ ] **T2.1** `pipeline/benchmark_v3.py`: implement target-family, Bemis–Murcko scaffold, and temporal splits;
      make scaffold exclusion real and log how many references it removes. **AC:** three splits produced; leakage
      guard asserts query analogues excluded at the configured threshold.
- [ ] **T2.2** Metrics: AUROC, BEDROC(α=20 & 80.5), EF@1%/5%, MRR, coverage — each with **bootstrap 95% CIs**.
      **AC:** every metric reported with CI, n, and split provenance.
- [ ] **T2.3** Property-matched decoys (DUD-E-style) integrated as negatives; keep existing cross-target decoys for
      the specificity margin. **AC:** decoy provenance recorded; enrichment computed against them.
- [ ] **T2.4** Compare 2D-only vs 3D-only vs fusion under each split. **AC:** a single table; honest reporting even
      if 3D does not help (constitution III).

## Phase 3 — External baseline
- [ ] **T3.1** `pipeline/baseline_external.py`: run one of SEA / SwissTargetPrediction / PIDGINv4 on the benchmark
      queries; cache results. **AC:** head-to-head table vs our method on identical queries; if the service is
      unavailable, mark **pending**, do not fabricate.

## Phase 4 — Calibration / combiner + sensitivity
- [ ] **T4.1** `pipeline/scoring_model.py`: learned combiner (regularised logistic or LGBMRanker) trained under the
      Phase-2 splits; features = the separate evidence fields; **Platt/isotonic calibration**; report reliability
      diagram + Brier score. **AC:** held-out metrics + calibration curve emitted; only calibrated outputs called
      "probability".
- [ ] **T4.2** `pipeline/sensitivity_analysis.py`: if labels too sparse for T4.1, or as a complement, perturb weights
      (±25–50%, leave-one-layer-out) and report rank stability (Kendall τ, RBO, bootstrap over references). **AC:**
      stability of top hypotheses quantified and plotted.

## Phase 5 — Applicability domain + reporting
- [ ] **T5.1** `pipeline/applicability_domain.py`: NN-distance AD flags (Tanimoto + USRCAT) per prediction; threshold
      from config. **AC:** every prediction row carries an AD flag; out-of-domain discounted in shortlist logic.
- [ ] **T5.2** Update `docs/` methods supplement: 3D method, benchmark protocol, calibration, AD, limitations, and
      the honest negative results. Regenerate figures. **AC:** supplement reproducible from a pinned run.
- [ ] **T5.3** End-to-end: `python run_pipeline.py --config config.yaml` reproduces all reported numbers from a
      pinned snapshot and writes the provenance manifest. **AC:** clean-clone reproduction succeeds.

## Cross-cutting acceptance (must all hold at the end)
- Golden v2 outputs unchanged; v3 is additive.
- Every reported metric carries uncertainty + split provenance.
- No fabricated data anywhere; gaps are logged as gaps.
- Rankings are either calibrated-and-validated or demonstrably stable.
- One command + one pinned snapshot = full reproduction.
