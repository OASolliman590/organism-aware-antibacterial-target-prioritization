# Feature Spec 001 — Publishable v3: 3D-Aware, Calibrated, Reproducible Prioritization

**Status:** Draft for implementation
**Owner:** (author) → executed by coding agent (Codex)
**Depends on:** current v2 engine (`pipeline/open_target_discovery_v2.py`), constitution (`memory/constitution.md`)

## 1. Problem statement (WHY)

The v2 pipeline is a transparent, organism-aware, multi-evidence *decision-support* framework, but it is
not yet publishable as a method because:

1. The chemical evidence layer is **2D-fingerprint similarity only** — a weak, promiscuity-prone signal for
   target identification that misses shape/pharmacophore complementarity and is vulnerable to 2D scaffold bias.
2. The scoring **weights are hand-tuned and uncalibrated**, with no held-out validation and no sensitivity analysis.
3. **Benchmarking is not leakage-controlled** across target/scaffold/time, and there is **no comparison to an
   external baseline**.
4. **Reproducibility depends on live APIs** with no pinned snapshots, no tests, no run provenance.

## 2. Objective (WHAT)

Elevate the framework to a defensible, reproducible method suitable for a methods-supported application paper by:

- **A. Adding a 3D shape + pharmacophore matching layer** as an additional, auditable evidence stream, fused
  with the existing 2D layer.
- **B. Replacing hand-tuned weights** with a calibrated, validated scoring model (or, minimally, a documented
  sensitivity/ablation analysis proving stability).
- **C. Establishing leakage-controlled retrospective benchmarking** with uncertainty and an external baseline.
- **D. Making every run reproducible** via pinned data snapshots, seeds, declarative config, tests, and provenance.

## 3. In scope
- New 3D-matching module and its evidence fields.
- Benchmarking harness with target/scaffold/time splits, enrichment metrics, bootstrap CIs, decoys, baseline.
- Score calibration/learning + sensitivity analysis.
- Reproducibility layer (config, snapshots, seeds, tests, provenance manifest).
- Documentation/figures updates for a methods supplement.

## 4. Out of scope
- Docking, MD, or any structure-based free-energy calculation (remain downstream, separate).
- Wet-lab validation.
- New organism additions beyond the current six (may be a later feature).
- A production web service.

## 5. Functional requirements

- **FR-1 (3D conformers):** Generate a reproducible conformer ensemble per query and per reference ligand
  (fixed seed, bounded ensemble size, energy window), cached to disk with provenance.
- **FR-2 (3D shape similarity):** Compute an alignment-free shape+pharmacophore descriptor similarity
  (USRCAT) AND an alignment-based shape+color overlay score (Open3DAlign/O3A) between each query and each
  target-class reference ligand. Emit max and top-k aggregates per (query, target_class).
- **FR-3 (pharmacophore):** Compute a 3D pharmacophore-feature similarity (RDKit feature factory / Gobbi
  pharmacophore fingerprint) per (query, target_class).
- **FR-4 (fusion):** Fuse 2D (ECFP4/MACCS), 3D shape, and pharmacophore evidence into a single
  `chemical_evidence_score_v3` via a documented rank/score-fusion rule, keeping each component as a separate field.
- **FR-5 (calibration):** Provide a calibrated scoring path: either a learned model (logistic regression or
  learning-to-rank) trained on benchmark labels with held-out evaluation and probability calibration, or the
  retained heuristic *plus* a mandatory sensitivity/ablation report. The active mode is config-selectable.
- **FR-6 (benchmark):** A harness that evaluates retrieval under (a) target-family holdout, (b) Bemis–Murcko
  scaffold holdout, (c) temporal split; reports AUROC, BEDROC(α=20 and 80.5), EF@1%/5%, coverage, MRR, all with
  bootstrap 95% CIs; and includes property-matched decoys.
- **FR-7 (baseline):** Run at least one external target-prediction baseline on the same benchmark queries and
  report head-to-head.
- **FR-8 (applicability domain):** Every prediction carries an AD flag from nearest-neighbour distance to the
  reference set, with a documented threshold.
- **FR-9 (reproducibility):** Declarative config; pinned, dated data snapshots with recorded versions; a
  provenance manifest per run (data versions, commit, config hash, seeds); deterministic outputs given the same
  snapshot+config+seed.
- **FR-10 (regression safety):** The v2 outputs remain reproducible; v3 adds fields/files rather than silently
  changing v2 semantics. A golden-file test guards this.

## 6. Non-functional requirements
- Runtime: full 3D layer for the current 2-compound × ~20-target-class reference set completes on a laptop in
  minutes; benchmark (hundreds of molecules) in <~1 hour with caching.
- Determinism: identical results for identical snapshot+config+seed.
- Test coverage on the new scoring/fusion/benchmark logic; CI runs the fast subset.

## 7. Success criteria (how we know it is publishable-ready)
1. 3D layer implemented, cached, and reported as separate fields; a documented case where 2D and 3D disagree
   (e.g. a scaffold-hop) is captured in the results.
2. Benchmark reports all three splits with bootstrap CIs and at least one external baseline; the scaffold- and
   time-split numbers are presented honestly (no leakage), and the framework meets or transparently underperforms
   the baseline with the reasons analysed.
3. Rankings are shown stable under weight perturbation (sensitivity analysis) OR replaced by a calibrated learned
   model whose calibration curve and held-out metrics are reported.
4. A single command reproduces all reported numbers from a pinned snapshot; the provenance manifest is emitted.
5. Methods supplement updated with the 3D method, benchmark protocol, and limitations.

## 8. Explicit non-goals / honesty guards
- If the benchmark cannot support a valid split (insufficient public metadata), report it as *unavailable*,
  not as a passing result (constitution III).
- If 3D matching does not improve retrieval, that negative result is reported, not hidden.
