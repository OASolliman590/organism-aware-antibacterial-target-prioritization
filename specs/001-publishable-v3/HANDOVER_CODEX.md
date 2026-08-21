# Handover Prompt — paste to Codex

You are taking over an existing cheminformatics repository to make it **publishable**. Work strictly from the
Spec-Kit documents already in the repo. Do not invent scientific methods or data.

## Repository
`organism-aware-antibacterial-target-prioritization` (local clone, branch `main`). It is an organism-aware,
multi-evidence antibacterial **target-prioritization** framework: for query compounds it ranks a public target
universe per ESKAPE organism by layering chemical similarity, target specificity, reference quality, species
sequence-transfer, biology, CARD resistance, RCSB pocket precedent, anti-target risk, and uncertainty. Current
core engine: `pipeline/open_target_discovery_v2.py`. It currently uses **2D fingerprint similarity only** and
**hand-tuned uncalibrated weights**, with no leakage-controlled benchmark and no reproducibility spine.

## Your mission
Execute Feature Spec 001 to add a **3D shape/pharmacophore matching layer**, replace/justify the scoring with
**calibration or a sensitivity analysis**, add a **leakage-controlled benchmark with uncertainty and an external
baseline**, and make every run **reproducible**. This elevates it from an internal heuristic to a method a
reviewer can trust.

## Read these first (in order), then follow them exactly
1. `memory/constitution.md` — non-negotiable principles. The most important: **never fabricate or impute data**;
   **no score is a probability unless calibrated**; **benchmark only under explicit target/scaffold/time splits
   with uncertainty**; **reproducibility (pinned snapshots, seeds, provenance) is a deliverable**; **keep every
   evidence layer as a separate auditable field**; **report negative results honestly**.
2. `specs/001-publishable-v3/spec.md` — WHAT and WHY, functional requirements FR-1…FR-10, success criteria.
3. `specs/001-publishable-v3/research.md` — the real, citable methods for each decision (ETKDGv3 conformers,
   USRCAT, Open3DAlign/O3A shape+color, RDKit pharmacophores, rank fusion, BEDROC/EF/AUROC + bootstrap CIs,
   Bemis–Murcko/temporal splits, DUD-E-style decoys, SEA/SwissTargetPrediction baseline, calibration, AD).
4. `specs/001-publishable-v3/plan.md` — architecture, new modules, config, phasing.
5. `specs/001-publishable-v3/tasks.md` — the ordered, agent-executable task list with acceptance criteria.

## Execution rules
- Work **phase by phase** in `tasks.md` order (Phase 0 reproducibility spine first — everything depends on
  determinism). Do not start a phase until the previous phase's acceptance criteria pass.
- **Phase 0 hard gate:** create the golden-file test that freezes current v2 outputs *before* changing anything,
  so v3 is provably additive.
- For each task: implement, add/extend tests, run them, and only then check the box. Keep commits small and
  scoped to one task; reference the task id in the commit message.
- If external data/APIs (ChEMBL, CARD, UniProt, RCSB, the baseline tool) are unavailable, **log the gap and mark
  the result pending — never simulate** (constitution I/III).
- If 3D matching or the learned combiner does **not** improve retrieval, report that honestly; it is a valid
  result, not a failure to hide.
- Pin RDKit and all deps; record versions in the run manifest. Fix seeds everywhere stochastic.

## Definition of done
All `spec.md` §7 success criteria met: 3D layer implemented as separate auditable fields with a 2D-vs-3D
disagreement report; benchmark reports target/scaffold/temporal splits with bootstrap CIs and an external
baseline; rankings are either calibrated+validated or shown stable via sensitivity analysis; one command +
one pinned snapshot reproduces every reported number and emits a provenance manifest; methods supplement and
tests updated and green.

## Context you should know (do not re-derive)
- Real query compounds used in the associated study are two novel scaffolds (xanthine-linked hybrids); against
  this pipeline they all score Low/Insufficient (expected for novel chemistry) and the β-lactamase target got
  **zero** chemical support — a useful honest negative control to keep in the benchmark/reporting.
- This tool is decision-support for *target hypotheses*; docking/MD/in-vivo are separate downstream evidence and
  are explicitly out of scope here.

Begin with Phase 0. Confirm your understanding of the constitution and the phase plan before writing code.
