# Project Constitution — Organism-Aware Antibacterial Target Prioritization

These are non-negotiable principles. Every spec, plan, and task must comply. When a change
conflicts with a principle, the principle wins or the principle is explicitly amended here first.

## I. Evidence integrity (no fabrication)
- Never synthesize, impute, or "fill in" experimental or database values. Missing data is recorded
  as missing (as the current code already does for the ChEMBL outage). A failed API is a logged gap,
  not a reason to simulate.
- Every numeric result must be traceable to a pinned data source (see Principle IV) and a code commit.

## II. Claims must be calibrated, not asserted
- No score is presented as a probability of binding or activity unless it has been calibrated against
  held-out labels (Platt/isotonic) and reported with uncertainty.
- Hand-tuned weights are provisional scaffolding only. Any weight that affects a published ranking must
  either be (a) learned from labelled data with held-out validation, or (b) accompanied by a sensitivity/
  ablation analysis demonstrating the conclusion is stable to reasonable perturbation.
- Decision-support language ("hypothesis", "prioritization") is retained; validation language
  ("validated target", "confirmed hit") is forbidden without wet-lab evidence.

## III. Honest, leakage-free benchmarking
- Retrospective performance is only credible under explicit splits: target-family holdout,
  Bemis–Murcko scaffold holdout, and temporal (time-split) holdout. Report each separately.
- Query molecules and near-analogues are excluded from the reference set at a stated similarity
  threshold; analog bias is measured, not hidden.
- Report AUROC, BEDROC, EF at 1%/5%, and coverage with bootstrap 95% CIs. Never report a single
  point metric without its uncertainty and its split provenance.
- Compare against at least one established external baseline (e.g. SEA, SwissTargetPrediction).

## IV. Reproducibility is a first-class deliverable
- All external data (ChEMBL, CARD, UniProt, RCSB) are pinned to a dated, versioned snapshot recorded
  in a run manifest. Live-API calls are for refresh only, never for a reported result.
- Deterministic seeds for every stochastic step (conformer generation, bootstrap, any learned model).
- Configuration is declarative (YAML/JSON), not environment-variable soup.
- Each run emits a provenance manifest: data versions, code commit, config hash, seeds, timestamps.

## V. Applicability domain is always stated
- Every prediction carries an applicability-domain flag (in/near/out of the reference chemical space),
  based on nearest-neighbour distance. Out-of-domain predictions are reported but explicitly discounted.

## VI. Separation of concerns preserved
- The existing design keeps each evidence layer (chemical, specificity, species, biology, resistance,
  structure, safety, uncertainty) as separate auditable fields. New evidence (3D shape/pharmacophore)
  is an additional layer with its own fields — never silently folded into an existing score.

## VII. Scientific scope discipline
- This tool ranks *testable target hypotheses*. It does not claim antibacterial activity, target
  engagement, or clinical relevance. Docking/MD/in-vivo remain downstream, separate evidence.
