# Feature Spec 002 — Honest species-mapping missingness in organism-aware priority

**Status:** Draft for implementation
**Owner:** (author) → Opus 5.5 spec, Haiku 5.5 implementation
**Depends on:** v2 engine (`pipeline/open_target_discovery_v2.py`), applicability-domain
module (`pipeline/applicability_domain.py`), constitution (`memory/constitution.md`)
**Origin:** `SCIENTIFIC_VALIDITY_REVIEW.md` §3.1

## 1. Problem statement (WHY)

In `apply_biology` (`pipeline/open_target_discovery_v2.py`), when an organism has **no
species-mapping record** for a target class (or the record carries a non-finite transfer
score), the code assigns:

```python
transfer = float(c.iloc[0].species_transfer_score) if len(c) and pd.notna(...) else 0.0
```

That `else 0.0` **imputes a measured value where the measurement is absent**. The imputed
zero then flows into the multiplicative priority as the transfer factor
`base + weight*transfer = 0.65 + 0.35*0 = 0.65`, so an *unmapped* target keeps 65 % of its
transfer factor and is indistinguishable, in both the stored evidence field and the ranking,
from a target that was genuinely mapped and found to have zero transferability.

This violates **Constitution Principle I** (missing data is recorded as missing, never
imputed) — the principle the chemical, fusion, and calibration layers already honour via
exclude-and-count / NaN discipline. It also lets "we do not know whether this organism even
has this protein" rank as `Moderate`.

## 2. Objective (WHAT)

Make species-mapping missingness **honest and auditable** without fabricating a value and
without introducing arbitrary new score magnitudes:

- **A.** A genuinely absent/non-finite species transfer is stored as **missing (`NaN`)** in
  the output, never as `0.0`. A real, mapped transfer of `0.0` remains `0.0` and stays
  distinguishable from missing.
- **B.** The multiplicative priority uses an **explicit, configurable** transfer factor for
  the missing case (`overall_factors.transfer.missing_factor`) instead of an implicit
  `base + weight*0`. The audit trail shows which factor was applied.
- **C.** Unmapped rows are **ordered after mapped rows** in organism shortlists, mirroring
  the existing applicability-domain `ordering_only_no_score_rescaling` convention
  (`pipeline/applicability_domain.py`). No score is rescaled. Applicability-domain
  eligibility remains the **primary** ordering key (it gates the reliability of the core
  chemical signal); species-mapping presence is a **secondary** tiebreaker applied within the
  same AD-eligibility tier.
- **D.** Every prediction row carries a boolean `species_mapping_present` flag and the
  existing `uncertainty_reasons` text continues to name the gap.

## 3. In scope

- `pipeline/open_target_discovery_v2.py` — `apply_biology` transfer handling, new output
  columns, and the two organism shortlist sorts (v2 `main()` and v3 `emit_v3_outputs`).
- `pipeline/applicability_domain.py` — a backward-compatible optional ordering hook so the
  v3 shortlist can prioritise mapped rows in addition to AD-eligible rows.
- `config.yaml` — add `V2_SCORING.overall_factors.transfer.missing_factor`.
- `docs_v2_evidence_schema.md` — document the NaN semantics and the new flag column.
- `tests/` — a new unit test for the missing-vs-zero-vs-mapped behaviour.

## 4. Out of scope / non-goals

- **No change to `score_query`** (the frozen v2 chemical path) — the golden test
  (`tests/test_v2_golden.py`) must remain byte-identical.
- No change to the chemical-evidence score, fusion, calibration, or benchmark code.
- No change to the numeric magnitude of the transfer factor for *mapped* rows, nor to the
  real-zero case.
- No new "harsher" discount beyond the configurable `missing_factor` (default preserves the
  current magnitude; tuning is the author's later decision, not this spec's).

## 5. Behavioural contract (acceptance criteria)

Let a prediction row be produced by `apply_biology` for an `(organism, target_class)`:

| Case | `species_transfer_score` (stored) | `species_mapping_present` | transfer factor used in `overall` | shortlist ordering |
|---|---|---|---|---|
| Mapped, finite transfer `t` | `t` (incl. `0.0`) | `True` | `base + weight*t` | with eligible rows |
| No mapping record | `NaN` | `False` | `missing_factor` | after mapped rows in same AD tier |
| Mapping record, non-finite transfer | `NaN` | `False` | `missing_factor` | after mapped rows in same AD tier |

Additional invariants:

1. **No imputation:** a missing transfer is never written as a finite number anywhere in any
   emitted table. (Constitution I)
2. **Distinguishability:** a mapped `transfer == 0.0` row and a missing-mapping row are
   distinguishable by `species_mapping_present` and by `species_transfer_score`
   (`0.0` vs `NaN`).
3. **Default magnitude preserved:** with `missing_factor` left at its default (`0.65`, equal
   to the mapped real-zero factor `base + weight*0`), the `overall_priority_score` of an
   unmapped row is numerically unchanged from today; only its stored transfer field and its
   shortlist position change. This keeps the change reviewable (score effects are opt-in via
   config, not silent).
4. **Confidence gate intact:** a missing-mapping row still cannot be `High`
   (`transfer_min=0.70` gate: `NaN >= 0.70` is `False`), can still be `Moderate`, and still
   carries the `species sequence mapping unresolved or weak` reason.
5. **Golden unchanged:** `tests/test_v2_golden.py` passes unmodified.
6. **AD shortlist unchanged by default:** `shortlist_with_applicability_domain` called without
   the new optional argument behaves exactly as before (existing
   `tests/test_applicability_domain.py` stays green).

## 6. Risks

- **Silent score drift.** Mitigated by criterion 3: default `missing_factor` preserves the
  current magnitude, so the only behavioural deltas are honest-missingness and ordering.
- **`NaN` propagation breaking the multiply.** The multiply must use a local numeric factor
  (`missing_factor`) for the missing case; `NaN` is stored in the *output field* only, never
  multiplied. Covered by criteria 1–3 and the new test.
- **Downstream readers assuming `species_transfer_score` is always numeric.** Checked: no
  pipeline/test consumer requires a finite value here; figures and exports treat it as
  optional. The new test and the full suite guard this.
