# Plan 002 — Honest species-mapping missingness (HOW)

Implements `spec.md` 002. Keep each change minimal and local; preserve all existing public
behaviour except where the behavioural contract requires a change.

## Change 1 — `config.yaml`

Under `V2_SCORING.overall_factors.transfer`, add `missing_factor`:

```yaml
  overall_factors:
    transfer: {base: 0.65, weight: 0.35, missing_factor: 0.65}
```

Rationale: `missing_factor` defaults to the mapped real-zero factor (`base + weight*0 =
0.65`) so the priority magnitude is unchanged by default (spec criterion 3). It is a separate
key so the author can later discount missing mappings by editing one value.

## Change 2 — `pipeline/open_target_discovery_v2.py`, `apply_biology`

### 2a. Transfer resolution (replace the `else 0.0` imputation)

Current:
```python
transfer=float(c.iloc[0].species_transfer_score) if len(c) and pd.notna(c.iloc[0].species_transfer_score) else 0.0
mapping_status=c.iloc[0].sequence_status if len(c) else 'no_mapping_record'
```

Replace with logic equivalent to:
```python
if len(c) and pd.notna(c.iloc[0].species_transfer_score):
    transfer = float(c.iloc[0].species_transfer_score)   # real value, including 0.0
    transfer_present = True
else:
    transfer = float("nan")                               # missing, never imputed
    transfer_present = False
mapping_status = c.iloc[0].sequence_status if len(c) else 'no_mapping_record'
```

`transfer` (possibly `NaN`) is what gets written to the `species_transfer_score` output
field (it already is — the `rr.update({... 'species_transfer_score':transfer ...})` line
stays, now carrying `NaN` when missing).

### 2b. Overall multiply — explicit missing factor

Where `overall` is computed, replace the transfer term with a precomputed factor:
```python
factors = V2_SCORING["overall_factors"]
if transfer_present:
    transfer_factor = float(factors["transfer"]["base"]) + float(factors["transfer"]["weight"]) * transfer
else:
    transfer_factor = float(factors["transfer"].get("missing_factor", factors["transfer"]["base"]))
overall = float(
    chemical_quality
    * transfer_factor
    * (float(factors["pocket"]["base"]) + float(factors["pocket"]["weight"]) * pocket)
    * (float(factors["biology"]["base"]) + float(factors["biology"]["weight"]) * biological)
    * (1 - float(factors["anti_target_penalty"]) * anti)
)
```
The `.get(..., base)` fallback keeps the function working if an older `config.yaml` lacks the
new key.

### 2c. New output column

In the `rr.update({...})` payload add:
```python
'species_mapping_present': bool(transfer_present),
```
Keep `species_transfer_score` and `sequence_mapping_status` as they are (the former now
carries `NaN` when missing).

Note: the confidence-class block already handles the missing case correctly —
`mapping_status != 'mapped'` fires the reason, and `NaN >= 0.70` / `NaN < threshold` are both
`False`, so no edit is needed there. Do **not** change the confidence thresholds.

## Change 3 — `pipeline/applicability_domain.py`, `shortlist_with_applicability_domain`

Add a backward-compatible optional parameter so callers can prioritise additional boolean
"eligible-first" columns ahead of the score, without disturbing AD ordering:

```python
def shortlist_with_applicability_domain(
    predictions, *, group_columns, score_column, top_n,
    additional_priority_columns: Sequence[str] = (),
):
    ...
    priority = [*additional_priority_columns, "ad_shortlist_eligible"]
    present = [col for col in priority if col in predictions.columns]
    ordered = predictions.sort_values(
        [*group_columns, *present, score_column],
        ascending=[*[True]*len(group_columns), *[False]*len(present), False],
        kind="mergesort",
    )
    return ordered.groupby(list(group_columns), sort=False).head(top_n)
```

Default `()` → identical behaviour to today (spec criterion 6). Only include columns that are
actually present so a missing column never raises.

## Change 4 — shortlist call sites order mapped rows first

- **v3** (`emit_v3_outputs`, `open_target_shortlist_by_organism_v3.csv`): pass
  `additional_priority_columns=["species_mapping_present"]` to
  `shortlist_with_applicability_domain`.
- **v2 `main()`** (`v2_open_target_shortlist_by_organism.csv`, the `ranked.sort_values([...])`
  line): add `species_mapping_present` as a sort key immediately after the group columns and
  before `overall_priority_score`, ascending=False, using a stable sort. Example:
  ```python
  ranked.sort_values(
      ['organism','query_id','species_mapping_present','overall_priority_score'],
      ascending=[True,True,False,False], kind='mergesort',
  ).groupby(['organism','query_id']).head(10).to_csv(...)
  ```
  Only do this for the `overall_priority_score` shortlist; leave the v2.1
  clinical-translation shortlist and the validation-priority CSV ordering unchanged (those
  are out of scope).

## Change 5 — `docs_v2_evidence_schema.md`

Document: `species_transfer_score` is `NaN` when no finite mapping exists (never `0.0`); add
the `species_mapping_present` boolean; state that unmapped rows are ordered after mapped rows
in organism shortlists (ordering only, no score rescaling), and that `missing_factor`
controls the missing-case transfer factor (default preserves the mapped real-zero magnitude).

## Change 6 — test

Add `tests/test_species_mapping_missingness.py` driving `apply_biology` on a tiny synthetic
frame + ontology + compatibility covering the three cases in the spec table. See tasks.md
T4 for the required assertions.

## Determinism / constitution check

- No new randomness, no network, no imputation (I ✓). No new probability claim (II ✓). No
  benchmark/leakage change (III ✓). No snapshot change (IV ✓). AD convention reused, not
  bypassed (V ✓). Transfer stays a separate auditable field; new flag is additive (VI ✓).
  No scope/claim change (VII ✓).
