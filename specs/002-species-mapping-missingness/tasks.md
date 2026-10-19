# Tasks — Spec 002 (ordered, agent-executable)

Rules: obey `memory/constitution.md`. Make only the changes described in `plan.md`. Do not
touch `score_query` or any chemical/fusion/benchmark/calibration code. Never impute a missing
value. Run the gates in T5 before reporting done. Preserve exact public behaviour except
where `spec.md` §5 requires a change.

- [ ] **T1 — config key.** Add `missing_factor: 0.65` to
      `V2_SCORING.overall_factors.transfer` in `config.yaml` (plan Change 1).
      **AC:** `load_config().value("v2_scoring.overall_factors.transfer.missing_factor")`
      returns `0.65` (note the config key is lowercase `v2_scoring`); YAML still parses.

- [ ] **T2 — honest missingness + explicit factor in `apply_biology`.** Apply plan Changes
      2a, 2b, 2c in `pipeline/open_target_discovery_v2.py`.
      **AC:**
      - missing/non-finite mapping ⇒ `species_transfer_score` is `NaN` (not `0.0`) and
        `species_mapping_present is False`;
      - mapped finite transfer (including `0.0`) ⇒ value preserved and
        `species_mapping_present is True`;
      - `overall` never evaluates to `NaN` solely because transfer is missing;
      - with default config the `overall_priority_score` for an unmapped row equals the
        pre-change value (within 1e-9).

- [ ] **T3 — shortlist ordering.** Apply plan Changes 3 and 4 (optional param on
      `shortlist_with_applicability_domain`; pass `species_mapping_present` from the v3 and
      v2 `overall_priority_score` shortlists).
      **AC:** within an `(organism, query_id)` group, every `species_mapping_present==True`
      row with a given AD-eligibility sorts ahead of any `False` row; calling the helper
      without the new argument is byte-identical to before.

- [ ] **T4 — unit test.** Add `tests/test_species_mapping_missingness.py` that calls
      `open_target_discovery_v2.apply_biology` on a minimal in-memory `raw` scored frame, a
      minimal ontology frame, and a `compat` frame that supplies a mapping for one
      `(organism, target_class)` with transfer `0.0`, and omits the mapping for another.
      **AC (assertions):**
      1. the omitted-mapping row has `pd.isna(species_transfer_score)` and
         `species_mapping_present == False`;
      2. the real-zero row has `species_transfer_score == 0.0` and
         `species_mapping_present == True`;
      3. both rows have a finite `overall_priority_score`;
      4. the omitted-mapping row's `confidence_class` is not `High`;
      5. the omitted-mapping row's `uncertainty_reasons` contains the species-mapping phrase;
      6. a sanity check that, after sorting a two-row group with
         `additional_priority_columns=["species_mapping_present"]`, the mapped row precedes
         the unmapped row.
      Build the `raw`/`ontology`/`compat` columns by inspecting the fields `apply_biology`
      reads (e.g. `target_class`, `chemical_evidence_score`, `chemical_quality_adjusted_score`,
      `reference_quality_grade`, `target_specificity_score`, and the ontology columns
      `organism_scope`, `clinical_status`, `essentiality_level`, `cellular_localization`,
      `resistance_relevance`, `target_role`, plus whatever `parent_class` needs). Keep the
      fixture tiny and offline. No RDKit required for this test.
      **Simplification:** `apply_biology` loops over every organism in the module-level
      `ORGANISMS` list for the one target class in `raw`. So a `compat` frame that supplies a
      mapping for exactly **one** `(organism, target_class)` pair (with
      `species_transfer_score == 0.0`, `sequence_status == 'mapped'`) and omits all other
      organisms yields, in a single `apply_biology` call, one mapped row (cases 2) and the
      remaining organisms as unmapped rows (cases 1). Pick a `target_class` with no
      `parent_target_class` in the ontology (so `parent_class` returns the class itself and
      the parent-fallback lookup cannot accidentally map the "unmapped" organisms). Import
      `ORGANISMS` from the module to assert on a specific organism name.

- [ ] **T5 — docs + gates.**
      - Update `docs_v2_evidence_schema.md` per plan Change 5.
      - Run the gates and paste results into the final report:
        1. `python -m pytest tests/test_v2_golden.py -q` → must pass (golden unchanged).
        2. `python -m pytest tests/test_applicability_domain.py tests/test_species_mapping_missingness.py -q` → must pass.
        3. `python -m pytest -q --basetemp=.pytest_tmp` → full suite must pass (baseline was
           134 passed; expect 135 with the new test).
      **AC:** all three gates green; report the exact counts.

## Notes for the implementer
- The virtualenv with pinned deps is at `/tmp/claude-0/venv`. Run pytest as
  `/tmp/claude-0/venv/bin/python -m pytest ...`.
- Do not reformat untouched lines; `apply_biology` is dense single-line-dict style — match it.
- If any gate fails, fix the cause; do not weaken an assertion or skip a test to get green
  (constitution III).
