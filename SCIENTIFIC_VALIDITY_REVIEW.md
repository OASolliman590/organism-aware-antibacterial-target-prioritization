# Scientific Validity Review

**Pipeline:** Organism-Aware Antibacterial Open Target Discovery (v3)
**Reviewed commit:** `d8d6954` (branch `ccr-8f3a6508-7s4zx1`)
**Date:** 2026-10-19
**Scope:** chemical-evidence engine, 3D/pharmacophore layer, rank fusion, organism-aware
priority, benchmarking/leakage control, calibration, and the stated interpretation limits.

---

## 1. Overall verdict

The pipeline is **methodologically honest and unusually well-guarded for a tool of this
class**. Its central scientific claim — that it produces *ranked, testable target
hypotheses*, not evidence of binding — is consistently enforced in code, documentation, and
output schema. The engineering around reproducibility (pinned snapshots, provenance
manifests, seeded determinism, leakage guards, "missing ≠ zero" discipline) is stronger than
most published cheminformatics target-prediction work.

The validity concerns below are therefore **not** "the method is wrong." They are places
where the headline numbers mean less than they appear to, where a stated principle is not
uniformly applied, or where a reader/reviewer could over-interpret the output. They are
ordered by how much they affect the scientific conclusions.

---

## 2. The load-bearing limitation: what the benchmark actually validates

**Finding (most important).** The entire chemical engine is nearest-neighbour similarity to
*known ligands of a target class* (ECFP4, MACCS, USRCAT, O3A shape/colour, 2D/3D
pharmacophore). By construction it measures "resembles a known binder of target X," which is
only a proxy for "binds target X."

The benchmark honestly separates three splits, but their scientific weight is very uneven:

- **Close-analogue split (≥0.85 ECFP4 guard):** this is the split that produces reportable
  retrieval numbers. It still permits references that are *just under* 0.85 Tanimoto to the
  query. For a similarity method, retrieval on this split is close to measuring
  memorisation of near-neighbours. It is the easy case.
- **Scaffold split (Bemis–Murcko):** reported transparently, but the docs already note that
  when the public benchmark contains no additional same-scaffold exclusions the numbers are
  identical to the analogue split — i.e. it is frequently *not actually testing* scaffold
  generalisation.
- **Target-family holdout:** the scientifically decisive test (can the method find the right
  target when it has *never seen a ligand of that family*?). The supplement correctly states
  a nearest-neighbour method "may consequently have zero positive-target coverage." That is
  not a data gap — it is the method's structural ceiling.

**Why it matters:** a reader scanning the top-1/top-3/MRR/EF numbers will attribute them to
the method's discriminative power. In reality they are dominated by the analogue split, which
is the split where leakage is *hardest to fully exclude* and generalisation is *least
tested*. **Recommendation:** in the abstract/README, state the headline retrieval numbers
*and the split they come from* in the same sentence, and label the analogue-split result as
"interpolation within known chemistry," not generalisation.

---

## 3. Specific validity issues

### 3.1 `species_transfer` conflates "missing mapping" with "zero transfer" in the priority score — contradicts the project's own missing≠zero rule

In `apply_biology` (`pipeline/open_target_discovery_v2.py`):

```python
transfer = float(c.iloc[0].species_transfer_score) if len(c) and pd.notna(...) else 0.0
mapping_status = c.iloc[0].sequence_status if len(c) else 'no_mapping_record'
```

A genuinely absent species mapping is assigned `transfer = 0.0` and then fed into the
multiplicative score as if it were a measured zero. Everywhere else the pipeline is strict
that *missing is not zero* (chem layers, fusion, calibration all exclude-and-count). Here the
principle is not applied: an unmapped target is treated as a maximally weak-but-real transfer
observation rather than as missing evidence.

The downstream effect is softened but not eliminated by the affine factor
(`config.yaml`: `transfer: {base: 0.65, weight: 0.35}`): an unmapped target keeps **65 %** of
its transfer factor. It cannot reach `High` confidence (gate requires `transfer ≥ 0.70`) but
it *can* rank highly and be labelled `Moderate`. So "we have no idea whether this organism
even has this target protein" costs a hypothesis only 35 % of one of four factors.

**Recommendation:** carry `no_mapping_record` as missing (NaN) through to the priority step
and either (a) exclude such rows from the ranked shortlist with a recorded reason, or (b)
apply a distinct, documented "unmapped" multiplier that is clearly harsher than
`base = 0.65`. Right now the `uncertainty_reasons` string flags it, but the *score* barely
does.

### 3.2 The organism-aware priority score is an unvalidated hand-weighted heuristic

```python
overall = chemical_quality
    * (0.65 + 0.35*transfer)
    * (0.75 + 0.25*pocket)
    * (0.50 + 0.50*biological)
    * (1 - 0.20*anti)
```

Every base/weight, the multiplicative (rather than additive or learned) composition, the
biology sub-weights, and the clinical/essentiality/accessibility ordinal maps are
investigator-chosen. There is no empirical fit for them because (as the supplement states) no
calibratable labelled set exists. This is **legitimate as declared decision-support**, and
the sensitivity analysis (`sensitivity_analysis.py`, Kendall-τ/RBO under ±25–50 % weight
perturbation and leave-one-layer-out) is the right mitigation and is a real strength.

But two things should be said more loudly:

- Sensitivity analysis shows the ranking is *robust to perturbing these weights*; it does
  **not** show the weights are *correct*. Robustness to perturbation of an arbitrary function
  is not evidence that the function tracks truth. The current wording risks being read as
  validation.
- Multiplicative composition means the factors interact non-linearly: a mediocre biological
  score and a mediocre pocket score compound. Whether that interaction is biologically
  intended (AND-like gating) or an artefact of the chosen form is not argued anywhere.

**Recommendation:** add one sentence stating explicitly that the priority weights are
expert-prior, not fitted, and that sensitivity analysis bounds *stability* not *accuracy*.

### 3.3 Early-recognition metrics (EF 1 %, BEDROC α=80.5) are near-degenerate at this universe size

The target universe is on the order of tens of target classes per query. In
`_tie_aware_enrichment`, `top_n = ceil(n_candidates * fraction)`. With ~30–60 candidates:

- **EF at 1 %** → `top_n = 1`. EF1% collapses to "is the single top-ranked candidate the
  active?" scaled by prevalence. It is not a 1 % enrichment in any meaningful sense.
- **BEDROC α=80.5** puts ~80 % of its weight in the top ~1.4 % of the list — again
  effectively the top one or two slots for these universe sizes, so it saturates and has
  very high variance.

The bootstrap CIs correctly *report* this variance, but a metric that is structurally pinned
to the top-1 slot should not be presented alongside EF5%/AUROC as if it adds independent
early-recognition information at this scale.

**Recommendation:** either report EF/BEDROC only when `n_candidates` is large enough for the
fraction to resolve more than ~1 slot, or annotate each EF/BEDROC cell with the effective
`top_n` so a reader sees EF1% == "top-1" for small queries. (The BEDROC formula itself is the
standard Truchon–Bayly implementation and is correct — this is about *applicability*, not
correctness.)

### 3.4 "Confidence class" is a label on an uncalibrated score, not statistical confidence

`High/Moderate/Low/Insufficient` are thresholds on `overall_priority_score`, which is
explicitly *not* a probability. Calling the top bin "High confidence" invites exactly the
over-interpretation the rest of the project works hard to prevent. A 0.50 `overall_min`
cutoff is a cut on an arbitrary heuristic, not a 50 % anything.

**Recommendation:** rename to something that cannot be read as calibrated certainty —
e.g. `priority_tier` with `Tier 1–4`, or `hypothesis_strength`. Keep the gating logic.

### 3.5 Assay heterogeneity in the reference sets

`target_scoring.py` and the reference loaders pool reference ligands and carry
`pchembl_value` / `standard_value` (IC50 nM) as "potency context" across what are, in ChEMBL,
heterogeneous assays, organisms, and measurement types. The `reference_quality_grade` layer
and pChEMBL filtering mitigate this, but the "best reference potency" fields shown next to a
hypothesis can imply a potency comparability that the underlying assays do not support.

**Recommendation:** label potency-context fields as assay-heterogeneous, or restrict the
displayed potency to a single assay-type/standard where possible.

### 3.6 Selection-conditioned 3D evidence (acknowledged, restated for completeness)

O3A alignment runs only on the USRCAT top-25 shortlist, so O3A-derived scores and their
bootstrap variance are conditional on a USRCAT pre-selection. The supplement discloses this.
It means the 3D-only and fused numbers are not independent of the shape prefilter — the "3D
adds information over 2D" comparison is partly "USRCAT+alignment over ECFP4," not a clean
orthogonal-evidence test. Fine as disclosed; worth keeping in the limitations of any figure
that contrasts 2D vs 3D vs fused.

### 3.7 Pharmacophore 3D normalisation and map-feature reuse (acknowledged)

The symmetric query↔reference coverage normalisation with FeatMap feature reuse is a
defensible, disclosed implementation choice. Reuse of a reference feature against multiple
probe features can inflate overlap for repeated pharmacophoric groups; the symmetric mean
dampens but does not remove this. No change required — the disclosure is adequate — but it
should remain in the method limitations, not only in code comments.

---

## 4. What is genuinely sound (and should be preserved)

- **Leakage control is real and asserted at runtime** (`>=0.85` analogue guard re-checked
  inside every split; scaffold removal counted). This is better than most published work.
- **Missing-data discipline** in the chem/fusion/calibration layers: `min_count`,
  exclude-and-count, no imputation, `*_is_probability = False` flags. Exemplary.
- **RRF fusion** divides by the *fixed* ideal contribution of all configured components, so a
  row with fewer available components cannot be rewarded for its missingness. Correct and
  subtle.
- **Calibration path refuses to emit a probability** unless disjoint train/calibration/test
  query roles exist (`assess_training_readiness`), and only the Platt output is ever named a
  probability. This is the single most important guard against the field's most common sin.
- **BEDROC/EF are tie-aware** (average ranks, fractional boundary handling) — many
  implementations silently mishandle ties and inflate early enrichment.
- **Negative-result policy**: empty/worse/no-qualifying-case outputs are authoritative and
  not back-filled. The disagreement and mode-comparison tables report "no gain" honestly.
- **Provenance**: snapshot hashing, create-only refresh writers, explicit `null` for
  unrecoverable legacy provenance rather than inferred values.

---

## 5. Prioritised recommendations

| # | Change | Type | Effort |
|---|--------|------|--------|
| 1 | State headline retrieval numbers *with their split* and label analogue-split results as interpolation, not generalisation (§2) | Framing | Low |
| 2 | Stop treating `no_mapping_record` as `transfer = 0.0` in the priority multiply; carry it as missing and downweight explicitly (§3.1) | **Code** | Low–Med |
| 3 | Rename `confidence_class` → a tier/strength label (§3.4) | Code/schema | Low |
| 4 | Annotate EF/BEDROC with effective `top_n`, or suppress sub-resolution fractions (§3.3) | Code | Low |
| 5 | Add one explicit sentence: priority weights are expert-prior, not fitted; sensitivity bounds stability not accuracy (§3.2) | Framing | Low |
| 6 | Flag potency-context fields as assay-heterogeneous (§3.5) | Framing/schema | Low |

Items 1, 3, 4, 5, 6 are documentation/labelling changes that raise scientific integrity with
near-zero risk. Item 2 is the one substantive code change I would make, because it is the
only place where a stated principle of the project is contradicted by the implementation in a
way that affects the ranking.

None of these undermine the pipeline's core validity as a **hypothesis-generation and
prioritisation framework**. They tighten the gap between what the numbers say and what a
reader will think they say.
