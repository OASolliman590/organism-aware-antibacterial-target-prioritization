# Implementation Plan — Spec 001

Maps requirements (spec.md) and method decisions (research.md) to concrete modules, data flow, and interfaces.
Preserves the existing v2 engine; v3 is additive.

## 1. Architecture overview

```
                         config.yaml  (paths, thresholds, seeds, fusion mode, splits)
                                 │
   data/ (pinned snapshots) ─────┼──────────────────────────────────────────────
   reference_ligands/*.json      │                                              │
   ontology/quality/compat/...   ▼                                              ▼
                         ┌───────────────────┐        ┌──────────────────────────────────┐
   queries SDF ─────────►│ chem2d (existing) │        │ chem3d (NEW)                     │
                         │ ECFP4/MACCS Tani  │        │ conformers→USRCAT→O3A shape/color│
                         └─────────┬─────────┘        │ + pharmacophore                  │
                                   │                  └───────────────┬──────────────────┘
                                   └──────────► fusion (NEW) ◄────────┘
                                                 │  chemical_evidence_score_v3 (+components)
                                                 ▼
                    apply_biology (existing v2) → species/biology/resistance/structure/anti-target
                                                 │
                                                 ▼
                          scoring/combiner  ──►  {heuristic+sensitivity}  OR  {learned+calibrated}
                                                 │        (config-selectable, FR-5)
                                                 ▼
             predictions_by_organism_v3.csv + applicability-domain flags + provenance manifest
                                                 │
                                                 ▼
                   benchmark harness (NEW): splits × metrics × CIs × decoys × external baseline
```

## 2. New / changed modules

| Module | Purpose | Key libs |
|---|---|---|
| `pipeline/chem3d_matching.py` (NEW) | conformers (ETKDGv3), USRCAT, O3A shape+color, pharmacophore; per-(query,class) 3D fields; disk cache | rdkit |
| `pipeline/evidence_fusion.py` (NEW) | rank/score fusion of 2D+3D+pharmacophore → `chemical_evidence_score_v3`, components retained | numpy/pandas |
| `pipeline/scoring_model.py` (NEW) | config-selectable combiner: (a) heuristic passthrough, (b) learned logistic/LTR + calibration | scikit-learn / lightgbm |
| `pipeline/sensitivity_analysis.py` (NEW) | weight perturbation + ablation; rank-stability metrics (Kendall τ, RBO), bootstrap | numpy/scipy |
| `pipeline/benchmark_v3.py` (NEW, extends benchmark_v2) | target/scaffold/time splits; AUROC/BEDROC/EF/MRR + bootstrap CIs; decoy handling | scikit-learn/scipy |
| `pipeline/baseline_external.py` (NEW) | run/parse an external tool (SEA or SwissTargetPrediction or PIDGIN) on benchmark queries | requests / local model |
| `pipeline/applicability_domain.py` (NEW) | NN-distance AD flags (Tanimoto + USRCAT) | rdkit/numpy |
| `pipeline/provenance.py` (NEW) | snapshot versions, commit, config hash, seeds → manifest.json | stdlib |
| `pipeline/config.py` (NEW) | load/validate declarative YAML; single source of thresholds & seeds | pydantic/pyyaml |
| `pipeline/open_target_discovery_v2.py` (EDIT, guarded) | consume fusion output when v3 enabled; keep v2 path intact behind config flag | — |
| `data/snapshots/` (NEW) | pinned, dated data snapshots + `SNAPSHOT_VERSIONS.json` | — |
| `tests/` (NEW) | unit + golden-file regression (v2 frozen) + tiny fixture set | pytest |

## 3. Data-flow contracts (see contracts/)
- `chem3d_matching` input: list of RDKit mols (queries) + reference dict {class: [mols]}; output: DataFrame with
  `usrcat_max, usrcat_top5_mean, o3a_shape_tanimoto_max, o3a_color_max,
  pharmacophore_2d_gobbi_sim_max, pharmacophore_3d_sim_max` per (query, class),
  plus a conformer-cache manifest. The legacy `pharmacophore_sim_max` alias
  remains the 2D Gobbi value for v2-compatible consumers.
- `evidence_fusion` input: 2D score frame + 3D frame; output: adds `chemical_evidence_score_v3` and keeps all
  component columns; deterministic given config.
- `benchmark_v3` output: one row per (split_type, metric) with point estimate + CI bounds + n + split provenance.

## 4. Config (single YAML) — illustrative keys
```yaml
run: { seed: 20240601, fusion_mode: rank_fusion, combiner: heuristic }  # or combiner: learned
chem3d: { n_confs: 30, prune_rms: 0.5, energy_window_kcal: 10, o3a_shortlist_top: 25 }
benchmark: { splits: [target_family, scaffold, temporal], time_cutoff: "2018-01-01",
             bootstrap_n: 1000, decoys: property_matched }
applicability_domain: { tanimoto_in: 0.4, tanimoto_out: 0.25 }
snapshots: { chembl: "ChEMBL_34", card: "2024-XX", uniprot: "2024_0X", rcsb_query_date: "2024-XX-XX" }
```

## 5. Phasing (delivery order)
1. **Reproducibility spine** — config, provenance, snapshot pinning, tests + golden v2 freeze. (Everything else
   depends on determinism; do this first.)
2. **3D matching layer** — chem3d + fusion, additive fields, cached, with the 2D-vs-3D disagreement report.
3. **Benchmark v3** — splits, metrics, CIs, decoys.
4. **External baseline** — head-to-head.
5. **Calibration/model + sensitivity** — learned combiner OR heuristic-stability report.
6. **Applicability domain** + docs/figures/methods-supplement update.

Each phase is independently reviewable and leaves the pipeline runnable.

## 6. Risks & mitigations
- **Sparse labels for a learned combiner** → fall back to constitution-compliant heuristic + sensitivity (FR-5b).
- **O3A cost on full reference set** → USRCAT pre-shortlist before O3A (research A3).
- **External baseline API limits / offline** → cache results; if unavailable, report baseline as pending, do not
  fabricate (constitution I/III).
- **Conformer nondeterminism** → fixed embed seed + pinned RDKit version in the manifest.
- **Scope creep** → docking/MD explicitly out of scope (spec §4).

## 7. Definition of done (per spec §7)
All success criteria met; `python run_pipeline.py --config config.yaml` reproduces every reported number from the
pinned snapshot and emits a provenance manifest; methods supplement updated; tests green in CI.
