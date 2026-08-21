# Research & Method Decisions — Spec 001

Grounding for every technical choice. Methods named here are established and citable; verify exact
reference formatting against the primary sources before manuscript submission. Where a decision is a
judgement call, the rationale and the alternative are given.

---

## A. Why add 3D matching (the scientific case)

2D fingerprint (ECFP4/Morgan) Tanimoto measures topological substructure overlap. It is fast and
interpretable but (i) penalises scaffold hops that preserve the pharmacophore, and (ii) rewards
topological similarity that does not imply shape/electrostatic complementarity with a binding site.
For **target identification**, the biologically relevant question is whether the query can occupy the
same pocket as reference actives — a shape + pharmacophore question. Adding a 3D layer directly targets
the weakest, load-bearing part of the current method.

### A1. Conformer generation
- **ETKDG (v3)** distance-geometry with experimental torsion knowledge — RDKit `EmbedMultipleConfs` with
  `ETKDGv3()`. Riniker & Landrum, *J. Chem. Inf. Model.* 2015 (ETKDG); v3 adds small-ring/macrocycle terms.
- Ensemble policy: fixed random seed; bounded N (e.g. 20–50) pruned by RMSD and an energy window
  (MMFF94 relative energy) to keep only low-energy, diverse conformers. Determinism required (constitution IV).
- Cache conformers on disk keyed by canonical SMILES + params; regeneration only on cache miss.

### A2. Alignment-free 3D similarity — **USRCAT**
- Ultrafast Shape Recognition with CREDO Atom Types (USRCAT), Schreyer & Blundell, *J. Cheminform.* 2012.
  RDKit `rdMolDescriptors.GetUSRCAT` + `GetUSRScore`.
- Encodes molecular shape *and* pharmacophoric atom categories (hydrophobic, aromatic, H-bond donor/acceptor)
  as a 60-D moment descriptor; similarity is alignment-free and extremely fast → suitable for scanning all
  reference ligands. Good first-pass 3D filter.

### A3. Alignment-based 3D similarity — **Open3DAlign (O3A)** + shape/color Tanimoto
- Open3DALIGN (O3A), Tosco, Balle, Shiri, *J. Comput. Aided Mol. Des.* 2011. RDKit
  `rdMolAlign.GetO3A` (MMFF or Crippen atom typing) aligns query conformers onto each reference, giving an
  alignment score; on the aligned pair compute **shape Tanimoto** (`rdShapeHelpers.ShapeTanimotoDist`) and a
  pharmacophore/"color" overlap. This is the closer analogue to ROCS-style shape+color scoring using only
  open-source tooling. Use it on the USRCAT-shortlisted references (expensive → shortlist first).

### A4. Pharmacophore layer
- 3D pharmacophore feature similarity via RDKit `ChemicalFeatures` (BaseFeatures fdef) + Gobbi/Poppe
  2D pharmacophore fingerprints (`Gobbi_Pharm2D`) as a complementary, alignment-free pharmacophore signal.
- Rationale: pharmacophore matching captures interaction-type complementarity independent of scaffold.

### A5. Fusion of 2D + 3D + pharmacophore
- Keep each component as a separate field (constitution VI). Fuse into `chemical_evidence_score_v3` via
  **rank-based data fusion** (e.g. Reciprocal Rank Fusion, Cormack et al. 2009) or normalized-score fusion.
  Rank fusion is robust to differing score scales and is the safer default over hand-weighted sums.
- Report per-component contributions so a reviewer can see whether 3D changed a call (e.g. a scaffold hop
  that 2D missed but shape/pharmacophore recovered) — this disagreement case is itself a result (FR success #1).

### A6. Honest expectation
- 3D methods do not universally beat 2D on retrieval; performance is dataset-dependent. The deliverable is a
  *fair* 2D-vs-3D-vs-fusion comparison under the benchmark, reported whichever way it falls (constitution III).

---

## B. Calibration & the weight problem

The current `overall_priority = chem_quality·(0.65+0.35·transfer)·(0.75+0.25·pocket)·(0.50+0.50·biology)·(1−0.20·anti)`
and the chemical normalisations `(ecfp4−0.10)/0.55` etc. are arbitrary. Two acceptable resolutions:

### B1. Preferred — learn the combiner from labelled benchmark data
- Frame as ranking/classification against the ESKAPE mechanism-labelled benchmark: features = the separate
  evidence fields (2D, 3D, specificity, reference quality, species transfer, biology, resistance, pocket,
  anti-target); label = correct mechanism target for each benchmark drug.
- Model: start with **regularised logistic regression** (interpretable coefficients → still auditable) or
  **learning-to-rank** (LambdaMART/`lightgbm` LGBMRanker) if per-query ranking is the target. Keep it simple
  and interpretable for a methods paper.
- **Probability calibration:** Platt scaling or isotonic regression on a held-out fold; report calibration
  curves (reliability diagram) and Brier score. Only then may outputs be called probabilities (constitution II).
- Train/evaluate strictly under the Section C splits — no leakage.

### B2. Minimum acceptable — retain heuristic + prove stability
- If labelled data are too sparse for a stable learned model, keep the transparent heuristic but add a
  **sensitivity/ablation analysis**: perturb each weight (e.g. ±25–50%, and one-at-a-time leave-one-layer-out),
  recompute rankings, and report the rank stability of the top hypotheses (e.g. rank-biased overlap, Kendall τ,
  bootstrap over reference ligands). Publish the stability, not just the point ranking.

---

## C. Leakage-controlled benchmarking (the credibility core)

### C1. Splits (report each separately; constitution III)
- **Target-family holdout:** hold out all references of a target family; can the method still retrieve it from
  chemistry alone? Tests generalisation across targets.
- **Bemis–Murcko scaffold holdout:** exclude references sharing the query's Murcko scaffold. Bemis & Murcko,
  *J. Med. Chem.* 1996. Tests scaffold-level generalisation (the current benchmark admits this may be trivial —
  fix it by making the exclusion real and reporting how many references it removes).
- **Temporal (time-split):** train on references before a cutoff date, test on later — the most realistic proxy
  for prospective use (Sheridan, *J. Chem. Inf. Model.* 2013).

### C2. Metrics with uncertainty
- **AUROC**, **BEDROC** (Truchon & Bayly, *J. Chem. Inf. Model.* 2007; report α=20 and α=80.5),
  **Enrichment Factor** at 1% and 5%, **MRR**, coverage. Early-recognition metrics (BEDROC/EF) matter more than
  AUROC for prioritisation. Report **bootstrap 95% CIs** for all.

### C3. Decoys / negative controls
- Use **property-matched decoys** (DUD-E methodology: Mysinger, Carchia, Irwin, Shoichet, *J. Med. Chem.* 2012)
  or DEKOIS 2.0 to avoid analog/artificial enrichment. Also retain the built-in cross-target decoys for the
  specificity margin. The observed negative control (KPC-2 → zero chemical support for OX-11/T2Z14) is a good
  real example to report.

### C4. External baseline (required)
- Compare on the same queries against at least one established tool: **SEA** (Keiser et al., *Nat. Biotechnol.*
  2007), **SwissTargetPrediction** (Daina, Michielin, Zoete, *Nucleic Acids Res.* 2019), or **PIDGINv4**
  (Mervin et al.). Report head-to-head; underperformance is acceptable if analysed honestly.

---

## D. Applicability domain
- Flag each prediction by nearest-neighbour Tanimoto (and USRCAT distance) to the reference set: in-domain /
  near / out-of-domain, with a documented threshold. Out-of-domain predictions (likely for these novel
  scaffolds) are reported but discounted. Standard QSAR AD practice (Netzeva et al., 2005; Sahigara et al., 2012).

---

## E. Reproducibility engineering
- **Data pinning:** record ChEMBL release (e.g. ChEMBL_XX), CARD version + download date, UniProt release,
  RCSB query date; store the fetched snapshots (or their hashes) and never overwrite silently.
- **Config:** single declarative YAML (paths, thresholds, ensemble size, fusion mode, split params, seeds).
- **Determinism:** global seed; conformer embedding seed; bootstrap seed; any model seed.
- **Provenance manifest** per run: data versions, git commit, config hash, seeds, timestamps, environment.
- **Tests:** unit tests for scoring/fusion/metrics; a golden-file regression test pinning v2 outputs (FR-10);
  a tiny fixture dataset for fast CI.
- **Packaging:** move from loose scripts + env-vars to an installable package with `pyproject.toml`, typed
  public functions, and structured logging.

---

## F. Framing for publishability
- Position the contribution honestly as **organism-aware, multi-evidence target-hypothesis prioritization with
  explicit uncertainty and applicability domain** — the organism-transfer + evidence-separation + honesty are the
  differentiators versus generic 2D target-prediction tools. The 3D layer + calibration + leakage-controlled
  benchmark are what move it from "internal heuristic" to "method a reviewer can trust."
- Report the negative/uncertain results (everything Low/Insufficient for the novel compounds; KPC-2 unsupported)
  as evidence of calibration, not weakness.
