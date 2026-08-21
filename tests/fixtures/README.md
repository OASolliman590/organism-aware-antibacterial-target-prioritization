# Tiny public test fixture

These records are exact copies from the pinned v2 public snapshot; they are not
synthetic examples. `queries.json` contains ciprofloxacin and trimethoprim rows from
`data/benchmark/eskape_benchmark_drugs.csv`. `references.json` contains the first two
tracked ligand records from each of `ref_ligands_GyrB.json`, `ref_ligands_DHFR.json`,
and `ref_ligands_FabI.json`. The fixture is intentionally too small for scientific
performance claims and is used only for fast deterministic interface tests.
