# Public ESKAPE antibacterial benchmark

This benchmark evaluates open chemical target discovery on **16 known antibacterial drugs** with PubChem structures and conservative mechanism-level labels. It is intended to test whether ligand-space evidence retrieves chemically compatible target families without using a fixed target panel.

## Drug and target coverage

The benchmark includes fluoroquinolones, trimethoprim, sulfamethoxazole, fosfomycin, triclosan, ceftaroline, rifampicin, ribosome-active drugs, vancomycin, daptomycin, meropenem, and colistin. The target labels are broad mechanism classes where the drug mechanism is naturally complex or multi-protein.

The current ChEMBL reference universe supports DHFR, DHPS, FabI, GyrB, TopoIV, MurA, FtsZ, LpxA/LpxC/LpxH, MurC, FabH, PBP2a, RpoB, LeuRS, D-Ala-D-Ala ligase, and 70S ribosome. Some benchmark mechanisms, including general membrane action and broad PBP inhibition, may have no suitable reference class in the current snapshot. These cases are reported as **coverage limitations**, not false negative predictions.

## Leakage control

For every benchmark query, reference molecules with ECFP4/Morgan Tanimoto similarity of at least 0.85 to the query are excluded before scoring. Exact or close analogues therefore cannot inflate the benchmark result. The query itself is never used as a positive reference.

## Metrics

The pipeline reports whether the known target family is represented, rank of the best accepted target-family alias, top-1/top-3/top-5 retrieval, reciprocal rank, number of candidate target classes, and top predicted target. The metrics are retrieval statistics rather than calibrated probabilities.

## Sources

Structures are fetched reproducibly from the [PubChem PUG REST service](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial). Target-family reference activities are retrieved from [ChEMBL](https://www.ebi.ac.uk/chembl/). Mechanism labels should be checked against current product labels and primary pharmacology literature before publication.
