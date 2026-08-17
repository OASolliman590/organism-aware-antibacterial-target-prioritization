# Target Prediction Pipeline — Methodology Notes

## Compounds (12)
- BI-1, BI-6: benzimidazole + thioacetamide phenylamide, OC(F)F (difluoromethoxy) substituent
- OX-11: methylxanthine (theobromine core) N6-thioacetate hydrazide + sulfonamide? no — theobromine linked via CH2 to 1,3,4-oxadiazole? Actually: Cn1cnc2c1c(=O)...n2C = theobromine; CSc1nnc(...)o1 = 1,3,4-oxadiazole-2-thiol linked to theobromine at position 6; acetylhydrazone to 2,5-dimethylphenyl
- T2Z5, T2Z6, T2Z9, T2Z14: 1,3-dimethylxanthine-7-acetylhydrazone of 1,2,4-triazolo[1,5-a]pyrimidine / 1,2,4-triazolo[4,3-a]pyrimidine; T2Z5/T2Z6 = 1,2,4-triazolo[1,5-a]pyrimidine with bromo/chlorophenyl; T2Z9/T2Z14 = triazolopyrimidine isomer with phenyl-substituted triazole, bromo/chlorophenyl
- X1V9, X1V11, X1V19: benzothiazole sulfonamide with aryl-hydrazone (X1V9 = o-tolyl, X1V11 = o-anisyl, X1V19 = p-cyanophenyl) — classic DHPS (sulfonamide) pharmacophore
- X1V20, X1V26: 4-amino-benzenesulfonamide + benzothiazole + hydrazone — PABA/DHPS mimics

## User target panels (may be wrong — pipeline must evaluate)
- K. pneumoniae: dhfr, fabl, topoisomerase iv, mura, lpxh (docking+MD done)
- B. cereus: fabl, dhfr, ftsz, gyrb, mura
- E. coli: mura, dhfr, fabl, gyrb, lpxc
- P. mirabilis: mura, dhfr, fabl, gyrb, lpxc
- A. baumannii: pbp2a?, fabl, gyrb, lpxa, murac (note: PBP2a is MRSA — user mislabeled?)
- MRSA: pbp2a, dhfr, fabl, ftsz, gyrb
- Prior CSV: compounds grouped by docking: BI-1/BI-6/T2Z5/T2Z9 → B. cereus beta-lactamase; OX-11/T2Z14 → K. pneumoniae KPC-2; T2Z6 → E. coli FtsZ/MurA; BI-6/X1V9 → AmpC; X1V9 → P. mirabilis class C beta-lactamase; X1V20 → OXA-23; X1V11 → LeuRS; X1V26 → FabH; X1V19 → MRSA BlaZ/DNA gyrase

## Prediction strategy (ligand-based, RDKit only)
1. **Reference ligand sets**: Collect known active ligands per target class from ChEMBL (via RDKit-free REST API: ChEMBL webservices REST API, public, no auth). Targets:
   - DHFR (CHEMBL226), FabI (CHEMBL1879), MurA (CHEMBL2093), MurC (CHEMBL4735), GyrB (CHEMBL216), FtsZ (CHEMBL5468/GTPase), LpxC (CHEMBL4217), LpxA (CHEMBL4186), LpxH (CHEMBL4277), Topo IV ParE/ParC (CHEMBL5910/5371), PBP2a/MecA (CHEMBL5861), FabH (CHEMBL2096), DHPS (CHEMBL211)
   - Pull confirmed actives (IC50/Ki < 1 uM, confidence >= 7)
2. **Scoring**: 
   - ECFP4/Tanimoto to each reference ligand → best/top-mean similarity per target
   - Morgan bit overlap / MACCS key similarity (justification figures)
   - Substructure/pharmacophore feature check: sulfonamide-aniline → DHPS; xanthine scaffold → no specific known; hydrazone/aroyl-hydrazone; 1,2,4-triazolo[1,5-a]pyrimidine
3. **Physicochemical/permeability filters**: Gram-negative (E. coli/K. pneumo/P. mirabilis/A. baumannii) need TPSA ≤ ~130 Å² for outer membrane permeation (Oprea rules / Veber; literature ~60-80 for good GN penetration, ≤100-120 tolerated). B. cereus/MRSA Gram-positive more permissive.
4. **Consensus ranking**: weighted composite of ligand similarity, pharmacophore evidence, permeability fit → shortlist top targets per organism.
5. **Validation**: compare predictions to prior docking/MD assignments (CSV) — agreement/discrepancy analysis.

## Figures required
- TMAP/MACS layout of 12 compounds + reference ligands colored by target class
- Heatmap: max Tanimoto per compound vs target class
- Scaffold/substructure panels (SAR)
- Physicochemical property plots (MW, TPSA, LogP) vs permeability thresholds
- Radar/justification summary per organism

## ChEMBL REST notes
- Base: https://www.ebi.ac.uk/chembl/api/data/
- target/{chembl_id}.json → uniProt, pref_name
- activity?target_chembl_id=X&standard_type=IC50&standard_value__lte=1000&pchembl_value__gte=6&limit=500
- molecule?molecule_chembl_id=X.json → canonical_smiles
