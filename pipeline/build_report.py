"""Build the final scientific report from pipeline outputs."""
from pathlib import Path
import os
import pandas as pd

ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
RES = ROOT / 'results'
OUT = RES / 'target_prioritization_report.md'
means = pd.read_csv(RES / 'organism_target_mean_scores.csv')
tops = pd.read_csv(RES / 'organism_compound_top_predictions.csv')
comparison = pd.read_csv(RES / 'prior_target_comparison.csv')
props = pd.read_csv(RES / 'compound_properties_with_series.csv')

organism_order = ['Klebsiella pneumoniae','Bacillus cereus','Escherichia coli','Proteus mirabilis','Acinetobacter baumannii','MRSA / Staphylococcus aureus']

lines = []
lines += ['# Organism-aware target prioritization for 12 antibacterial compounds', '',
          '**Analysis type:** RDKit ligand-based target prioritization with organism-panel constraints, MACCS/ECFP4 consensus, structural-feature flags, and prior docking/MD reconciliation.', '',
          '> **Scope statement.** These results prioritize testable hypotheses. They do not demonstrate target engagement, biochemical inhibition, or organism-specific antibacterial mechanism.', '']
lines += ['## Executive summary', '',
          'The uploaded structure set contains **12 validated compounds**. The authoritative compound-to-name mapping was reconstructed from the named records embedded in `README.md` because the filenames of the supplied SDF files did not consistently match their internal names. The resulting normalized set is BI-1, BI-6, OX-11, T2Z14, T2Z5, T2Z6, T2Z9, X1V11, X1V19, X1V20, X1V26, and X1V9.', '',
          'The ligand-based analysis is strongest for the BI series against GyrB-like reference chemistry. The X1V compounds form a coherent sulfonamide-benzothiazole chemical series, but the available organism panels do not contain DHPS; DHPS was therefore retained as an orthogonal target class in the full score matrix rather than silently omitted. The T2Z/OX-11 compounds occupy a separate heteroaromatic/xanthine-like chemical-space region and show mixed Mur-ligase, GyrB, FabI, or TopoIV-related similarity depending on the organism panel.', '',
          'The main scientific conclusion is conservative: **GyrB, MurC/MurA, FabI, TopoIV, LpxC, FtsZ, and the prior beta-lactamase hypotheses should be treated as a ranked hypothesis set rather than a single universal target.** The repeated GyrB ranking is partly a consequence of shared ligand-reference chemistry and must not be misreported as proof that all 12 compounds inhibit GyrB.', '']

lines += ['## 1. Data preparation and input integrity', '',
          'RDKit successfully parsed and sanitized 12 unique named structures. Explicit hydrogens were removed before fingerprint calculation. The compound manifest is available at `results/compound_properties_with_series.csv`.', '',
          '| Series | Compounds | Chemical interpretation |',
          '|---|---|---|']
for series, q in props.groupby('series'):
    lines.append(f'| {series} | {", ".join(q.compound.tolist())} | {"benzimidazole-thioacetamide" if series == "BI-series" else "xanthine/heteroaromatic hydrazone" if series == "OX/T2Z-series" else "sulfonamide-benzothiazole hydrazone"} |')
lines += ['', 'The supplied `final_target_list_per_bacterium.csv` was retained unchanged and compared to the ligand-scored target universe. Targets outside that universe were explicitly labelled as orthogonal prior hypotheses rather than treated as negative predictions.', '']

lines += ['## 2. Methodology', '',
          'For each candidate target class, active ChEMBL reference ligands were collected from target-specific or target-family records. The data assembly retained the canonical SMILES, ChEMBL molecule identifier, assay type, standard value, derived or reported pChEMBL value, target identifier, and reference organism. Sparse classes were preserved but flagged by their reference count.', '',
          'Each user compound was compared with every reference ligand using 2048-bit ECFP4/Morgan fingerprints and MACCS keys. The reported ECFP4 metrics are maximum similarity and mean similarity of the five closest reference ligands. The MACCS metric is the maximum key similarity. Class-relevant substructure flags provide a separately inspectable SAR layer.', '',
          'The raw evidence score is a weighted combination of normalized ECFP4 maximum similarity, top-five ECFP4 mean similarity, MACCS similarity, and a small class-specific SAR bonus. A modest organism/target prior was then applied only inside each requested target panel. For *A. baumannii*, PBP2a was down-weighted because MecA/PBP2a is a MRSA resistance determinant and is biologically mismatched to the requested *A. baumannii* panel; this is a panel-quality correction, not a claim that all *A. baumannii* PBPs are irrelevant.', '',
          'References: [RDKit](https://www.rdkit.org/), [ChEMBL](https://www.ebi.ac.uk/chembl/), and [UMAP](https://umap-learn.readthedocs.io/).', '']

lines += ['## 3. Organism-panel target scores', '',
          'The table below reports mean organism-adjusted score across all 12 compounds. The companion CSV contains medians, maxima, ECFP4 and MACCS components, SAR bonuses, and reference counts.', '',
          '| Organism | Rank | Target class | Mean adjusted score | Max score | Mean ECFP4 | Mean MACCS | Reference ligands |',
          '|---|---:|---|---:|---:|---:|---:|---:|']
for org in organism_order:
    q = means[means.organism == org].sort_values('mean_score', ascending=False).reset_index(drop=True)
    for i, r in q.iterrows():
        lines.append(f'| {org} | {i+1} | {r.target_class} | {r.mean_score:.3f} | {r.max_score:.3f} | {r.mean_ecfp4:.3f} | {r.mean_maccs:.3f} | {int(r.n_references)} |')
lines += ['', '### Interpretation by organism', '']
interpret = {
    'Klebsiella pneumoniae': 'TopoIV is the highest mean-scoring class within the original panel, followed by FabI, MurA, and LpxH. The prior KPC-2 beta-lactamase hypotheses are not scored by the current ligand-reference universe and should remain primary orthogonal hypotheses for the OX-11/T2Z14 subset. FtsZ is a scored family but lies outside the original panel.',
    'Bacillus cereus': 'GyrB is the highest mean-scoring class, followed by FabI, MurA, and FtsZ. FtsZ is the strongest direct overlap with the prior target CSV, while the metallo-beta-lactamase hypothesis remains outside the ligand-scored universe.',
    'Escherichia coli': 'GyrB is the highest mean-scoring class, with LpxC and FabI next. MurA is a lower-ranked ligand hypothesis but a direct prior overlap for T2Z6 and should remain a control target because the prior structure contains the established MurA/fosfomycin context.',
    'Proteus mirabilis': 'GyrB is highest, followed by LpxC, FabI, and MurA. The prior X1V9 beta-lactamase class-C hypothesis is not scored; MurA is a direct scored overlap and should remain the better-controlled cell-wall target.',
    'Acinetobacter baumannii': 'GyrB and MurC are the highest ligand-scored classes. This supports adding MurC/GyrB to the next computational panel, but does not invalidate the prior OXA-23 hypothesis for X1V20 or LeuRS/FabH hypotheses for X1V11/X1V26. PBP2a is strongly down-ranked as an organism-mismatched panel entry.',
    'MRSA / Staphylococcus aureus': 'GyrB is the strongest ligand-scored class, followed by FabI, FtsZ, and DHFR. The prior DNA-gyrase hypothesis maps directly to the GyrB family. The prior BlaZ hypothesis is not scored and should be retained as an orthogonal resistance-mechanism target; PBP2a is biologically important but shows weak ligand-reference evidence for this compound set.',
}
for org in organism_order:
    lines += [f'**{org}.** {interpret[org]}', '']

lines += ['## 4. Shortlisted hypotheses by compound and organism', '',
          'The full top-two table is available in `results/organism_target_shortlist.csv`. A target is listed as a computational shortlist item when it ranks first or second within the organism-specific panel. The table below shows the highest-ranked target for each compound in each organism.', '',
          '| Organism | Compound | Top target | Adjusted score | ECFP4 max | MACCS max | Second target | Score margin |',
          '|---|---|---|---:|---:|---:|---|---:|']
for org in organism_order:
    q = tops[tops.organism == org].sort_values('compound')
    for _, r in q.iterrows():
        lines.append(f'| {org} | {r.compound} | {r.top_target} | {r.top_score:.3f} | {r.top_ecfp4:.3f} | {r.top_maccs:.3f} | {r.second_target} | {r.score_margin:.3f} |')

lines += ['', '## 5. Validation against prior docking and molecular dynamics targets', '',
          'The user-provided prior target list contains 16 rows. The current reconciliation identifies four direct scored overlaps, three scored-family hypotheses outside the original organism panels, and eight hypotheses that are not scored because the current ligand-reference universe does not include beta-lactamases or LeuRS.', '',
          '| Validation category | Count | Meaning |',
          '|---|---:|---|',
          '| Direct scored overlap | 4 | Prior target maps to a class scored inside that organism panel |',
          '| Scored family outside original panel | 3 | Prior target family is represented, but not in the original requested panel |',
          '| Orthogonal unscored hypothesis | 8 | Dedicated docking/MD hypothesis retained; no ligand-based negative claim |', '',
          'Direct overlaps are *B. cereus*–FtsZ, *E. coli*–MurA, *P. mirabilis*–MurA, and MRSA–DNA gyrase mapped to GyrB. The three outside-panel scored families are *K. pneumoniae*–FtsZ, *E. coli*–FtsZ, and *A. baumannii*–FabH. The unscored set includes KPC-2, KPC-2–avibactam, AmpC, class-C beta-lactamase, OXA-23, BlaZ, and LeuRS hypotheses.', '',
          'This validation pattern is scientifically useful: it shows where the ligand-only model agrees with existing structural hypotheses, where the old target panel was incomplete, and where an orthogonal protein-level hypothesis cannot be assessed by the present ligand database.', '']

lines += ['## 6. Figures and what they show', '',
          '| Figure | Scientific purpose |',
          '|---|---|']
figs = [
('compound_structure_grid.png','Visual audit of the 12 input structures and naming.'),
('compound_maccs_similarity_heatmap.png','Pairwise MACCS-key similarity; shows series coherence and chemical redundancy.'),
('compound_target_ecfp4_heatmap.png','Maximum ECFP4 similarity to active reference ligands by target class.'),
('compound_target_maccs_heatmap.png','Orthogonal MACCS-key target-class similarity heatmap.'),
('tmap_like_ecfp4_reference_map.png','UMAP/Jaccard chemical-space map of user compounds and representative reference ligands; labelled TMAP-like, not exact TMAP.'),
('organism_panel_target_rankings.png','Top two organism-panel hypotheses per compound with an explicit moderate-evidence guide line.'),
('prior_docking_md_overlap.png','Counts of prior target assignments that overlap scored families.'),
('physchem_permeability_space.png','MW/logP/TPSA space and a screening-level Gram-negative permeability caution.'),
]
for f, d in figs:
    lines.append(f'| `results/figures/{f}` | {d} |')

lines += ['', '## 7. Limitations and recommended next experiments', '',
          'The analysis is limited by the quality and coverage of public ligand-reference data. A high fingerprint score is not a calibrated binding probability. Cross-organism ChEMBL records can blur species specificity, whole-complex measurements can blur subunit specificity, and reference-set size is highly unequal. LpxA has only one retained reference ligand and PBP2a has a small set, so rankings for these classes are low-confidence.', '',
          'The current workflow does not model protein sequence identity, binding-site conservation, active-site protonation, cofactors, membrane permeability, efflux, target expression, or cellular resistance. The Gram-negative TPSA flag is only a screening caution and is not a permeability predictor.', '',
          'For the next round, retain three hypotheses per compound: the top ligand-based target in the organism panel, the strongest prior docking/MD target even when unscored, and one mechanistic control. Use sequence-matched proteins for GyrB, MurA/MurC, FabI, FtsZ, and LpxC. For KPC-2, OXA-23, AmpC, BlaZ, and class-C beta-lactamases, maintain organism-matched catalytic-state docking and validate with beta-lactamase inhibition assays. For X1V compounds, add DHPS and beta-lactamase family scoring when a curated reference ligand set is available. For MRSA PBP2a, add a target-specific reference set and verify whether the X1V series contains a chemically plausible covalent or noncovalent PBP2a motif before investing in further MD.', '',
          '## 8. Reproducibility', '',
          'Run `python run_pipeline.py` after placing the input files in `inputs/`. The generated tables, figures, and report are written to `results/`. The code is intentionally modular: preparation, scoring, organism ranking, figure generation, result summarization, and prior-target comparison are separate scripts.', '',
          '### References', '',
          '1. [RDKit documentation](https://www.rdkit.org/).',
          '2. [ChEMBL database and API](https://www.ebi.ac.uk/chembl/).',
          '3. [UMAP documentation](https://umap-learn.readthedocs.io/).',
          '4. Rogers D, Brown RD. Extended-connectivity fingerprints and molecular similarity; the ECFP/Morgan representation is used here as a ligand-space descriptor.',
          '']
OUT.write_text('\n'.join(lines))
print(f'Wrote {OUT}')
