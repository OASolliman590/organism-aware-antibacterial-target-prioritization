"""Organism-aware target ranking and publication-quality figures.

Important interpretation rule: this is a ligand-based prioritization pipeline, not a
validated target-identification experiment. Similarity, MACCS agreement, and SAR flags
are reported separately so high rankings remain auditable.
"""
from pathlib import Path
import os
import json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Draw
from sklearn.decomposition import PCA

try:
    import umap
except Exception:
    umap = None

WORK = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parents[1]))
RES = WORK / 'results'
FIG = RES / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

TARGET_PANEL = {
    'Klebsiella pneumoniae': ['DHFR', 'FabI', 'TopoIV', 'MurA', 'LpxH'],
    'Bacillus cereus': ['FabI', 'DHFR', 'FtsZ', 'GyrB', 'MurA'],
    'Escherichia coli': ['MurA', 'DHFR', 'FabI', 'GyrB', 'LpxC'],
    'Proteus mirabilis': ['MurA', 'DHFR', 'FabI', 'GyrB', 'LpxC'],
    'Acinetobacter baumannii': ['PBP2a', 'FabI', 'GyrB', 'LpxA', 'MurC'],
    'MRSA / Staphylococcus aureus': ['PBP2a', 'DHFR', 'FabI', 'FtsZ', 'GyrB'],
}

GRAM_NEGATIVE = {'Klebsiella pneumoniae', 'Escherichia coli', 'Proteus mirabilis', 'Acinetobacter baumannii'}

TARGET_RATIONALE = {
    'DHFR': 'FolA/DHFR is an essential folate-cycle enzyme; evidence is strongest when diaminopyrimidine-like or folate-mimetic chemistry is present and the nearest active reference ligands are similar.',
    'DHPS': 'Sulfonamide-aniline/PABA mimicry is mechanistically relevant to DHPS, so this class is an important orthogonal rescue candidate for sulfonamide-bearing X1V compounds.',
    'FabI': 'FabI/FabL is the enoyl-ACP reductase of type-II fatty-acid synthesis; ranking is ligand-similarity based because the series lacks a canonical triclosan-like motif.',
    'FtsZ': 'FtsZ is an essential cytokinetic GTPase; prioritization is credible only when similarity to known FtsZ ligands is accompanied by a compatible heteroaromatic pharmacophore.',
    'GyrB': 'GyrB is the ATPase subunit of DNA gyrase; target evidence is based on similarity to GyrB reference ligands and should be distinguished from whole-gyrase/topoisomerase-IV predictions.',
    'Gyr': 'Whole DNA gyrase reference data provide a broader quinolone/topoisomerase control space; agreement with GyrB strengthens, but does not prove, a gyrase hypothesis.',
    'TopoIV': 'Topoisomerase IV is a DNA-decatenation complex; because reference chemistry is often quinolone-like, low similarity should be treated as weak evidence.',
    'MurA': 'MurA catalyzes the first committed enolpyruvyl-transfer step of peptidoglycan synthesis; the best-established covalent benchmark is fosfomycin, so non-electrophilic similarity is provisional.',
    'MurC': 'MurC and related Mur ligases are ATP-dependent peptidoglycan precursor enzymes; the reference set includes related Mur-ligase chemistry and is lower-confidence than FabI/GyrB.',
    'LpxC': 'LpxC is a Zn-dependent deacetylase in lipid-A biosynthesis; hydroxamate-like chelation is a strong mechanistic flag, while similarity alone is insufficient.',
    'LpxH': 'LpxH is a metal-dependent phosphodiesterase in lipid-A biosynthesis; target ranking is exploratory because the reference set is structurally heterogeneous.',
    'LpxA': 'LpxA is the acyltransferase initiating lipid-A biosynthesis; very sparse public ligand data make this a low-confidence computational control.',
    'PBP2a': 'PBP2a/MecA is a validated MRSA resistance determinant, but it is not an appropriate organism-specific target for A. baumannii; the A. baumannii panel entry is therefore flagged as biologically mismatched.',
    'FabH': 'FabH is a beta-ketoacyl-ACP synthase; sparse reference data make this a family-level exploratory target rather than a primary assignment.',
}

# Prior assignments from the user-provided CSV / previous work, used only for validation.
PRIOR = {
    'Klebsiella pneumoniae': ['KPC-2 beta-lactamase', 'KPC-2 beta-lactamase-avibactam', 'FtsZ'],
    'Bacillus cereus': ['Metallo-beta-lactamase / beta-lactamase', 'FtsZ'],
    'Escherichia coli': ['FtsZ', 'MurA', 'AmpC beta-lactamase'],
    'Proteus mirabilis': ['Beta-lactamase class-C family', 'MurA'],
    'Acinetobacter baumannii': ['OXA-23 beta-lactamase', 'LeuRS', 'FabH / 3-oxoacyl-ACP synthase III'],
    'MRSA / Staphylococcus aureus': ['BlaZ beta-lactamase', 'DNA gyrase'],
}


def load_compounds():
    out = []
    sdf = Chem.SDMolSupplier(str(WORK / 'data/compounds/compounds_normalized.sdf'), removeHs=True)
    for mol in sdf:
        if mol is not None:
            out.append(mol)
    return out


def fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def bv_array(bv):
    arr = np.zeros((bv.GetNumBits(),), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(bv, arr)
    return arr


def compute_flags(mol):
    patterns = {
        'sulfonamide': '[SX3](=O)(=O)[NX3]',
        'beta_lactam': 'C1C(=O)N1',
        'hydroxamate': 'C(=O)N[OX2H]',
        'diaminopyrimidine': '[nX2]1c[nX2]c(N)n1',
        'xanthine': 'O=[CX3]1[NX3][CX3](=O)[NX3][CX3]1',
        'hydrazone': 'C(=O)NN=C',
        'benzimidazole': 'c1c[nH]c2ccccc12',
        'benzothiazole': 'c1c2ccccc2nc1',
        'aryl_halide': '[a][F,Cl,Br]',
    }
    out = {}
    for k, smarts in patterns.items():
        query = Chem.MolFromSmarts(smarts)
        out[k] = int(query is not None and mol.HasSubstructMatch(query))
    return out


def organism_target_prior(organism, target):
    # Transparent, modest priors encode target-context biology without overpowering
    # ligand evidence. Values are not probabilities and are reported explicitly.
    priors = {
        'Klebsiella pneumoniae': {'TopoIV': 1.05, 'LpxH': 1.05},
        'Bacillus cereus': {'FtsZ': 1.05, 'GyrB': 1.03},
        'Escherichia coli': {'LpxC': 1.05, 'GyrB': 1.03},
        'Proteus mirabilis': {'LpxC': 1.05, 'GyrB': 1.03},
        'Acinetobacter baumannii': {'MurC': 1.05, 'GyrB': 1.03, 'PBP2a': 0.20, 'LpxA': 0.85},
        'MRSA / Staphylococcus aureus': {'PBP2a': 1.10, 'GyrB': 1.05, 'FtsZ': 1.03},
    }
    return priors.get(organism, {}).get(target, 1.0)


def add_pharmacophore_evidence(row, flags):
    cls = row['target_class']
    bonus = 0.0
    if cls == 'DHPS' and flags['sulfonamide']:
        bonus += 0.20
    if cls == 'PBP2a' and flags['beta_lactam']:
        bonus += 0.20
    if cls == 'LpxC' and flags['hydroxamate']:
        bonus += 0.20
    if cls == 'DHFR' and flags['diaminopyrimidine']:
        bonus += 0.10
    if cls in {'Gyr', 'GyrB', 'TopoIV'} and flags['xanthine']:
        bonus += 0.03
    if cls in {'FabI', 'FtsZ', 'Gyr', 'GyrB', 'TopoIV', 'MurA', 'MurC', 'LpxH'} and flags['hydrazone']:
        bonus += 0.02
    return bonus


def make_rankings():
    scores = pd.read_csv(RES / 'scores.csv')
    props = pd.read_csv(RES / 'compound_properties.csv')
    prop_map = props.set_index('compound').to_dict('index')
    sdf_mols = {m.GetProp('_Name'): m for m in load_compounds() if m.HasProp('_Name')}
    flags_map = {name: compute_flags(mol) for name, mol in sdf_mols.items()}
    rows = []
    for _, r in scores.iterrows():
        c = r['compound']
        flags = flags_map[c]
        sim = max(0.0, min(1.0, (r['ecfp4_max_tanimoto'] - 0.10) / 0.55))
        top5 = max(0.0, min(1.0, (r['ecfp4_top5_mean'] - 0.08) / 0.45))
        maccs = max(0.0, min(1.0, (r['maccs_max_tanimoto'] - 0.10) / 0.70))
        sar = add_pharmacophore_evidence(r, flags)
        evidence_score = 0.50 * sim + 0.25 * top5 + 0.15 * maccs + sar
        evidence_score = min(1.0, evidence_score)
        # Applied after the raw evidence components, and intentionally modest.
        # The PBP2a penalty for A. baumannii flags a biological mismatch in the
        # user-specified panel rather than claiming absence of all PBPs.
        rows.append({**r.to_dict(), **prop_map[c], **flags,
                     'similarity_component': sim, 'top5_component': top5,
                     'maccs_component': maccs, 'sar_bonus': sar,
                     'evidence_score': evidence_score,
                     'organism_prior': organism_target_prior('unassigned', r['target_class']),
                     'evidence_tier': ('high' if evidence_score >= 0.55 else 'moderate' if evidence_score >= 0.35 else 'low')})
    ranked = pd.DataFrame(rows)

    panel_rows = []
    for organism, classes in TARGET_PANEL.items():
        sub = ranked[ranked.target_class.isin(classes)].copy()
        for compound in ranked.compound.unique():
            q = sub[sub.compound == compound].copy()
            q['organism'] = organism
            q['organism_prior'] = q['target_class'].map(lambda t: organism_target_prior(organism, t))
            q['organism_adjusted_score'] = q['evidence_score'] * q['organism_prior']
            q = q.sort_values('organism_adjusted_score', ascending=False)
            q['panel_rank'] = np.arange(1, len(q) + 1)
            q['panel_target_count'] = len(q)
            q['organism_target_panel'] = '; '.join(classes)
            panel_rows.append(q)
    panel = pd.concat(panel_rows, ignore_index=True)
    panel.to_csv(RES / 'organism_panel_rankings.csv', index=False)

    # Shortlist = top 2 target classes per organism and compound, with explicit cutoff.
    ranked.to_csv(RES / 'ranked_target_evidence.csv', index=False)
    shortlist = panel[panel.panel_rank <= 2].copy()
    shortlist['justification'] = shortlist.apply(lambda x: TARGET_RATIONALE.get(x['target_class'], ''), axis=1)
    shortlist.to_csv(RES / 'organism_target_shortlist.csv', index=False)
    json.dump({'target_rationale': TARGET_RATIONALE, 'target_panels': TARGET_PANEL, 'prior_assignments': PRIOR},
              open(RES / 'target_context.json', 'w'), indent=2)
    return ranked, panel, shortlist, flags_map


def plot_compound_grid(mols):
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(350, 250),
                               legends=[m.GetProp('_Name') for m in mols], useSVG=False)
    img.save(FIG / 'compound_structure_grid.png')


def plot_similarity_heatmaps(mols, ranked):
    names = [m.GetProp('_Name') for m in mols]
    maccs = [MACCSkeys.GenMACCSKeys(m) for m in mols]
    mat = np.array([[DataStructs.TanimotoSimilarity(a, b) for b in maccs] for a in maccs])
    plt.figure(figsize=(10, 8))
    sns.heatmap(mat, xticklabels=names, yticklabels=names, cmap='viridis', vmin=0, vmax=1, square=True)
    plt.title('Pairwise MACCS-key similarity of the 12 compounds')
    plt.tight_layout(); plt.savefig(FIG / 'compound_maccs_similarity_heatmap.png', dpi=300); plt.close()

    pivot = ranked.pivot(index='compound', columns='target_class', values='ecfp4_max_tanimoto')
    pivot = pivot.loc[names]
    plt.figure(figsize=(15, 8))
    sns.heatmap(pivot, cmap='magma', vmin=0, vmax=0.8, annot=True, fmt='.2f', linewidths=.25)
    plt.title('Maximum ECFP4 similarity to active ChEMBL reference ligands')
    plt.xlabel('Candidate target class'); plt.ylabel('User compound')
    plt.tight_layout(); plt.savefig(FIG / 'compound_target_ecfp4_heatmap.png', dpi=300); plt.close()

    pivot = ranked.pivot(index='compound', columns='target_class', values='maccs_max_tanimoto').loc[names]
    plt.figure(figsize=(15, 8))
    sns.heatmap(pivot, cmap='crest', vmin=0, vmax=1.0, annot=True, fmt='.2f', linewidths=.25)
    plt.title('Maximum MACCS-key similarity to active ChEMBL reference ligands')
    plt.xlabel('Candidate target class'); plt.ylabel('User compound')
    plt.tight_layout(); plt.savefig(FIG / 'compound_target_maccs_heatmap.png', dpi=300); plt.close()


def plot_tmap_like(mols, ranked):
    # Combine user compounds and the top-3 reference ligand per target class.
    records = []
    labels = []
    groups = []
    for m in mols:
        records.append(bv_array(fp(m))); labels.append(m.GetProp('_Name')); groups.append('user')
    top_refs = json.load(open(RES / 'top_refs.json'))
    for cls, refs in top_refs.items():
        for i, ref in enumerate(refs[:3]):
            mol = Chem.MolFromSmiles(ref['canonical_smiles'])
            if mol is not None:
                records.append(bv_array(fp(mol))); labels.append(f'{cls}:{i+1}'); groups.append(cls)
    X = np.vstack(records)
    if umap is not None:
        coords = umap.UMAP(n_neighbors=min(12, len(X)-1), min_dist=0.25, metric='jaccard', random_state=42).fit_transform(X)
        method = 'UMAP on ECFP4/Jaccard distance'
    else:
        coords = PCA(n_components=2, random_state=42).fit_transform(X)
        method = 'PCA on ECFP4 bits'
    plt.figure(figsize=(14, 10))
    unique = ['user'] + sorted(set(groups) - {'user'})
    palette = sns.color_palette('tab20', n_colors=max(3, len(unique)))
    cmap = dict(zip(unique, palette))
    for group in unique:
        idx = [i for i, g in enumerate(groups) if g == group]
        plt.scatter(coords[idx, 0], coords[idx, 1], s=55 if group == 'user' else 24,
                    alpha=0.9 if group == 'user' else 0.45, label=group, c=[cmap[group]])
    for i, label in enumerate(labels[:len(mols)]):
        plt.text(coords[i, 0] + 0.04, coords[i, 1] + 0.04, label, fontsize=8, weight='bold')
    plt.title(f'TMAP-like target-reference chemical space ({method})')
    plt.xlabel('embedding dimension 1'); plt.ylabel('embedding dimension 2')
    plt.legend(ncol=3, fontsize=7, loc='best', frameon=True)
    plt.tight_layout(); plt.savefig(FIG / 'tmap_like_ecfp4_reference_map.png', dpi=300); plt.close()


def plot_property_space(props):
    props = props.copy()
    props['permeability_flag'] = np.where((props.tpsa > 130) & (props.compound.notna()), 'high TPSA / Gram-negative caution', 'not flagged')
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=props, x='logp', y='tpsa', size='mw', hue='permeability_flag', sizes=(70, 450), palette={'high TPSA / Gram-negative caution': '#d62728', 'not flagged': '#1f77b4'})
    for _, r in props.iterrows():
        plt.text(r.logp + .03, r.tpsa + 1.0, r.compound, fontsize=8)
    plt.axhline(130, color='grey', linestyle='--', linewidth=1, label='TPSA 130 Å² screening line')
    plt.title('Physicochemical space and Gram-negative permeability caution')
    plt.tight_layout(); plt.savefig(FIG / 'physchem_permeability_space.png', dpi=300); plt.close()


def plot_panel_rankings(panel):
    top = panel[panel.panel_rank <= 2].copy()
    top['label'] = top['compound'] + ' → ' + top['target_class']
    g = sns.catplot(data=top, x='organism_adjusted_score', y='compound', col='organism', col_wrap=2,
                    hue='target_class', kind='bar', height=4, aspect=1.4, sharex=True, legend='brief')
    g.set_titles('{col_name}')
    for ax in g.axes.flat:
        ax.axvline(.35, color='grey', linestyle='--', linewidth=.7)
        ax.set_xlabel('organism-adjusted target score')
        ax.set_ylabel('compound')
    if g._legend is not None:
        g._legend.set_title('target class')
    g.fig.suptitle('Top two organism-panel target hypotheses per compound', y=1.02)
    g.fig.tight_layout(); g.fig.savefig(FIG / 'organism_panel_target_rankings.png', dpi=300); plt.close()


def plot_prior_validation(panel):
    # Summarize direct prior overlaps at target-family level.
    mapping = {'FtsZ': 'FtsZ', 'MurA': 'MurA', 'DNA gyrase': 'GyrB', 'FabH / 3-oxoacyl-ACP synthase III': 'FabH'}
    rows = []
    for org, priors in PRIOR.items():
        for p in priors:
            matched = [k for k, v in mapping.items() if k in p]
            rows.append({'organism': org, 'prior_target': p, 'mapped_candidate': matched[0] if matched else 'outside scored panel', 'direct_overlap': int(bool(matched))})
    val = pd.DataFrame(rows)
    val.to_csv(RES / 'prior_validation_assignments.csv', index=False)
    summary = val.groupby('organism', as_index=False)['direct_overlap'].sum()
    plt.figure(figsize=(11, 5))
    sns.barplot(data=summary, x='organism', y='direct_overlap', color='#4c78a8')
    plt.xticks(rotation=35, ha='right'); plt.ylabel('number of prior assignments overlapping scored families')
    plt.title('Validation against prior docking/MD target families')
    plt.tight_layout(); plt.savefig(FIG / 'prior_docking_md_overlap.png', dpi=300); plt.close()


def main():
    ranked, panel, shortlist, flags = make_rankings()
    mols = load_compounds()
    props = pd.read_csv(RES / 'compound_properties.csv')
    plot_compound_grid(mols)
    plot_similarity_heatmaps(mols, ranked)
    plot_tmap_like(mols, ranked)
    plot_property_space(props)
    plot_panel_rankings(panel)
    plot_prior_validation(panel)
    print('ranked rows:', len(ranked), 'panel rows:', len(panel), 'shortlist rows:', len(shortlist))
    print('figures:', ', '.join(p.name for p in sorted(FIG.glob('*.png'))))

if __name__ == '__main__':
    main()
