"""V2 open target discovery for protected local compounds and public benchmarks.

The output keeps chemical evidence, reference quality, species transfer, biological
priority, anti-target risk, and overall prioritization as separate auditable fields.
Cross-target molecules are decoys for specificity analysis only; they are not called
experimentally inactive.
"""
import json
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys

try:
    from pipeline.applicability_domain import (
        assign_applicability_domain,
        shortlist_with_applicability_domain,
    )
    from pipeline.chem3d_matching import (
        aggregate_reference_evidence,
        score_reference_evidence_by_target,
    )
    from pipeline.config import load_config
    from pipeline.disagreement_report import build_disagreement_report
    from pipeline.evidence_fusion import fuse_evidence
    from pipeline.execution import ordered_process_map
    from pipeline.snapshots import verify_snapshot
except ModuleNotFoundError:  # direct ``python pipeline/<script>.py`` execution
    from applicability_domain import (
        assign_applicability_domain,
        shortlist_with_applicability_domain,
    )
    from chem3d_matching import (
        aggregate_reference_evidence,
        score_reference_evidence_by_target,
    )
    from config import load_config
    from disagreement_report import build_disagreement_report
    from evidence_fusion import fuse_evidence
    from execution import ordered_process_map
    from snapshots import verify_snapshot


CONFIG = load_config()
ROOT = CONFIG.root
REF = CONFIG.path_for("reference_ligands")
ONTO = CONFIG.path_for("target_ontology")
ONTO_SUB = CONFIG.path_for("target_subtype_ontology")
QUALITY = CONFIG.path_for("reference_quality")
COMPAT = CONFIG.path_for("species_compatibility")
CARD_SUM = CONFIG.path_for("card_resistance_summary")
CARD_SNP = CONFIG.path_for("card_snp_summary")
CARD_SNP_ORG = CONFIG.path_for("card_snp_organism_summary")
STRUCT_SUM = CONFIG.path_for("structure_summary")
ANN_OLD=ROOT/'data'/'target_annotations.csv'
BENCH = CONFIG.path_for("benchmark")
BENCH_STRUCTURES = CONFIG.path_for("benchmark_structures")
PRIVATE_COMPOUNDS = CONFIG.path_for("private_compounds")
RES = CONFIG.path_for("results"); RES.mkdir(parents=True, exist_ok=True)

ORGANISMS = list(CONFIG.value("organisms.names"))
ORGANISM_ALIASES = dict(CONFIG.value("organisms.aliases"))
GRAM_NEG = set(CONFIG.value("organisms.gram_negative"))
GRAM_POS = set(CONFIG.value("organisms.gram_positive"))
CHEM2D = CONFIG.value("chem2d")
V2_SCORING = CONFIG.value("v2_scoring")


def fp(m):
    return AllChem.GetMorganGenerator(
        radius=int(CHEM2D["fingerprint_radius"]),
        fpSize=int(CHEM2D["fingerprint_bits"]),
    ).GetFingerprint(m)
def maccs(m): return MACCSkeys.GenMACCSKeys(m)
def sim(a,b): return float(DataStructs.TanimotoSimilarity(a,b))
def mol(s):
    x=Chem.MolFromSmiles(str(s)) if s and str(s)!='nan' else None
    return Chem.RemoveHs(x) if x else None


def load_refs():
    refs={}
    for p in sorted(REF.glob('ref_ligands_*.json')):
        if p.stem.endswith('summary'): continue
        cls=p.stem.replace('ref_ligands_',''); cls='GyrB' if cls=='Gyr' else cls
        for r in json.loads(p.read_text()):
            m=mol(r.get('canonical_smiles_standardized') or r.get('canonical_smiles'))
            if m is None: continue
            smi=Chem.MolToSmiles(m)
            refs.setdefault(cls,[]).append({**r,'_mol':m,'_fp':fp(m),'_maccs':maccs(m),'_smi':smi})
    out={}
    for cls,rows in refs.items():
        seen=set(); out[cls]=[]
        for r in rows:
            if r['_smi'] in seen: continue
            seen.add(r['_smi']); out[cls].append(r)
    return out


def load_queries(path):
    out=[]
    if not path.exists(): return out
    for m in Chem.SDMolSupplier(str(path),removeHs=True):
        if m is None: continue
        name=m.GetProp('_Name') if m.HasProp('_Name') else 'private_compound'
        out.append({'query_id':name,'query_name':name,'mol':m,'fp':fp(m),'maccs':maccs(m),'source':'private_local'})
    return out


def load_ontology():
    base=pd.read_csv(ONTO)
    if ONTO_SUB.exists():
        base=pd.concat([base,pd.read_csv(ONTO_SUB)],ignore_index=True).drop_duplicates('target_class',keep='first')
    return base


def parent_class(target, ann):
    value=ann.get('parent_target_class') if hasattr(ann,'get') else None
    if value is not None and pd.notna(value) and str(value) not in {'','nan'}:
        return str(value)
    return target


def quality_score(grade):
    return float(V2_SCORING["reference_quality"].get(str(grade), 0.0))

def clinical_value(x):
    x=str(x).lower()
    values = V2_SCORING["clinical"]
    if 'approved' in x: return float(values["approved"])
    if 'clinical' in x: return float(values["clinical"])
    if 'validated' in x: return float(values["validated"])
    if 'preclinical' in x: return float(values["preclinical"])
    return float(values["default"])

def ordinal(x):
    values = V2_SCORING["ordinal"]
    return float(values.get(str(x).lower(), values["default"]))

def access_value(x):
    values = V2_SCORING["accessibility"]
    return float(values.get(str(x).lower(), values["default"]))

def scope_score(scope,org,target):
    values = V2_SCORING["organism_scope"]
    s=str(scope).lower()
    if target=='PBP2a' and org!='Staphylococcus aureus': return float(values["incompatible_pbp2a"])
    if 'mrsa' in s and org=='Staphylococcus aureus': return float(values["exact_or_group_match"])
    if 'gram-negative' in s and org in GRAM_NEG: return float(values["exact_or_group_match"])
    if 'gram-positive' in s and org in GRAM_POS: return float(values["exact_or_group_match"])
    if 'broad' in s or 'bacteria' in s: return float(values["broad_bacterial"])
    if 'staphylococci' in s and org=='Staphylococcus aureus': return float(values["exact_or_group_match"])
    if org in s: return float(values["exact_or_group_match"])
    return float(values["default"])

def anti_target_annotation(target):
    # Annotation-only risk prior; no claim of measured human off-target activity.
    values = V2_SCORING["anti_target"]
    high={'DHFR':('human DHFR homolog; inspect selectivity','high'),'LeuRS':('human LARS homolog; inspect selectivity','high'),'70S_ribosome':('mitochondrial translation selectivity','medium'),'30S_ribosome':('mitochondrial translation selectivity','medium'),'50S_ribosome':('mitochondrial translation selectivity','medium')}
    if target.startswith('30S_') or target.startswith('50S_'):
        return (float(values["medium"]),'mitochondrial translation selectivity; annotation-only risk','annotation_only')
    if target in high: note,r=high[target]; return (float(values[r]),note,'annotation_only')
    return (float(values["default"]),'no direct human orthologue risk assigned in ontology; safety remains untested','annotation_only')

def validation_plan(target,role):
    direct={'GyrB','TopoIV','FtsZ','DHFR','FabI','FabH','MurA','MurC','MurE','LpxA','LpxC','LpxH','RpoB','LeuRS','PBP','PBP2a','PBP1A','PBP1B','PBP2','PBP2B','PBP2X','PBP3','PBP4'}
    if target in direct:
        return 'purified-protein inhibition or binding; species-orthologue assay; resistant-mutant or complementation test'
    if target.startswith('Beta-lactamase_'):
        return 'enzyme inhibition plus antibiotic-rescue assay; resistance-isolate panel; target-dependence control'
    if 'ribosome' in target or target.startswith('30S_') or target.startswith('50S_') or target=='D-Ala-D-Ala':
        return 'target-complex binding or biochemical translation assay; cell-based target-dependence experiment'
    return 'phenotypic assay followed by mechanism-specific orthogonal validation'

def chemical_evidence_score(ecfp4_max, ecfp4_top5_mean, maccs_max):
    """Compute the legacy v2 ranking score from declarative parameters."""
    chemical = V2_SCORING["chemical"]
    normalization = chemical["normalization"]
    weights = chemical["weights"]
    values = {
        "ecfp4_max": ecfp4_max,
        "ecfp4_top5_mean": ecfp4_top5_mean,
        "maccs_max": maccs_max,
    }
    components = {
        name: np.clip(
            (float(values[name]) - float(params["offset"]))
            / float(params["scale"]),
            0,
            1,
        )
        for name, params in normalization.items()
    }
    return float(
        np.clip(
            sum(float(weights[name]) * components[name] for name in weights),
            0,
            1,
        )
    )


def score_query(q,refs,quality,compat,ontology,exclude_close=False,cutoff=None):
    cutoff = (
        float(CHEM2D["close_analogue_cutoff"]) if cutoff is None else float(cutoff)
    )
    rows=[]
    for cls,rs in refs.items():
        kept=[]; excluded=0
        for r in rs:
            if exclude_close and sim(q['fp'],r['_fp'])>=cutoff: excluded+=1
            else: kept.append(r)
        if not kept: continue
        e=np.array([sim(q['fp'],r['_fp']) for r in kept]); k=np.array([sim(q['maccs'],r['_maccs']) for r in kept])
        order=np.argsort(-e); topn=min(int(CHEM2D["top_k"]),len(e)); best=kept[int(order[0])]
        # Cross-target molecules are decoys for specificity only.
        other=[]
        for other_cls,other_rs in refs.items():
            if other_cls==cls: continue
            other.extend(other_rs)
        decoy_max=max([sim(q['fp'],r['_fp']) for r in other],default=0.0)
        margin=float(e[order[0]]-decoy_max)
        qrow=quality[quality.target_class==cls]
        if len(qrow): qr=qrow.iloc[0]; qscore=quality_score(qr.quality_grade); nref=int(qr.n_valid_ligands); nscaf=int(qr.n_unique_scaffolds); grade=qr.quality_grade
        else: qscore=0.0; nref=len(kept); nscaf=0; grade='insufficient'
        chem=chemical_evidence_score(float(e[order[0]]),float(e[order[:topn]].mean()),float(k.max()))
        specificity_config = V2_SCORING["specificity"]
        specificity=float(np.clip(
            float(specificity_config["baseline"])
            + float(specificity_config["margin_weight"])
            * margin
            / float(specificity_config["margin_scale"]),
            0,
            1,
        ))
        chem_quality=float(chem*qscore*specificity)
        rows.append({'query_id':q['query_id'],'source':q['source'],'target_class':cls,
                     'ecfp4_max':float(e[order[0]]),'ecfp4_top5_mean':float(e[order[:topn]].mean()),'maccs_max':float(k.max()),
                     'cross_target_decoy_max':decoy_max,'target_specificity_margin':margin,'target_specificity_score':specificity,
                     'chemical_evidence_score':chem,'reference_quality_score':qscore,'chemical_quality_adjusted_score':chem_quality,
                     'reference_quality_grade':grade,'n_references_after_exclusion':len(kept),'n_close_references_excluded':excluded,
                     'n_unique_scaffolds':nscaf,'best_reference_molecule':best.get('molecule_chembl_id',''),'best_reference_organism':best.get('organism',''),
                     'query_target_label':q.get('target_class',''),'query_mechanism_class':q.get('mechanism_class',''),'query_organisms':q.get('organisms','')})
    return pd.DataFrame(rows)


def _references_for_v3(q, refs, *, exclude_close=False, cutoff=None):
    """Apply the exact v2 analogue-exclusion rule before any 3D comparison."""

    threshold = (
        float(CHEM2D["close_analogue_cutoff"])
        if cutoff is None
        else float(cutoff)
    )
    selected = {}
    for target_class, records in refs.items():
        kept = [
            record
            for record in records
            if not (
                exclude_close and sim(q["fp"], record["_fp"]) >= threshold
            )
        ]
        if kept:
            selected[target_class] = kept
    return selected


def score_query_v3(
    q,
    refs,
    quality,
    compat,
    ontology,
    *,
    config=CONFIG,
    exclude_close=False,
    cutoff=None,
    cache_dir=None,
    return_reference_evidence=False,
):
    """Add 3D/pharmacophore evidence and deterministic fusion to one v2 frame."""

    chem2d = score_query(
        q,
        refs,
        quality,
        compat,
        ontology,
        exclude_close=exclude_close,
        cutoff=cutoff,
    )
    selected_refs = _references_for_v3(
        q, refs, exclude_close=exclude_close, cutoff=cutoff
    )
    if chem2d.empty or not selected_refs:
        if return_reference_evidence:
            return chem2d, pd.DataFrame()
        return chem2d
    reference_evidence = score_reference_evidence_by_target(
        q["query_id"], q["mol"], selected_refs, config, cache_dir=cache_dir
    )
    chem3d = aggregate_reference_evidence(reference_evidence, config)
    fused = fuse_evidence(chem2d, chem3d, config)
    fused["chemical_quality_adjusted_score_v3"] = (
        fused["chemical_evidence_score_v3"]
        * fused["reference_quality_score"]
        * fused["target_specificity_score"]
    )
    fused["chemical_quality_adjusted_score_v3_is_probability"] = False
    fused = assign_applicability_domain(fused, config)
    return (fused, reference_evidence) if return_reference_evidence else fused


_EMIT_V3_WORKER_CONTEXT = None


def _initialize_emit_v3_worker(
    refs, quality, compat, ontology, config, cache_dir
):
    global _EMIT_V3_WORKER_CONTEXT
    _EMIT_V3_WORKER_CONTEXT = {
        "refs": refs,
        "quality": quality,
        "compat": compat,
        "ontology": ontology,
        "config": config,
        "cache_dir": cache_dir,
    }


def _score_emit_v3_worker(task):
    """Score one independent query in an initialized spawned process."""

    if _EMIT_V3_WORKER_CONTEXT is None:
        raise RuntimeError("v3 scoring worker context was not initialized")
    dataset_scope, query = task
    context = _EMIT_V3_WORKER_CONTEXT
    return score_query_v3(
        query,
        context["refs"],
        context["quality"],
        context["compat"],
        context["ontology"],
        config=context["config"],
        exclude_close=dataset_scope == "benchmark",
        cache_dir=context["cache_dir"],
    )

def add_unscored_classes(scored,q,ontology):
    """Keep the broad ontology visible: missing references are explicit, not absent targets."""
    present=set(scored.target_class) if not scored.empty else set()
    missing=[]
    for cls in ontology.target_class.astype(str):
        if cls in present: continue
        missing.append({'query_id':q['query_id'],'source':q['source'],'target_class':cls,
            'ecfp4_max':0.0,'ecfp4_top5_mean':0.0,'maccs_max':0.0,
            'cross_target_decoy_max':0.0,'target_specificity_margin':0.0,'target_specificity_score':0.0,
            'chemical_evidence_score':0.0,'reference_quality_score':0.0,'chemical_quality_adjusted_score':0.0,
            'reference_quality_grade':'insufficient','n_references_after_exclusion':0,'n_close_references_excluded':0,
            'n_unique_scaffolds':0,'best_reference_molecule':'','best_reference_organism':'',
            'query_target_label':q.get('target_class',''),'query_mechanism_class':q.get('mechanism_class',''),
            'query_organisms':q.get('organisms','')})
    return pd.concat([scored,pd.DataFrame(missing)],ignore_index=True) if missing else scored


def apply_biology(
    raw,
    ontology,
    compat,
    card_summary=None,
    snp_summary=None,
    snp_org=None,
    struct_summary=None,
    *,
    chemical_evidence_column="chemical_evidence_score",
    chemical_quality_column="chemical_quality_adjusted_score",
):
    if raw.empty: return raw
    ann=ontology.set_index('target_class'); rows=[]
    for _,r in raw.iterrows():
        if r.target_class not in ann.index: continue
        a=ann.loc[r.target_class]
        chemical_evidence = float(r[chemical_evidence_column])
        chemical_quality = float(r[chemical_quality_column])
        parent=parent_class(r.target_class,a)
        for org in ORGANISMS:
            c=compat[(compat.organism==org)&(compat.target_class==r.target_class)] if not compat.empty else pd.DataFrame()
            transfer_source=r.target_class
            if c.empty and parent != r.target_class and not compat.empty:
                c=compat[(compat.organism==org)&(compat.target_class==parent)]
                transfer_source=parent
            transfer=float(c.iloc[0].species_transfer_score) if len(c) and pd.notna(c.iloc[0].species_transfer_score) else 0.0
            mapping_status=c.iloc[0].sequence_status if len(c) else 'no_mapping_record'
            scope=scope_score(a.organism_scope,org,r.target_class)
            clinical=clinical_value(a.clinical_status); essential=ordinal(a.essentiality_level); access=access_value(a.cellular_localization); resistance=ordinal(a.resistance_relevance)
            resistance_family={'GyrB':'GyrB/TopoIV resistance','TopoIV':'GyrB/TopoIV resistance','RpoB':'RpoB resistance','PBP2a':'PBP2a','DHFR':'DHFR resistance','DHPS':'DHFR/DHPS resistance','D-Ala-D-Ala':'D-Ala-D-Ala resistance','70S_ribosome':'Ribosome resistance','30S_ribosome':'Ribosome resistance','50S_ribosome':'Ribosome resistance','LpxA':'Lipid-A/envelope resistance','LpxC':'Lipid-A/envelope resistance','LpxH':'Lipid-A/envelope resistance','MurA':'MurA-pathway resistance','Beta-lactamase':'Beta-lactamase'}.get(r.target_class,'')
            if not resistance_family:
                resistance_family={'PBP1A':'PBP','PBP1B':'PBP','PBP2':'PBP','PBP2B':'PBP','PBP2X':'PBP','PBP3':'PBP','PBP4':'PBP','Beta-lactamase_class_A':'Beta-lactamase','Beta-lactamase_class_B':'Beta-lactamase','Beta-lactamase_class_C':'Beta-lactamase','Beta-lactamase_class_D':'Beta-lactamase'}.get(r.target_class,'')
            card_models=0; snp_rows=0; org_snp_rows=0
            if card_summary is not None and resistance_family:
                z=card_summary[card_summary.target_resistance_family==resistance_family]
                card_models=int(z.n_models.iloc[0]) if len(z) else 0
            if snp_summary is not None and resistance_family:
                z=snp_summary[snp_summary.resistance_family==resistance_family]
                snp_rows=int(z.n_snp_rows.iloc[0]) if len(z) else 0
            if snp_org is not None and resistance_family:
                z=snp_org[(snp_org.resistance_family==resistance_family)&(snp_org.organism.str.contains(org.split()[0],case=False,na=False))]
                org_snp_rows=int(z.n_snp_rows.sum()) if len(z) else 0
            context = V2_SCORING["context"]
            card_context=float(np.clip(
                float(context["card_model_weight"])*(1 if card_models else 0)
                +float(context["card_snp_weight"])*(1 if snp_rows else 0),0,1))
            if struct_summary is not None:
                z=struct_summary[struct_summary.target_class==r.target_class]
                if z.empty and parent != r.target_class: z=struct_summary[struct_summary.target_class==parent]
                struct_candidates=int(z.n_search_candidates.iloc[0]) if len(z) else 0; co_crystal=int(z.n_with_co_crystal_ligand.iloc[0]) if len(z) else 0
            else: struct_candidates=0; co_crystal=0
            pocket=float(
                context["pocket_with_cocrystal"] if co_crystal
                else context["pocket_with_structure"] if struct_candidates
                else context["pocket_missing"]
            )
            biology_weights = V2_SCORING["biology_weights"]
            biological=(
                float(biology_weights["organism_scope"])*scope
                +float(biology_weights["clinical"])*clinical
                +float(biology_weights["essentiality"])*essential
                +float(biology_weights["accessibility"])*access
                +float(biology_weights["resistance"])*resistance
                +float(biology_weights["card_context"])*card_context
            )
            resistance_burden=float(np.clip(
                float(context["card_model_weight"])*(1 if card_models else 0)
                +float(context["card_snp_weight"])*(1 if snp_rows else 0),0,1))
            translation_weights = V2_SCORING["clinical_translation_weights"]
            if r.target_class.startswith('Beta-lactamase_'):
                clinical_translation=float(np.clip(
                    float(translation_weights["clinical"])*clinical
                    +float(translation_weights["organism_scope"])*scope
                    +float(translation_weights["essentiality"])*essential
                    +float(translation_weights["accessibility"])*access
                    +float(translation_weights["pocket"])*pocket
                    +float(translation_weights["resistance_or_card"])*card_context,
                    0,1))
            else:
                clinical_translation=float(np.clip(
                    float(translation_weights["clinical"])*clinical
                    +float(translation_weights["organism_scope"])*scope
                    +float(translation_weights["essentiality"])*essential
                    +float(translation_weights["accessibility"])*access
                    +float(translation_weights["pocket"])*pocket
                    +float(translation_weights["resistance_or_card"])*(1-float(context["resistance_burden_discount"])*resistance_burden),
                    0,1))
            anti,anti_note,anti_status=anti_target_annotation(r.target_class)
            factors = V2_SCORING["overall_factors"]
            overall=float(
                chemical_quality
                *(float(factors["transfer"]["base"])+float(factors["transfer"]["weight"])*transfer)
                *(float(factors["pocket"]["base"])+float(factors["pocket"]["weight"])*pocket)
                *(float(factors["biology"]["base"])+float(factors["biology"]["weight"])*biological)
                *(1-float(factors["anti_target_penalty"])*anti)
            )
            reasons=[]
            if r.reference_quality_grade in {'low','insufficient'}: reasons.append('reference coverage limited')
            reason_thresholds = V2_SCORING["uncertainty_reasons"]
            if mapping_status!='mapped' or transfer<float(reason_thresholds["weak_transfer_below"]): reasons.append('species sequence mapping unresolved or weak')
            if r.target_specificity_score<float(reason_thresholds["weak_specificity_below"]): reasons.append('similarity is not target-specific versus cross-target decoys')
            if anti>=float(reason_thresholds["anti_target_at_least"]): reasons.append('human-homologue or mitochondrial selectivity risk is annotation-only')
            if pocket==0: reasons.append('no RCSB co-crystal/pocket evidence in bounded public catalog')
            confidence = V2_SCORING["confidence"]
            if (overall>=float(confidence["high"]["overall_min"])
                    and chemical_quality>=float(confidence["high"]["chemical_quality_min"])
                    and transfer>=float(confidence["high"]["transfer_min"])): conf='High'
            elif (overall>=float(confidence["moderate"]["overall_min"])
                    and chemical_quality>=float(confidence["moderate"]["chemical_quality_min"])): conf='Moderate'
            elif chemical_evidence>=float(confidence["low_chemical_min"]): conf='Low'
            else: conf='Insufficient'
            rr=r.to_dict(); rr.update({'organism':org,'parent_target_class':parent,'target_subtype':a.get('target_subtype',r.target_class),'binding_site_or_mechanism':a.get('binding_site_or_mechanism',a.get('mechanism_granularity','')),'organism_transfer_source':transfer_source,'species_transfer_score':transfer,'sequence_mapping_status':mapping_status,
                'organism_scope_score':scope,'clinical_priority_score':clinical,'essentiality_score':essential,'cellular_access_score':access,
                'resistance_relevance_score':resistance,'card_resistance_context_score':card_context,'resistance_burden_score':resistance_burden,'card_model_count':card_models,'card_snp_row_count':snp_rows,'organism_specific_snp_row_count':org_snp_rows,
                'rcsb_structure_candidate_count':struct_candidates,'rcsb_co_crystal_ligand_count':co_crystal,'pocket_evidence_score':pocket,'biological_priority_score':biological,'clinical_translation_score':clinical_translation,'chemical_hypothesis_score':chemical_quality,'anti_target_risk_score':anti,
                'anti_target_evidence_status':anti_status,'anti_target_note':anti_note,'overall_priority_score':overall,
                'confidence_class':conf,'uncertainty_reasons':'; '.join(reasons) if reasons else 'none',
                'recommended_validation':validation_plan(r.target_class,a.target_role),'clinical_status':a.clinical_status,
                'target_role':a.target_role,'organism_scope':a.organism_scope,'cellular_localization':a.cellular_localization,
                'resistance_relevance':a.resistance_relevance})
            if chemical_quality_column != "chemical_quality_adjusted_score":
                rr["chemical_hypothesis_score_source"] = chemical_quality_column
            rows.append(rr)
    return pd.DataFrame(rows)


def _complete_v3_missingness(frame, config=CONFIG):
    """Make unscored ontology rows explicit without fabricating component values."""

    completed = frame.copy()
    components = list(config.value("fusion.components"))
    if "chemical_evidence_score_v3" not in completed:
        completed["chemical_evidence_score_v3"] = np.nan
    for component in components:
        if component not in completed:
            completed[component] = np.nan
    missing = ";".join(components)
    unscored = completed["chemical_evidence_score_v3"].isna()
    completed.loc[unscored, components] = np.nan
    completed.loc[unscored, "fusion_component_count"] = 0
    completed.loc[unscored, "fusion_missing_components"] = missing
    completed.loc[unscored, "fusion_method"] = "unavailable_no_reference"
    completed.loc[unscored, "chemical_evidence_score_v3_is_probability"] = False
    completed.loc[
        unscored, "chemical_quality_adjusted_score_v3_is_probability"
    ] = False
    component_count = completed["fusion_component_count"].fillna(0).astype(int)
    completed["v3_evidence_status"] = np.select(
        [component_count == 0, component_count < len(components)],
        ["unavailable_no_reference", "partial_component_coverage"],
        default="complete_component_coverage",
    )
    return completed


def emit_v3_outputs(
    private,
    bench,
    refs,
    quality,
    ontology,
    compat,
    card_summary=None,
    snp_summary=None,
    snp_org=None,
    struct_summary=None,
    *,
    config=CONFIG,
    result_dir=RES,
    cache_dir=None,
):
    """Run the additive v3 path and return every output path written."""

    if str(config.value("run.fusion_mode")) != "rank_fusion":
        raise ValueError("v3 output currently requires run.fusion_mode=rank_fusion")
    result_dir.mkdir(parents=True, exist_ok=True)
    written = []
    disagreement_reports = []

    private_scores = []
    scoring_workers = int(config.value("chem3d.scoring_workers"))

    scoring_tasks = [
        *(("private", query) for query in private),
        *(("benchmark", query) for query in bench),
    ]
    scored_queries = ordered_process_map(
        _score_emit_v3_worker,
        scoring_tasks,
        workers=scoring_workers,
        initializer=_initialize_emit_v3_worker,
        initargs=(refs, quality, compat, ontology, config, cache_dir),
    )
    scored_private = scored_queries[: len(private)]
    for query, score in zip(private, scored_private, strict=True):
        score = add_unscored_classes(score, query, ontology)
        score = _complete_v3_missingness(score, config)
        score = assign_applicability_domain(score, config)
        if not score.empty:
            private_scores.append(score)
    if private_scores:
        raw = pd.concat(private_scores, ignore_index=True)
        path = result_dir / "open_target_scores_private_v3.csv"
        raw.to_csv(path, index=False)
        written.append(path)
        private_for_report = raw.copy()
        private_for_report["dataset_scope"] = "private_local"
        disagreement_reports.append(
            build_disagreement_report(
                private_for_report,
                component_columns=list(config.value("fusion.components")),
                minimum_absolute_rank_shift=float(
                    config.value("fusion.disagreement_min_absolute_rank_shift")
                ),
            )
        )
        if str(config.value("run.combiner")) == "heuristic":
            ranked = apply_biology(
                raw,
                ontology,
                compat,
                card_summary,
                snp_summary,
                snp_org,
                struct_summary,
                chemical_evidence_column="chemical_evidence_score_v3",
                chemical_quality_column="chemical_quality_adjusted_score_v3",
            )
            path = result_dir / "open_target_predictions_by_organism_v3.csv"
            ranked.to_csv(path, index=False)
            written.append(path)
            path = result_dir / "open_target_shortlist_by_organism_v3.csv"
            shortlist_with_applicability_domain(
                ranked,
                group_columns=["organism", "query_id"],
                score_column="overall_priority_score",
                top_n=10,
            ).to_csv(path, index=False)
            written.append(path)
            path = result_dir / "validation_priority_candidates_v3.csv"
            (
                ranked[ranked.confidence_class.isin(["High", "Moderate"])]
                .sort_values("overall_priority_score", ascending=False)
                .to_csv(path, index=False)
            )
            written.append(path)

    benchmark_scores = []
    scored_benchmark = scored_queries[len(private) :]
    for score in scored_benchmark:
        if not score.empty:
            benchmark_scores.append(score)
    if benchmark_scores:
        benchmark = pd.concat(benchmark_scores, ignore_index=True)
        benchmark = _complete_v3_missingness(benchmark, config)
        benchmark = assign_applicability_domain(benchmark, config)
        path = result_dir / "benchmark_open_target_scores_v3.csv"
        benchmark.to_csv(path, index=False)
        written.append(path)
        benchmark_for_report = benchmark.copy()
        benchmark_for_report["dataset_scope"] = "public_benchmark"
        disagreement_reports.append(
            build_disagreement_report(
                benchmark_for_report,
                component_columns=list(config.value("fusion.components")),
                minimum_absolute_rank_shift=float(
                    config.value("fusion.disagreement_min_absolute_rank_shift")
                ),
            )
        )

    if disagreement_reports:
        disagreements = pd.concat(disagreement_reports, ignore_index=True)
    else:
        empty_columns = [
            "query_id",
            "target_class",
            "dataset_scope",
            "chemical_evidence_score",
            "chemical_evidence_score_v3",
            *list(config.value("fusion.components")),
            "fusion_component_count",
            "fusion_missing_components",
        ]
        disagreements = build_disagreement_report(
            pd.DataFrame(columns=empty_columns),
            component_columns=list(config.value("fusion.components")),
            minimum_absolute_rank_shift=float(
                config.value("fusion.disagreement_min_absolute_rank_shift")
            ),
        )
    path = result_dir / "chemical_evidence_disagreements_v3.csv"
    disagreements.to_csv(path, index=False)
    written.append(path)

    coverage = pd.DataFrame(
        {
            "target_class": sorted(refs),
            "n_reference_ligands": [len(refs[key]) for key in sorted(refs)],
            "chem3d_method": (
                "ETKDGv3_USRCAT_O3A_BaseFeatures3D_Gobbi_Pharm2D"
            ),
        }
    )
    path = result_dir / "open_target_reference_coverage_v3.csv"
    coverage.to_csv(path, index=False)
    written.append(path)
    return written

def main():
    if bool(CONFIG.value("snapshots.verify_on_load")):
        verify_snapshot(CONFIG)
    refs=load_refs(); quality=pd.read_csv(QUALITY) if QUALITY.exists() else pd.DataFrame()
    ontology=load_ontology(); compat=pd.read_csv(COMPAT) if COMPAT.exists() else pd.DataFrame()
    card_summary=pd.read_csv(CARD_SUM) if CARD_SUM.exists() else pd.DataFrame()
    snp_summary=pd.read_csv(CARD_SNP) if CARD_SNP.exists() else pd.DataFrame()
    snp_org=pd.read_csv(CARD_SNP_ORG) if CARD_SNP_ORG.exists() else pd.DataFrame(columns=['organism','resistance_family','n_snp_rows'])
    struct_summary=pd.read_csv(STRUCT_SUM) if STRUCT_SUM.exists() else pd.DataFrame()
    private=load_queries(PRIVATE_COMPOUNDS)
    bench=load_queries(BENCH_STRUCTURES)
    if not bench and BENCH.exists():
        bdf=pd.read_csv(BENCH)
        for _,b in bdf.iterrows():
            m=mol(b.canonical_smiles)
            if m: bench.append({'query_id':b.drug,'query_name':b.drug,'mol':m,'fp':fp(m),'maccs':maccs(m),'source':'eskape_benchmark',**b.to_dict()})
    private_scores=[]
    for q in private:
        s=score_query(q,refs,quality,compat,ontology,exclude_close=False)
        s=add_unscored_classes(s,q,ontology)
        if not s.empty: private_scores.append(s)
    if private_scores:
        raw=pd.concat(private_scores,ignore_index=True); raw.to_csv(RES/'v2_open_target_scores_private.csv',index=False)
        ranked=apply_biology(raw,ontology,compat,card_summary,snp_summary,snp_org,struct_summary); ranked.to_csv(RES/'v2_open_target_predictions_by_organism.csv',index=False)
        ranked.sort_values(['organism','query_id','overall_priority_score'],ascending=[True,True,False]).groupby(['organism','query_id']).head(10).to_csv(RES/'v2_open_target_shortlist_by_organism.csv',index=False)
        ranked[ranked.chemical_hypothesis_score>0].sort_values(['organism','query_id','clinical_translation_score'],ascending=[True,True,False]).groupby(['organism','query_id']).head(10).to_csv(RES/'v21_clinical_translation_shortlist_by_organism.csv',index=False)
        ranked[ranked.confidence_class.isin(['High','Moderate'])].sort_values('overall_priority_score',ascending=False).to_csv(RES/'v2_validation_priority_candidates.csv',index=False)
        print('private queries',len(private),'raw rows',len(raw),'ranked rows',len(ranked))
    else: print('No protected private structures available; public-only run completed.')
    bench_scores=[]
    for q in bench:
        s=score_query(q,refs,quality,compat,ontology,exclude_close=True)
        if not s.empty: bench_scores.append(s)
    if bench_scores:
        bp=pd.concat(bench_scores,ignore_index=True); bp.to_csv(RES/'v2_benchmark_open_target_scores.csv',index=False)
        # Benchmarks are summarized in the dedicated v2 benchmark script.
        print('benchmark queries',len(bench),'rows',len(bp))
    pd.DataFrame({'target_class':sorted(refs),'n_reference_ligands':[len(refs[k]) for k in sorted(refs)]}).to_csv(RES/'v2_open_target_reference_coverage.csv',index=False)
    if str(CONFIG.value("run.fusion_mode")) == "rank_fusion":
        written = emit_v3_outputs(
            private,
            bench,
            refs,
            quality,
            ontology,
            compat,
            card_summary,
            snp_summary,
            snp_org,
            struct_summary,
        )
        print("v3 outputs", len(written))

if __name__=='__main__': main()
