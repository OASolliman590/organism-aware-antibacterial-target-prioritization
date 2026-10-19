from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline import open_target_discovery_v2 as v2
from pipeline.applicability_domain import shortlist_with_applicability_domain

TARGET = "GyrB"
MAPPED_ORG, UNMAPPED_ORG, NONFINITE_ORG = v2.ORGANISMS[0], v2.ORGANISMS[1], v2.ORGANISMS[2]


def _fixtures():
    raw = pd.DataFrame(
        [
            {
                "query_id": "q1",
                "query_name": "q1",
                "target_class": TARGET,
                "chemical_evidence_score": 0.8,
                "chemical_quality_adjusted_score": 0.9,
                "reference_quality_grade": "high",
                "target_specificity_score": 0.9,
            }
        ]
    )
    ontology = pd.DataFrame(
        [
            {
                "target_class": TARGET,
                "parent_target_class": np.nan,
                "organism_scope": "broad bacteria",
                "clinical_status": "clinical",
                "essentiality_level": "essential",
                "cellular_localization": "cytoplasm",
                "resistance_relevance": "high",
                "target_role": "direct viability target",
                "mechanism_granularity": "enzyme",
            }
        ]
    )
    compat = pd.DataFrame(
        [
            {
                "organism": MAPPED_ORG,
                "target_class": TARGET,
                "species_transfer_score": 0.0,
                "sequence_status": "mapped",
            },
            {
                "organism": NONFINITE_ORG,
                "target_class": TARGET,
                "species_transfer_score": np.nan,
                "sequence_status": "unmapped",
            },
        ]
    )
    return raw, ontology, compat


def test_missing_mapping_is_nan_not_zero_and_ordered_after_mapped() -> None:
    raw, ontology, compat = _fixtures()
    ranked = v2.apply_biology(raw, ontology, compat)
    by_org = ranked.set_index("organism")
    mapped = by_org.loc[MAPPED_ORG]
    unmapped = by_org.loc[UNMAPPED_ORG]

    # 1. omitted mapping: NaN transfer, flagged as absent
    assert pd.isna(unmapped.species_transfer_score)
    assert unmapped.species_mapping_present == False  # noqa: E712
    # 2. real mapped zero: stays 0.0 and is flagged as present
    assert mapped.species_transfer_score == 0.0
    assert mapped.species_mapping_present == True  # noqa: E712
    # 3. both rows have a finite overall score
    assert np.isfinite(mapped.overall_priority_score)
    assert np.isfinite(unmapped.overall_priority_score)
    # 4. missing mapping cannot be High
    assert unmapped.confidence_class != "High"
    # 5. the species-mapping gap is named in the uncertainty reasons
    assert "species sequence mapping unresolved or weak" in unmapped.uncertainty_reasons
    # Non-finite mapping record behaves like an omitted mapping
    nonfinite = by_org.loc[NONFINITE_ORG]
    assert pd.isna(nonfinite.species_transfer_score)
    assert nonfinite.species_mapping_present == False  # noqa: E712
    assert nonfinite.confidence_class != "High"
    # Default config: missing_factor equals the mapped real-zero factor, so the score is unchanged.
    assert np.isclose(
        unmapped.overall_priority_score, mapped.overall_priority_score, rtol=0, atol=1e-9
    )

    # 6. in a two-row group, the mapped row precedes the unmapped row
    group = ranked[ranked.organism.isin([MAPPED_ORG, UNMAPPED_ORG])].copy()
    group["ad_shortlist_eligible"] = True
    # Make the unmapped row score higher, so only the mapping priority can put mapped first.
    group.loc[group.organism == UNMAPPED_ORG, "overall_priority_score"] = (
        float(mapped.overall_priority_score) + 1.0
    )
    shortlist = shortlist_with_applicability_domain(
        group,
        group_columns=["query_id"],
        score_column="overall_priority_score",
        top_n=2,
        additional_priority_columns=["species_mapping_present"],
    )
    assert list(shortlist.organism) == [MAPPED_ORG, UNMAPPED_ORG]
