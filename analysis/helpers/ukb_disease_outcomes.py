"""UKB disease outcome definitions — config-driven.

Single source of truth for case/control classification, control exclusion
codes, and timing classification. All disease-specific values come from
configs/disease.yaml via disease_config.load_disease_config().

To replicate on a different disease cohort, edit configs/disease.yaml.
Do not modify this file.

Usage:
    from analysis.helpers.ukb_disease_outcomes import (
        DISEASE_ICD, DISEASE_EXCLUDE_CONTROL_CODES,
        classify_disease_timing,
    )
"""
from __future__ import annotations

import math
from typing import FrozenSet

from analysis.helpers.disease_config import load_disease_config

_cfg = load_disease_config()

# ============================================================================
# ICD-10 Codes (loaded from configs/disease.yaml)
# ============================================================================

DISEASE_ICD: FrozenSet[str] = frozenset(_cfg.icd_codes)

DISEASE_FIRST_OCCURRENCE_FIELD: str = _cfg.first_occurrence_field

DEMYELINATING_CODES: FrozenSet[str] = frozenset(
    _cfg.control_exclusion_codes.get("demyelinating", [])
)
EPILEPSY_CODES: FrozenSet[str] = frozenset(
    _cfg.control_exclusion_codes.get("epilepsy", [])
)
MOVEMENT_CODES: FrozenSet[str] = frozenset(
    _cfg.control_exclusion_codes.get("movement", [])
)
DEMENTIA_CODES: FrozenSet[str] = frozenset(
    _cfg.control_exclusion_codes.get("dementia", [])
)

DISEASE_EXCLUDE_CONTROL_CODES: FrozenSet[str] = _cfg.all_exclusion_codes


# ============================================================================
# Classification Functions
# ============================================================================

def is_disease_case(age_at_diagnosis: float | None) -> bool:
    """Return True if participant has a recorded diagnosis."""
    return age_at_diagnosis is not None and not math.isnan(float(age_at_diagnosis))


def is_valid_disease_control(icd_codes: FrozenSet[str] | set) -> bool:
    """Return True if participant has no exclusion codes."""
    return len(set(icd_codes) & DISEASE_EXCLUDE_CONTROL_CODES) == 0


def classify_disease_timing(
    age_at_sampling: float,
    age_at_diagnosis: float | None,
) -> str:
    """Classify a participant's disease status relative to blood sampling.

    Returns one of cfg.status_values:
        'pre_onset'  — blood drawn BEFORE diagnosis (presymptomatic)
        'post_onset' — blood drawn AFTER diagnosis
        'hc'         — no diagnosis (healthy control)
    """
    sv = _cfg.status_values
    if age_at_diagnosis is None or math.isnan(float(age_at_diagnosis)):
        return sv["control"]
    years_to_diagnosis = float(age_at_diagnosis) - float(age_at_sampling)
    return sv["pre_onset"] if years_to_diagnosis > 0 else sv["post_onset"]


def years_to_disease_diagnosis(
    age_at_sampling: float,
    age_at_diagnosis: float | None,
) -> float | None:
    """Return years from blood draw to diagnosis.

    Positive = presymptomatic. Negative = blood after diagnosis. None = no Dx.
    """
    if age_at_diagnosis is None or math.isnan(float(age_at_diagnosis)):
        return None
    return float(age_at_diagnosis) - float(age_at_sampling)


# ============================================================================
# Key biomarkers (from configs/disease.yaml)
# ============================================================================

CNS_INJURY_MARKERS = tuple(_cfg.cns_injury_markers)
TOP_PRESYMPTOMATIC_DEPS = tuple(_cfg.top_presymptomatic_deps)


if __name__ == "__main__":
    sv = _cfg.status_values
    assert classify_disease_timing(40.0, 45.0) == sv["pre_onset"]
    assert classify_disease_timing(50.0, 45.0) == sv["post_onset"]
    assert classify_disease_timing(50.0, None) == sv["control"]
    assert is_valid_disease_control(set(_cfg.icd_codes)) is False
    assert is_valid_disease_control({"I10"}) is True
    assert years_to_disease_diagnosis(40.0, 45.0) == 5.0
    print(f"ukb_disease_outcomes.py: all checks passed (cohort={_cfg.cohort_short})")
