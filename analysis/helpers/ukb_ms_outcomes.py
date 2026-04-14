"""UKB Multiple Sclerosis outcome definitions.

Single source of truth for:
  - G35 ICD-10 case/control classification
  - Control exclusion codes (neurological)
  - Timing classification (pre-onset / post-onset)

Mirror of pattern from CADASIL repo: analysis/helpers/ukb_outcomes.py

Usage:
    from helpers.ukb_ms_outcomes import (
        MS_ICD, MS_EXCLUDE_CONTROL_CODES, classify_ms_timing
    )
"""
from __future__ import annotations

import math
from typing import FrozenSet

# ============================================================================
# ICD-10 Codes
# ============================================================================

MS_ICD: FrozenSet[str] = frozenset({"G35"})

# UKB first-occurrence field for G35.
# Field 131042 = "Date G35 first reported (all sources)"
# In the first-occurrence CSV: column is "first_date_icd10_p131042"
# NOTE: p131040 is G32 (other degenerative), NOT G35.
MS_FIRST_OCCURRENCE_FIELD: str = "p131042"

# Control exclusion — participants with any of these codes are excluded
# from the control group to ensure neurologically clean controls.
DEMYELINATING_CODES: FrozenSet[str] = frozenset({"G35", "G36", "G37"})
EPILEPSY_CODES: FrozenSet[str] = frozenset({f"G{n}" for n in range(40, 42)})
MOVEMENT_CODES: FrozenSet[str] = frozenset({f"G{n}" for n in range(20, 27)})
DEMENTIA_CODES: FrozenSet[str] = frozenset({f"F{n:02d}" for n in range(0, 10)})

MS_EXCLUDE_CONTROL_CODES: FrozenSet[str] = (
    DEMYELINATING_CODES | EPILEPSY_CODES | MOVEMENT_CODES | DEMENTIA_CODES
)


# ============================================================================
# Classification Functions
# ============================================================================

def is_ms_case(age_at_diagnosis: float | None) -> bool:
    """Return True if participant has a recorded MS diagnosis."""
    return age_at_diagnosis is not None and not math.isnan(float(age_at_diagnosis))


def is_valid_ms_control(icd_codes: FrozenSet[str] | set) -> bool:
    """Return True if participant has no neurological exclusion codes."""
    return len(set(icd_codes) & MS_EXCLUDE_CONTROL_CODES) == 0


def classify_ms_timing(
    age_at_sampling: float,
    age_at_diagnosis: float | None,
) -> str:
    """Classify a participant's MS status relative to blood sampling.

    Returns:
        'pre_onset'  — blood drawn BEFORE MS diagnosis (presymptomatic)
        'post_onset' — blood drawn AFTER MS diagnosis
        'hc'         — no MS diagnosis (healthy control)
    """
    if age_at_diagnosis is None or math.isnan(float(age_at_diagnosis)):
        return "hc"
    years_to_diagnosis = float(age_at_diagnosis) - float(age_at_sampling)
    return "pre_onset" if years_to_diagnosis > 0 else "post_onset"


def years_to_ms_diagnosis(
    age_at_sampling: float,
    age_at_diagnosis: float | None,
) -> float | None:
    """Return years from blood draw to MS diagnosis.

    Positive = presymptomatic (blood before diagnosis).
    Negative = blood after diagnosis.
    None = no diagnosis.
    """
    if age_at_diagnosis is None or math.isnan(float(age_at_diagnosis)):
        return None
    return float(age_at_diagnosis) - float(age_at_sampling)


# ============================================================================
# Key MS proteins (Abdelhak et al. 2026)
# ============================================================================

CNS_INJURY_MARKERS = ("nefl", "mog", "gfap")   # Fig. 3 proteins
TOP_PRESYMPTOMATIC_DEPS = (                       # Top pre-onset DEPs from paper
    "spag8", "cdin1", "lctl", "dock10", "vwa1",
    "insm2", "kcnd2", "igsf11", "kif26b", "ccdc40",
    "nefl", "mog", "gfap", "nrgn", "il3", "aif1",
    "sirt6",
)


if __name__ == "__main__":
    assert classify_ms_timing(40.0, 45.0) == "pre_onset"
    assert classify_ms_timing(50.0, 45.0) == "post_onset"
    assert classify_ms_timing(50.0, None) == "hc"
    assert is_valid_ms_control({"G35"}) is False
    assert is_valid_ms_control({"I10"}) is True
    assert years_to_ms_diagnosis(40.0, 45.0) == 5.0
    print("ukb_ms_outcomes.py: all checks passed OK")
