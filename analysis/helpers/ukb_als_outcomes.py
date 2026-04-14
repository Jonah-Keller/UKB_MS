"""UKB ALS / Motor Neuron Disease outcome definitions.

Single source of truth for:
  - G12.2 ICD-10 case/control classification
  - Control exclusion codes
  - Timing classification (pre-onset / post-onset)
  - Death confirmation logic (Chia et al. 2025 required G12.2 at death)

Usage:
    from helpers.ukb_als_outcomes import (
        ALS_ICD, ALS_EXCLUDE_CONTROL_CODES, classify_als_timing
    )
"""
from __future__ import annotations

import math
from typing import FrozenSet

# ============================================================================
# ICD-10 Codes
# ============================================================================

# UKB stores ICD codes without dots: G12.2 → G122
ALS_ICD: FrozenSet[str] = frozenset({"G12.2", "G122"})

# UKB first-occurrence field for G12 (Motor neuron disease)
# Field 131036 = "Date G12 first reported (all sources)"
# TODO: Verify exact field for G12.2 subtype vs G12 broadly
ALS_FIRST_OCCURRENCE_FIELD: str = "p131036"

# Death cause fields: 40001 (primary), 40002 (contributory)
ALS_DEATH_FIELDS: tuple[str, ...] = ("p40001_i0", "p40002_i0",
                                      "p40001_i1", "p40002_i1")

# Control exclusion codes (Chia et al. 2025 Methods)
NEUROPATHY_CODES: FrozenSet[str] = frozenset({f"G{n}" for n in range(60, 65)})
MYOPATHY_CODES: FrozenSet[str] = frozenset({f"G{n}" for n in range(70, 74)})
MS_CODES: FrozenSet[str] = frozenset({"G35"})

ALS_EXCLUDE_CONTROL_CODES: FrozenSet[str] = (
    NEUROPATHY_CODES | MYOPATHY_CODES | MS_CODES | ALS_ICD
)

# ICD prefixes that confirm ALS in death records
ALS_DEATH_PREFIXES: tuple[str, ...] = ("G12", "G122")


# ============================================================================
# Classification Functions
# ============================================================================

def is_als_case(age_at_diagnosis: float | None) -> bool:
    """Return True if participant has a recorded ALS/MND diagnosis."""
    return age_at_diagnosis is not None and not math.isnan(float(age_at_diagnosis))


def is_valid_als_control(icd_codes: FrozenSet[str] | set) -> bool:
    """Return True if participant has no ALS/neuropathy/myopathy codes."""
    return len(set(icd_codes) & ALS_EXCLUDE_CONTROL_CODES) == 0


def is_als_confirmed_at_death(death_codes: list[str]) -> bool:
    """Return True if any death cause code matches ALS.

    Chia et al. 2025: UKB ALS cases required G12.2 confirmed at death.
    """
    for code in death_codes:
        if code is None:
            continue
        code_clean = str(code).upper().replace(".", "")
        if any(code_clean.startswith(p.replace(".", "")) for p in ALS_DEATH_PREFIXES):
            return True
    return False


def classify_als_timing(
    age_at_sampling: float,
    age_at_diagnosis: float | None,
) -> str:
    """Classify participant's ALS status relative to blood sampling.

    Returns:
        'pre_onset'  — blood drawn BEFORE ALS diagnosis (presymptomatic)
        'post_onset' — blood drawn AFTER ALS diagnosis
        'hc'         — no ALS diagnosis
    """
    if age_at_diagnosis is None or math.isnan(float(age_at_diagnosis)):
        return "hc"
    years_to_diagnosis = float(age_at_diagnosis) - float(age_at_sampling)
    return "pre_onset" if years_to_diagnosis > 0 else "post_onset"


def years_to_als_diagnosis(
    age_at_sampling: float,
    age_at_diagnosis: float | None,
) -> float | None:
    """Return years from blood draw to ALS diagnosis (positive = presymptomatic)."""
    if age_at_diagnosis is None or math.isnan(float(age_at_diagnosis)):
        return None
    return float(age_at_diagnosis) - float(age_at_sampling)


# ============================================================================
# Key ALS proteins (Chia et al. 2025, 33 discovery DEPs)
# ============================================================================

CHIA_33_DEPS: tuple[str, ...] = (
    "nefl", "chit1", "gpnmb", "lgals3", "s100a12",
    "chi3l1", "fcgr2b", "ctsh", "spp1", "havcr1",
    "trem2", "cxcl10", "mmp9", "mmp7", "pigr",
    "fcgr3a", "il6", "ccl18", "lcn2", "olfm4",
    "s100a8", "s100a9", "vim", "lgals1", "cd163",
    "msr1", "clec7a", "mrc1", "mki67", "top2a",
    "hmgb1", "thbs1", "fbn1",
)


if __name__ == "__main__":
    assert classify_als_timing(55.0, 60.0) == "pre_onset"
    assert classify_als_timing(65.0, 60.0) == "post_onset"
    assert classify_als_timing(55.0, None) == "hc"
    assert is_valid_als_control({"G60"}) is False
    assert is_valid_als_control({"I10"}) is True
    assert is_als_confirmed_at_death(["G122", "I25"]) is True
    assert is_als_confirmed_at_death(["I21", "E11"]) is False
    print("ukb_als_outcomes.py: all checks passed OK")
