"""Build disease cohort for UKB proteomics analysis (config-driven).

All disease-specific constants (ICD codes, first-occurrence field, status
column name, output filename prefix) come from configs/disease.yaml. To
replicate this study on a different disease cohort, edit that file.

Replicating Abdelhak et al. 2026 (Nature Medicine) cohort approach in UKB
for the configured disease (default: MS / G35):
  - Cases:    cfg.icd_codes first-occurrence (field cfg.first_occurrence_field)
  - Controls: no codes from cfg.control_exclusion_codes
  - Timing:   classify each Olink participant as pre_onset, post_onset, control
  - Output:   data/ukb/cohort/{cfg.cohort_short}_cohort.csv

Data sources:
  - First-occurrence dates: CADASIL repo diagnoses/ukb_first_occurrence_FINAL_*.csv
  - Olink instance 0: CADASIL repo olink/i0/olink_instance_0_extracted_data.csv
  - Age at instance: CADASIL repo misc/ukb_age_of_instance.csv

Field ID reference (UKB Data Dictionary, verified Apr 2026):
  G12=p131016  G20=p131022  G21=p131024  G30=p131036  G31=p131038
  G32=p131040  G35=p131042  G36=p131044  G37=p131046  G40=p131048

Usage:
    python 01_build_ms_cohort.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CADASIL_ROOT = REPO_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from helpers.disease_config import load_disease_config
from helpers.ukb_ms_outcomes import classify_ms_timing, years_to_ms_diagnosis

cfg = load_disease_config()
STATUS_COL = cfg.cohort_status_col
STATUS_CONTROL = cfg.status_values["control"]
STATUS_PRE = cfg.status_values["pre_onset"]
STATUS_POST = cfg.status_values["post_onset"]

# --- PATHS ---
FIRST_OCC_PATH = CADASIL_ROOT / "data" / "ukb" / "diagnoses" / "ukb_first_occurrence_FINAL_20251031_144424.csv"
OLINK_PATH = CADASIL_ROOT / "data" / "ukb" / "olink" / "i0" / "olink_instance_0_extracted_data.csv"
AGE_INST_PATH = CADASIL_ROOT / "data" / "ukb" / "misc" / "ukb_age_of_instance.csv"
OUTPUT_FILE = REPO_ROOT / "data" / "ukb" / "cohort" / f"{cfg.cohort_short}_cohort.csv"
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# --- FIELD MAPPING (loaded from configs/disease.yaml) ---
# e.g. cfg.first_occurrence_field == "p131042" → "first_date_icd10_p131042" (G35 MS)
MS_FIELD = f"first_date_icd10_{cfg.first_occurrence_field}"

# Control-exclusion fields: the case field plus every other disqualifying
# first-occurrence column derived from configs/disease.yaml. Static mapping
# of UKB ICD-10 root → first-occurrence field id (extend here as new diseases
# are configured in disease.yaml).
_ICD_TO_FIELD = {
    "G12": "p131016",  # Motor neuron disease
    "G20": "p131022",  # Parkinson's disease
    "G21": "p131024",  # Secondary parkinsonism
    "G30": "p131036",  # Alzheimer's disease
    "G31": "p131038",  # Other degenerative diseases of nervous system
    "G32": "p131040",  # Other degenerative disorders of nervous system
    "G35": "p131042",  # Multiple sclerosis
    "G36": "p131044",  # Other acute disseminated demyelination
    "G37": "p131046",  # Other demyelinating diseases of CNS
    "G40": "p131048",  # Epilepsy
}

EXCLUSION_FIELDS = sorted({
    f"first_date_icd10_{_ICD_TO_FIELD[code]}"
    for code in cfg.all_exclusion_codes
    if code in _ICD_TO_FIELD
} | {MS_FIELD})

# Only read these columns from the 1582-column CSV to keep load time fast
FIRST_OCC_COLS = ["participant_id", "sex", "year_birth", "age_recruitment",
                  MS_FIELD] + [f for f in EXCLUSION_FIELDS if f != MS_FIELD]


def build_ms_cohort() -> pd.DataFrame:
    """Build disease cohort by cross-referencing diagnoses with Olink participants."""

    # ---- 1. Load first-occurrence dates (only needed columns) ----
    print(f"Loading first-occurrence dates: {FIRST_OCC_PATH.name}")
    avail_cols = pd.read_csv(FIRST_OCC_PATH, nrows=0).columns.tolist()
    use_cols = [c for c in FIRST_OCC_COLS if c in avail_cols]
    first_occ = pd.read_csv(FIRST_OCC_PATH, usecols=use_cols, low_memory=False)
    first_occ = first_occ.rename(columns={"participant_id": "eid"})
    first_occ["eid"] = first_occ["eid"].astype(int)
    print(f"  {len(first_occ):,} participants loaded ({len(use_cols)} columns)")

    # ---- 2. Load Olink instance 0 EIDs ----
    print(f"Loading Olink EIDs: {OLINK_PATH.name}")
    olink_eid_col = "olink_instance_0.eid"
    olink_eids = pd.read_csv(OLINK_PATH, usecols=[olink_eid_col], low_memory=False)
    olink_eid_set = set(olink_eids[olink_eid_col].astype(int))
    print(f"  {len(olink_eid_set):,} participants with Olink data")

    # ---- 3. Load age at instance 0 ----
    print(f"Loading age at instance: {AGE_INST_PATH.name}")
    age_inst = pd.read_csv(AGE_INST_PATH, low_memory=False)
    if "participant.eid" in age_inst.columns:
        age_inst = age_inst.rename(columns={"participant.eid": "eid"})
    age_inst["eid"] = age_inst["eid"].astype(int)
    age_i0_col = next(
        (c for c in age_inst.columns if "21003" in c and ("i0" in c or "_0" in c)), None
    )
    if age_i0_col:
        age_map = age_inst.set_index("eid")[age_i0_col]
        print(f"  Using age column: {age_i0_col}")
    else:
        age_map = None
        print("  WARNING: age-at-instance-0 column not found; falling back to age_recruitment")

    # ---- 4. Restrict to Olink participants ----
    df = first_occ[first_occ["eid"].isin(olink_eid_set)].copy()
    print(f"  {len(df):,} participants with both Olink + diagnosis data")

    # ---- 5. Age at diagnosis (year of first case-ICD record - year of birth) ----
    if MS_FIELD in df.columns:
        ms_year = pd.to_datetime(df[MS_FIELD], errors="coerce").dt.year
        df["age_at_diagnosis"] = (ms_year - df["year_birth"]).where(ms_year.notna())
    else:
        print(f"  WARNING: {MS_FIELD} not found!")
        df["age_at_diagnosis"] = np.nan

    # ---- 6. Age at Olink blood draw ----
    if age_map is not None:
        df["age_at_sampling"] = df["eid"].map(age_map).astype(float)
    elif "age_recruitment" in df.columns:
        df["age_at_sampling"] = df["age_recruitment"].astype(float)
    else:
        df["age_at_sampling"] = np.nan

    # ---- 7. Exclusion mask for controls ----
    exclusion_mask = pd.Series(False, index=df.index)
    for field in EXCLUSION_FIELDS:
        if field in df.columns and field != MS_FIELD:
            exclusion_mask |= pd.to_datetime(df[field], errors="coerce").notna()

    # ---- 8. Classify (vectorized) ----
    has_case = df["age_at_diagnosis"].notna()
    has_sampling = df["age_at_sampling"].notna()
    is_excluded = exclusion_mask & ~has_case

    keep = has_sampling & ~is_excluded
    df = df[keep].copy()
    has_case = has_case[keep]

    # Controls
    df[STATUS_COL] = STATUS_CONTROL
    df["years_to_diagnosis"] = np.nan

    # Cases
    case_idx = df.index[has_case[keep]]
    age_samp_case = df.loc[case_idx, "age_at_sampling"]
    age_dx_case = df.loc[case_idx, "age_at_diagnosis"]
    ytd = age_dx_case - age_samp_case  # positive = presymptomatic
    df.loc[case_idx, STATUS_COL] = ytd.apply(
        lambda y: STATUS_PRE if y > 0 else STATUS_POST
    )
    df.loc[case_idx, "years_to_diagnosis"] = ytd.round(2)

    # Clear age_at_diagnosis for controls
    df.loc[~has_case[keep], "age_at_diagnosis"] = np.nan

    out = df[["eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
              "years_to_diagnosis", "sex"]].copy()
    out["age_at_sampling"] = out["age_at_sampling"].round(2)
    out["age_at_diagnosis"] = out["age_at_diagnosis"].round(2)
    out["olink_instance"] = 0
    return out.reset_index(drop=True)


def print_summary(cohort: pd.DataFrame) -> None:
    print(f"\n=== {cfg.disease_short_caps} Cohort Summary ===")
    counts = cohort[STATUS_COL].value_counts()
    n_pre = counts.get(STATUS_PRE, 0)
    n_post = counts.get(STATUS_POST, 0)
    n_hc = counts.get(STATUS_CONTROL, 0)
    print(f"Total participants:  {len(cohort):,}")
    print(f"{cfg.disease_short_caps} cases:           {n_pre + n_post:,}")
    if n_pre > 0:
        med_ytd = cohort.loc[cohort[STATUS_COL] == STATUS_PRE,
                             "years_to_diagnosis"].median()
        print(f"  Pre-onset:        {n_pre:,}  (median {med_ytd:.1f} yrs before Dx)")
    print(f"  Post-onset:       {n_post:,}")
    print(f"Controls:           {n_hc:,}")
    if "sex" in cohort.columns:
        print(f"Sex distribution:   {dict(cohort['sex'].value_counts())}")
    print(f"Output:             {OUTPUT_FILE}")


def main() -> None:
    print(f"Building {cfg.disease_short_caps} cohort from real UKB data...")
    cohort = build_ms_cohort()
    cohort.to_csv(OUTPUT_FILE, index=False)
    print_summary(cohort)


if __name__ == "__main__":
    main()
