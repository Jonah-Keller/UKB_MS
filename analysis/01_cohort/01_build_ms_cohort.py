"""Build MS cohort for UKB proteomics analysis.

Replicating Abdelhak et al. 2026 (Nature Medicine) cohort approach in UKB:
  - Cases:    G35 first-occurrence (field p131042, 2,619 UKB participants)
  - Controls: no demyelinating / epilepsy / movement / dementia codes
  - Timing:   classify each Olink participant as pre_onset, post_onset, or control
  - Output:   data/ukb/cohort/ms_cohort.csv

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

from helpers.ukb_ms_outcomes import classify_ms_timing, years_to_ms_diagnosis

# --- PATHS ---
FIRST_OCC_PATH = CADASIL_ROOT / "data" / "ukb" / "diagnoses" / "ukb_first_occurrence_FINAL_20251031_144424.csv"
OLINK_PATH = CADASIL_ROOT / "data" / "ukb" / "olink" / "i0" / "olink_instance_0_extracted_data.csv"
AGE_INST_PATH = CADASIL_ROOT / "data" / "ukb" / "misc" / "ukb_age_of_instance.csv"
OUTPUT_FILE = REPO_ROOT / "data" / "ukb" / "cohort" / "ms_cohort.csv"
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# --- FIELD MAPPING (verified Apr 2026) ---
MS_FIELD = "first_date_icd10_p131042"   # G35 MS (2,619 participants)

EXCLUSION_FIELDS = [
    "first_date_icd10_p131042",  # G35 MS
    "first_date_icd10_p131044",  # G36 Other acute disseminated demyelination
    "first_date_icd10_p131046",  # G37 Other demyelinating diseases of CNS
    "first_date_icd10_p131022",  # G20 Parkinson's disease
    "first_date_icd10_p131024",  # G21 Secondary parkinsonism
    "first_date_icd10_p131036",  # G30 Alzheimer's disease
    "first_date_icd10_p131038",  # G31 Other degenerative diseases of nervous system
    "first_date_icd10_p131048",  # G40 Epilepsy
]

# Only read these columns from the 1582-column CSV to keep load time fast
FIRST_OCC_COLS = ["participant_id", "sex", "year_birth", "age_recruitment",
                  MS_FIELD] + [f for f in EXCLUSION_FIELDS if f != MS_FIELD]


def build_ms_cohort() -> pd.DataFrame:
    """Build MS cohort by cross-referencing diagnoses with Olink participants."""

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

    # ---- 5. Age at MS diagnosis (year of first G35 record - year of birth) ----
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
    has_ms = df["age_at_diagnosis"].notna()
    has_sampling = df["age_at_sampling"].notna()
    is_excluded = exclusion_mask & ~has_ms

    keep = has_sampling & ~is_excluded
    df = df[keep].copy()
    has_ms = has_ms[keep]

    # Controls
    df["ms_status"] = "control"
    df["years_to_diagnosis"] = np.nan

    # MS cases
    ms_idx = df.index[has_ms[keep]]
    age_samp_ms = df.loc[ms_idx, "age_at_sampling"]
    age_dx_ms = df.loc[ms_idx, "age_at_diagnosis"]
    ytd = age_dx_ms - age_samp_ms  # positive = presymptomatic
    df.loc[ms_idx, "ms_status"] = ytd.apply(
        lambda y: "pre_onset" if y > 0 else "post_onset"
    )
    df.loc[ms_idx, "years_to_diagnosis"] = ytd.round(2)

    # Clear age_at_diagnosis for controls
    df.loc[~has_ms[keep], "age_at_diagnosis"] = np.nan

    out = df[["eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
              "years_to_diagnosis", "sex"]].copy()
    out["age_at_sampling"] = out["age_at_sampling"].round(2)
    out["age_at_diagnosis"] = out["age_at_diagnosis"].round(2)
    out["olink_instance"] = 0
    return out.reset_index(drop=True)


def print_summary(cohort: pd.DataFrame) -> None:
    print("\n=== MS Cohort Summary ===")
    counts = cohort["ms_status"].value_counts()
    n_pre = counts.get("pre_onset", 0)
    n_post = counts.get("post_onset", 0)
    n_hc = counts.get("control", 0)
    print(f"Total participants:  {len(cohort):,}")
    print(f"MS cases:           {n_pre + n_post:,}")
    if n_pre > 0:
        med_ytd = cohort.loc[cohort["ms_status"] == "pre_onset",
                             "years_to_diagnosis"].median()
        print(f"  Pre-onset:        {n_pre:,}  (median {med_ytd:.1f} yrs before Dx)")
    print(f"  Post-onset:       {n_post:,}")
    print(f"Controls:           {n_hc:,}")
    if "sex" in cohort.columns:
        print(f"Sex distribution:   {dict(cohort['sex'].value_counts())}")
    print(f"Output:             {OUTPUT_FILE}")


def main() -> None:
    print("Building MS cohort from real UKB data...")
    cohort = build_ms_cohort()
    cohort.to_csv(OUTPUT_FILE, index=False)
    print_summary(cohort)


if __name__ == "__main__":
    main()
