"""Build ALS cohort for UKB proteomics analysis.

Replicating Chia et al. 2025 (Nature Medicine) UKB cohort approach:
  - Cases:    G12.2 first-occurrence (field p131016 = G12 Motor neuron disease, 918 UKB)
  - Controls: no neuropathy (G60-G64), myopathy (G70-73), or MS (G35)
  - Timing:   pre-onset (blood before Dx) vs. post-onset
  - Death:    Check if ALS confirmed at death via ado_participant p42028
  - Output:   data/ukb/cohort/als_cohort.csv

Note from Chia paper: "13 ALS, 23,601 controls; ICD-10 G12.2 confirmed at
death for cases" — UKB ALS n is very small (~15-50 with Olink data).

Data sources:
  - First-occurrence dates: CADASIL repo diagnoses/ukb_first_occurrence_FINAL_*.csv
  - Olink instance 0: CADASIL repo olink/i0/olink_instance_0_extracted_data.csv
  - Age at instance: CADASIL repo misc/ukb_age_of_instance.csv
  - Death dates (MND): ado_participant.csv field p42028

Field ID reference (UKB Data Dictionary, verified Apr 2026):
  G12=p131016  G35=p131042  G60=p131082  G70=p131092

Usage:
    python 02_build_als_cohort.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
CADASIL_ROOT = REPO_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from helpers.ukb_als_outcomes import classify_als_timing, years_to_als_diagnosis

# --- PATHS ---
FIRST_OCC_PATH = CADASIL_ROOT / "data" / "ukb" / "diagnoses" / "ukb_first_occurrence_FINAL_20251031_144424.csv"
OLINK_PATH = CADASIL_ROOT / "data" / "ukb" / "olink" / "i0" / "olink_instance_0_extracted_data.csv"
AGE_INST_PATH = CADASIL_ROOT / "data" / "ukb" / "misc" / "ukb_age_of_instance.csv"
ADO_PATH = CADASIL_ROOT / "data" / "ukb" / "misc" / "ado_participant.csv"
OUTPUT_FILE = REPO_ROOT / "data" / "ukb" / "cohort" / "als_cohort.csv"
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# --- FIELD MAPPING (verified Apr 2026) ---
ALS_FIELD = "first_date_icd10_p131016"   # G12 Motor neuron disease (918 participants)
MND_DEATH_FIELD = "p42028"               # ado_participant: Date of MND report

EXCLUSION_FIELDS = [
    "first_date_icd10_p131016",  # G12 ALS/MND (cases handled via has_als flag)
    "first_date_icd10_p131042",  # G35 MS
    "first_date_icd10_p131082",  # G60 Hereditary/idiopathic neuropathy
    "first_date_icd10_p131092",  # G70 Myasthenia gravis / NMJ disorders
]

# Only read these columns from the 1582-column CSV
FIRST_OCC_COLS = ["participant_id", "sex", "year_birth", "age_recruitment",
                  ALS_FIELD] + [f for f in EXCLUSION_FIELDS if f != ALS_FIELD]


def build_als_cohort() -> pd.DataFrame:
    """Build ALS cohort by cross-referencing diagnoses with Olink participants."""

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
        print("  WARNING: age-at-instance-0 column not found")

    # ---- 4. Load MND death confirmation from ado_participant ----
    mnd_death_eids: set[int] = set()
    if ADO_PATH.exists():
        print(f"Loading MND death dates: {ADO_PATH.name}")
        ado = pd.read_csv(ADO_PATH, low_memory=False)
        if "participant_id" in ado.columns:
            ado = ado.rename(columns={"participant_id": "eid"})
        if "eid" in ado.columns:
            ado["eid"] = ado["eid"].astype(int)
        if MND_DEATH_FIELD in ado.columns:
            mnd_mask = pd.to_datetime(ado[MND_DEATH_FIELD], errors="coerce").notna()
            mnd_death_eids = set(ado.loc[mnd_mask, "eid"])
            print(f"  {len(mnd_death_eids):,} participants with MND death confirmation")
        else:
            print(f"  WARNING: {MND_DEATH_FIELD} not in ado_participant")
    else:
        print(f"  WARNING: ado_participant not found at {ADO_PATH}")

    # ---- 5. Restrict to Olink participants ----
    df = first_occ[first_occ["eid"].isin(olink_eid_set)].copy()
    print(f"  {len(df):,} participants with both Olink + diagnosis data")

    # ---- 6. Age at G12 diagnosis ----
    if ALS_FIELD in df.columns:
        als_year = pd.to_datetime(df[ALS_FIELD], errors="coerce").dt.year
        df["age_at_diagnosis"] = (als_year - df["year_birth"]).where(als_year.notna())
    else:
        print(f"  WARNING: {ALS_FIELD} not found!")
        df["age_at_diagnosis"] = np.nan

    # ---- 7. Age at Olink blood draw ----
    if age_map is not None:
        df["age_at_sampling"] = df["eid"].map(age_map).astype(float)
    elif "age_recruitment" in df.columns:
        df["age_at_sampling"] = df["age_recruitment"].astype(float)
    else:
        df["age_at_sampling"] = np.nan

    # ---- 8. Exclusion mask for controls ----
    exclusion_mask = pd.Series(False, index=df.index)
    for field in EXCLUSION_FIELDS:
        if field in df.columns and field != ALS_FIELD:
            exclusion_mask |= pd.to_datetime(df[field], errors="coerce").notna()

    # ---- 9. Classify (vectorized) ----
    has_als = df["age_at_diagnosis"].notna()
    has_sampling = df["age_at_sampling"].notna()
    is_excluded = exclusion_mask & ~has_als

    keep = has_sampling & ~is_excluded
    df = df[keep].copy()
    has_als = has_als[keep]

    df["als_status"] = "control"
    df["years_to_diagnosis"] = np.nan
    df["als_confirmed_at_death"] = np.nan

    als_idx = df.index[has_als[keep]]
    age_samp_als = df.loc[als_idx, "age_at_sampling"]
    age_dx_als = df.loc[als_idx, "age_at_diagnosis"]
    ytd = age_dx_als - age_samp_als
    df.loc[als_idx, "als_status"] = ytd.apply(
        lambda y: "pre_onset" if y > 0 else "post_onset"
    )
    df.loc[als_idx, "years_to_diagnosis"] = ytd.round(2)
    df.loc[als_idx, "als_confirmed_at_death"] = df.loc[als_idx, "eid"].isin(
        mnd_death_eids
    ).astype(float)

    df.loc[~has_als[keep], "age_at_diagnosis"] = np.nan

    out = df[["eid", "als_status", "age_at_sampling", "age_at_diagnosis",
              "years_to_diagnosis", "sex", "als_confirmed_at_death"]].copy()
    out["age_at_sampling"] = out["age_at_sampling"].round(2)
    out["age_at_diagnosis"] = out["age_at_diagnosis"].round(2)
    out["olink_instance"] = 0
    return out.reset_index(drop=True)


def print_summary(cohort: pd.DataFrame) -> None:
    print("\n=== ALS Cohort Summary ===")
    counts = cohort["als_status"].value_counts()
    n_pre = counts.get("pre_onset", 0)
    n_post = counts.get("post_onset", 0)
    n_hc = counts.get("control", 0)
    n_conf = int(cohort["als_confirmed_at_death"].sum()) if "als_confirmed_at_death" in cohort.columns else 0
    print(f"Total participants:      {len(cohort):,}")
    print(f"ALS cases:               {n_pre + n_post:,}")
    print(f"  Pre-onset:             {n_pre:,}")
    print(f"  Post-onset:            {n_post:,}")
    print(f"  Confirmed at death:    {n_conf:,}")
    print(f"Controls:                {n_hc:,}")
    if "sex" in cohort.columns:
        print(f"Sex distribution:        {dict(cohort['sex'].value_counts())}")
    print(f"Output:                  {OUTPUT_FILE}")


def main() -> None:
    print("Building ALS cohort from real UKB data...")
    cohort = build_als_cohort()
    cohort.to_csv(OUTPUT_FILE, index=False)
    print_summary(cohort)


if __name__ == "__main__":
    main()
