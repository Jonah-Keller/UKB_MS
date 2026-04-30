"""Download RAP extraction outputs and integrate into local data.

Run AFTER the RAP notebook (03_rap_extraction_notebook.py) has completed.

Usage:
    python 04_download_and_integrate.py [--download-only] [--integrate-only]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
RAP_OUTPUT_DIR = REPO_ROOT / "data" / "ukb" / "rap_extraction"
RAP_PROJECT = "project-GxX43xjJpjp2G7XK6ffz0qb1"
RAP_REMOTE_DIR = "/data/ms_als_extraction/"

# Output directories
GENETICS_DIR = REPO_ROOT / "data" / "ukb" / "genetics"
COVARIATES_DIR = REPO_ROOT / "data" / "ukb" / "covariates"
HES_DIR = REPO_ROOT / "data" / "ukb" / "hes"
COHORT_DIR = REPO_ROOT / "data" / "ukb" / "cohort"

for d in [RAP_OUTPUT_DIR, GENETICS_DIR, COVARIATES_DIR, HES_DIR, COHORT_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def download_from_rap():
    """Download all extraction CSVs from RAP."""
    print("=" * 60)
    print("Downloading from RAP...")
    print("=" * 60)

    # List files on RAP
    result = subprocess.run(
        ["dx", "ls", f"{RAP_REMOTE_DIR}"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"Error listing RAP files: {result.stderr}")
        return False

    files = [f.strip() for f in result.stdout.strip().split("\n") if f.strip().endswith(".csv")]
    print(f"Found {len(files)} CSV files on RAP")

    for f in files:
        remote_path = f"{RAP_REMOTE_DIR}{f}"
        local_path = RAP_OUTPUT_DIR / f
        print(f"  Downloading {f}...")
        dl_result = subprocess.run(
            ["dx", "download", f"{RAP_PROJECT}:{remote_path}", "-o", str(local_path), "-f"],
            capture_output=True, text=True,
        )
        if dl_result.returncode != 0:
            print(f"    ERROR: {dl_result.stderr}")
        else:
            size_mb = local_path.stat().st_size / 1e6
            print(f"    OK ({size_mb:.1f} MB)")

    # Download MRI data from parent project extracted directory
    print("  Downloading Brain MRI data...")
    mri_dest = REPO_ROOT / "data" / "ukb" / "brain_mri" / "i2"
    mri_dest.mkdir(parents=True, exist_ok=True)
    mri_path = mri_dest / "brain_mri_i2.csv"
    
    dl_mri = subprocess.run(
        ["dx", "download", f"{RAP_PROJECT}:/data/brain_mri/i2/brain_mri_i2.csv", "-o", str(mri_path), "-f"],
        capture_output=True, text=True,
    )
    if dl_mri.returncode != 0:
        print(f"    MRI ERROR: {dl_mri.stderr}")
    else:
        print(f"    MRI OK ({(mri_path.stat().st_size / 1e6):.1f} MB)")

    return True


def integrate_genetic_pcs():
    """Process genetic PCs into a clean table."""
    pc_path = RAP_OUTPUT_DIR / "genetic_pcs.csv"
    if not pc_path.exists():
        print("  genetic_pcs.csv not found, skipping")
        return

    print("\nProcessing Genetic PCs...")
    df = pd.read_csv(pc_path)
    
    # Rename columns: p22009_a1 → PC1, etc.
    rename = {"eid": "eid"}
    for col in df.columns:
        if col.startswith("p22009_a"):
            pc_num = col.split("_a")[-1]
            rename[col] = f"PC{pc_num}"
    
    df = df.rename(columns=rename)
    
    # Keep first 10 PCs (standard for GWAS covariates)
    pc_cols = ["eid"] + [f"PC{i}" for i in range(1, 11) if f"PC{i}" in df.columns]
    df_out = df[pc_cols].dropna(subset=["PC1"])
    
    out_path = GENETICS_DIR / "genetic_pcs_top10.csv"
    df_out.to_csv(out_path, index=False)
    print(f"  Saved: {out_path} ({len(df_out)} participants)")


def integrate_hla():
    """Extract HLA-DRB1*15:01 dosage."""
    hla_path = RAP_OUTPUT_DIR / "hla_imputation.csv"
    if not hla_path.exists():
        print("  hla_imputation.csv not found, skipping")
        return

    print("\nProcessing HLA imputation...")
    df = pd.read_csv(hla_path, low_memory=False)
    
    # Find the DRB1*15:01 column
    # HLA field 22182 has 362 array indices, each for a different allele
    # We need to identify which array index = DRB1*15:01
    # From UKB coding 2966: DRB1*1501 often at specific indices
    # For now, save the full HLA table and we'll map later
    
    out_path = GENETICS_DIR / "hla_imputation.csv"
    df.to_csv(out_path, index=False)
    print(f"  Saved full HLA: {out_path} ({len(df)} participants, {len(df.columns)} columns)")
    
    # Try to identify DRB1*15:01 by checking field descriptions
    print("  Note: DRB1*15:01 allele index needs to be identified from UKB coding 2966")


def integrate_hes_icd10():
    """Process HES ICD-10 to extract G12.2 (ALS) cases specifically."""
    hes_path = RAP_OUTPUT_DIR / "hes_icd10_diagnoses.csv"
    if not hes_path.exists():
        print("  hes_icd10_diagnoses.csv not found, skipping")
        return

    print("\nProcessing HES ICD-10 diagnoses...")
    df = pd.read_csv(hes_path, low_memory=False)
    
    # Find all G12.2 cases in any diagnosis column
    diag_cols = [c for c in df.columns if c.startswith("p41270") or c.startswith("p41202") or c.startswith("p41204")]
    
    g12_2_eids = set()
    g35_eids = set()
    g12_all_eids = set()
    
    for col in diag_cols:
        vals = df[col].dropna().astype(str)
        # G12.2 = Motor neuron disease (ALS specifically)
        g12_2_eids.update(df.loc[vals.str.startswith("G122"), "eid"].values)
        g12_2_eids.update(df.loc[vals == "G12.2", "eid"].values)
        # G35 = Multiple sclerosis
        g35_eids.update(df.loc[vals.str.startswith("G35"), "eid"].values)
        # All G12
        g12_all_eids.update(df.loc[vals.str.startswith("G12"), "eid"].values)
    
    print(f"  G12 (any): {len(g12_all_eids)} participants")
    print(f"  G12.2 (ALS specific): {len(g12_2_eids)} participants")
    print(f"  G35 (MS): {len(g35_eids)} participants")
    
    # Save G12.2 EID list
    g12_2_df = pd.DataFrame({"eid": sorted(g12_2_eids), "has_G12_2": True})
    g12_2_df.to_csv(HES_DIR / "g12_2_als_eids.csv", index=False)
    
    g35_df = pd.DataFrame({"eid": sorted(g35_eids), "has_G35_hes": True})
    g35_df.to_csv(HES_DIR / "g35_ms_eids.csv", index=False)
    
    # Save full HES for reference
    df.to_csv(HES_DIR / "hes_icd10_full.csv", index=False)
    print(f"  Saved G12.2 EIDs: {HES_DIR / 'g12_2_als_eids.csv'}")


def integrate_death_causes():
    """Process death cause ICD-10 codes."""
    death_path = RAP_OUTPUT_DIR / "death_causes.csv"
    if not death_path.exists():
        print("  death_causes.csv not found, skipping")
        return

    print("\nProcessing death causes...")
    df = pd.read_csv(death_path, low_memory=False)
    
    # Find cases where G12.2 is in any death cause field
    death_cols = [c for c in df.columns if c.startswith("p40001") or c.startswith("p40002")]
    
    g12_2_death = set()
    for col in death_cols:
        vals = df[col].dropna().astype(str)
        g12_2_death.update(df.loc[vals.str.contains("G12.2|G122", regex=True, na=False), "eid"].values)
    
    print(f"  G12.2 as cause of death: {len(g12_2_death)} participants")
    
    death_df = pd.DataFrame({"eid": sorted(g12_2_death), "g12_2_death": True})
    death_df.to_csv(HES_DIR / "g12_2_death_confirmed.csv", index=False)
    
    # Save full death data
    df.to_csv(HES_DIR / "death_causes_full.csv", index=False)


def integrate_genetic_qc():
    """Process genetic QC fields."""
    qc_path = RAP_OUTPUT_DIR / "genetic_qc.csv"
    if not qc_path.exists():
        print("  genetic_qc.csv not found, skipping")
        return

    print("\nProcessing genetic QC...")
    df = pd.read_csv(qc_path, low_memory=False)
    
    # Rename for clarity
    rename_map = {}
    for col in df.columns:
        if "22006" in col: rename_map[col] = "genetic_ethnic_grouping"
        elif "22019" in col: rename_map[col] = "sex_chrom_aneuploidy"
        elif "22021" in col: rename_map[col] = "genetic_kinship"
        elif "22027" in col: rename_map[col] = "het_missing_outlier"
        elif "22001" in col: rename_map[col] = "genetic_sex"
        elif "22000" in col: rename_map[col] = "genotype_batch"
        elif "22020" in col: rename_map[col] = "used_in_pca"
    
    df = df.rename(columns=rename_map)
    
    # Apply standard genetic QC filters
    n_total = len(df)
    passers = df.copy()
    
    if "sex_chrom_aneuploidy" in passers.columns:
        passers = passers[passers["sex_chrom_aneuploidy"].isna()]
    if "het_missing_outlier" in passers.columns:
        passers = passers[passers["het_missing_outlier"].isna()]
    if "used_in_pca" in passers.columns:
        passers = passers[passers["used_in_pca"] == 1]
    
    print(f"  Total: {n_total}")
    print(f"  Pass genetic QC: {len(passers)}")
    
    out_path = GENETICS_DIR / "genetic_qc.csv"
    df.to_csv(out_path, index=False)
    
    # Save pass list
    passers[["eid"]].to_csv(GENETICS_DIR / "genetic_qc_pass_eids.csv", index=False)
    print(f"  Saved: {out_path}")


def integrate_covariates():
    """Process covariates."""
    cov_path = RAP_OUTPUT_DIR / "covariates_demographics.csv"
    if not cov_path.exists():
        print("  covariates_demographics.csv not found, skipping")
        return

    print("\nProcessing covariates...")
    df = pd.read_csv(cov_path, low_memory=False)
    
    out_path = COVARIATES_DIR / "rap_covariates.csv"
    df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path} ({len(df)} participants)")


def build_refined_als_cohort():
    """Build the most refined ALS cohort using G12.2 + death confirmation."""
    als_cohort_path = COHORT_DIR / "als_cohort.csv"
    g12_2_path = HES_DIR / "g12_2_als_eids.csv"
    death_path = HES_DIR / "g12_2_death_confirmed.csv"
    algo_path = COHORT_DIR / "als_cohort_algo_refined.csv"
    
    if not als_cohort_path.exists():
        print("  als_cohort.csv not found, skipping")
        return
    
    print("\nBuilding refined ALS cohort...")
    als = pd.read_csv(als_cohort_path)
    als_cases = als[als["als_status"] != "control"]
    als_controls = als[als["als_status"] == "control"]
    
    print(f"  Original ALS cases: {len(als_cases)}")
    
    # Apply G12.2 HES filter if available
    if g12_2_path.exists():
        g12_2 = pd.read_csv(g12_2_path)
        g12_2_eids = set(g12_2["eid"])
        als_g12_2 = als_cases[als_cases["eid"].isin(g12_2_eids)]
        print(f"  With G12.2 in HES: {len(als_g12_2)}")
    else:
        als_g12_2 = als_cases
        print("  G12.2 HES data not available, using all G12")
    
    # Apply death confirmation if available
    if death_path.exists():
        death = pd.read_csv(death_path)
        death_eids = set(death["eid"])
        als_death = als_g12_2[als_g12_2["eid"].isin(death_eids)]
        print(f"  G12.2 + death confirmed: {len(als_death)}")
    
    # Build Chia-style cohort: G12.2 + death confirmed
    if g12_2_path.exists() and death_path.exists():
        chia_cases = als_cases[als_cases["eid"].isin(g12_2_eids & death_eids)]
        chia_cohort = pd.concat([chia_cases, als_controls])
        chia_path = COHORT_DIR / "als_cohort_chia_refined.csv"
        chia_cohort.to_csv(chia_path, index=False)
        print(f"  Chia-style refined: {len(chia_cases)} cases, {len(als_controls)} controls")
        print(f"  Saved: {chia_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--integrate-only", action="store_true")
    args = parser.parse_args()

    if not args.integrate_only:
        success = download_from_rap()
        if not success:
            print("Download failed. Run --integrate-only after manually downloading.")
    
    if not args.download_only:
        print("\n" + "=" * 60)
        print("INTEGRATING EXTRACTED DATA")
        print("=" * 60)
        
        integrate_genetic_pcs()
        integrate_hla()
        integrate_hes_icd10()
        integrate_death_causes()
        integrate_genetic_qc()
        integrate_covariates()
        build_refined_als_cohort()
        
        print("\n" + "=" * 60)
        print("INTEGRATION COMPLETE")
        print("=" * 60)
        print("\nFiles created:")
        for d in [GENETICS_DIR, COVARIATES_DIR, HES_DIR]:
            for f in sorted(d.glob("*.csv")):
                size = f.stat().st_size / 1e6
                print(f"  {f.relative_to(REPO_ROOT)} ({size:.1f} MB)")


if __name__ == "__main__":
    main()
