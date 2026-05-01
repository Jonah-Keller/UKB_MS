"""Build shared covariate table for UKB MS/ALS analyses.

Gathers: age at each instance, sex, BMI, Townsend deprivation, assessment
center, smoking status, 40 genetic PCs, UMAP1/2 of top-10 PCs (for ALS
limma covariate per Chia 2025), and HLA-DRB1*15:01 carrier status (IC
proxy for MS subgroup analysis per Abdelhak 2026).

Output: data/ukb/covariates/covariate_table.csv

Usage:
    python 01_build_covariate_table.py
    python 01_build_covariate_table.py --debug
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))

# --- CONFIG ---
SEED = 42
N_SUBSET = 500

DATA_DIR = REPO_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal" / "data" / "ukb"
OUTPUT_FILE = REPO_ROOT / "data" / "ukb" / "covariates" / "covariate_table.csv"
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# UKB field IDs — verify against data/ukb/field_id_dictionary.csv
FIELD_IDS = {
    "age_i0":            "p21003_i0",
    "age_i1":            "p21003_i1",
    "age_i2":            "p21003_i2",
    "age_i3":            "p21003_i3",
    "sex":               "p31",
    "bmi":               "p21001_i0",
    "townsend":          "p189",
    "assessment_center": "p54_i0",
    "smoking":           "p20116_i0",
    "age_recruitment":   "p21022",
}
PC_FIELD_PREFIX = "p22009_a"   # 40 PCs: p22009_a1 .. p22009_a40
HLA_DRB1_1501_FIELD = "p22182" # HLA imputation dosages


# ============================================================================
# File discovery
# ============================================================================

def _find_file(root: Path, *patterns: str) -> Path | None:
    for pattern in patterns:
        hits = sorted(root.rglob(pattern))
        if hits:
            return hits[0]
    return None


def _col(df: pd.DataFrame, field_id: str) -> str | None:
    """Find column in df matching field_id, handling 'participant.' prefix."""
    for prefix in ("", "participant."):
        name = f"{prefix}{field_id}"
        if name in df.columns:
            return name
    # Partial match
    matches = [c for c in df.columns if field_id in c]
    return matches[0] if matches else None


def _extract_field(df: pd.DataFrame, field_id: str) -> pd.Series:
    col = _col(df, field_id)
    if col is None:
        return pd.Series(np.nan, index=df.index, name=field_id)
    return df[col].astype(float)


# ============================================================================
# Loaders
# ============================================================================

def load_misc_file(label: str, *patterns: str) -> pd.DataFrame | None:
    path = _find_file(DATA_DIR / "misc", *patterns)
    if path is None:
        path = _find_file(DATA_DIR, *patterns)
    if path is None:
        print(f"  {label}: not found")
        return None
    print(f"  {label}: {path.name}")
    df = pd.read_csv(path, low_memory=False)
    if "participant.eid" in df.columns:
        df.rename(columns={"participant.eid": "eid"}, inplace=True)
    elif "participant_id" in df.columns:
        df.rename(columns={"participant_id": "eid"}, inplace=True)
    if "eid" in df.columns:
        df["eid"] = df["eid"].astype(int)
    return df


def load_demographics() -> pd.DataFrame | None:
    """Try to load master participant table."""
    return load_misc_file(
        "demographics",
        "data_participant.csv",
        "ukb_demographics*.csv",
        "participant*.csv",
    )


def load_bmi_age() -> pd.DataFrame | None:
    return load_misc_file("BMI/age", "ukb_bmi_age*.csv", "*bmi*.csv")


def load_age_of_instance() -> pd.DataFrame | None:
    return load_misc_file("age-at-instance", "ukb_age_of_instance*.csv",
                          "age_of_instance*.csv", "age_instance*.csv")


def load_genetic_pcs() -> pd.DataFrame | None:
    path = _find_file(DATA_DIR / "genetics",
                      "*22009*.csv", "*genetic_pc*.csv", "*pcs*.csv")
    if path is None:
        print("  genetic PCs: not found")
        return None
    print(f"  genetic PCs: {path.name}")
    df = pd.read_csv(path, low_memory=False)
    if "participant.eid" in df.columns:
        df.rename(columns={"participant.eid": "eid"}, inplace=True)
    elif "participant_id" in df.columns:
        df.rename(columns={"participant_id": "eid"}, inplace=True)
    if "eid" in df.columns:
        df["eid"] = df["eid"].astype(int)
    return df


def load_hla() -> pd.DataFrame | None:
    path = _find_file(DATA_DIR / "genetics",
                      "*22182*.csv", "*hla*.csv", "*HLA*.csv")
    if path is None:
        print("  HLA imputation: not found")
        return None
    print(f"  HLA imputation: {path.name}")
    df = pd.read_csv(path, low_memory=False)
    if "participant.eid" in df.columns:
        df.rename(columns={"participant.eid": "eid"}, inplace=True)
    elif "participant_id" in df.columns:
        df.rename(columns={"participant_id": "eid"}, inplace=True)
    if "eid" in df.columns:
        df["eid"] = df["eid"].astype(int)
    return df


# ============================================================================
# PC standardisation and UMAP
# ============================================================================

def extract_pcs(pc_df: pd.DataFrame) -> pd.DataFrame:
    """Detect and rename PC columns to PC1..PC40."""
    # Patterns: p22009_a1, PC1, PC_1, etc.
    pc_candidates = sorted([
        c for c in pc_df.columns
        if any(pat in c for pat in ["22009", "PC", "_a"])
        and c != "eid"
    ])
    # Map to PC1..PC40 by order
    rename = {}
    for i, col in enumerate(pc_candidates[:40], start=1):
        rename[col] = f"PC{i}"
    result = pc_df[["eid"] + list(rename.keys())].rename(columns=rename)
    print(f"  Found {len(rename)} genetic PCs")
    return result


def compute_umap(cov_df: pd.DataFrame, n_pcs: int = 10) -> pd.DataFrame:
    """Compute UMAP of top genetic PCs and add UMAP1/UMAP2 columns."""
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1) if f"PC{i}" in cov_df.columns]
    if len(pc_cols) < 2:
        print(f"  UMAP: only {len(pc_cols)} PCs available, skipping")
        cov_df["UMAP1"] = np.nan
        cov_df["UMAP2"] = np.nan
        return cov_df

    try:
        from umap import UMAP
    except ImportError:
        print("  UMAP: umap-learn not installed — pip install umap-learn")
        cov_df["UMAP1"] = np.nan
        cov_df["UMAP2"] = np.nan
        return cov_df

    pc_data = cov_df[pc_cols].dropna()
    if len(pc_data) < 10:
        cov_df["UMAP1"] = np.nan
        cov_df["UMAP2"] = np.nan
        return cov_df

    print(f"  Computing UMAP from {len(pc_cols)} PCs on {len(pc_data)} participants...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = UMAP(n_components=2, random_state=SEED,
                      n_neighbors=15, min_dist=0.1).fit_transform(pc_data.values)

    cov_df.loc[pc_data.index, "UMAP1"] = coords[:, 0]
    cov_df.loc[pc_data.index, "UMAP2"] = coords[:, 1]
    print("  UMAP done")
    return cov_df


# ============================================================================
# HLA-DRB1*15:01 extraction
# ============================================================================

def extract_hla_drb1_1501(hla_df: pd.DataFrame) -> pd.Series:
    """Extract HLA-DRB1*15:01 carrier status (dosage > 0.5).

    UKB field 22182 stores imputed HLA allele dosages. DRB1*15:01 is the
    primary MS risk allele (strongest GWAS hit, r2~0.9 with rs3135388).

    Column name varies — search for '15:01', '1501', or 'DRB1' in column names.
    """
    drb1_col = next(
        (c for c in hla_df.columns
         if any(x in c for x in ["15:01", "1501", "DRB1*15", "DRB1_1501"])),
        None,
    )
    if drb1_col is None:
        print("  HLA-DRB1*15:01 column not found in HLA file.")
        print("  TODO: Map UKB field 22182 allele codes to DRB1*15:01")
        return pd.Series(pd.NA, index=hla_df.index, dtype="boolean")

    dosage = pd.to_numeric(hla_df[drb1_col], errors="coerce")
    carrier = (dosage > 0.5).astype("boolean")
    n_carriers = carrier.sum()
    print(f"  HLA-DRB1*15:01 carriers: {n_carriers} "
          f"({100*n_carriers/len(carrier):.1f}%)")
    return carrier


# ============================================================================
# Main
# ============================================================================

def build_covariate_table(debug: bool = False) -> pd.DataFrame:
    print("Loading covariate files...")

    # --- Demographics / misc fields ---
    cov = None
    for loader in [load_demographics, load_bmi_age, load_age_of_instance]:
        df = loader()
        if df is None:
            continue
        # Extract fields we want
        fields_found = {}
        for name, fid in FIELD_IDS.items():
            s = _extract_field(df, fid)
            if s.notna().any():
                fields_found[name] = s
        if fields_found and cov is None:
            cov = df[["eid"]].copy()
        if fields_found and cov is not None:
            for name, s in fields_found.items():
                cov[name] = s.values[:len(cov)] if len(s) == len(cov) else s

    if cov is None:
        raise FileNotFoundError(
            "No covariate source files found. Extract covariates from UKB-RAP "
            "and place under data/ukb/covariates/ before running this stage."
        )

    if debug:
        cov = cov.sample(n=min(N_SUBSET, len(cov)), random_state=SEED).reset_index(drop=True)

    # --- Genetic PCs ---
    pc_df = load_genetic_pcs()
    if pc_df is not None:
        pc_clean = extract_pcs(pc_df)
        cov = cov.merge(pc_clean, on="eid", how="left")
    else:
        for i in range(1, 41):
            cov[f"PC{i}"] = np.nan

    # --- UMAP ---
    cov = compute_umap(cov)

    # --- HLA-DRB1*15:01 ---
    hla_df = load_hla()
    if hla_df is not None:
        carrier = extract_hla_drb1_1501(hla_df)
        cov["hla_drb1_1501_carrier"] = carrier.values[:len(cov)] if len(carrier) == len(cov) else np.nan
    else:
        cov["hla_drb1_1501_carrier"] = pd.NA
        print("  HLA-DRB1*15:01: not available (set to NA)")
        print("  TODO: Extract from UKB field 22182 when genetics data available")

    return cov


def print_summary(cov: pd.DataFrame) -> None:
    print("\n=== Covariate Table Summary ===")
    print(f"Participants: {len(cov)}")
    for col in ["age_i0", "sex", "bmi", "townsend"]:
        if col in cov.columns and cov[col].notna().any():
            print(f"  {col}: mean={cov[col].mean():.1f}, n_valid={cov[col].notna().sum()}")
    pc_avail = sum(1 for i in range(1, 41) if f"PC{i}" in cov.columns
                   and cov[f"PC{i}"].notna().any())
    print(f"  Genetic PCs available: {pc_avail}/40")
    umap_ok = cov["UMAP1"].notna().any() if "UMAP1" in cov.columns else False
    print(f"  UMAP computed: {'yes' if umap_ok else 'no'}")
    hla = cov.get("hla_drb1_1501_carrier", pd.Series([]))
    n_carriers = hla.sum() if hasattr(hla, "sum") else 0
    print(f"  HLA-DRB1*15:01 carriers: {n_carriers}")
    print(f"Output: {OUTPUT_FILE}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    print(f"Building covariate table (debug={args.debug})...")
    cov = build_covariate_table(debug=args.debug)
    cov.to_csv(OUTPUT_FILE, index=False)
    print_summary(cov)


if __name__ == "__main__":
    main()
