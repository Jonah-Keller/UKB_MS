"""Olink QC and Cohort Merge for MS and ALS analyses.

Steps:
  1. Load Olink instance 0 NPX data from CADASIL repo
  2. Merge with MS and ALS cohort tables
  3. Perform IQR-based global outlier detection (3*IQR from median NPX)
  4. Output cleaned, merged datasets for downstream limma/ML

Outputs:
  data/ukb/olink/processed/ms_olink_qc.csv
  data/ukb/olink/processed/als_olink_qc.csv
  results/qc/ms_npx_distribution.png
  results/qc/als_npx_distribution.png

Usage:
    python 01_olink_qc.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[2]
CADASIL_ROOT = REPO_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"

# --- PATHS ---
OLINK_PATH = CADASIL_ROOT / "data" / "ukb" / "olink" / "i0" / "olink_instance_0_extracted_data.csv"
MS_COHORT_PATH = REPO_ROOT / "data" / "ukb" / "cohort" / "ms_cohort.csv"
ALS_COHORT_PATH = REPO_ROOT / "data" / "ukb" / "cohort" / "als_cohort.csv"

OUT_DIR = REPO_ROOT / "data" / "ukb" / "olink" / "processed"
QC_DIR = REPO_ROOT / "results" / "qc"
OUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)


def load_olink() -> pd.DataFrame:
    """Load Olink NPX data, strip column prefix."""
    print(f"Loading Olink NPX: {OLINK_PATH.name}")
    df = pd.read_csv(OLINK_PATH, low_memory=False)
    # Strip 'olink_instance_0.' prefix
    df.columns = [c.replace("olink_instance_0.", "") for c in df.columns]
    df["eid"] = df["eid"].astype(int)
    print(f"  {len(df)} participants, {len(df.columns) - 1} proteins")
    return df


def run_qc(
    disease: str,
    cohort: pd.DataFrame,
    olink: pd.DataFrame,
    status_col: str,
) -> pd.DataFrame | None:
    """Merge cohort with Olink, run QC, save output."""
    print(f"\n=== {disease} Olink QC ===")

    # Merge
    merged = cohort.merge(olink, on="eid", how="inner")
    print(f"  Merged: {len(merged)} participants")
    if len(merged) == 0:
        print(f"  WARNING: No overlap between {disease} cohort and Olink data!")
        return None

    # Identify protein columns (everything not in cohort)
    cohort_cols = set(cohort.columns)
    protein_cols = [c for c in merged.columns if c not in cohort_cols and c != "eid"]
    print(f"  Protein columns: {len(protein_cols)}")

    # Global median NPX per participant
    protein_matrix = merged[protein_cols].astype(float)
    global_median = protein_matrix.median(axis=1)

    # IQR outlier detection (3*IQR)
    Q1 = global_median.quantile(0.25)
    Q3 = global_median.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 3 * IQR
    upper = Q3 + 3 * IQR
    is_outlier = (global_median < lower) | (global_median > upper)
    n_outlier = is_outlier.sum()
    print(f"  IQR outliers: {n_outlier} (bounds: [{lower:.3f}, {upper:.3f}])")

    # Missingness per-protein summary
    pct_missing = protein_matrix.isna().mean() * 100
    n_high_miss = (pct_missing > 50).sum()
    print(f"  Proteins with >50% missing: {n_high_miss}")

    # Plot distribution by status
    merged["global_median_npx"] = global_median
    merged["qc_outlier"] = is_outlier

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Violin + box by status
    ax = axes[0]
    status_order = sorted(merged[status_col].unique())
    colors = {"control": "#4A90D9", "pre_onset": "#E67E22", "post_onset": "#E74C3C", "hc": "#4A90D9"}
    palette = {s: colors.get(s, "#888") for s in status_order}
    sns.violinplot(data=merged, x=status_col, y="global_median_npx",
                   order=status_order, palette=palette, alpha=0.6, ax=ax)
    sns.boxplot(data=merged, x=status_col, y="global_median_npx",
                order=status_order, width=0.15, fliersize=0,
                boxprops=dict(alpha=0.8), ax=ax)
    ax.set_title(f"{disease} — Global Median NPX by Status")
    ax.set_ylabel("Median NPX")
    ax.set_xlabel("")

    # Histogram showing outlier cutoffs
    ax = axes[1]
    ax.hist(global_median, bins=80, color="#4A90D9", alpha=0.7, edgecolor="white")
    ax.axvline(lower, color="red", ls="--", label=f"Lower: {lower:.2f}")
    ax.axvline(upper, color="red", ls="--", label=f"Upper: {upper:.2f}")
    ax.set_title(f"{disease} — QC Outlier Detection ({n_outlier} flagged)")
    ax.set_xlabel("Global Median NPX")
    ax.legend()

    plt.tight_layout()
    fig_path = QC_DIR / f"{disease.lower()}_npx_distribution.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved QC plot: {fig_path}")

    # Filter and save
    clean = merged[~is_outlier].drop(columns=["global_median_npx", "qc_outlier"])
    out_path = OUT_DIR / f"{disease.lower()}_olink_qc.csv"
    clean.to_csv(out_path, index=False)
    print(f"  Saved QC'd data: {out_path} ({len(clean)} participants)")

    # Status breakdown after QC
    status_counts = clean[status_col].value_counts()
    for status, count in status_counts.items():
        print(f"    {status}: {count}")

    return clean


def main() -> None:
    olink = load_olink()

    # MS
    ms_cohort = pd.read_csv(MS_COHORT_PATH)
    ms_cohort["eid"] = ms_cohort["eid"].astype(int)
    run_qc("MS", ms_cohort, olink, "ms_status")

    # ALS
    als_cohort = pd.read_csv(ALS_COHORT_PATH)
    als_cohort["eid"] = als_cohort["eid"].astype(int)
    run_qc("ALS", als_cohort, olink, "als_status")


if __name__ == "__main__":
    main()
