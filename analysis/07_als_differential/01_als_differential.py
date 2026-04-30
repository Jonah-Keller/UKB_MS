"""DEPRECATED — ALS-only legacy implementation. Use the cohort-agnostic R version:
   analysis/04_differential/01_limma_ms_vs_hc.R
   This file is retained for reproducing Chia et al. 2025 in the original
   paper context. For new disease replications DO NOT run this; the templated
   pipeline writes the same {cohort_short}_*_vs_hc.csv outputs config-style.

ALS Differential Protein Abundance Analysis — replicating Chia et al. 2025.

Three comparisons:
  1. All ALS (ALGO_MND refined) vs HC
  2. Pre-onset ALS vs HC
  3. Post-onset ALS vs HC

Covariates (matching Chia):
  - Age at sampling
  - Sex
  Note: tube type and UMAP not available in UKB Olink extract

Expected top hits from Chia:
  NEFL (log2FC=2.34), LIF, CSRP3, HSPB6, MEGF10, TPM3

Outputs:
  results/differential/als_all_vs_hc.csv
  results/differential/als_pre_vs_hc.csv
  results/differential/als_post_vs_hc.csv
  results/differential/als_volcano_plots.png

Usage:
    python 01_als_differential.py [--use-refined]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from helpers.limma_runner import run_limma_analysis, make_volcano_data, print_summary

# --- PATHS ---
ALS_OLINK_QC = REPO_ROOT / "data" / "ukb" / "olink" / "processed" / "als_olink_qc.csv"
ALS_REFINED = REPO_ROOT / "data" / "ukb" / "cohort" / "als_cohort_algo_refined.csv"
OLINK_PATH = Path(__file__).resolve().parents[2].parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal" / "data" / "ukb" / "olink" / "i0" / "olink_instance_0_extracted_data.csv"
OUT_DIR = REPO_ROOT / "results" / "differential"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Chia's 33 DEPs + ML features
ALS_KEY_PROTEINS = {
    "nefl", "lif", "csrp3", "hspb6", "megf10", "tpm3", "rnase3",
    "dusp13", "gfap", "igfbp2", "mmp9", "chi3l1", "gdf15",
    "sod1", "fus", "eif2s2", "hpcal1", "pdap1",
}


def load_data(use_refined: bool = True) -> tuple[pd.DataFrame, list[str]]:
    """Load QC'd ALS+Olink data.

    If use_refined=True, restrict cases to ALGO_MND-confirmed only.
    """
    print("Loading ALS Olink QC data...")
    data = pd.read_csv(ALS_OLINK_QC, low_memory=False)
    print(f"  Full dataset: {len(data)} participants")

    if use_refined and ALS_REFINED.exists():
        # Load refined cohort EIDs and filter
        refined = pd.read_csv(ALS_REFINED)
        refined_case_eids = set(refined[refined["als_status"] != "control"]["eid"])
        # Keep all controls + only ALGO_MND confirmed cases
        mask = (data["als_status"] == "control") | data["eid"].isin(refined_case_eids)
        data = data[mask].copy()
        n_cases = (data["als_status"] != "control").sum()
        print(f"  After ALGO_MND refinement: {len(data)} ({n_cases} cases)")

    cohort_cols = {"eid", "als_status", "age_at_sampling", "age_at_diagnosis",
                   "years_to_diagnosis", "sex", "als_confirmed_at_death", "olink_instance"}
    protein_cols = [c for c in data.columns if c not in cohort_cols]
    print(f"  {len(protein_cols)} proteins")

    if "sex" in data.columns:
        data["sex"] = pd.to_numeric(data["sex"], errors="coerce")

    return data, protein_cols


def plot_volcanoes(results_dict: dict[str, pd.DataFrame], output_path: Path) -> None:
    """Create a 3-panel volcano plot."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    titles = {
        "all": "All ALS vs HC",
        "pre": "Pre-onset ALS vs HC",
        "post": "Post-onset ALS vs HC",
    }

    for ax, (key, results) in zip(axes, results_dict.items()):
        if results.empty:
            ax.set_title(f"{titles.get(key, key)} — No results")
            continue

        vdata = make_volcano_data(results)

        sig5 = vdata["fdr"] < 0.05
        sig10 = vdata["fdr"] < 0.10
        not_sig = ~sig10

        ax.scatter(vdata.loc[not_sig, "log2fc"], vdata.loc[not_sig, "-log10(p)"],
                   c="#CCCCCC", s=8, alpha=0.5, zorder=1)
        ax.scatter(vdata.loc[sig10 & ~sig5, "log2fc"], vdata.loc[sig10 & ~sig5, "-log10(p)"],
                   c="#F4A460", s=12, alpha=0.7, zorder=2)

        up = sig5 & (vdata["log2fc"] > 0)
        down = sig5 & (vdata["log2fc"] < 0)
        ax.scatter(vdata.loc[up, "log2fc"], vdata.loc[up, "-log10(p)"],
                   c="#E74C3C", s=16, alpha=0.8, zorder=3)
        ax.scatter(vdata.loc[down, "log2fc"], vdata.loc[down, "-log10(p)"],
                   c="#3498DB", s=16, alpha=0.8, zorder=3)

        # Label key ALS proteins
        for _, row in vdata.iterrows():
            prot_name = row["protein"].split(".")[-1] if "." in row["protein"] else row["protein"]
            if prot_name.lower() in ALS_KEY_PROTEINS and row["fdr"] < 0.10:
                ax.annotate(prot_name.upper(), (row["log2fc"], row["-log10(p)"]),
                            fontsize=7, fontweight="bold", ha="center", va="bottom",
                            xytext=(0, 5), textcoords="offset points")

        ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.5, alpha=0.5)
        ax.axvline(0, color="gray", ls="-", lw=0.5, alpha=0.3)

        n_sig = sig5.sum()
        ax.set_title(f"{titles.get(key, key)}\n({n_sig} DEPs at FDR<5%)")
        ax.set_xlabel("log₂ fold change (ALS vs HC)")
        ax.set_ylabel("-log₁₀(p-value)")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved volcano plots: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--use-refined", action="store_true", default=True,
                        help="Use ALGO_MND-refined cohort (default: True)")
    parser.add_argument("--no-refined", dest="use_refined", action="store_false")
    args = parser.parse_args()

    data, protein_cols = load_data(use_refined=args.use_refined)
    covariates = ["age_at_sampling", "sex"]
    results_dict = {}

    # --- 1. All ALS vs HC ---
    print("\n[1/3] All ALS vs HC...")
    data["_als_binary"] = data["als_status"].apply(
        lambda x: "als" if x in ("pre_onset", "post_onset") else x
    )
    res_all = run_limma_analysis(
        data, protein_cols,
        group_col="_als_binary", case_value="als", control_value="control",
        covariates=covariates,
    )
    res_all.to_csv(OUT_DIR / "als_all_vs_hc.csv", index=False)
    print_summary(res_all, "All ALS vs HC")
    results_dict["all"] = res_all

    # --- 2. Pre-onset vs HC ---
    print("\n[2/3] Pre-onset ALS vs HC...")
    res_pre = run_limma_analysis(
        data, protein_cols,
        group_col="als_status", case_value="pre_onset", control_value="control",
        covariates=covariates,
    )
    res_pre.to_csv(OUT_DIR / "als_pre_vs_hc.csv", index=False)
    print_summary(res_pre, "Pre-onset ALS vs HC")
    results_dict["pre"] = res_pre

    # --- 3. Post-onset vs HC ---
    print("\n[3/3] Post-onset ALS vs HC...")
    res_post = run_limma_analysis(
        data, protein_cols,
        group_col="als_status", case_value="post_onset", control_value="control",
        covariates=covariates,
    )
    res_post.to_csv(OUT_DIR / "als_post_vs_hc.csv", index=False)
    print_summary(res_post, "Post-onset ALS vs HC")
    results_dict["post"] = res_post

    # --- Volcano plots ---
    plot_volcanoes(results_dict, OUT_DIR / "als_volcano_plots.png")

    # --- Chia replication check: how many of their 33 DEPs replicate? ---
    if not res_all.empty:
        chia_33 = {"nefl", "lif", "csrp3", "hspb6", "megf10", "tpm3", "rnase3",
                   "igfbp2", "mmp9", "chi3l1", "gdf15", "sod1", "dusp13",
                   "fhl1", "myom2", "tnnt3", "neb", "mb", "ampd1",
                   "actn2", "casq1", "ldb3", "pdlim3", "pdlim5",
                   "eif2s2", "hpcal1", "jpt2", "mtif3", "pdap1", "smad3",
                   "kynu", "cst3", "stom"}
        our_sig = set()
        for _, r in res_all.iterrows():
            pname = r["protein"].split(".")[-1] if "." in r["protein"] else r["protein"]
            if pname.lower() in chia_33 and r["fdr"] < 0.10:
                our_sig.add(pname.lower())
        print(f"\nChia replication: {len(our_sig)}/{len(chia_33)} of their DEPs replicate at FDR<10%")
        if our_sig:
            print(f"  Replicated: {', '.join(sorted(our_sig))}")

    print("\nDone! Results saved to:", OUT_DIR)


if __name__ == "__main__":
    main()
