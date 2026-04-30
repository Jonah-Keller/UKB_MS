"""DEPRECATED — MS-only legacy implementation. Use the cohort-agnostic R version:
   analysis/04_differential/01_limma_ms_vs_hc.R
   This file is retained for reproducing Abdelhak et al. 2026 in the original
   paper context. For new disease replications (stroke, DM2, PD, ...) DO NOT
   run this; the templated pipeline writes the same {cohort_short}_*_vs_hc.csv
   outputs config-style.

MS Differential Protein Abundance Analysis — replicating Abdelhak et al. 2026.

Three comparisons, each using limma-style empirical Bayes:
  1. All MS vs HC (combined pre+post onset)
  2. Pre-onset MS vs HC
  3. Post-onset MS vs HC

Covariates (matching Abdelhak):
  - Age at sampling
  - Sex
  - Combined NPX (average across all proteins per sample)

Outputs:
  results/differential/ms_all_vs_hc.csv
  results/differential/ms_pre_vs_hc.csv
  results/differential/ms_post_vs_hc.csv
  results/differential/ms_volcano_plots.png

Usage:
    python 01_ms_differential.py
"""
from __future__ import annotations

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
MS_OLINK_QC = REPO_ROOT / "data" / "ukb" / "olink" / "processed" / "ms_olink_qc.csv"
OUT_DIR = REPO_ROOT / "results" / "differential"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Key proteins to highlight on volcano plots (from Abdelhak et al.)
MS_KEY_PROTEINS = {
    "nefl", "mog", "gfap", "il3", "aif1", "sirt6", "dock10",
    "mtf1", "sanbr", "kcnd2", "tbx22", "nrgn",
}


def load_data() -> tuple[pd.DataFrame, list[str]]:
    """Load QC'd MS+Olink data and identify protein columns."""
    print("Loading MS Olink QC data...")
    data = pd.read_csv(MS_OLINK_QC, low_memory=False)
    print(f"  {len(data)} participants")

    # Identify protein columns (everything not in cohort metadata)
    cohort_cols = {"eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
                   "years_to_diagnosis", "sex", "olink_instance"}
    protein_cols = [c for c in data.columns if c not in cohort_cols]
    print(f"  {len(protein_cols)} proteins")

    # Convert sex to numeric if needed
    if "sex" in data.columns:
        data["sex"] = pd.to_numeric(data["sex"], errors="coerce")

    return data, protein_cols


def plot_volcanoes(
    results_dict: dict[str, pd.DataFrame],
    output_path: Path,
) -> None:
    """Create a 3-panel volcano plot."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    titles = {
        "all": "All MS vs HC",
        "pre": "Pre-onset MS vs HC",
        "post": "Post-onset MS vs HC",
    }

    for ax, (key, results) in zip(axes, results_dict.items()):
        if results.empty:
            ax.set_title(f"{titles.get(key, key)} — No results")
            continue

        vdata = make_volcano_data(results)

        # Color by significance
        sig5 = vdata["fdr"] < 0.05
        sig10 = vdata["fdr"] < 0.10
        not_sig = ~sig10

        ax.scatter(vdata.loc[not_sig, "log2fc"], vdata.loc[not_sig, "-log10(p)"],
                   c="#CCCCCC", s=8, alpha=0.5, zorder=1)
        ax.scatter(vdata.loc[sig10 & ~sig5, "log2fc"], vdata.loc[sig10 & ~sig5, "-log10(p)"],
                   c="#F4A460", s=12, alpha=0.7, zorder=2, label=f"FDR<10% (n={sig10.sum()})")
        
        # FDR<5%: color by direction
        up = sig5 & (vdata["log2fc"] > 0)
        down = sig5 & (vdata["log2fc"] < 0)
        ax.scatter(vdata.loc[up, "log2fc"], vdata.loc[up, "-log10(p)"],
                   c="#E74C3C", s=16, alpha=0.8, zorder=3)
        ax.scatter(vdata.loc[down, "log2fc"], vdata.loc[down, "-log10(p)"],
                   c="#3498DB", s=16, alpha=0.8, zorder=3)

        # Label key MS proteins
        for _, row in vdata.iterrows():
            prot_name = row["protein"].split(".")[-1] if "." in row["protein"] else row["protein"]
            if prot_name.lower() in MS_KEY_PROTEINS and row["fdr"] < 0.10:
                ax.annotate(prot_name.upper(), (row["log2fc"], row["-log10(p)"]),
                            fontsize=7, fontweight="bold", ha="center", va="bottom",
                            xytext=(0, 5), textcoords="offset points")

        # Thresholds
        ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.5, alpha=0.5)
        ax.axvline(0, color="gray", ls="-", lw=0.5, alpha=0.3)

        n_sig = sig5.sum()
        ax.set_title(f"{titles.get(key, key)}\n({n_sig} DEPs at FDR<5%)")
        ax.set_xlabel("log₂ fold change (MS vs HC)")
        ax.set_ylabel("-log₁₀(p-value)")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved volcano plots: {output_path}")


def main() -> None:
    data, protein_cols = load_data()

    covariates = ["age_at_sampling", "sex"]
    results_dict = {}

    # --- 1. All MS vs HC ---
    print("\n[1/3] All MS vs HC...")
    # Create a combined 'case' label for all MS
    data["_ms_binary"] = data["ms_status"].apply(
        lambda x: "ms" if x in ("pre_onset", "post_onset") else x
    )
    res_all = run_limma_analysis(
        data, protein_cols,
        group_col="_ms_binary", case_value="ms", control_value="control",
        covariates=covariates, add_combined_npx=True,
    )
    res_all.to_csv(OUT_DIR / "ms_all_vs_hc.csv", index=False)
    print_summary(res_all, "All MS vs HC")
    results_dict["all"] = res_all

    # --- 2. Pre-onset vs HC ---
    print("\n[2/3] Pre-onset MS vs HC...")
    res_pre = run_limma_analysis(
        data, protein_cols,
        group_col="ms_status", case_value="pre_onset", control_value="control",
        covariates=covariates, add_combined_npx=True,
    )
    res_pre.to_csv(OUT_DIR / "ms_pre_vs_hc.csv", index=False)
    print_summary(res_pre, "Pre-onset MS vs HC")
    results_dict["pre"] = res_pre

    # --- 3. Post-onset vs HC ---
    print("\n[3/3] Post-onset MS vs HC...")
    res_post = run_limma_analysis(
        data, protein_cols,
        group_col="ms_status", case_value="post_onset", control_value="control",
        covariates=covariates, add_combined_npx=True,
    )
    res_post.to_csv(OUT_DIR / "ms_post_vs_hc.csv", index=False)
    print_summary(res_post, "Post-onset MS vs HC")
    results_dict["post"] = res_post

    # --- Volcano plots ---
    plot_volcanoes(results_dict, OUT_DIR / "ms_volcano_plots.png")

    # --- Cross-comparison: overlap between pre and all ---
    if not res_all.empty and not res_pre.empty:
        sig_all = set(res_all.loc[res_all["fdr"] < 0.10, "protein"])
        sig_pre = set(res_pre.loc[res_pre["fdr"] < 0.10, "protein"])
        overlap = sig_all & sig_pre
        print(f"\nDEPs (FDR<10%): all={len(sig_all)}, pre={len(sig_pre)}, overlap={len(overlap)}")
        if sig_all:
            print(f"Overlap fraction: {len(overlap)/len(sig_all)*100:.1f}% of 'all' DEPs shared with pre-onset")

    print("\nDone! Results saved to:", OUT_DIR)


if __name__ == "__main__":
    main()
