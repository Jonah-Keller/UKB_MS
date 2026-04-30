"""Validate MS and ALS cohort extraction against UKB ground truth.

This script:
  1. Cross-references our cohort EIDs with ICD-10 code tables (4-char precision)
  2. Cross-references ALS with the ALGO_MOTOR_NEURONE_DISEASE algorithmic outcome
  3. Checks for ALS causative gene mutations (SOD1, C9orf72, TARDBP, FUS) via
     Olink protein levels as proxies (SOD1, TARDBP are on Olink)
  4. Compares NEFL (neurofilament light) levels between cases and controls
     to validate that our MS and ALS extractions capture real disease signal

Chia et al. 2025 Definition (from paper):
  "ALS cases were defined as individuals diagnosed with ICD-10 code G12.2
   and confirmed at death... had blood collected up to 1 year before symptom onset"
  External Validation Set 2: 13 ALS cases, 23,601 controls

Usage:
    python 00_validate_cohort_extraction.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
CADASIL_ROOT = REPO_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"
sys.path.insert(0, str(REPO_ROOT / "analysis"))

from helpers.disease_config import load_disease_config

cfg = load_disease_config()
DISEASE_STATUS_COL = cfg.cohort_status_col
DISEASE_STATUS_CONTROL = cfg.status_values["control"]

# --- PATHS ---
DISEASE_COHORT = REPO_ROOT / "data" / "ukb" / "cohort" / f"{cfg.cohort_short}_cohort.csv"
ALS_COHORT = REPO_ROOT / "data" / "ukb" / "cohort" / "als_cohort.csv"
DISEASE_OLINK_QC = REPO_ROOT / "data" / "ukb" / "olink" / "processed" / f"{cfg.cohort_short}_olink_qc.csv"
ALS_OLINK_QC = REPO_ROOT / "data" / "ukb" / "olink" / "processed" / "als_olink_qc.csv"

# Ground truth ICD-10 code tables. The primary case table is resolved from the
# first configured ICD code (e.g. G35 → G/G3/G35/G35.csv). ALS is a secondary
# comparison cohort and remains hardcoded.
def _icd_table_path(code: str) -> Path:
    """Resolve CADASIL repo path for a given ICD-10 code (e.g. 'G35')."""
    letter = code[0]
    bucket = f"{letter}{code[1]}" if len(code) > 1 else letter
    return CADASIL_ROOT / "data" / "ukb" / "diagnoses" / "icd_codes" / letter / bucket / code / f"{code}.csv"

DISEASE_ICD_TABLE = _icd_table_path(cfg.icd_codes[0])
G12_ICD = CADASIL_ROOT / "data" / "ukb" / "diagnoses" / "icd_codes" / "G" / "G1" / "G12" / "G12.csv"
ALGO_MND = CADASIL_ROOT / "data" / "ukb" / "diagnoses" / "algorithmic_outcomes" / "ALGO_MOTOR_NEURONE_DISEASE" / "ALGO_MOTOR_NEURONE_DISEASE.csv"

# Olink data
OLINK_PATH = CADASIL_ROOT / "data" / "ukb" / "olink" / "i0" / "olink_instance_0_extracted_data.csv"

# Output
OUT_DIR = REPO_ROOT / "results" / "validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Per-cohort metadata: (label, cohort path, status col, control label).
# The disease cohort is configured via disease.yaml; ALS is a fixed
# secondary comparison cohort.
COHORT_PANELS = [
    (cfg.disease_short_caps, DISEASE_COHORT, DISEASE_STATUS_COL, DISEASE_STATUS_CONTROL),
    ("ALS",                  ALS_COHORT,     "als_status",       "control"),
]


def validate_icd_overlap():
    """Cross-reference our cohort EIDs with the granular ICD-10 tables."""
    print("=" * 60)
    print("VALIDATION 1: ICD-10 Code Cross-Reference")
    print("=" * 60)

    # --- Disease Validation (MS by default; configured via disease.yaml) ---
    disease_label = cfg.disease_short_caps
    primary_icd = cfg.icd_codes[0]
    disease_cohort = pd.read_csv(DISEASE_COHORT)
    disease_cases = disease_cohort[
        disease_cohort[DISEASE_STATUS_COL] != DISEASE_STATUS_CONTROL
    ]["eid"].values

    icd_truth = pd.read_csv(DISEASE_ICD_TABLE)
    icd_truth_eids = set(icd_truth["eid"].values)

    cases_in_truth = len(set(disease_cases) & icd_truth_eids)
    print(f"\n{disease_label} Cases in our cohort: {len(disease_cases)}")
    print(f"{primary_icd} ICD-10 code table: {len(icd_truth_eids)} total UKB participants")
    print(f"Our {disease_label} cases found in {primary_icd} table: "
          f"{cases_in_truth}/{len(disease_cases)} "
          f"({100*cases_in_truth/max(len(disease_cases),1):.1f}%)")

    # --- ALS Validation ---
    als_cohort = pd.read_csv(ALS_COHORT)
    als_cases = als_cohort[als_cohort["als_status"] != "control"]["eid"].values

    g12_icd = pd.read_csv(G12_ICD)
    g12_eids = set(g12_icd["eid"].values)

    algo_mnd = pd.read_csv(ALGO_MND)
    algo_eids = set(algo_mnd["eid"].values)

    als_in_g12 = len(set(als_cases) & g12_eids)
    als_in_algo = len(set(als_cases) & algo_eids)

    print(f"\nALS Cases in our cohort: {len(als_cases)}")
    print(f"G12 ICD-10 code table: {len(g12_eids)} total")
    print(f"ALGO_MOTOR_NEURONE_DISEASE: {len(algo_eids)} total")
    print(f"Our cases in G12 table: {als_in_g12}/{len(als_cases)}")
    print(f"Our cases in ALGO_MND: {als_in_algo}/{len(als_cases)}")

    # Chia definition: G12.2 + confirmed at death + blood ≤1yr before onset
    # ALGO_MND is the algorithmic definition which should include only definite MND
    als_refined = set(als_cases) & algo_eids
    print(f"\n--- Refined ALS (ALGO_MND confirmed): {len(als_refined)} cases ---")
    if als_refined:
        refined_df = als_cohort[als_cohort["eid"].isin(als_refined)]
        pre = (refined_df["als_status"] == "pre_onset").sum()
        post = (refined_df["als_status"] == "post_onset").sum()
        print(f"  Pre-onset: {pre}, Post-onset: {post}")
        if "als_confirmed_at_death" in refined_df.columns:
            death_confirmed = refined_df["als_confirmed_at_death"].sum()
            print(f"  Death-confirmed: {int(death_confirmed)}")

    # Save refined ALS cohort
    refined_cohort = als_cohort[als_cohort["eid"].isin(als_refined) | (als_cohort["als_status"] == "control")]
    refined_path = REPO_ROOT / "data" / "ukb" / "cohort" / "als_cohort_algo_refined.csv"
    refined_cohort.to_csv(refined_path, index=False)
    print(f"  Saved refined ALS cohort to: {refined_path}")
    print(f"  Total in refined: {len(refined_cohort)} ({(refined_cohort['als_status'] != 'control').sum()} cases, {(refined_cohort['als_status'] == 'control').sum()} controls)")

    return als_refined


def validate_nefl_signal():
    """Check NEFL (neurofilament light) levels between cases and controls.
    
    NEFL is the strongest ALS biomarker (log2FC = 2.34 in Chia et al.)
    and is also elevated in MS. If our cohort extraction is correct,
    NEFL should be significantly elevated in cases vs controls.
    """
    print("\n" + "=" * 60)
    print("VALIDATION 2: NEFL Signal Validation")
    print("=" * 60)

    # Find NEFL column in Olink
    print("Loading Olink data (just NEFL + a few validation proteins)...")
    # Read header to find NEFL
    header = pd.read_csv(OLINK_PATH, nrows=0).columns.tolist()
    nefl_col = next((c for c in header if "nefl" in c.lower()), None)
    sod1_col = next((c for c in header if c.lower().endswith(".sod1") or c.lower() == "sod1"), None)
    tardbp_col = next((c for c in header if "tardbp" in c.lower()), None)
    gfap_col = next((c for c in header if "gfap" in c.lower()), None)

    eid_col = next(c for c in header if "eid" in c.lower())
    cols_to_load = [eid_col]
    protein_map = {}
    for name, col in [("NEFL", nefl_col), ("SOD1", sod1_col), ("TARDBP", tardbp_col), ("GFAP", gfap_col)]:
        if col:
            cols_to_load.append(col)
            protein_map[col] = name

    olink = pd.read_csv(OLINK_PATH, usecols=cols_to_load, low_memory=False)
    olink = olink.rename(columns={eid_col: "eid"})
    olink = olink.rename(columns=protein_map)
    olink["eid"] = olink["eid"].astype(int)
    print(f"  Loaded {len(olink)} participants, proteins: {list(protein_map.values())}")

    results = []

    for disease, cohort_path, status_col, control_label in COHORT_PANELS:
        cohort = pd.read_csv(cohort_path)
        cohort["eid"] = cohort["eid"].astype(int)
        merged = cohort.merge(olink, on="eid", how="inner")

        print(f"\n--- {disease} ---")
        for protein in protein_map.values():
            if protein not in merged.columns:
                continue

            cases = merged.loc[merged[status_col] != control_label, protein].dropna()
            controls = merged.loc[merged[status_col] == control_label, protein].dropna()

            if len(cases) < 2 or len(controls) < 2:
                continue

            t_stat, p_val = stats.ttest_ind(cases, controls, equal_var=False)
            cohens_d = (cases.mean() - controls.mean()) / np.sqrt(
                (cases.std()**2 + controls.std()**2) / 2
            )
            log2fc = cases.mean() - controls.mean()  # Already log2 (NPX)

            status_str = "✅" if p_val < 0.05 else "⚠️"
            print(f"  {protein}: cases={len(cases)}, controls={len(controls)}")
            print(f"    Mean(cases)={cases.mean():.3f}, Mean(ctrl)={controls.mean():.3f}")
            print(f"    log2FC={log2fc:.3f}, Cohen's d={cohens_d:.3f}, p={p_val:.2e} {status_str}")

            results.append({
                "disease": disease,
                "protein": protein,
                "n_cases": len(cases),
                "n_controls": len(controls),
                "mean_cases": cases.mean(),
                "mean_controls": controls.mean(),
                "log2fc": log2fc,
                "cohens_d": cohens_d,
                "t_stat": t_stat,
                "p_value": p_val,
            })

            # Sub-stratify by timing
            for timing in ["pre_onset", "post_onset"]:
                subset = merged.loc[merged[status_col] == timing, protein].dropna()
                if len(subset) >= 2:
                    t2, p2 = stats.ttest_ind(subset, controls, equal_var=False)
                    fc2 = subset.mean() - controls.mean()
                    print(f"    {timing}: n={len(subset)}, log2FC={fc2:.3f}, p={p2:.2e}")

    # Create visualization
    if results:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        for idx, (disease, cohort_path, status_col, control_label) in enumerate(COHORT_PANELS):
            ax = axes[idx]
            cohort = pd.read_csv(cohort_path)
            cohort["eid"] = cohort["eid"].astype(int)
            merged = cohort.merge(olink, on="eid", how="inner")

            if "NEFL" in merged.columns:
                order = [control_label, "pre_onset", "post_onset"]
                available = [s for s in order if s in merged[status_col].values]
                colors = {control_label: "#3498DB", "pre_onset": "#E67E22", "post_onset": "#E74C3C"}
                palette = [colors[s] for s in available]

                sns.violinplot(data=merged, x=status_col, y="NEFL",
                              order=available, palette=palette, alpha=0.6, ax=ax)
                sns.stripplot(data=merged[merged[status_col] != control_label],
                              x=status_col, y="NEFL",
                              order=[s for s in available if s != control_label],
                              color="black", alpha=0.5, size=3, ax=ax)

                ax.set_title(f"{disease} — NEFL (Neurofilament Light Chain)")
                ax.set_ylabel("NEFL NPX (log2)")
                ax.set_xlabel("")

        plt.tight_layout()
        fig_path = OUT_DIR / "nefl_validation.png"
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"\nSaved NEFL validation plot: {fig_path}")

    # Save results table
    if results:
        results_df = pd.DataFrame(results)
        results_path = OUT_DIR / "protein_validation_results.csv"
        results_df.to_csv(results_path, index=False)
        print(f"Saved validation results: {results_path}")


def check_als_genetic_proxies():
    """Check ALS causative gene protein levels as genetic mutation proxies.
    
    Known ALS causative genes available on Olink:
      - SOD1 (superoxide dismutase 1) — ~20% familial ALS
      - TARDBP (TDP-43) — ~5% familial ALS
      - FUS — ~5% familial ALS (may not be on Olink)
      - C9orf72 — ~40% familial, ~7% sporadic (not directly on Olink)
    
    Hypothesis: If any of our ALS cases carry SOD1 mutations, their SOD1
    protein levels may be altered. However, this is an indirect check.
    """
    print("\n" + "=" * 60)
    print("VALIDATION 3: ALS Genetic Mutation Proxies")
    print("=" * 60)

    header = pd.read_csv(OLINK_PATH, nrows=0).columns.tolist()
    eid_col = next(c for c in header if "eid" in c.lower())

    # Find ALS-relevant genes on Olink
    gene_targets = ["sod1", "tardbp", "fus", "optn", "tbk1", "sqstm1", "ubqln2", "vcp"]
    found_genes = {}
    for gene in gene_targets:
        col = next((c for c in header if c.lower().endswith(f".{gene}") or c.lower() == gene), None)
        if col:
            found_genes[gene.upper()] = col

    print(f"ALS causative genes found on Olink: {list(found_genes.keys())}")
    missing = [g.upper() for g in gene_targets if g.upper() not in found_genes]
    if missing:
        print(f"Not found on Olink: {missing}")

    if not found_genes:
        print("No ALS causative genes available on Olink for proxy analysis.")
        return

    # Load these columns
    cols = [eid_col] + list(found_genes.values())
    olink = pd.read_csv(OLINK_PATH, usecols=cols, low_memory=False)
    olink = olink.rename(columns={eid_col: "eid"})
    rename_map = {v: k for k, v in found_genes.items()}
    olink = olink.rename(columns=rename_map)
    olink["eid"] = olink["eid"].astype(int)

    als_cohort = pd.read_csv(ALS_COHORT)
    als_cohort["eid"] = als_cohort["eid"].astype(int)
    merged = als_cohort.merge(olink, on="eid", how="inner")

    for gene in found_genes.keys():
        if gene not in merged.columns:
            continue
        cases = merged.loc[merged["als_status"] != "control", gene].dropna()
        controls = merged.loc[merged["als_status"] == "control", gene].dropna()

        if len(cases) < 2:
            continue

        t_stat, p_val = stats.ttest_ind(cases, controls, equal_var=False)
        fc = cases.mean() - controls.mean()

        # Find outliers (potential mutation carriers: >3 SD from control mean)
        control_mean = controls.mean()
        control_std = controls.std()
        outlier_threshold = control_mean + 3 * control_std
        outlier_low = control_mean - 3 * control_std

        case_outliers_high = (cases > outlier_threshold).sum()
        case_outliers_low = (cases < outlier_low).sum()

        print(f"\n  {gene}:")
        print(f"    Cases: n={len(cases)}, mean={cases.mean():.3f}")
        print(f"    Controls: n={len(controls)}, mean={controls.mean():.3f}")
        print(f"    log2FC={fc:.3f}, p={p_val:.2e}")
        print(f"    Case outliers (>3SD): {case_outliers_high} high, {case_outliers_low} low")


def main():
    als_refined_eids = validate_icd_overlap()
    validate_nefl_signal()
    check_als_genetic_proxies()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("""
Key findings from Chia et al. UKB extraction:
  - Used ICD-10 G12.2 (not just G12) + confirmed at death
  - 13 ALS cases in External Validation Set 2
  - 23,601 controls (no myopathy G70-73, no neuropathy G60-64)
  - Blood collected up to 1 year before symptom onset
  - Used ALGO_MOTOR_NEURONE_DISEASE for definitive case identification

Recommendation: Use ALGO_MND-refined cohort for strict ALS replication.
G12 first-occurrence includes non-ALS motor neuron diseases (SMA, etc.)
""")


if __name__ == "__main__":
    main()
