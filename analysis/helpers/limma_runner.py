"""Limma-equivalent differential protein abundance analysis in Python.

Replicates the limma pipeline used in:
  - Abdelhak et al. 2026 (MS): limma with age, sex, combined NPX
  - Chia et al. 2025 (ALS): limma with age, sex, tube type, UMAP1/2

Implementation:
  - OLS regression per protein: NPX ~ disease_status + age + sex [+ combined_npx]
  - Empirical Bayes variance shrinkage (limma's eBayes)
  - FDR correction via Benjamini-Hochberg

Usage:
    from helpers.limma_runner import run_limma_analysis

Author: Jonah Keller (AI-assisted)
Date: 2026-04-14
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


@dataclass
class LimmaResult:
    """Container for differential protein abundance results."""
    protein: str
    log2fc: float
    se: float
    t_stat: float
    p_value: float
    fdr: float  # populated after all proteins are tested
    n_cases: int
    n_controls: int
    mean_cases: float
    mean_controls: float


def _empirical_bayes_shrinkage(
    coefficients: np.ndarray,
    standard_errors: np.ndarray,
    df_residuals: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply limma-style empirical Bayes moderation to t-statistics.

    This shrinks protein-specific variances toward a common prior,
    reducing false positives from low-variance proteins and increasing
    power for high-variance proteins.

    Implements the Smyth (2004) eBayes procedure:
      s2_post = (d0*s0^2 + df*s2) / (d0 + df)
      t_mod = coef / sqrt(s2_post / n)

    Returns moderated t-stats, moderated p-values, and posterior variances.
    """
    # Protein-specific sample variances
    s2 = standard_errors ** 2

    # Estimate prior parameters via method of moments
    # (simplified — full implementation uses trigamma/digamma fitting)
    valid = np.isfinite(s2) & (s2 > 0) & np.isfinite(df_residuals) & (df_residuals > 0)
    if valid.sum() < 3:
        # Not enough data for shrinkage, return unmodified
        t_mod = coefficients / standard_errors
        p_mod = 2 * stats.t.sf(np.abs(t_mod), df_residuals)
        return t_mod, p_mod, s2

    log_s2 = np.log(s2[valid])
    df_valid = df_residuals[valid]

    # Estimate prior degrees of freedom (d0) and prior variance (s0^2)
    # Using the method from Smyth (2004) — simplified trigamma approach
    mean_log_s2 = np.mean(log_s2)
    var_log_s2 = np.var(log_s2, ddof=1)
    mean_df = np.mean(df_valid)

    # trigamma(df/2) ≈ 2/df for large df
    trigamma_approx = 2.0 / mean_df
    d0 = max(2.0 * trigamma_approx / max(var_log_s2 - trigamma_approx, 1e-10), 0.5)
    d0 = min(d0, 1000)  # Cap prior df

    s0_sq = np.exp(mean_log_s2 - (np.log(d0 / 2.0) - np.log(mean_df / 2.0) if d0 > 0 else 0))
    s0_sq = max(s0_sq, 1e-10)

    # Posterior variance: weighted average of prior and sample variance
    s2_post = np.full_like(s2, np.nan)
    s2_post[valid] = (d0 * s0_sq + df_valid * s2[valid]) / (d0 + df_valid)

    # Moderated t-statistic
    t_mod = np.full_like(coefficients, np.nan)
    se_post = np.sqrt(np.where(s2_post > 0, s2_post, np.nan))
    mask = np.isfinite(se_post) & (se_post > 0) & np.isfinite(coefficients)
    t_mod[mask] = coefficients[mask] / se_post[mask]

    # Moderated p-value with increased df
    df_total = np.full_like(df_residuals, np.nan, dtype=float)
    df_total[valid] = d0 + df_valid
    p_mod = np.full_like(t_mod, np.nan)
    finite_mask = np.isfinite(t_mod) & np.isfinite(df_total) & (df_total > 0)
    p_mod[finite_mask] = 2 * stats.t.sf(
        np.abs(t_mod[finite_mask]),
        df_total[finite_mask],
    )

    return t_mod, p_mod, s2_post


def run_limma_analysis(
    data: pd.DataFrame,
    protein_cols: list[str],
    group_col: str,
    case_value: str,
    control_value: str = "control",
    covariates: Optional[list[str]] = None,
    add_combined_npx: bool = False,
    fdr_threshold: float = 0.05,
    min_nonmissing: int = 5,
) -> pd.DataFrame:
    """Run limma-style differential protein abundance analysis.

    Parameters
    ----------
    data : pd.DataFrame
        Merged cohort + Olink data.
    protein_cols : list[str]
        Column names for proteins (NPX values).
    group_col : str
        Column containing disease status labels.
    case_value : str
        Value in group_col indicating disease cases.
    control_value : str
        Value in group_col indicating controls.
    covariates : list[str], optional
        Additional covariates to adjust for (e.g., ["age_at_sampling", "sex"]).
    add_combined_npx : bool
        If True, add mean NPX across all proteins as a covariate
        (Abdelhak correction for sample quality).
    fdr_threshold : float
        FDR threshold for significance reporting.
    min_nonmissing : int
        Minimum non-missing values per group to test a protein.

    Returns
    -------
    pd.DataFrame
        Results table sorted by p-value, with columns:
        protein, log2fc, se, t_stat, p_value, fdr, n_cases, n_controls,
        mean_cases, mean_controls, significant
    """
    if covariates is None:
        covariates = []

    # Filter to cases and controls only
    mask = data[group_col].isin([case_value, control_value])
    df = data.loc[mask].copy()

    # Create binary disease indicator
    df["_disease"] = (df[group_col] == case_value).astype(float)

    # Optionally compute combined NPX (average across all proteins per sample)
    # Use nanmean to handle missing values — avoids dropna removing rows
    if add_combined_npx:
        prot_matrix = df[protein_cols].values
        df["_combined_npx"] = np.nanmean(prot_matrix, axis=1)
        design_cols = ["_disease"] + covariates + ["_combined_npx"]
    else:
        design_cols = ["_disease"] + covariates

    # Build design matrix
    design_available = [c for c in design_cols if c in df.columns]

    # Storage for per-protein results
    coefficients = []
    standard_errors = []
    df_residuals = []
    protein_names = []
    n_cases_list = []
    n_controls_list = []
    mean_cases_list = []
    mean_controls_list = []

    n_proteins = len(protein_cols)
    for i, prot in enumerate(protein_cols):
        if (i + 1) % 500 == 0 or i == 0:
            print(f"  Processing protein {i + 1}/{n_proteins}...")

        # Get complete cases for this protein + design
        cols_needed = [prot] + design_available
        sub = df[cols_needed].dropna()

        n_case = (sub["_disease"] == 1).sum()
        n_ctrl = (sub["_disease"] == 0).sum()

        if n_case < min_nonmissing or n_ctrl < min_nonmissing:
            continue

        y = sub[prot].values.astype(np.float64)
        X = np.column_stack(
            [np.ones(len(sub), dtype=np.float64)]
            + [sub[c].values.astype(np.float64) for c in design_available]
        )

        try:
            # OLS: y = Xβ + ε
            beta, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
            y_hat = X @ beta
            resid = y - y_hat
            n = len(y)
            p = X.shape[1]
            df_resid = n - p

            if df_resid <= 0:
                continue

            mse = np.sum(resid**2) / df_resid
            # Standard error of the disease coefficient (index 1)
            try:
                cov_matrix = mse * np.linalg.inv(X.T @ X)
                se_disease = np.sqrt(cov_matrix[1, 1])
            except np.linalg.LinAlgError:
                continue

            coef_disease = beta[1]  # log2FC for disease
            coefficients.append(coef_disease)
            standard_errors.append(se_disease)
            df_residuals.append(df_resid)
            protein_names.append(prot)
            n_cases_list.append(n_case)
            n_controls_list.append(n_ctrl)

            case_vals = sub.loc[sub["_disease"] == 1, prot]
            ctrl_vals = sub.loc[sub["_disease"] == 0, prot]
            mean_cases_list.append(case_vals.mean())
            mean_controls_list.append(ctrl_vals.mean())

        except Exception:
            continue

    if not coefficients:
        return pd.DataFrame()

    # Convert to arrays
    coefs = np.array(coefficients)
    ses = np.array(standard_errors)
    dfs = np.array(df_residuals)

    # Apply empirical Bayes shrinkage
    t_mod, p_mod, s2_post = _empirical_bayes_shrinkage(coefs, ses, dfs)

    # FDR correction
    valid_p = np.isfinite(p_mod)
    fdr_vals = np.full_like(p_mod, np.nan)
    if valid_p.sum() > 0:
        _, fdr_corrected, _, _ = multipletests(
            p_mod[valid_p], method="fdr_bh", alpha=fdr_threshold,
        )
        fdr_vals[valid_p] = fdr_corrected

    # Build results table
    results = pd.DataFrame({
        "protein": protein_names,
        "log2fc": coefs,
        "se": ses,
        "t_stat": t_mod,
        "p_value": p_mod,
        "fdr": fdr_vals,
        "n_cases": n_cases_list,
        "n_controls": n_controls_list,
        "mean_cases": mean_cases_list,
        "mean_controls": mean_controls_list,
    })

    results["significant_fdr05"] = results["fdr"] < 0.05
    results["significant_fdr10"] = results["fdr"] < 0.10
    results = results.sort_values("p_value").reset_index(drop=True)

    return results


def make_volcano_data(results: pd.DataFrame) -> pd.DataFrame:
    """Add volcano plot columns to results."""
    df = results.copy()
    df["-log10(p)"] = -np.log10(df["p_value"].clip(lower=1e-300))
    df["-log10(fdr)"] = -np.log10(df["fdr"].clip(lower=1e-300))
    return df


def print_summary(
    results: pd.DataFrame,
    label: str = "",
    fdr_threshold: float = 0.05,
) -> None:
    """Print a summary of differential analysis results."""
    if results.empty:
        print(f"\n{'=' * 50}")
        if label:
            print(f"  {label}")
        print("  No results (empty DataFrame)")
        print(f"{'=' * 50}")
        return

    n_sig_05 = (results["fdr"] < 0.05).sum()
    n_sig_10 = (results["fdr"] < 0.10).sum()
    n_tested = len(results)

    print(f"\n{'=' * 50}")
    if label:
        print(f"  {label}")
    print(f"  Proteins tested:          {n_tested}")
    print(f"  Significant (FDR < 5%):   {n_sig_05}")
    print(f"  Significant (FDR < 10%):  {n_sig_10}")

    if n_sig_05 > 0:
        top = results[results["fdr"] < 0.05].head(10)
        print(f"\n  Top hits (FDR < 5%):")
        for _, r in top.iterrows():
            direction = "↑" if r["log2fc"] > 0 else "↓"
            prot = r["protein"].split(".")[-1] if "." in r["protein"] else r["protein"]
            print(f"    {direction} {prot:20s} log2FC={r['log2fc']:+.3f}  FDR={r['fdr']:.2e}")
    print(f"{'=' * 50}")
