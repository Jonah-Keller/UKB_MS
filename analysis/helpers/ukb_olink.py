"""UKB Olink Explore 3072 data loading utilities.

Handles the wide-format CSV that UKB exports from the RAP, strips the
instance prefix from column names, and provides helpers for identifying
protein vs. metadata columns.

Usage:
    from helpers.ukb_olink import load_olink_wide, get_protein_cols
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Columns that are metadata, not protein NPX values
_METADATA_PATTERNS = (
    "eid", "ms_status", "als_status", "age_at_sampling", "age_at_diagnosis",
    "years_to_diagnosis", "olink_instance", "sex", "combined_npx",
    "age_i", "bmi", "townsend", "assessment_center", "smoking",
    "als_confirmed_at_death", "hla_drb1", "UMAP", "PC",
)


def load_olink_wide(path: Path | str) -> pd.DataFrame:
    """Load UKB Olink instance-0 CSV and strip the column prefix.

    UKB exports columns named 'olink_instance_0.nefl', etc.
    This function strips the prefix to give plain gene-symbol columns.

    Args:
        path: Path to olink_instance_0_extracted_data.csv

    Returns:
        Wide DataFrame with columns: eid, <gene_symbol>, ...
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Olink file not found: {path}\n"
            "Expected: data/ukb/olink/i0/olink_instance_0_extracted_data.csv\n"
            "Extract from UKB RAP or use make_mock_olink() for testing."
        )
    df = pd.read_csv(path, low_memory=False)
    df.columns = [c.removeprefix("olink_instance_0.") for c in df.columns]
    # Ensure eid is int
    if "eid" in df.columns:
        df["eid"] = df["eid"].astype(int)
    return df


def get_protein_cols(df: pd.DataFrame) -> list[str]:
    """Return list of protein (NPX) column names, excluding metadata."""
    def _is_meta(col: str) -> bool:
        return any(col.lower().startswith(p.lower()) for p in _METADATA_PATTERNS)
    return [c for c in df.columns if not _is_meta(c)]


def compute_combined_npx(df: pd.DataFrame, protein_cols: list[str]) -> pd.Series:
    """Compute per-sample mean NPX across all proteins.

    Used as a covariate in Abdelhak et al. combined MS analysis to
    account for sample-quality / storage-time variation.
    """
    return df[protein_cols].mean(axis=1)


def flag_npx_outliers(
    df: pd.DataFrame,
    protein_cols: list[str],
    sd_threshold: float = 3.0,
) -> pd.Series:
    """Flag samples whose median NPX is >sd_threshold SDs from cohort median.

    Returns:
        Boolean Series — True = outlier, indexed like df.
    """
    medians = df[protein_cols].median(axis=1)
    mu = medians.mean()
    sigma = medians.std()
    return (medians - mu).abs() > sd_threshold * sigma


def make_mock_olink(
    n: int = 300,
    n_proteins: int = 50,
    seed: int = 42,
    eids: Optional[list[int]] = None,
) -> pd.DataFrame:
    """Create a mock Olink wide DataFrame for pipeline testing.

    Includes key proteins from both papers so downstream scripts can
    find them by name.
    """
    rng = np.random.default_rng(seed)
    named_proteins = [
        "nefl", "mog", "gfap", "nrgn", "il3", "aif1", "sirt6", "dock10",
        "kcnd2", "sanbr", "vwa1", "spag8", "cdin1", "lctl", "insm2",
        "chit1", "gpnmb", "lgals3", "s100a12", "chi3l1", "ctsh", "spp1",
    ]
    filler = [f"protein_{i}" for i in range(n_proteins - len(named_proteins))]
    all_proteins = named_proteins + filler

    data = rng.standard_normal((n, len(all_proteins)))
    df = pd.DataFrame(data, columns=all_proteins)
    df.insert(0, "eid", eids if eids is not None else range(1000, 1000 + n))
    return df
