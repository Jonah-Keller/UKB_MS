#!/usr/bin/env python3
"""Parse HES R/S/T/U/V/Y/Z ICD-10 chapters into per-code icd_codes/ files.

Run after downloading icd_hes_rstz_first_occurrence.csv from RAP:
  python analysis/00_rap_extraction/10_parse_hes_rstz.py

Input:
  data/ukb/rap_extraction/icd_hes_rstz_first_occurrence.csv

Output:
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/R/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/S/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/T/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/U/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/V/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/Y/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/Z/...
  Same eid, diagnosis, age_at_diagnosis format as existing G-Q files.

Summary:
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/parse_rstz_summary.csv
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJ_ROOT   = Path(__file__).resolve().parents[2]
CADASIL_DIR = PROJ_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"
ICD_DIR     = CADASIL_DIR / "data" / "ukb" / "diagnoses" / "icd_codes"
RAP_DIR     = PROJ_ROOT / "data" / "ukb" / "rap_extraction"

TARGET_CHAPTERS = set("RSTUVYZ")

# Minimum cases to write a file (avoid empty/noise files)
MIN_CASES = 5


def icd_output_path(icd_code: str, base: Path) -> Path:
    clean = icd_code.replace(".", "")
    l1 = clean[0]
    l2 = clean[:2] if len(clean) >= 2 else l1
    l3 = clean[:3] if len(clean) >= 3 else l2
    if len(clean) <= 3:
        return base / l1 / l2 / l3 / f"{icd_code}.csv"
    return base / l1 / l2 / l3 / icd_code / f"{icd_code}.csv"


def parse_hes_rstz(csv_path: Path, out_dir: Path, min_cases: int = MIN_CASES) -> pd.DataFrame:
    log.info("Loading %s ...", csv_path.name)
    df = pd.read_csv(csv_path, low_memory=False)
    log.info("  %d rows × %d cols", *df.shape)

    # Expected columns: eid, icd3, chapter, first_occurrence_date, age_at_diagnosis
    if "icd3" not in df.columns:
        raise ValueError("Expected 'icd3' column — is this the right file?")

    # Filter to target chapters (should already be filtered, but be safe)
    df = df[df["chapter"].isin(TARGET_CHAPTERS)].copy()
    log.info("  After chapter filter: %d rows", len(df))

    summary_rows = []
    for icd3, grp in df.groupby("icd3"):
        if icd3[0] not in TARGET_CHAPTERS:
            continue

        # One row per eid — keep the earliest date (already done in RAP script, but enforce)
        grp_dedup = grp.sort_values("first_occurrence_date").drop_duplicates("eid")

        n = len(grp_dedup)
        if n < min_cases:
            log.debug("  %s: n=%d < %d, skipping", icd3, n, min_cases)
            continue

        out = pd.DataFrame({
            "eid":              grp_dedup["eid"].values,
            "diagnosis":        1,
            "age_at_diagnosis": grp_dedup["age_at_diagnosis"].values,
        })

        out_path = icd_output_path(icd3, out_dir)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_path, index=False)

        median_age = grp_dedup["age_at_diagnosis"].dropna().median()
        log.info("  %-6s  %6d cases  median_age=%.1f", icd3, n, median_age if not np.isnan(median_age) else -1)
        summary_rows.append({
            "icd_code":    icd3,
            "chapter":     icd3[0],
            "n_events":    n,
            "median_age":  median_age,
            "output_path": str(out_path),
        })

    return pd.DataFrame(summary_rows)


def main():
    ap = argparse.ArgumentParser(description="Parse HES R/S/T/U/V/Y/Z chapters")
    ap.add_argument(
        "--input", type=Path,
        default=RAP_DIR / "icd_hes_rstz_first_occurrence.csv",
    )
    ap.add_argument("--out-dir", type=Path, default=ICD_DIR)
    ap.add_argument("--min-cases", type=int, default=MIN_CASES,
                    help="Minimum cases to write a file (default: 5)")
    args = ap.parse_args()

    if not args.input.exists():
        log.error("Input not found: %s", args.input)
        log.error("Download from RAP first:")
        log.error('  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_hes_rstz_first_occurrence.csv" \\')
        log.error('      --output data/ukb/rap_extraction/icd_hes_rstz_first_occurrence.csv')
        raise SystemExit(1)

    summary = parse_hes_rstz(args.input, args.out_dir, min_cases=args.min_cases)

    summary_path = args.out_dir.parent / "parse_rstz_summary.csv"
    summary.to_csv(summary_path, index=False)
    log.info("\nParsed %d R/S/T/U/V/Y/Z codes (min_cases=%d)", len(summary), args.min_cases)
    log.info("Summary: %s", summary_path)

    if len(summary) > 0:
        log.info("\nPer-chapter breakdown:")
        for ch, grp in summary.groupby("chapter"):
            top = grp.loc[grp["n_events"].idxmax(), "icd_code"]
            log.info("  Chapter %s: %d codes, total=%d, top=%s (n=%d)",
                     ch, len(grp), grp["n_events"].sum(), top, grp["n_events"].max())

        log.info("\nTop 20 R/S/T/Z codes by n_events:")
        top20 = summary.sort_values("n_events", ascending=False).head(20)
        log.info("\n%s", top20[["icd_code", "chapter", "n_events", "median_age"]].to_string(index=False))


if __name__ == "__main__":
    main()
