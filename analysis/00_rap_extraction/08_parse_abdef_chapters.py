#!/usr/bin/env python3
"""Parse missing ICD-10 chapters A, B, D, E, F from RAP first-occurrence extract.

Run after downloading the two files from RAP:
  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_abdef_first_occurrence.csv"
  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/abdef_field_mapping.json"

Place both in:  data/ukb/rap_extraction/

Then run:
  python analysis/00_rap_extraction/08_parse_abdef_chapters.py

Output:
  Per-code CSV files added to the shared icd_codes/ directory:
    CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/A/
    CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/B/
    CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/D/  (supplements cancer data)
    CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/E/
    CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/F/

  Same format as existing G-Q files:  eid, diagnosis, age_at_diagnosis

Summary CSV:
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/parse_abdef_summary.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJ_ROOT   = Path(__file__).resolve().parents[2]
CADASIL_DIR = PROJ_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"
ICD_DIR     = CADASIL_DIR / "data" / "ukb" / "diagnoses" / "icd_codes"
RAP_DIR     = PROJ_ROOT / "data" / "ukb" / "rap_extraction"

CHAPTERS_TO_PARSE = set("ABDEF")


def icd_output_path(icd_code: str, base: Path) -> Path:
    clean = icd_code.replace(".", "")
    l1 = clean[0]
    l2 = clean[:2] if len(clean) >= 2 else l1
    l3 = clean[:3] if len(clean) >= 3 else l2
    if len(clean) <= 3:
        return base / l1 / l2 / l3 / f"{icd_code}.csv"
    return base / l1 / l2 / l3 / icd_code / f"{icd_code}.csv"


def load_field_mapping(json_path: Path) -> dict[str, tuple[str, str]]:
    """Load field_id (as str) -> (icd10_code, description) from JSON."""
    if not json_path.exists():
        raise FileNotFoundError(f"Field mapping not found: {json_path}")
    with open(json_path) as f:
        raw = json.load(f)
    result = {}
    for k, v in raw.items():
        col = f"p{k}" if not str(k).startswith("p") else str(k)
        result[col] = (str(v[0]), str(v[1]))
    log.info("Loaded %d field mappings from %s", len(result), json_path)
    return result


def parse_chapters(
    csv_path: Path,
    field_map: dict[str, tuple[str, str]],
    out_dir: Path,
    chapters: set[str] = CHAPTERS_TO_PARSE,
) -> pd.DataFrame:
    log.info("Loading %s ...", csv_path.name)
    df = pd.read_csv(csv_path, low_memory=False)
    log.info("  %d rows × %d cols", *df.shape)

    # Resolve eid and year_birth columns
    eid_col = next((c for c in df.columns if c in ("eid", "participant_id")), None)
    yob_col = next((c for c in df.columns if c in ("year_birth", "p34")), None)
    if eid_col is None:
        raise ValueError("No eid or participant_id column found")
    if yob_col is None:
        raise ValueError("No year_birth or p34 column found — needed for age computation")

    log.info("  eid_col=%s  yob_col=%s", eid_col, yob_col)

    summary_rows = []
    processed = 0
    skipped_chapter = 0
    skipped_empty = 0

    for col, (icd_code, desc) in field_map.items():
        if icd_code[0] not in chapters:
            skipped_chapter += 1
            continue
        if col not in df.columns:
            log.debug("Column %s (%s) not in CSV — skipping", col, icd_code)
            skipped_empty += 1
            continue

        sub = df[[eid_col, yob_col, col]].dropna(subset=[col]).copy()
        if len(sub) == 0:
            skipped_empty += 1
            continue

        sub["_date"] = pd.to_datetime(sub[col], errors="coerce")
        sub = sub.dropna(subset=["_date"])
        if len(sub) == 0:
            skipped_empty += 1
            continue

        sub["age_at_diagnosis"] = sub["_date"].dt.year - sub[yob_col].fillna(0).astype(int)
        sub.loc[sub["age_at_diagnosis"] < 0, "age_at_diagnosis"] = np.nan

        out = pd.DataFrame({
            "eid":              sub[eid_col].values,
            "diagnosis":        1,
            "age_at_diagnosis": sub["age_at_diagnosis"].values,
        })

        out_path = icd_output_path(icd_code, out_dir)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_path, index=False)
        processed += 1

        log.info("  %-6s  %-55s  %6d events", icd_code, desc[:55], len(out))
        summary_rows.append({
            "icd_code":    icd_code,
            "chapter":     icd_code[0],
            "description": desc,
            "n_events":    len(out),
            "median_age":  out["age_at_diagnosis"].median(),
            "output_path": str(out_path),
        })

    log.info(
        "\nParsed: %d codes | Skipped empty: %d | Skipped wrong chapter: %d",
        processed, skipped_empty, skipped_chapter,
    )
    return pd.DataFrame(summary_rows)


def main():
    ap = argparse.ArgumentParser(description="Parse A/B/D/E/F ICD chapters from RAP extract")
    ap.add_argument(
        "--input", type=Path,
        default=RAP_DIR / "icd_abdef_first_occurrence.csv",
        help="Path to icd_abdef_first_occurrence.csv  [default: data/ukb/rap_extraction/]",
    )
    ap.add_argument(
        "--mapping", type=Path,
        default=RAP_DIR / "abdef_field_mapping.json",
        help="Path to abdef_field_mapping.json",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=ICD_DIR,
        help="Root icd_codes/ directory",
    )
    ap.add_argument(
        "--chapters", type=str, default="ABDEF",
        help="Which chapters to parse  [default: ABDEF]",
    )
    args = ap.parse_args()

    if not args.input.exists():
        log.error("Input CSV not found: %s", args.input)
        log.error("Download from RAP first:")
        log.error('  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_abdef_first_occurrence.csv"')
        log.error("  Place in: %s", RAP_DIR)
        sys.exit(1)

    field_map = load_field_mapping(args.mapping)
    summary = parse_chapters(
        args.input, field_map, args.out_dir, chapters=set(args.chapters.upper())
    )

    summary_path = args.out_dir.parent / "parse_abdef_summary.csv"
    summary.to_csv(summary_path, index=False)
    log.info("Summary: %s (%d rows)", summary_path, len(summary))

    # Per-chapter breakdown
    if len(summary) > 0:
        log.info("\nPer-chapter summary:")
        for ch, grp in summary.groupby("chapter"):
            log.info(
                "  Chapter %s: %d codes, total events=%d, "
                "top: %s (n=%d)",
                ch, len(grp), grp["n_events"].sum(),
                grp.loc[grp["n_events"].idxmax(), "icd_code"],
                grp["n_events"].max(),
            )

        log.info("\nTop 20 codes by n_events:")
        top = summary.sort_values("n_events", ascending=False).head(20)
        log.info("\n%s", top[["icd_code","description","n_events","median_age"]].to_string(index=False))


if __name__ == "__main__":
    main()
