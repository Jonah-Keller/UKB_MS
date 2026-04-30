#!/usr/bin/env python3
"""Parse E/F ICD-10 first-occurrence chapters into per-code icd_codes files.

Run after downloading icd_ef_first_occurrence.csv from RAP:
  python analysis/00_rap_extraction/06_parse_ef_chapters.py

Input:
  data/ukb/rap_extraction/icd_ef_first_occurrence.csv  (or path from --input)

Output:
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/E/...
  CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/F/...
  Same eid, diagnosis, age_at_diagnosis format as existing G-Q files.

Usage:
  python 06_parse_ef_chapters.py
  python 06_parse_ef_chapters.py --input /path/to/icd_ef_first_occurrence.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

PROJ_ROOT   = Path(__file__).resolve().parents[2]
CADASIL_DIR = PROJ_ROOT.parent / "CADASIL_Proteome_ML_Keller_2024_Rebuttal"
ICD_DIR     = CADASIL_DIR / "data" / "ukb" / "diagnoses" / "icd_codes"
RAP_DIR     = PROJ_ROOT / "data" / "ukb" / "rap_extraction"

# Embedded field map — used when ef_field_mapping.json is not present
EF_FIELD_MAP: dict[str, tuple[str, str]] = {
    "p130692": ("E01", "iodine-deficiency-related thyroid disorders"),
    "p130694": ("E02", "subclinical iodine-deficiency hypothyroidism"),
    "p130696": ("E03", "other hypothyroidism"),
    "p130698": ("E04", "other non-toxic goitre"),
    "p130700": ("E05", "thyrotoxicosis [hyperthyroidism]"),
    "p130702": ("E06", "thyroiditis"),
    "p130704": ("E07", "other disorders of thyroid"),
    "p130706": ("E10", "insulin-dependent diabetes mellitus"),
    "p130708": ("E11", "non-insulin-dependent diabetes mellitus"),
    "p130712": ("E13", "other specified diabetes mellitus"),
    "p130714": ("E14", "unspecified diabetes mellitus"),
    "p130720": ("E20", "hypoparathyroidism"),
    "p130722": ("E21", "hyperparathyroidism"),
    "p130724": ("E22", "hyperfunction of pituitary gland"),
    "p130726": ("E23", "hypofunction and other disorders of pituitary gland"),
    "p130728": ("E24", "cushings syndrome"),
    "p130732": ("E26", "hyperaldosteronism"),
    "p130734": ("E27", "other disorders of adrenal gland"),
    "p130736": ("E28", "ovarian dysfunction"),
    "p130742": ("E31", "polyglandular dysfunction"),
    "p130746": ("E34", "other endocrine disorders"),
    "p130774": ("E55", "vitamin d deficiency"),
    "p130792": ("E66", "obesity"),
    "p130814": ("E78", "disorders of lipoprotein metabolism and other lipidaemias"),
    "p130822": ("E84", "cystic fibrosis"),
    "p130824": ("E85", "amyloidosis"),
    "p130828": ("E87", "other disorders of fluid electrolyte and acid-base balance"),
    "p130830": ("E88", "other metabolic disorders"),
    "p130836": ("F00", "dementia in alzheimers disease"),
    "p130838": ("F01", "vascular dementia"),
    "p130840": ("F02", "dementia in other diseases classified elsewhere"),
    "p130842": ("F03", "unspecified dementia"),
    "p130846": ("F05", "delirium not induced by alcohol"),
    "p130854": ("F10", "mental and behavioural disorders due to use of alcohol"),
    "p130868": ("F17", "mental and behavioural disorders due to use of tobacco"),
    "p130874": ("F20", "schizophrenia"),
    "p130884": ("F25", "schizoaffective disorders"),
    "p130890": ("F30", "manic episode"),
    "p130892": ("F31", "bipolar affective disorder"),
    "p130894": ("F32", "depressive episode"),
    "p130896": ("F33", "recurrent depressive disorder"),
    "p130898": ("F34", "persistent mood affective disorders"),
    "p130904": ("F40", "phobic anxiety disorders"),
    "p130906": ("F41", "other anxiety disorders"),
    "p130908": ("F42", "obsessive-compulsive disorder"),
    "p130910": ("F43", "reaction to severe stress and adjustment disorders"),
    "p130912": ("F44", "dissociative conversion disorders"),
    "p130914": ("F45", "somatoform disorders"),
    "p130918": ("F50", "eating disorders"),
    "p130932": ("F60", "specific personality disorders"),
    "p130970": ("F84", "pervasive developmental disorders"),
    "p130976": ("F90", "hyperkinetic disorders"),
}


def load_field_map(mapping_json: Optional[Path]) -> dict[str, tuple[str, str]]:
    if mapping_json and mapping_json.exists():
        with open(mapping_json) as f:
            raw = json.load(f)
        # Keys may be ints or strings; normalise to "p<id>"
        result = {}
        for k, v in raw.items():
            key = f"p{k}" if not str(k).startswith("p") else str(k)
            result[key] = (v[0], v[1])
        log.info("Loaded %d fields from %s", len(result), mapping_json)
        return result
    log.warning("No ef_field_mapping.json found, using embedded subset (%d fields)", len(EF_FIELD_MAP))
    return EF_FIELD_MAP


def icd_output_path(icd_code: str, base: Path) -> Path:
    clean = icd_code.replace(".", "")
    l1 = clean[0]
    l2 = clean[:2] if len(clean) >= 2 else l1
    l3 = clean[:3] if len(clean) >= 3 else l2
    if len(clean) <= 3:
        return base / l1 / l2 / l3 / f"{icd_code}.csv"
    return base / l1 / l2 / l3 / icd_code / f"{icd_code}.csv"


def parse_ef(csv_path: Path, field_map: dict, out_dir: Path) -> pd.DataFrame:
    log.info("Loading %s ...", csv_path.name)
    df = pd.read_csv(csv_path, low_memory=False)
    log.info("  %d rows x %d cols", *df.shape)

    # Resolve year_birth column
    yob_col = next((c for c in df.columns if "year_birth" in c or c == "p34"), None)
    if yob_col is None:
        raise ValueError("year_birth or p34 column not found in CSV")

    eid_col = next((c for c in df.columns if c in ("eid", "participant_id")), None)
    if eid_col is None:
        raise ValueError("eid or participant_id column not found")

    summary_rows = []
    for col, (icd_code, desc) in field_map.items():
        if col not in df.columns:
            log.debug("Column %s not in CSV, skipping %s", col, icd_code)
            continue

        sub = df[[eid_col, yob_col, col]].dropna(subset=[col]).copy()
        if len(sub) == 0:
            continue

        sub["date_parsed"] = pd.to_datetime(sub[col], errors="coerce")
        sub = sub.dropna(subset=["date_parsed"])
        if len(sub) == 0:
            continue

        sub["age_at_diagnosis"] = sub["date_parsed"].dt.year - sub[yob_col].astype(int)
        sub.loc[sub["age_at_diagnosis"] < 0, "age_at_diagnosis"] = np.nan

        out = pd.DataFrame({
            "eid": sub[eid_col].values,
            "diagnosis": 1,
            "age_at_diagnosis": sub["age_at_diagnosis"].values,
        })

        out_path = icd_output_path(icd_code, out_dir)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_path, index=False)

        log.info("  %s (%s): %d events -> %s", icd_code, desc[:40], len(out), out_path)
        summary_rows.append({
            "icd_code": icd_code,
            "description": desc,
            "n_events": len(out),
            "median_age": out["age_at_diagnosis"].median(),
            "output_path": str(out_path),
        })

    return pd.DataFrame(summary_rows)


def main():
    parser = argparse.ArgumentParser(description="Parse E/F ICD chapters from RAP extract")
    parser.add_argument("--input", type=Path,
                        default=RAP_DIR / "icd_ef_first_occurrence.csv",
                        help="Path to icd_ef_first_occurrence.csv")
    parser.add_argument("--mapping", type=Path,
                        default=RAP_DIR / "ef_field_mapping.json",
                        help="Optional: ef_field_mapping.json from RAP extraction")
    parser.add_argument("--out-dir", type=Path, default=ICD_DIR,
                        help="icd_codes output directory")
    args = parser.parse_args()

    if not args.input.exists():
        log.error("Input file not found: %s", args.input)
        log.error("Download from RAP first:")
        log.error('  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_ef_first_occurrence.csv"')
        raise SystemExit(1)

    field_map = load_field_map(args.mapping)
    summary = parse_ef(args.input, field_map, args.out_dir)

    summary_path = args.out_dir.parent / "parse_ef_chapters_summary.csv"
    summary.to_csv(summary_path, index=False)
    log.info("\nParsed %d E/F ICD codes", len(summary))
    log.info("Summary: %s", summary_path)

    if len(summary) > 0:
        log.info("\nTop codes by n_events:")
        log.info(summary.sort_values("n_events", ascending=False).head(10).to_string(index=False))


if __name__ == "__main__":
    main()
