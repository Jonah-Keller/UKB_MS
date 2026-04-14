#!/usr/bin/env python3
"""Extract MS/ALS-relevant fields from UKB RAP using dx CLI table-exporter.

This script runs LOCALLY (not on RAP JupyterLab). It uses the dx CLI
table-exporter app to extract specific fields from the dispensed UKB dataset
and download the results.

Extracts:
  1. Genetic PCs 1-40 (field 22009) — population stratification covariates
  2. HLA imputed alleles (field 22182) — DRB1*15:01 MS risk allele
  3. Death cause codes (fields 40001, 40002) — ALS death confirmation
  4. Genetic QC fields (22001, 22006, 22019, 22027) — sample QC

Prerequisites:
  - dxpy installed: pip install dxpy
  - Authenticated: dx login
  - Project selected: dx select project-GxX43xjJpjp2G7XK6ffz0qb1

Pattern from: CADASIL repo analysis/ukb/00_data_preparation/

Usage:
    python 02_extract_via_dx_cli.py --dry-run    # print commands only
    python 02_extract_via_dx_cli.py              # run extractions
    python 02_extract_via_dx_cli.py --download   # download completed results

Author: Jonah Keller (AI-assisted)
Date: 2026-04-14
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# --- Configuration ---
PROJECT_ID = "project-GxX43xjJpjp2G7XK6ffz0qb1"
DATASET_ID = "record-GxX62v8J88vBkPBk3pXXv0J8"
RESULTS_FOLDER = "/ms_als_extraction/"

# dx CLI path (use the one from the CADASIL venv)
DX = str(Path(__file__).resolve().parents[2].parent
         / "CADASIL_Proteome_ML_Keller_2024_Rebuttal" / ".venv" / "bin" / "dx")

# Local output directory
REPO_ROOT = Path(__file__).resolve().parents[2]
LOCAL_GENETICS_DIR = REPO_ROOT / "data" / "ukb" / "genetics"
LOCAL_MISC_DIR = REPO_ROOT / "data" / "ukb" / "misc"

# --- Field definitions ---
# Each extraction is a dict with: name, table, fields, output_name
EXTRACTIONS = [
    {
        "name": "Genetic PCs (field 22009)",
        "table": "participant",
        "fields": ["eid"] + [f"p22009_a{i}" for i in range(1, 41)],
        "output_name": "ukb_genetic_pcs_40",
        "local_dir": LOCAL_GENETICS_DIR,
    },
    {
        "name": "HLA Imputed Alleles (field 22182)",
        "table": "participant",
        # HLA field has many array indices — extract broadly
        "fields": ["eid"] + [f"p22182_i{i}" for i in range(0, 362)],
        "output_name": "ukb_hla_imputed",
        "local_dir": LOCAL_GENETICS_DIR,
    },
    {
        "name": "Death Causes (fields 40001, 40002)",
        "table": "participant",
        "fields": (
            ["eid"]
            + [f"p40001_i{i}" for i in range(2)]
            + [f"p40002_i{i}" for i in range(15)]
        ),
        "output_name": "ukb_death_causes",
        "local_dir": LOCAL_MISC_DIR,
    },
    {
        "name": "Genetic QC fields",
        "table": "participant",
        "fields": [
            "eid", "p22001", "p22006", "p22019", "p22027",
            "p21022", "p31",
            "p53_i0", "p53_i1", "p53_i2", "p53_i3",
        ],
        "output_name": "ukb_genetic_qc",
        "local_dir": LOCAL_MISC_DIR,
    },
]


def run_dx(cmd: list[str], capture: bool = True) -> Optional[str]:
    """Run a dx CLI command."""
    try:
        result = subprocess.run(cmd, capture_output=capture, text=True, check=True)
        return result.stdout.strip() if capture else None
    except subprocess.CalledProcessError as e:
        log.error("Command failed: %s", " ".join(cmd))
        if e.stderr:
            log.error("stderr: %s", e.stderr.strip()[:500])
        return None
    except FileNotFoundError:
        log.error("dx CLI not found at: %s", cmd[0])
        sys.exit(1)


def submit_extraction(extraction: dict, dry_run: bool = False) -> Optional[str]:
    """Submit a table-exporter job for the given field set."""
    name = extraction["name"]
    table = extraction["table"]
    fields = extraction["fields"]
    output_name = extraction["output_name"]

    log.info("Submitting: %s (%d fields)", name, len(fields))

    # Build the field list as comma-separated for dx
    field_str = ",".join(fields)

    # Use the swiss-army-knife approach: run a SQL query via the dataset
    # Actually, table-exporter is better for this use case
    cmd = [
        DX, "run", "app-table-exporter",
        f"-idataset_or_cohort_or_dashboard={DATASET_ID}",
        f"-ioutput={output_name}",
        f"-ioutput_format=CSV",
        f"-iheader_style=FIELD-NAME",
        f"-ifield_names={field_str}",
        f"-ientity=participant",
        "--destination", RESULTS_FOLDER,
        "--name", f"extract_{output_name}",
        "--brief",
        "-y",
    ]

    if dry_run:
        log.info("[DRY RUN] %s", " ".join(cmd))
        return None

    job_id = run_dx(cmd)
    if job_id:
        log.info("  -> Job ID: %s", job_id)
    else:
        log.error("  -> Failed to submit %s", name)

        # Fallback: try using the record ID format
        log.info("  Trying alternative approach...")
        cmd_alt = [
            DX, "run", "app-table-exporter",
            f"-idataset_or_cohort_or_dashboard=record-{DATASET_ID.split('_')[0]}",
            f"-ioutput={output_name}",
            f"-ioutput_format=CSV",
            f"-ifield_names={field_str}",
            "--destination", RESULTS_FOLDER,
            "--name", f"extract_{output_name}",
            "--brief",
            "-y",
        ]
        job_id = run_dx(cmd_alt)
        if job_id:
            log.info("  -> Job ID (alt): %s", job_id)

    return job_id


def check_status() -> None:
    """Check status of running extraction jobs."""
    log.info("Checking running jobs...")
    output = run_dx([
        DX, "find", "jobs",
        "--project", PROJECT_ID,
        "--state", "running",
        "--name", "extract_*",
        "--brief",
    ])
    if not output:
        log.info("No running extraction jobs.")
        # Also check recently completed
        output = run_dx([
            DX, "find", "jobs",
            "--project", PROJECT_ID,
            "--state", "done",
            "--name", "extract_*",
            "--brief",
            "-n", "10",
        ])
        if output:
            log.info("Recently completed jobs:")
            for line in output.splitlines():
                desc = run_dx([DX, "describe", line.strip(), "--json"])
                if desc:
                    info = json.loads(desc)
                    log.info("  %s: %s (%s)", info.get("name"), info.get("state"), line.strip())
        return

    for line in output.splitlines():
        job_id = line.strip()
        desc = run_dx([DX, "describe", job_id, "--json"])
        if desc:
            info = json.loads(desc)
            log.info("  %s: %s", info.get("name"), info.get("state"))


def download_results() -> None:
    """Download completed extraction results."""
    log.info("Downloading results from %s...", RESULTS_FOLDER)

    # Ensure output dirs exist
    LOCAL_GENETICS_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_MISC_DIR.mkdir(parents=True, exist_ok=True)

    # List files in results folder
    output = run_dx([DX, "ls", RESULTS_FOLDER])
    if not output:
        log.error("No files found in %s", RESULTS_FOLDER)
        return

    for line in output.splitlines():
        fname = line.strip()
        if not fname or fname.endswith("/"):
            continue

        # Determine local destination
        if "genetic_pc" in fname.lower() or "hla" in fname.lower():
            local_dir = LOCAL_GENETICS_DIR
        else:
            local_dir = LOCAL_MISC_DIR

        local_path = local_dir / fname
        if local_path.exists():
            log.info("  Skipping (exists): %s", fname)
            continue

        log.info("  Downloading: %s -> %s", fname, local_dir)
        run_dx([
            DX, "download",
            f"{RESULTS_FOLDER}{fname}",
            "-o", str(local_path),
        ], capture=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract UKB fields via dx CLI")
    parser.add_argument("--dry-run", "-n", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--status", "-s", action="store_true",
                        help="Check job status")
    parser.add_argument("--download", "-d", action="store_true",
                        help="Download completed results")
    parser.add_argument("--extract", "-e", type=str, default=None,
                        help="Extract specific job by name (e.g. 'genetic_pcs')")
    args = parser.parse_args()

    if args.status:
        check_status()
        return

    if args.download:
        download_results()
        return

    # Ensure results folder exists
    run_dx([DX, "mkdir", "-p", RESULTS_FOLDER])

    # Submit extractions
    submitted = []
    for ext in EXTRACTIONS:
        if args.extract and args.extract not in ext["output_name"]:
            continue
        job_id = submit_extraction(ext, dry_run=args.dry_run)
        if job_id:
            submitted.append((ext["name"], job_id))

    if args.dry_run:
        log.info("Dry run complete. %d jobs would be submitted.", len(EXTRACTIONS))
    else:
        log.info("=== Submitted %d jobs ===", len(submitted))
        for name, jid in submitted:
            log.info("  %s -> %s", name, jid)
        log.info("")
        log.info("Monitor with: python %s --status", Path(__file__).name)
        log.info("Download with: python %s --download", Path(__file__).name)


if __name__ == "__main__":
    main()
