"""Extract MS/ALS-relevant fields from UKB RAP via dxdata.

Run this script inside a UKB RAP JupyterLab session (Python 3, ≥16 GB RAM).
It extracts:
  1. Genetic PCs (field 22009, 40 PCs) — for population stratification control
  2. HLA-DRB1*15:01 imputed dosage (field 22182) — MS immune cluster proxy
  3. Death cause codes (fields 40001, 40002) — ALS death confirmation
  4. Assessment dates (field 53) — for precise time-to-diagnosis calculation

Outputs are written to the JupyterLab working directory as TSV files.
Download them via `dx download` and place in data/ukb/genetics/ and data/ukb/misc/.

Pattern mirrors: CADASIL repo analysis/ukb/00_data_preparation/01_prepare_pc2_phenotype.py

Usage (RAP JupyterLab terminal):
    python 01_extract_rap_fields.py

Author: Jonah Keller (AI-assisted)
Date:   2026-04-14
"""

from __future__ import annotations

import logging
import sys
from typing import Dict, List

import pandas as pd

# ---------------------------------------------------------------------------
# RAP-specific imports — only available inside a RAP JupyterLab session
# ---------------------------------------------------------------------------
try:
    import dxdata
    import dxpy
except ImportError:
    print(
        "ERROR: dxdata / dxpy not found. This script must be run inside a "
        "UKB RAP JupyterLab session.",
        file=sys.stderr,
    )
    sys.exit(1)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════

# RAP project ID (same as CADASIL repo)
PROJECT_ID = "project-GxX43xjJpjp2G7XK6ffz0qb1"

# Number of genetic PCs to extract (field 22009)
N_GENETIC_PCS = 40

# HLA field 22182: Imputed HLA allele dosages
# DRB1*15:01 is the primary MS risk allele
HLA_FIELD_ID = 22182

# Death cause fields
DEATH_FIELDS: Dict[str, str] = {
    "p40001_i0": "primary_death_cause_i0",
    "p40001_i1": "primary_death_cause_i1",
    "p40002_i0": "contributory_death_cause_i0",
    "p40002_i1": "contributory_death_cause_i1",
    "p40002_i2": "contributory_death_cause_i2",
    "p40002_i3": "contributory_death_cause_i3",
    "p40002_i4": "contributory_death_cause_i4",
    "p40002_i5": "contributory_death_cause_i5",
    "p40002_i6": "contributory_death_cause_i6",
    "p40002_i7": "contributory_death_cause_i7",
    "p40002_i8": "contributory_death_cause_i8",
    "p40002_i9": "contributory_death_cause_i9",
    "p40002_i10": "contributory_death_cause_i10",
    "p40002_i11": "contributory_death_cause_i11",
    "p40002_i12": "contributory_death_cause_i12",
    "p40002_i13": "contributory_death_cause_i13",
}

# Assessment dates (field 53) and genetic QC fields
PARTICIPANT_CORE: Dict[str, str] = {
    "eid": "eid",
    "p21022": "age_at_recruitment",
    "p31": "sex",
    "p22001": "genetic_sex",
    "p22006": "genetic_ethnic_grouping",
    "p22019": "sex_chromosome_aneuploidy",
    "p22027": "het_missing_outliers",
    "p53_i0": "assessment_date_i0",
    "p53_i1": "assessment_date_i1",
    "p53_i2": "assessment_date_i2",
    "p53_i3": "assessment_date_i3",
}

# Output file names
GENETIC_PCS_OUT = "ukb_genetic_pcs_40.tsv"
HLA_OUT = "ukb_hla_drb1_1501.tsv"
DEATH_CAUSES_OUT = "ukb_death_causes.tsv"
CORE_FIELDS_OUT = "ukb_core_participant_fields.tsv"


# ═══════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════

def load_dispensed_dataset() -> "dxdata.Dataset":
    """Load the dispensed UKB dataset on the RAP."""
    dispensed_id = dxpy.find_one_data_object(
        typename="Dataset",
        name="app*",
        folder="/",
        name_mode="glob",
        describe=True,
    )
    dataset_id = dispensed_id["describe"]["id"]
    log.info("Found dispensed dataset: %s", dataset_id)
    return dxdata.load_dataset(id=dataset_id)


def safe_retrieve(participant, field_names: List[str], engine) -> pd.DataFrame:
    """Retrieve fields, skipping any that don't exist in the dataset."""
    fields = []
    missing = []
    for fname in field_names:
        try:
            fields.append(participant[fname])
        except (LookupError, KeyError):
            missing.append(fname)

    if missing:
        log.warning("%d fields not found: %s", len(missing), missing[:10])

    if not fields:
        raise RuntimeError("No fields could be retrieved!")

    spark_df = participant.retrieve_fields(fields=fields, engine=engine)
    df = spark_df.toPandas()
    # Strip table prefix from column names
    df.columns = [c.split(".")[-1] if "." in c else c for c in df.columns]
    return df


# ═══════════════════════════════════════════════════════════════════════════
# EXTRACTION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def extract_genetic_pcs(dataset, engine, n_pcs: int = N_GENETIC_PCS) -> None:
    """Extract genetic principal components (field 22009)."""
    log.info("=== Extracting %d Genetic PCs (field 22009) ===", n_pcs)
    participant = dataset["participant"]

    # Field 22009 array indices are 1-indexed on RAP: p22009_a1, p22009_a2, ...
    field_names = ["eid"] + [f"p22009_a{i}" for i in range(1, n_pcs + 1)]

    df = safe_retrieve(participant, field_names, engine)

    # Rename to PC1..PCn
    rename = {f"p22009_a{i}": f"PC{i}" for i in range(1, n_pcs + 1)}
    df = df.rename(columns=rename)

    # Drop rows with all-NaN PCs (no genetic data)
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1) if f"PC{i}" in df.columns]
    n_before = len(df)
    df = df.dropna(subset=pc_cols, how="all")
    log.info("Dropped %d rows with no genetic data; %d remain", n_before - len(df), len(df))

    df.to_csv(GENETIC_PCS_OUT, sep="\t", index=False)
    log.info("Wrote %s (%d rows, %d cols)", GENETIC_PCS_OUT, df.shape[0], df.shape[1])


def extract_hla(dataset, engine) -> None:
    """Extract HLA-DRB1*15:01 dosage (field 22182).

    UKB field 22182 stores imputed classical HLA alleles.
    The array indices correspond to specific alleles — we need to identify
    which index maps to DRB1*15:01.

    Strategy: extract all 22182 array values and identify DRB1*15:01 post-hoc,
    or use the known allele coding from UKB documentation.
    """
    log.info("=== Extracting HLA Imputation (field 22182) ===")
    participant = dataset["participant"]

    # Field 22182 has multiple array indices for different HLA alleles
    # We extract all available indices. Typical: p22182_i0 through p22182_iN
    # The exact number depends on the dataset version
    hla_fields = ["eid"]
    for i in range(0, 400):  # Conservative upper bound
        hla_fields.append(f"p22182_i{i}")

    # Also try the array notation (p22182_a{i})
    for i in range(0, 400):
        hla_fields.append(f"p22182_a{i}")

    df = safe_retrieve(participant, hla_fields, engine)

    # The result will only contain columns that actually exist
    log.info("Retrieved %d HLA columns", len(df.columns) - 1)

    # Drop rows with all-NaN HLA (no imputed data)
    hla_cols = [c for c in df.columns if c != "eid"]
    if hla_cols:
        n_before = len(df)
        df = df.dropna(subset=hla_cols, how="all")
        log.info("Dropped %d rows with no HLA data; %d remain", n_before - len(df), len(df))

    df.to_csv(HLA_OUT, sep="\t", index=False)
    log.info("Wrote %s (%d rows, %d cols)", HLA_OUT, df.shape[0], df.shape[1])
    log.info("NOTE: Post-hoc parsing needed to identify DRB1*15:01 from allele codes.")
    log.info("See: https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=22182")


def extract_death_causes(dataset, engine) -> None:
    """Extract death cause ICD-10 codes (fields 40001, 40002)."""
    log.info("=== Extracting Death Cause Codes (fields 40001, 40002) ===")
    participant = dataset["participant"]

    field_names = ["eid"] + list(DEATH_FIELDS.keys())
    df = safe_retrieve(participant, field_names, engine)

    rename = {k: v for k, v in DEATH_FIELDS.items() if k in df.columns}
    df = df.rename(columns=rename)

    # Drop rows where ALL death cause columns are NaN (alive)
    death_cols = [c for c in df.columns if c != "eid"]
    if death_cols:
        n_before = len(df)
        df = df.dropna(subset=death_cols, how="all")
        log.info("Dropped %d alive participants; %d with death data", n_before - len(df), len(df))

    df.to_csv(DEATH_CAUSES_OUT, sep="\t", index=False)
    log.info("Wrote %s (%d rows, %d cols)", DEATH_CAUSES_OUT, df.shape[0], df.shape[1])


def extract_core_fields(dataset, engine) -> None:
    """Extract core participant fields (demographics, genetic QC, dates)."""
    log.info("=== Extracting Core Participant Fields ===")
    participant = dataset["participant"]

    field_names = list(PARTICIPANT_CORE.keys())
    df = safe_retrieve(participant, field_names, engine)

    rename = {k: v for k, v in PARTICIPANT_CORE.items() if k in df.columns and k != "eid"}
    df = df.rename(columns=rename)

    df.to_csv(CORE_FIELDS_OUT, sep="\t", index=False)
    log.info("Wrote %s (%d rows, %d cols)", CORE_FIELDS_OUT, df.shape[0], df.shape[1])


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    log.info("=" * 60)
    log.info("UKB RAP Field Extraction for MS/ALS Proteomics")
    log.info("Project: %s", PROJECT_ID)
    log.info("=" * 60)

    # 1. Load dispensed dataset and create engine
    log.info("[STEP 1] Loading dispensed dataset...")
    dataset = load_dispensed_dataset()
    engine = dxdata.connect()
    log.info("[STEP 1] Dataset and Spark engine ready.")

    # 2. Extract genetic PCs
    log.info("[STEP 2] Genetic PCs...")
    extract_genetic_pcs(dataset, engine)

    # 3. Extract HLA imputation
    log.info("[STEP 3] HLA imputation...")
    extract_hla(dataset, engine)

    # 4. Extract death causes
    log.info("[STEP 4] Death causes...")
    extract_death_causes(dataset, engine)

    # 5. Extract core participant fields
    log.info("[STEP 5] Core participant fields...")
    extract_core_fields(dataset, engine)

    log.info("")
    log.info("=" * 60)
    log.info("EXTRACTION COMPLETE")
    log.info("Output files:")
    log.info("  %s  (genetic PCs 1-40)", GENETIC_PCS_OUT)
    log.info("  %s  (HLA-DRB1*15:01 dosage)", HLA_OUT)
    log.info("  %s  (death cause ICD codes)", DEATH_CAUSES_OUT)
    log.info("  %s  (demographics, QC, dates)", CORE_FIELDS_OUT)
    log.info("")
    log.info("NEXT STEPS:")
    log.info("  1. dx download %s %s %s %s", GENETIC_PCS_OUT, HLA_OUT, DEATH_CAUSES_OUT, CORE_FIELDS_OUT)
    log.info("  2. Place genetic PCs in:  data/ukb/genetics/")
    log.info("  3. Place HLA in:          data/ukb/genetics/")
    log.info("  4. Place death causes in: data/ukb/misc/")
    log.info("  5. Place core fields in:  data/ukb/misc/")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
