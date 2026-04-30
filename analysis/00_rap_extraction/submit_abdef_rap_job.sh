#!/usr/bin/env bash
# submit_abdef_rap_job.sh
#
# Submits the ICD chapter A/B/D/E/F extraction to UKB RAP via dx CLI.
#
# Prerequisites:
#   dx login   (if not already authenticated)
#   dx select  (select the UKB project)
#
# Usage:
#   bash analysis/00_rap_extraction/submit_abdef_rap_job.sh
#
# The script:
#   1. Uploads 07_rap_extract_abdef_chapters.py to the RAP project
#   2. Launches a JupyterLab instance job with sufficient RAM
#   3. Prints the job ID and a command to monitor progress
#   4. After job completes, run the download block at the bottom
#
# Estimated runtime: 15-25 minutes (341 fields, 500k participants)
# Recommended instance: mem1_ssd1_v2_x16  (16 CPUs, 120 GB RAM)

set -euo pipefail

SCRIPT_PATH="analysis/00_rap_extraction/07_rap_extract_abdef_chapters.py"
RAP_PROJECT="project-J2P0fqjJ1Gg4gByJ6BJ51pYj"
RAP_SCRIPT_DEST="/code/07_rap_extract_abdef_chapters.py"
RAP_OUTPUT_DIR="/data/ms_als_extraction"
INSTANCE_TYPE="mem1_ssd1_v2_x16"   # 16 vCPUs, 120 GB RAM — good for Spark

echo "========================================================"
echo " UKB RAP: ICD-10 A/B/D/E/F Chapter Extraction"
echo "========================================================"
echo ""

# ── 1. Check dx login ────────────────────────────────────────────────────────
echo "[1/4] Checking dx login..."
if ! dx whoami &>/dev/null; then
    echo "ERROR: Not logged into dx. Run:  dx login"
    exit 1
fi
DX_USER=$(dx whoami --id)
echo "  Logged in as: ${DX_USER}"

# ── 2. Upload extraction script ──────────────────────────────────────────────
echo ""
echo "[2/4] Uploading extraction script to RAP..."
dx mkdir -p "${RAP_OUTPUT_DIR}" 2>/dev/null || true
dx upload "${SCRIPT_PATH}" --destination "${RAP_PROJECT}:${RAP_SCRIPT_DEST}" --brief
echo "  Uploaded: ${RAP_SCRIPT_DEST}"

# ── 3. Launch JupyterLab job (non-interactive, runs script directly) ─────────
echo ""
echo "[3/4] Launching extraction job on ${INSTANCE_TYPE}..."
JOB_ID=$(dx run dxjupyterlab \
    --instance-type "${INSTANCE_TYPE}" \
    --project "${RAP_PROJECT}" \
    --destination "${RAP_OUTPUT_DIR}" \
    -iscript="${RAP_PROJECT}:${RAP_SCRIPT_DEST}" \
    --name "icd_abdef_extraction" \
    --brief \
    --yes)

echo "  Job submitted: ${JOB_ID}"

# ── 4. Monitor ───────────────────────────────────────────────────────────────
echo ""
echo "[4/4] Monitor progress:"
echo "  dx watch ${JOB_ID}"
echo "  dx describe ${JOB_ID}"
echo ""
echo "========================================================"
echo " After job completes, download results:"
echo "========================================================"
cat << 'DOWNLOAD_BLOCK'

mkdir -p data/ukb/rap_extraction

dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_abdef_first_occurrence.csv" \
    --output data/ukb/rap_extraction/icd_abdef_first_occurrence.csv

dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/abdef_field_mapping.json" \
    --output data/ukb/rap_extraction/abdef_field_mapping.json

# Parse into per-code icd_codes/ files
python analysis/00_rap_extraction/08_parse_abdef_chapters.py

# Verify new chapters appeared
ls /Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/data/ukb/diagnoses/icd_codes/

DOWNLOAD_BLOCK

echo ""
echo "Job ID: ${JOB_ID}"
echo "Done."
