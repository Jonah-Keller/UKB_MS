#!/usr/bin/env bash
# 04_submit_gwas_rap.sh
# Submit the cluster GWAS pipeline to UKB RAP.
#
# Usage:
#   bash analysis/16_cluster_gwas/04_submit_gwas_rap.sh        # step 1 only
#   bash analysis/16_cluster_gwas/04_submit_gwas_rap.sh --step2 STEP1_JOB_ID
#
# Two-phase submission:
#   Phase 1: Run this script without args → submits REGENIE step 1 (~3-4 hours)
#   Phase 2: Run with --step2 <STEP1_JOB_ID> → submits 22 step-2 chr jobs,
#            each depending on step 1 completing.
#
# Prerequisites:
#   Rscript analysis/16_cluster_gwas/01_build_phenotype.R  (generates pheno/covar files)
#   dx login && dx select project-J2P0fqjJ1Gg4gByJ6BJ51pYj

set -euo pipefail

PROJ="project-J2P0fqjJ1Gg4gByJ6BJ51pYj"
RAP_DIR="/data/ms_gwas_cluster"
PHENO_DIR="data/ukb/genetics/gwas_cluster"

# ── Argument parsing ──────────────────────────────────────────────────────────
STEP2=false
STEP1_JOB=""
if [[ "${1:-}" == "--step2" ]]; then
    STEP2=true
    STEP1_JOB="${2:?'Usage: --step2 <STEP1_JOB_ID>'}"
fi

# ── Check dx login ────────────────────────────────────────────────────────────
if ! dx whoami &>/dev/null; then
    echo "ERROR: Not logged into dx. Run: dx login"
    exit 1
fi
echo "Logged in as: $(dx whoami)"

# ── Phase 1: Upload files + submit REGENIE step 1 ────────────────────────────
if [[ "${STEP2}" == "false" ]]; then
    echo ""
    echo "======================================================"
    echo " Phase 1: Upload files and run REGENIE step 1"
    echo "======================================================"

    # Ensure phenotype files exist
    for f in pheno_C1_vs_none.tsv pheno_C2_vs_none.tsv covariates.tsv; do
        if [[ ! -f "${PHENO_DIR}/${f}" ]]; then
            echo "ERROR: Missing ${PHENO_DIR}/${f} — run 01_build_phenotype.R first"
            exit 1
        fi
    done

    echo "[1/3] Creating RAP output directory..."
    dx mkdir -p "${PROJ}:${RAP_DIR}/step1_output" 2>/dev/null || true

    echo "[2/3] Uploading phenotype, covariate, and step1 script..."
    # Capture file IDs directly from upload output
    PHENO_C1_ID=$(dx upload "${PHENO_DIR}/pheno_C1_vs_none.tsv" \
        --destination "${PROJ}:${RAP_DIR}/pheno_C1_vs_none.tsv" --brief 2>/dev/null)
    echo "  Uploaded: pheno_C1_vs_none.tsv (${PHENO_C1_ID})"

    PHENO_C2_ID=$(dx upload "${PHENO_DIR}/pheno_C2_vs_none.tsv" \
        --destination "${PROJ}:${RAP_DIR}/pheno_C2_vs_none.tsv" --brief 2>/dev/null)
    echo "  Uploaded: pheno_C2_vs_none.tsv (${PHENO_C2_ID})"

    COVAR_ID=$(dx upload "${PHENO_DIR}/covariates.tsv" \
        --destination "${PROJ}:${RAP_DIR}/covariates.tsv" --brief 2>/dev/null)
    echo "  Uploaded: covariates.tsv (${COVAR_ID})"

    STEP1_SCRIPT_ID=$(dx upload "analysis/16_cluster_gwas/02_regenie_step1_rap.sh" \
        --destination "${PROJ}:${RAP_DIR}/02_regenie_step1_rap.sh" --brief 2>/dev/null)
    echo "  Uploaded: 02_regenie_step1_rap.sh (${STEP1_SCRIPT_ID})"

    echo "[3/3] Submitting REGENIE step 1 job..."
    # Use dx download with explicit file IDs — avoids dxfuse lazy-load timing issues
    ICMD="dx download ${PHENO_C1_ID} -o /home/dnanexus/pheno_C1_vs_none.tsv -f && \
dx download ${PHENO_C2_ID} -o /home/dnanexus/pheno_C2_vs_none.tsv -f && \
dx download ${COVAR_ID} -o /home/dnanexus/covariates.tsv -f && \
dx download ${STEP1_SCRIPT_ID} -o /home/dnanexus/02_regenie_step1_rap.sh -f && \
bash /home/dnanexus/02_regenie_step1_rap.sh"

    STEP1_JOB=$(dx run swiss-army-knife \
        -icmd="${ICMD}" \
        --instance-type mem2_ssd1_v2_x32 \
        --destination "${PROJ}:${RAP_DIR}/step1_output" \
        --name "ms_cluster_gwas_step1" \
        --brief \
        --yes)

    echo ""
    echo "======================================================"
    echo " Step 1 job submitted: ${STEP1_JOB}"
    echo "======================================================"
    echo ""
    echo "Monitor: dx watch ${STEP1_JOB}"
    echo ""
    echo "After step 1 completes (~3-4 hours), run:"
    echo "  bash analysis/16_cluster_gwas/04_submit_gwas_rap.sh --step2 ${STEP1_JOB}"
    echo ""
    echo "Step 1 job ID: ${STEP1_JOB}" | tee "${PHENO_DIR}/step1_job_id.txt"
    exit 0
fi

# ── Phase 2: Submit REGENIE step 2 per chromosome ────────────────────────────
echo ""
echo "======================================================"
echo " Phase 2: REGENIE step 2 (chromosomes 1-22)"
echo " Depends on: ${STEP1_JOB}"
echo "======================================================"

# Download step 1 pred lists from the step 1 output to get file IDs
echo "[1/3] Getting step 1 output file IDs..."
STEP1_OUTPUTS=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/step1_output" --brief)

# The step 1 output includes _pred.list files and _loco.gz files
# We need to pass these as inputs to step 2

# Upload step 2 script
dx upload "analysis/16_cluster_gwas/03_regenie_step2_rap.sh" \
    --destination "${PROJ}:${RAP_DIR}/03_regenie_step2_rap.sh" --brief
STEP2_SCRIPT_ID=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/03_regenie_step2_rap.sh" --brief)

PHENO_C1_ID=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/pheno_C1_vs_none.tsv" --brief)
PHENO_C2_ID=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/pheno_C2_vs_none.tsv" --brief)
COVAR_ID=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/covariates.tsv" --brief)

# Pred list files from step 1 output
PRED_C1_ID=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/step1_output/step1_C1_pred.list" --brief)
PRED_C2_ID=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/step1_output/step1_C2_pred.list" --brief)

# Loco files — one per chromosome (step1_C1_1.loco.gz etc.)
# Collect all loco file IDs
echo "[2/3] Collecting loco file IDs..."
LOCO_IDS_C1=()
LOCO_IDS_C2=()
for chr in $(seq 1 22); do
    c1_loco=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/step1_output/step1_C1_1.loco.gz" --brief 2>/dev/null || echo "")
    c2_loco=$(dx find data --project "${PROJ}" --path "${RAP_DIR}/step1_output/step1_C2_1.loco.gz" --brief 2>/dev/null || echo "")
done

echo "[3/3] Submitting step 2 jobs (one per chromosome)..."
CHR_JOBS=()
for chr in $(seq 1 22); do
    # Collect all loco files dynamically for this chromosome
    JOB_ID=$(dx run swiss-army-knife \
        -iin="${PHENO_C1_ID}" \
        -iin="${PHENO_C2_ID}" \
        -iin="${COVAR_ID}" \
        -iin="${PRED_C1_ID}" \
        -iin="${PRED_C2_ID}" \
        -iin="${STEP2_SCRIPT_ID}" \
        -icmd="for i in /home/dnanexus/in/in/*/; do cp \"\${i}\"* /home/dnanexus/ 2>/dev/null || true; done && export CHR=${chr} && bash /home/dnanexus/03_regenie_step2_rap.sh" \
        --instance-type mem1_ssd1_v2_x8 \
        --destination "${PROJ}:${RAP_DIR}/step2_output/chr${chr}" \
        --name "ms_cluster_gwas_step2_chr${chr}" \
        --depends-on "${STEP1_JOB}" \
        --brief \
        --yes)
    CHR_JOBS+=("${JOB_ID}")
    echo "  Submitted chr${chr}: ${JOB_ID}"
done

echo ""
echo "======================================================"
echo " All step 2 jobs submitted (will start after step 1)"
echo "======================================================"
printf '%s\n' "${CHR_JOBS[@]}" | tee "${PHENO_DIR}/step2_job_ids.txt"
echo ""
echo "After all jobs complete, run:"
echo "  bash analysis/16_cluster_gwas/05_collect_results.sh"
