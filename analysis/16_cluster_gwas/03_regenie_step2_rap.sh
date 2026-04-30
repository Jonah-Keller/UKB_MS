#!/usr/bin/env bash
# 03_regenie_step2_rap.sh
# REGENIE step 2: per-chromosome GWAS on imputed variants
#
# Called once per chromosome. CHR environment variable must be set.
# Runs on RAP inside Swiss Army Knife (mem1_ssd1_v2_x8).
# Step 1 loco files are copied in by dx run.
#
# Runtime estimate: 1-2 hours per chromosome

set -euo pipefail

CHR="${CHR:-1}"   # chromosome number (1-22), passed as env var
ln -sfn "/mnt/project/Bulk/Imputation/UKB imputation from genotype" /tmp/imputed
BGEN="/tmp/imputed/ukb22828_c${CHR}_b0_v3.bgen"
SAMPLE="/tmp/imputed/ukb22828_c${CHR}_b0_v3.sample"
OUT_DIR="/home/dnanexus/out/out"
mkdir -p "${OUT_DIR}"

echo "=== REGENIE step 2: chromosome ${CHR} ==="

# ── C1 GWAS ───────────────────────────────────────────────────────────────────
regenie \
    --step 2 \
    --bgen "${BGEN}" \
    --sample "${SAMPLE}" \
    --phenoFile /home/dnanexus/pheno_C1_vs_none.tsv \
    --covarFile /home/dnanexus/covariates.tsv \
    --pred /home/dnanexus/step1_C1_pred.list \
    --bt \
    --firth --approx \
    --pThresh 0.05 \
    --bsize 400 \
    --minINFO 0.3 \
    --af-cc \
    --threads 8 \
    --out "${OUT_DIR}/chr${CHR}_C1" \
    --chr "${CHR}"

# ── C2 GWAS ───────────────────────────────────────────────────────────────────
regenie \
    --step 2 \
    --bgen "${BGEN}" \
    --sample "${SAMPLE}" \
    --phenoFile /home/dnanexus/pheno_C2_vs_none.tsv \
    --covarFile /home/dnanexus/covariates.tsv \
    --pred /home/dnanexus/step1_C2_pred.list \
    --bt \
    --firth --approx \
    --pThresh 0.05 \
    --bsize 400 \
    --minINFO 0.3 \
    --af-cc \
    --threads 8 \
    --out "${OUT_DIR}/chr${CHR}_C2" \
    --chr "${CHR}"

echo "=== Step 2 chromosome ${CHR} complete ==="
ls -lh "${OUT_DIR}/"
