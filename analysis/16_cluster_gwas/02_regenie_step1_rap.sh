#!/usr/bin/env bash
# 02_regenie_step1_rap.sh
# REGENIE step 1: whole-genome regression null model
#
# Runs on RAP inside Swiss Army Knife (mem2_ssd1_v2_x32 recommended).
# Input files are copied into /home/dnanexus/ by dx run.
# Genotype PLINK files are accessed via /mnt/project/.
#
# Runtime estimate: 2-4 hours on mem2_ssd1_v2_x32

set -euo pipefail

# Symlink dxfuse paths that contain spaces so plink2 can parse them
ln -sfn "/mnt/project/Bulk/Genotype Results/Genotype calls" /tmp/geno_calls
ln -sfn "/mnt/project/Bulk/Imputation/UKB imputation from genotype" /tmp/imputed

GENO_BASE="/tmp/geno_calls"
OUT_DIR="/home/dnanexus/out/out"
mkdir -p "${OUT_DIR}"

# ── 1. Merge + QC genotyped variants (autosomes 1-22) ─────────────────────────
echo "=== Merging genotyped PLINK files ==="
# Build merge list (chromosomes 2-22; chr1 is the base)
MERGE_LIST="/tmp/merge_list.txt"
> "${MERGE_LIST}"
for chr in $(seq 2 22); do
    echo "${GENO_BASE}/ukb22418_c${chr}_b0_v2" >> "${MERGE_LIST}"
done

plink2 \
    --bfile "${GENO_BASE}/ukb22418_c1_b0_v2" \
    --pmerge-list "${MERGE_LIST}" bfile \
    --make-pgen \
    --out /tmp/ukb_merged \
    --threads 24 \
    --memory 180000

echo "Merge complete. Running QC..."

# QC: MAF>1%, HWE p>1e-15, missingness<1%
# Note: no HWE filter because these are a mix of cases/controls
plink2 \
    --pfile /tmp/ukb_merged \
    --maf 0.01 \
    --geno 0.01 \
    --hwe 1e-15 midp \
    --make-pgen \
    --out /tmp/ukb_qc \
    --threads 24 \
    --memory 180000

echo "QC complete."

# LD pruning for step 1 (REGENIE recommends LD-pruned set)
plink2 \
    --pfile /tmp/ukb_qc \
    --indep-pairwise 1000 100 0.1 \
    --out /tmp/ukb_prune \
    --threads 24 \
    --memory 180000

plink2 \
    --pfile /tmp/ukb_qc \
    --extract /tmp/ukb_prune.prune.in \
    --make-pgen \
    --out /tmp/ukb_step1 \
    --threads 24 \
    --memory 180000

echo "LD pruning complete. Variants retained for step 1:"
wc -l /tmp/ukb_step1.pvar

# ── 2. REGENIE step 1 ──────────────────────────────────────────────────────────
echo "=== Running REGENIE step 1 ==="

regenie \
    --step 1 \
    --pgen /tmp/ukb_step1 \
    --phenoFile /home/dnanexus/pheno_C1_vs_none.tsv \
    --covarFile /home/dnanexus/covariates.tsv \
    --bt \
    --bsize 1000 \
    --loocv \
    --lowmem \
    --lowmem-prefix /tmp/regenie_tmp \
    --threads 24 \
    --out "${OUT_DIR}/step1_C1" \
    --force-step1

regenie \
    --step 1 \
    --pgen /tmp/ukb_step1 \
    --phenoFile /home/dnanexus/pheno_C2_vs_none.tsv \
    --covarFile /home/dnanexus/covariates.tsv \
    --bt \
    --bsize 1000 \
    --loocv \
    --lowmem \
    --lowmem-prefix /tmp/regenie_tmp2 \
    --threads 24 \
    --out "${OUT_DIR}/step1_C2" \
    --force-step1

echo "=== REGENIE step 1 complete ==="
ls -lh "${OUT_DIR}/"
