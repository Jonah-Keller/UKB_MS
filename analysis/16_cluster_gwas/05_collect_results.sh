#!/usr/bin/env bash
# 05_collect_results.sh
# Download REGENIE GWAS results from RAP and compile.
#
# Run after all step 2 jobs have completed.
# Usage: bash analysis/16_cluster_gwas/05_collect_results.sh

set -euo pipefail

PROJ="project-J2P0fqjJ1Gg4gByJ6BJ51pYj"
RAP_DIR="/data/ms_gwas_cluster/step2_output"
LOCAL_DIR="results/gwas_cluster"
mkdir -p "${LOCAL_DIR}/raw"

echo "======================================================"
echo " Downloading GWAS results from RAP"
echo "======================================================"

for chr in $(seq 1 22); do
    for phenotype in C1 C2; do
        REMOTE="${RAP_DIR}/chr${chr}/chr${chr}_${phenotype}_${phenotype}.regenie"
        LOCAL="${LOCAL_DIR}/raw/chr${chr}_${phenotype}.regenie"
        if [[ ! -f "${LOCAL}" ]]; then
            echo "  Downloading chr${chr} ${phenotype}..."
            dx download "${PROJ}:${REMOTE}" --output "${LOCAL}" --overwrite 2>/dev/null || \
                echo "    WARN: ${REMOTE} not found (job may still be running)"
        else
            echo "  Already exists: chr${chr}_${phenotype}"
        fi
    done
done

echo ""
echo "Compiling chromosome files..."

for phenotype in C1 C2; do
    OUT="${LOCAL_DIR}/${phenotype}_gwas_all_chr.tsv.gz"
    if ls "${LOCAL_DIR}/raw/chr"*"_${phenotype}.regenie" 1>/dev/null 2>&1; then
        # Add header from chr1, then append all chromosomes without header
        head -1 "${LOCAL_DIR}/raw/chr1_${phenotype}.regenie" > /tmp/gwas_${phenotype}.tsv
        for chr in $(seq 1 22); do
            f="${LOCAL_DIR}/raw/chr${chr}_${phenotype}.regenie"
            [[ -f "${f}" ]] && tail -n +2 "${f}" >> /tmp/gwas_${phenotype}.tsv
        done
        gzip -c /tmp/gwas_${phenotype}.tsv > "${OUT}"
        echo "  ${phenotype}: $(wc -l < /tmp/gwas_${phenotype}.tsv) variants → ${OUT}"
    else
        echo "  WARN: No files found for ${phenotype}"
    fi
done

echo ""
echo "Done. Next step:"
echo "  Rscript analysis/16_cluster_gwas/06_ms_gwas_lookup.R"
