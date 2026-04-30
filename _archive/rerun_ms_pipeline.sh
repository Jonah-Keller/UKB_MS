#!/usr/bin/env bash
# rerun_ms_pipeline.sh — sequential re-run of MS-labeled analysis scripts after
# the UMAP1/UMAP2 -> PC1/PC2 covariate switch in 03_cns_trajectories.R.
# Each script is run with Rscript, logged to logs/<script_basename>.log, and
# the wall-clock / exit code appended to logs/rerun_summary.tsv.
#
# Usage:  bash rerun_ms_pipeline.sh
#
# Scripts run in approximate dependency order.  A failure in one script does
# NOT block subsequent scripts — failures are captured in the summary TSV.

set -u
cd "$(dirname "$0")"
mkdir -p logs
SUMMARY=logs/rerun_summary.tsv
: > "$SUMMARY"
printf "script\tstatus\texit_code\twall_sec\n" >> "$SUMMARY"

SCRIPTS=(
    analysis/04_differential/01_limma_ms_vs_hc.R
    analysis/04_differential/04_proteome_wide_ms.R
    analysis/04_differential/03_cns_trajectories.R
    analysis/06_ml/01_ms_ml_classifier.R
    analysis/06_ml/04_ms_preonset_ml.R
    analysis/06_ml/05_ms_sex_stratified_ml.R
    analysis/07_hla_subgroup/01_hla_ms_subgroup.R
    analysis/07_hla_subgroup/01_hla_interaction_limma.R
    analysis/08_mri/02_brain_mri_dti_ms.R
    analysis/08_mri/03_cognitive_ms.R
    analysis/09_sex_stratified/01_sex_stratified_ms_deps.R
    analysis/10_prs/01_ms_prs.R
    analysis/10_prs/01_ms_prs_protein_interaction.R
    analysis/11_cross_disease/01_ms_als_comparison.R
    analysis/12_clustering/01_ms_patient_clustering.R
    analysis/12_clustering/02_preonset_endophenotype.R
    analysis/13_survival/01_ms_protein_cox.R
    analysis/14_comorbidity/02a_ms_prems_phewas.R
    analysis/15_endophenotype/02_ms_prems_comorbidity_endophenotype.R
    analysis/15_endophenotype/03_ms_genetic_stratification.R
    analysis/15_endophenotype/04_ms_cluster_proteomics.R
    analysis/15_endophenotype/05_ms_cluster_pathways.R
    analysis/15_endophenotype/06_ms_cluster_ppi.R
    analysis/15_endophenotype/07_ms_cluster_celltypes.R
    analysis/15_endophenotype/08_postms_cluster_proteomics.R
    analysis/15_endophenotype/09_postms_cluster_pathways.R
    analysis/15_endophenotype/10_postms_cluster_celltypes.R
    analysis/15_endophenotype/11_ms_cluster_prems_ml.R
    analysis/15_endophenotype/12_ms_gmm_copresentation.R
    analysis/16_cluster_gwas/06_ms_gwas_lookup.R
)

for script in "${SCRIPTS[@]}"; do
    if [[ ! -f "$script" ]]; then
        printf "%s\tMISSING\t-\t-\n" "$script" >> "$SUMMARY"
        continue
    fi
    log="logs/$(echo "$script" | tr '/' '_').log"
    echo "=== $(date '+%H:%M:%S') START  $script ==="
    t0=$(date +%s)
    if Rscript "$script" >"$log" 2>&1; then
        status=OK; rc=0
    else
        rc=$?; status=FAIL
    fi
    t1=$(date +%s)
    printf "%s\t%s\t%d\t%d\n" "$script" "$status" "$rc" "$((t1 - t0))" >> "$SUMMARY"
    echo "=== $(date '+%H:%M:%S') $status  ($((t1 - t0))s, rc=$rc)  $script ==="
done

echo "=== pipeline complete.  See $SUMMARY ==="
