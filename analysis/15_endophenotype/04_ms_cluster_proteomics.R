#!/usr/bin/env Rscript
# 04_ms_cluster_proteomics.R — Cluster-stratified differential proteomics
#
# For each pre-MS endophenotype cluster (C0, C1, C2), identifies differentially
# expressed proteins (DEPs) vs the MS-None reference (MS patients with no
# pre-MS features). Time-to-diagnosis is included as a covariate to separate
# cluster-specific protein signatures from MS-stage effects.
#
# Model (limma, within MS only):
#   NPX ~ cluster_group + age_at_sampling + sex + years_to_diagnosis
#
# Comparisons:
#   C0-None, C1-None, C2-None   (each cluster vs no pre-MS comorbidity)
#   C1-C0, C2-C0, C2-C1         (cross-cluster pairwise)
#
# Outputs (results/endophenotype/cluster_proteomics/):
#   cluster_deps_<contrast>.csv  — full DEP table per contrast
#   cluster_deps_summary.csv     — top hits across all contrasts
#   panel_g_<c>.pdf              — volcano per cluster vs None (→ figure5.R)

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(limma)
    library(MatchIt)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))    # CLUST_COLS
source(file.path(PROJ_DIR, "analysis", "helpers", "limma_utils.R"))  # run_cluster_psm_limma, make_volcano_plot

QC_FILE    <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
CLUST_FILE <- file.path(PROJ_DIR, "results", "endophenotype", "ms_prems_cluster_assignments.csv")
FIG_DIR    <- file.path(PROJ_DIR, "results", "figures", "5")
OUT_DIR    <- file.path(PROJ_DIR, "results", "endophenotype", "cluster_proteomics")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_THR <- 0.05

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc    <- fread(QC_FILE, showProgress = FALSE)
clust <- fread(CLUST_FILE)

# Protein columns: exclude metadata
meta_cols <- c("eid","ms_status","age_at_sampling","sex","age_at_diagnosis",
               "years_to_diagnosis","qc_outlier","plate_id","well_position",
               "hla_carrier","prs_score","comorbidity_cluster")
prot_cols <- setdiff(names(qc), meta_cols)
prot_cols <- prot_cols[vapply(prot_cols, function(x) is.numeric(qc[[x]]), logical(1))]
cat(sprintf("  Protein columns detected: %d\n", length(prot_cols)))

# Merge cluster assignments onto QC table
ms_only <- qc[qc_outlier == FALSE &
              ms_status %in% c("pre_onset", "post_onset") &
              !is.na(age_at_sampling) & !is.na(sex)]
ms_only <- merge(ms_only, clust[, .(eid, cluster, n_enr_blocks,
                                      years_to_diagnosis)],
                 by = "eid", all.x = TRUE,
                 suffixes = c("", "_clust"))

# Use cluster-assigned years_to_diagnosis (should be same as qc's if column exists)
if ("years_to_diagnosis_clust" %in% names(ms_only)) {
    ms_only[is.na(years_to_diagnosis), years_to_diagnosis := years_to_diagnosis_clust]
    ms_only[, years_to_diagnosis_clust := NULL]
}

ms_only[, cluster_f := factor(
    ifelse(is.na(cluster), "None", paste0("C", cluster)),
    levels = c("None", "C0", "C1", "C2")
)]

n_by_group <- ms_only[, .N, by = cluster_f]
cat("  MS cohort by cluster:\n"); print(n_by_group)

# Filter to proteins with <30% missingness within this cohort
miss_pct <- ms_only[, lapply(.SD, function(x) mean(is.na(x))), .SDcols = prot_cols]
keep_prot <- names(miss_pct)[as.numeric(miss_pct) < 0.30]
cat(sprintf("  Proteins retained (<30%% missing): %d / %d\n",
            length(keep_prot), length(prot_cols)))

# ── 2. PSM + limma (via helper) ───────────────────────────────────────────────
cat("\nRunning PSM + limma...\n")
res     <- run_cluster_psm_limma(ms_only, keep_prot, time_var = "years_to_diagnosis",
                                  psm_ratio = 10L, psm_cal = 0.2, fdr_thr = FDR_THR)
all_dep <- res$all_dep
bonf    <- res$bonf

# ── 3. Write DEP tables ────────────────────────────────────────────────────────
cat("Extracting DEP tables...\n")
for (cname in names(all_dep)) {
    fwrite(all_dep[[cname]], file.path(OUT_DIR, sprintf("cluster_deps_%s.csv", cname)))
    n_up   <- sum(all_dep[[cname]]$sig & all_dep[[cname]]$logFC > 0, na.rm = TRUE)
    n_down <- sum(all_dep[[cname]]$sig & all_dep[[cname]]$logFC < 0, na.rm = TRUE)
    cat(sprintf("  %-12s  FDR<0.05: %d up, %d down\n", cname, n_up, n_down))
}
fwrite(rbindlist(all_dep), file.path(OUT_DIR, "cluster_deps_all_contrasts.csv"))

# ── 4. Volcano panels (g, h, i) ───────────────────────────────────────────────
cat("\nBuilding volcano panels...\n")

volcano_specs <- list(
    list(cname = "C0_vs_None", label = "g", title = "C0 vs MS-None  \u00b7  spine / connective tissue"),
    list(cname = "C1_vs_None", label = "h", title = "C1 vs MS-None  \u00b7  cranial nerve / demyelinating"),
    list(cname = "C2_vs_None", label = "i", title = "C2 vs MS-None  \u00b7  neurological sx + EBV")
)

for (spec in volcano_specs) {
    p   <- make_volcano_plot(all_dep[[spec$cname]], spec$title, spec$label,
                             clust_col = CLUST_COLS[sub("_vs_None", "", spec$cname)],
                             bonf      = bonf)
    out <- file.path(FIG_DIR, sprintf("panel_%s.pdf", spec$label))
    ggsave(out, p, width = 3.5, height = 3.8, device = cairo_pdf)
    cat(sprintf("  panel_%s.pdf\n", spec$label))
}

cat("\n04_ms_cluster_proteomics.R complete.\n")
cat("Next: 05_ms_cluster_pathways.R (pathway enrichment on DEPs)\n")
cat("      06_ms_cluster_ppi.R       (STRING PPI network)\n")
cat("      07_ms_cluster_celltypes.R (scRNA-seq cell type mapping)\n")
