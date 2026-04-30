#!/usr/bin/env Rscript
# 04_ms_cluster_proteomics.R — Cluster-stratified differential proteomics
#
# For each pre-onset endophenotype cluster (C0, C1, C2), identifies
# differentially expressed proteins (DEPs) vs the cohort-None reference
# (cases with no pre-onset features). Time-to-diagnosis is included as a
# covariate to separate cluster-specific protein signatures from disease-stage
# effects.
#
# Model (limma, within cases only):
#   NPX ~ cluster_group + age_at_sampling + sex + years_to_diagnosis
#
# Comparisons:
#   C0-None, C1-None, C2-None   (each cluster vs no pre-onset comorbidity)
#   C1-C0, C2-C0, C2-C1         (cross-cluster pairwise)

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(here)
    library(glue)
    library(limma)
    library(MatchIt)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))    # CLUST_COLS
source(here::here("analysis", "helpers", "limma_utils.R"))  # run_cluster_psm_limma, make_volcano_plot

cfg <- load_disease_config()

QC_FILE    <- here::here("data", "ukb", "olink", "processed",
                         glue("{cfg$cohort_short}_olink_qc.csv"))
CLUST_FILE <- here::here("results", "endophenotype",
                         glue("{cfg$cohort_short}_prems_cluster_assignments.csv"))
FIG_DIR    <- here::here("results", "figures", "5")
OUT_DIR    <- here::here("results", "endophenotype", "cluster_proteomics")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_THR     <- 0.05
NONE_LABEL  <- glue("{cfg$disease_short_caps}-None")
status_col  <- cfg$cohort_status_col
case_states <- c(cfg$status_values$pre_onset, cfg$status_values$post_onset)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc    <- fread(QC_FILE, showProgress = FALSE)
clust <- fread(CLUST_FILE)

# Protein columns: exclude metadata
meta_cols <- c("eid", status_col, "age_at_sampling","sex","age_at_diagnosis",
               "years_to_diagnosis","qc_outlier","plate_id","well_position",
               "hla_carrier","prs_score","comorbidity_cluster")
prot_cols <- setdiff(names(qc), meta_cols)
prot_cols <- prot_cols[vapply(prot_cols, function(x) is.numeric(qc[[x]]), logical(1))]
cat(sprintf("  Protein columns detected: %d\n", length(prot_cols)))

# Merge cluster assignments onto QC table
case_only <- qc[qc_outlier == FALSE &
                get(status_col) %in% case_states &
                !is.na(age_at_sampling) & !is.na(sex)]
case_only <- merge(case_only, clust[, .(eid, cluster, n_enr_blocks,
                                        years_to_diagnosis)],
                   by = "eid", all.x = TRUE,
                   suffixes = c("", "_clust"))

if ("years_to_diagnosis_clust" %in% names(case_only)) {
    case_only[is.na(years_to_diagnosis), years_to_diagnosis := years_to_diagnosis_clust]
    case_only[, years_to_diagnosis_clust := NULL]
}

case_only[, cluster_f := factor(
    ifelse(is.na(cluster), "None", paste0("C", cluster)),
    levels = c("None", "C0", "C1", "C2")
)]

n_by_group <- case_only[, .N, by = cluster_f]
cat(glue("  {cfg$disease_short_caps} cohort by cluster:\n")); print(n_by_group)

# Filter to proteins with <30% missingness within this cohort
miss_pct  <- case_only[, lapply(.SD, function(x) mean(is.na(x))), .SDcols = prot_cols]
keep_prot <- names(miss_pct)[as.numeric(miss_pct) < 0.30]
cat(sprintf("  Proteins retained (<30%% missing): %d / %d\n",
            length(keep_prot), length(prot_cols)))

# ── 2. PSM + limma (via helper) ───────────────────────────────────────────────
cat("\nRunning PSM + limma...\n")
res     <- run_cluster_psm_limma(case_only, keep_prot, time_var = "years_to_diagnosis",
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
    list(cname = "C0_vs_None", label = "g", title = glue("C0 vs {NONE_LABEL}  ·  spine / connective tissue")),
    list(cname = "C1_vs_None", label = "h", title = glue("C1 vs {NONE_LABEL}  ·  cranial nerve / demyelinating")),
    list(cname = "C2_vs_None", label = "i", title = glue("C2 vs {NONE_LABEL}  ·  neurological sx + EBV"))
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
