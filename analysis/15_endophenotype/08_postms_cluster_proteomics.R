#!/usr/bin/env Rscript
# 08_postms_cluster_proteomics.R — Post-onset cluster differential proteomics
#
# Mirror of 04_ms_cluster_proteomics.R restricted to post-onset cases.
# Used for Figure 5 Supplement (5S) to show that endophenotype cluster protein
# signatures persist / differ in established disease.
#
# Model (limma, within post-onset only):
#   NPX ~ cluster_group + age_at_sampling + sex + years_since_diagnosis
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
FIG_DIR    <- here::here("results", "figures", "5S")
OUT_DIR    <- here::here("results", "endophenotype",
                         glue("{cfg$cohort_short}_postms_cluster_proteomics"))
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_THR    <- 0.05
NONE_LABEL <- glue("{cfg$disease_short_caps}-None")
POST_TAG   <- glue("post-{cfg$disease_short_caps}")
status_col <- cfg$cohort_status_col
post_state <- cfg$status_values$post_onset

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc    <- fread(QC_FILE, showProgress = FALSE)
clust <- fread(CLUST_FILE)

meta_cols <- c("eid", status_col, "age_at_sampling","sex","age_at_diagnosis",
               "years_to_diagnosis","qc_outlier","plate_id","well_position",
               "hla_carrier","prs_score","comorbidity_cluster")
prot_cols <- setdiff(names(qc), meta_cols)
prot_cols <- prot_cols[vapply(prot_cols, function(x) is.numeric(qc[[x]]), logical(1))]
cat(sprintf("  Protein columns detected: %d\n", length(prot_cols)))

# Restrict to post-onset cases only
post_only <- qc[qc_outlier == FALSE &
                get(status_col) == post_state &
                !is.na(age_at_sampling) & !is.na(sex)]
post_only <- merge(post_only, clust[, .(eid, cluster, n_enr_blocks,
                                        years_to_diagnosis)],
                   by = "eid", all.x = TRUE,
                   suffixes = c("", "_clust"))

if ("years_to_diagnosis_clust" %in% names(post_only)) {
    post_only[is.na(years_to_diagnosis), years_to_diagnosis := years_to_diagnosis_clust]
    post_only[, years_to_diagnosis_clust := NULL]
}

# years_since_diagnosis: positive = time elapsed since Dx
post_only[, years_since_diagnosis := -years_to_diagnosis]

.unique_clusters <- sort(unique(post_only$cluster[!is.na(post_only$cluster)]))
.cluster_labels  <- paste0("C", .unique_clusters)
post_only[, cluster_f := factor(
    ifelse(is.na(cluster), "None", paste0("C", cluster)),
    levels = c("None", .cluster_labels)
)]

n_by_group <- post_only[, .N, by = cluster_f]
cat("  Post-onset cohort by cluster:\n"); print(n_by_group)

miss_pct  <- post_only[, lapply(.SD, function(x) mean(is.na(x))), .SDcols = prot_cols]
keep_prot <- names(miss_pct)[as.numeric(miss_pct) < 0.30]
cat(sprintf("  Proteins retained (<30%% missing): %d / %d\n",
            length(keep_prot), length(prot_cols)))

# ── 2. PSM + limma (via helper) ───────────────────────────────────────────────
cat("\nRunning PSM + limma...\n")
res     <- run_cluster_psm_limma(post_only, keep_prot, time_var = "years_since_diagnosis",
                                  psm_ratio = 10L, psm_cal = 0.2, fdr_thr = FDR_THR)
all_dep <- res$all_dep
bonf    <- res$bonf

# ── 3. Write DEP tables ────────────────────────────────────────────────────────
cat("Extracting DEP tables...\n")
for (cname in names(all_dep)) {
    fwrite(all_dep[[cname]], file.path(OUT_DIR, sprintf("postms_deps_%s.csv", cname)))
    n_up   <- sum(all_dep[[cname]]$sig & all_dep[[cname]]$logFC > 0, na.rm = TRUE)
    n_down <- sum(all_dep[[cname]]$sig & all_dep[[cname]]$logFC < 0, na.rm = TRUE)
    cat(sprintf("  %-12s  FDR<0.05: %d up, %d down\n", cname, n_up, n_down))
}
fwrite(rbindlist(all_dep), file.path(OUT_DIR, "postms_deps_all_contrasts.csv"))

# ── 4. Volcano panels (a, b, c) ───────────────────────────────────────────────
cat("\nBuilding volcano panels...\n")

.volcano_palette <- cluster_palette(length(.cluster_labels))
volcano_specs <- lapply(seq_along(.cluster_labels), function(i) {
    cname <- paste0(.cluster_labels[i], "_vs_None")
    list(cname = cname,
         label = letters[i],   # a, b, c, d, ...
         title = glue("{.cluster_labels[i]} vs {NONE_LABEL}  (post-onset)"),
         color = .volcano_palette[[.cluster_labels[i]]])
})

for (spec in volcano_specs) {
    if (is.null(all_dep[[spec$cname]])) {
        cat(sprintf("  Skipping panel_%s.pdf: no DEPs for %s\n",
                    spec$label, spec$cname))
        next
    }
    p   <- make_volcano_plot(all_dep[[spec$cname]], spec$title, spec$label,
                             clust_col       = spec$color,
                             bonf            = bonf,
                             subtitle_suffix = glue("post-onset {cfg$disease_short_caps} only"))
    out <- file.path(FIG_DIR, sprintf("panel_%s.pdf", spec$label))
    ggsave(out, p, width = 3.5, height = 3.8, device = cairo_pdf)
    cat(sprintf("  panel_%s.pdf  (%s)\n", spec$label, spec$cname))
}

cat("\n08_postms_cluster_proteomics.R complete.\n")
