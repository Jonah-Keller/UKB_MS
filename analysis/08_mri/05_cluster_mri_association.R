#!/usr/bin/env Rscript
# 05_cluster_mri_association.R
# Test whether pre-onset disease clusters differ in brain MRI metrics.
#
# Design note: MRI was acquired at UKB instance 2 (~2014-2020); Olink blood draw
# at instance 0 (~2006-2010). By MRI time many "pre-onset" patients (defined at
# blood draw) will have been diagnosed — so results are exploratory / hypothesis
# generating only.
#
# Tests:
#   1. WMH volume: Wilcoxon rank-sum C1 vs C2
#   2. DTI FA (27 tracts): Wilcoxon per tract, FDR correction
#   3. DTI MD (27 tracts): same
#   4. Scatter plot: WMH by cluster
#   5. Heatmap: FA z-scores by cluster
#
# Output: results/mri/
#   cluster_wmh_boxplot.pdf
#   cluster_dti_heatmap.pdf
#   cluster_dti_tests.csv

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
    library(glue)
    library(here)
})

source(here::here("analysis", "helpers", "ukb_theme.R"))
source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

STATUS_COL <- cfg$cohort_status_col
PRE        <- cfg$status_values$pre_onset

CLUSTER_FILE <- here::here("results", "clustering",
                           glue::glue("{cfg$cohort_short}_cluster_assignments.csv"))
DTI_FILE     <- file.path(dirname(here::here()),
                          "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                          "data", "ukb", "brain_mri", "dti", "dti_tract_protein_data.tsv")
OUT_DIR      <- here::here("results", "mri")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load ───────────────────────────────────────────────────────────────────
cat("Loading data...\n")
clusters <- fread(CLUSTER_FILE)[get(STATUS_COL) == PRE]
dti      <- fread(DTI_FILE)

dt <- merge(clusters[, .(eid, cluster)], dti, by = "eid")
dt[, cluster_f := factor(cluster, levels = c(1, 2), labels = c("C1", "C2"))]
cat(sprintf("  Pre-onset patients with MRI: %d (C1=%d, C2=%d)\n",
            nrow(dt), sum(dt$cluster == 1), sum(dt$cluster == 2)))

# ── 2. WMH test ───────────────────────────────────────────────────────────────
wmh_test <- wilcox.test(wmh_volume ~ cluster, data = dt, exact = FALSE)
cat(sprintf("  WMH Wilcoxon p = %.3f\n", wmh_test$p.value))

p_wmh <- ggplot(dt, aes(cluster_f, log1p(wmh_volume), fill = cluster_f)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c(C1 = unname(CLUST_COLS["C1"]), C2 = unname(CLUST_COLS["C2"])),
                      guide = "none") +
    labs(
        x     = "Pre-onset cluster",
        y     = "WMH volume (log mm³)",
        title = "White matter hyperintensity volume by cluster",
        subtitle = sprintf("Wilcoxon p = %.2f  |  n = %d (C1=%d, C2=%d)",
                           wmh_test$p.value, nrow(dt),
                           sum(dt$cluster == 1), sum(dt$cluster == 2))
    ) +
    theme_ukb()

# ── 3. DTI tract tests ────────────────────────────────────────────────────────
fa_cols <- grep("^fa_", colnames(dt), value = TRUE)
md_cols <- grep("^md_", colnames(dt), value = TRUE)

run_wilcox_tracts <- function(cols) {
    rbindlist(lapply(cols, function(col) {
        vals <- dt[[col]]
        if (sum(!is.na(vals)) < 4) return(NULL)
        res <- wilcox.test(vals ~ dt$cluster, exact = FALSE)
        list(
            tract   = sub("^(fa|md|icvf)_", "", col),
            metric  = sub("_.*", "", col),
            p_value = res$p.value,
            median_c1 = median(vals[dt$cluster == 1], na.rm = TRUE),
            median_c2 = median(vals[dt$cluster == 2], na.rm = TRUE)
        )
    }))
}

fa_res <- run_wilcox_tracts(fa_cols)
md_res <- run_wilcox_tracts(md_cols)
all_res <- rbind(fa_res, md_res)
all_res[, fdr := p.adjust(p_value, method = "BH")]
all_res[, direction := ifelse(median_c1 > median_c2, "C1>C2", "C2>C1")]
setorder(all_res, p_value)

cat(sprintf("  DTI tracts tested: %d FA + %d MD\n", nrow(fa_res), nrow(md_res)))
cat(sprintf("  FDR < 0.05: %d tracts\n", sum(all_res$fdr < 0.05, na.rm = TRUE)))
if (any(all_res$p_value < 0.05)) {
    cat("  Nominal p < 0.05:\n")
    print(all_res[p_value < 0.05, .(metric, tract, p_value, fdr, direction)])
}

# ── 4. FA heatmap ─────────────────────────────────────────────────────────────
fa_long <- melt(dt[, c("cluster_f", fa_cols), with = FALSE],
                id.vars = "cluster_f",
                variable.name = "tract", value.name = "fa")
fa_long[, tract := sub("^fa_", "", tract)]

# z-score within each tract
fa_long[, fa_z := scale(fa)[,1], by = tract]

# Order tracts by median difference
tract_order <- fa_res[order(abs(median_c1 - median_c2), decreasing = TRUE), tract]
fa_long[, tract := factor(tract, levels = rev(tract_order))]

p_heatmap <- ggplot(fa_long, aes(cluster_f, tract, fill = fa_z)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = COL_DOWN, mid = "white", high = COL_PRE_UP,
                         midpoint = 0, name = "FA\n(z-score)") +
    labs(
        x     = "Pre-onset cluster",
        y     = NULL,
        title = "DTI fractional anisotropy by pre-onset cluster",
        subtitle = sprintf("n = %d pre-onset %s patients with MRI data",
                           nrow(dt), cfg$disease_short_caps)
    ) +
    theme_ukb(base_size = 9) +
    theme(axis.text.y = element_text(size = 7))

# ── 5. Top FA volcano-style plot ──────────────────────────────────────────────
fa_plot <- all_res[metric == "fa"]
fa_plot[, neg_log_p := -log10(p_value)]
fa_plot[, delta_fa  := median_c1 - median_c2]
fa_plot[, label := fifelse(p_value < 0.1,
                           gsub("_", " ", tract),
                           NA_character_)]

p_fa_scatter <- ggplot(fa_plot, aes(delta_fa, neg_log_p, label = label)) +
    geom_point(aes(colour = direction), size = 2.5, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
    geom_text_repel(size = 2.8, max.overlaps = 20, na.rm = TRUE) +
    scale_colour_manual(values = c("C1>C2" = COL_PRE_UP, "C2>C1" = COL_DOWN),
                        name = "Higher FA in") +
    labs(
        x     = "Median FA difference (C1 - C2)",
        y     = "-log10 p",
        title = "DTI FA: cluster differences"
    ) +
    theme_ukb()

# ── 6. Save outputs ───────────────────────────────────────────────────────────
fwrite(all_res, file.path(OUT_DIR, "cluster_dti_tests.csv"))

pdf(file.path(OUT_DIR, "cluster_wmh_boxplot.pdf"), width = 4, height = 4.5)
print(p_wmh)
dev.off()

pdf(file.path(OUT_DIR, "cluster_dti_heatmap.pdf"), width = 5.5, height = 8)
print(p_heatmap)
dev.off()

pdf(file.path(OUT_DIR, "cluster_dti_scatter.pdf"), width = 5.5, height = 4.5)
print(p_fa_scatter)
dev.off()

cat("\nDone. Outputs in", OUT_DIR, "\n")
