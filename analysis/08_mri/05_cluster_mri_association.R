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
#   1. WMH volume: Wilcoxon rank-sum on the two largest clusters
#   2. DTI FA (27 tracts): Wilcoxon per tract, FDR correction
#   3. DTI MD (27 tracts): same
#   4. Scatter plot: WMH by cluster
#   5. Heatmap: FA z-scores by cluster
#
# Output: results/mri/<cohort>_cluster_*.{pdf,csv}

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

# The cluster-assignments CSV writes its case-status column as a generic
# `status` column (see 15_endophenotype/02 line ~245), not the cohort-specific
# `<cohort>_status` name from the disease YAML.  Status values are still
# read from cfg.
STATUS_COL <- "status"
PRE        <- cfg$status_values$pre_onset

# Cluster assignments live in results/endophenotype/<cohort>_prems_cluster_assignments.csv
# (per stage 15_endophenotype/02 conventions). The legacy results/clustering/...
# path was the MS-only template name; tolerate both.
CLUSTER_FILE_NEW <- here::here("results", "endophenotype",
                               glue::glue("{cfg$cohort_short}_prems_cluster_assignments.csv"))
CLUSTER_FILE_OLD <- here::here("results", "clustering",
                               glue::glue("{cfg$cohort_short}_cluster_assignments.csv"))
CLUSTER_FILE <- if (file.exists(CLUSTER_FILE_NEW)) CLUSTER_FILE_NEW else CLUSTER_FILE_OLD

DTI_FILE     <- file.path(dirname(here::here()),
                          "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                          "data", "ukb", "brain_mri", "dti", "dti_tract_protein_data.tsv")
OUT_DIR      <- here::here("results", "mri")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(CLUSTER_FILE)) {
    cat(sprintf("05_cluster_mri_association.R skipped: %s not found.\n", basename(CLUSTER_FILE)))
    cat("  Run analysis/15_endophenotype/02_*.R first.\n")
    quit(status = 0)
}

# ── 1. Load ───────────────────────────────────────────────────────────────────
cat("Loading data...\n")
clusters <- fread(CLUSTER_FILE)[get(STATUS_COL) == PRE]
dti      <- fread(DTI_FILE)

dt <- merge(clusters[, .(eid, cluster)], dti, by = "eid")
# Discover the cluster set at runtime so this works for any k.  Drop the
# "None" / NA / empty assignments which mean "no cluster".  Pairwise tests
# use the two largest clusters; if k>2 the remainder are kept in plots but
# excluded from the pairwise inference.
dt <- dt[!is.na(cluster) & !(as.character(cluster) %in% c("None", ""))]
clust_levels <- sort(unique(as.character(dt$cluster)))
if (length(clust_levels) < 2L) {
    cat(sprintf("05_cluster_mri_association.R skipped: only %d non-None cluster(s) on disk.\n",
                length(clust_levels)))
    quit(status = 0)
}
dt <- dt[as.character(cluster) %in% clust_levels]
dt[, cluster_f := factor(as.character(cluster), levels = clust_levels)]
n_per_cluster <- dt[, .N, by = cluster_f][order(-N)]
# Two reference clusters for pairwise testing (largest by sample size)
ref_levels  <- as.character(n_per_cluster$cluster_f[1:min(2, nrow(n_per_cluster))])
c1_lab      <- paste0("C", ref_levels[1])
c2_lab      <- if (length(ref_levels) >= 2L) paste0("C", ref_levels[2]) else NA_character_
clust_label <- function(x) paste0("C", x)
n_label_str <- paste(sprintf("%s=%d", clust_label(n_per_cluster$cluster_f),
                              n_per_cluster$N), collapse = ", ")
cat(sprintf("  Pre-onset patients with MRI: %d (%s)\n", nrow(dt), n_label_str))

# Resolve cluster colours: prefer CLUST_COLS["C<id>"] if defined, else
# fall back to a default palette to support arbitrary cluster IDs.
.fill_for <- function(lab) {
    v <- tryCatch(unname(CLUST_COLS[lab]), error = function(e) NA_character_)
    if (is.null(v) || length(v) == 0L || is.na(v)) NA_character_ else v
}
default_pal <- c("#CC0066","#2B4C7E","#E6A817","#56B4E9","#3FA34D","#7A5195","#BC5090")
fill_vals <- setNames(
    vapply(seq_along(clust_levels), function(i) {
        lab <- clust_label(clust_levels[i])
        v   <- .fill_for(lab)
        if (is.na(v)) default_pal[((i - 1L) %% length(default_pal)) + 1L] else v
    }, character(1)),
    clust_levels
)

# ── 2. WMH test (pairwise on two largest clusters) ────────────────────────────
dt_pair <- dt[as.character(cluster) %in% ref_levels]
if (!is.na(c2_lab) && nrow(dt_pair) >= 4L) {
    wmh_test <- wilcox.test(wmh_volume ~ cluster, data = dt_pair, exact = FALSE)
    wmh_p <- wmh_test$p.value
    cat(sprintf("  WMH Wilcoxon (%s vs %s) p = %.3f\n", c1_lab, c2_lab, wmh_p))
} else {
    wmh_p <- NA_real_
    cat("  WMH test skipped (need ≥2 clusters with sufficient n).\n")
}

p_wmh <- ggplot(dt, aes(cluster_f, log1p(wmh_volume), fill = cluster_f)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    labs(
        x     = "Pre-onset cluster",
        y     = "WMH volume (log mm³)",
        title = "White matter hyperintensity volume by cluster",
        subtitle = sprintf("Wilcoxon (%s vs %s) p = %s  |  n = %d (%s)",
                           c1_lab,
                           if (is.na(c2_lab)) "—" else c2_lab,
                           if (is.na(wmh_p)) "NA" else sprintf("%.2f", wmh_p),
                           nrow(dt), n_label_str)
    ) +
    theme_ukb()

# ── 3. DTI tract tests ────────────────────────────────────────────────────────
fa_cols <- grep("^fa_", colnames(dt), value = TRUE)
md_cols <- grep("^md_", colnames(dt), value = TRUE)

run_wilcox_tracts <- function(cols) {
    if (is.na(c2_lab)) return(data.table())
    rbindlist(lapply(cols, function(col) {
        sub <- dt_pair[!is.na(get(col))]
        if (nrow(sub) < 4) return(NULL)
        res <- wilcox.test(sub[[col]] ~ sub$cluster, exact = FALSE)
        list(
            tract   = sub("^(fa|md|icvf)_", "", col),
            metric  = sub("_.*", "", col),
            p_value = res$p.value,
            median_c1 = median(sub[[col]][as.character(sub$cluster) == ref_levels[1]], na.rm = TRUE),
            median_c2 = median(sub[[col]][as.character(sub$cluster) == ref_levels[2]], na.rm = TRUE)
        )
    }))
}

fa_res <- run_wilcox_tracts(fa_cols)
md_res <- run_wilcox_tracts(md_cols)
all_res <- rbind(fa_res, md_res)
all_res[, fdr := p.adjust(p_value, method = "BH")]
all_res[, direction := ifelse(median_c1 > median_c2,
                              paste0(c1_lab, ">", c2_lab),
                              paste0(c2_lab, ">", c1_lab))]
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

.dir_levels <- c(paste0(c1_lab, ">", c2_lab), paste0(c2_lab, ">", c1_lab))
.dir_cols   <- setNames(c(COL_PRE_UP, COL_DOWN), .dir_levels)

p_fa_scatter <- ggplot(fa_plot, aes(delta_fa, neg_log_p, label = label)) +
    geom_point(aes(colour = direction), size = 2.5, alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
    geom_text_repel(size = 2.8, max.overlaps = 20, na.rm = TRUE) +
    scale_colour_manual(values = .dir_cols, name = "Higher FA in") +
    labs(
        x     = sprintf("Median FA difference (%s - %s)", c1_lab, c2_lab),
        y     = "-log10 p",
        title = "DTI FA: cluster differences"
    ) +
    theme_ukb()

# ── 6. Save outputs ───────────────────────────────────────────────────────────
.prefix <- function(name) file.path(OUT_DIR, sprintf("%s_%s", cfg$cohort_short, name))

fwrite(all_res, .prefix("cluster_dti_tests.csv"))

pdf(.prefix("cluster_wmh_boxplot.pdf"), width = 4, height = 4.5)
print(p_wmh)
dev.off()

pdf(.prefix("cluster_dti_heatmap.pdf"), width = 5.5, height = 8)
print(p_heatmap)
dev.off()

pdf(.prefix("cluster_dti_scatter.pdf"), width = 5.5, height = 4.5)
print(p_fa_scatter)
dev.off()

cat("\nDone. Outputs in", OUT_DIR, "\n")
