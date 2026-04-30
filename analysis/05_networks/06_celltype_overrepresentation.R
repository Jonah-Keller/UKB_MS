#!/usr/bin/env Rscript
# 06_celltype_overrepresentation.R
# Cell-type marker enrichment analysis for disease DEPs — refactored pipeline.
#
# Previous pipeline: Fisher's exact on binary up/down DEP lists vs. HPA
# nTPM>=10 marker sets. This failed because (a) marker sets covered 35-50%
# of the Olink universe, (b) overlap across cell types was massive, and
# (c) direction and effect size were collapsed into a hard threshold.
#
# Current pipeline:
#   Markers  : specificity-filtered (tau + per-cell-type z on log1p(nTPM))
#              via load_celltype_markers() in helpers/celltype_overrep_plot.R.
#   Test     : preranked fgsea on the full limma moderated-t statistic
#              (signed, rank-based, uses all proteins in the universe).
#              Optional confirmatory cameraPR run (inter-gene corr aware).
#   Legacy   : for side-by-side reproducibility the old Fisher test is run
#              with the old nTPM>=10 markers and emitted as a separate CSV.
#
# Cohort prefix is read from configs/disease.yaml. Analyses (unchanged):
#   combined  — {cohort_short}_combined_vs_hc.csv
#   pre       — {cohort_short}_pre_vs_hc.csv
#   post      — {cohort_short}_post_vs_hc.csv
#   preVsPost — {cohort_short}_pre_vs_post.csv
#
# Outputs (results/networks/):
#   {cohort_short}_celltype_gsea.csv              NEW primary file; columns:
#                                     cell_type, analysis, analysis_label,
#                                     size, NES, ES, pval, padj, leading_edge
#   {cohort_short}_celltype_camerapr.csv          NEW secondary (inter-gene corr confirmatory)
#   {cohort_short}_celltype_markers_summary.csv   NEW marker-set sizes for reporting
#   {cohort_short}_celltype_overrep.csv           LEGACY Fisher results
#   {cohort_short}_celltype_gsea_plot.pdf         NEW summary dot plot

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(glue)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "celltype_overrep_plot.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config(file.path(PROJ_DIR, "configs", "disease.yaml"))
COHORT  <- cfg$cohort_short
DISEASE <- cfg$disease_short_caps

DIFF_DIR    <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR     <- file.path(PROJ_DIR, "results", "networks")
HPA_FILE    <- file.path(PROJ_DIR, "data", "reference", "hpa_rna_single_cell_type.tsv.zip")
WALCHLI_CSV <- file.path(PROJ_DIR, "data", "reference", "walchli_avg_expr_by_celltype.csv")
CACHE_DIR   <- file.path(PROJ_DIR, "data", "reference")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(HPA_FILE) && !file.exists(sub("\\.zip$", ".tsv", HPA_FILE)))
    stop("HPA file missing. Download: https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip")

# ---------------------------------------------------------------------------
# 1. Load specificity-based marker sets (primary) and legacy sets (reference)
# ---------------------------------------------------------------------------
cat("Loading cell-type markers (specificity-filtered tau+z)...\n")
marker_sets <- load_celltype_markers(
    walchli_csv = WALCHLI_CSV,
    hpa_file    = HPA_FILE,
    method      = "specificity",
    hpa_tau_thr = 0.60, hpa_z_thr = 1.50, hpa_min_ntpm = 1.0,
    walchli_tau_thr = 0.45, walchli_z_thr = 1.00, walchli_min_ntpm = 0.5,
    verbose     = TRUE
)

cat("\nLoading legacy nTPM>=10 markers (for side-by-side reference)...\n")
marker_sets_legacy <- load_celltype_markers(
    walchli_csv = WALCHLI_CSV,
    hpa_file    = HPA_FILE,
    method      = "legacy_ntpm",
    hpa_ntpm_thr   = 10,
    walchli_pctile = 0.65, walchli_min_lp = 1.0,
    verbose     = TRUE
)

# Marker-set size summary CSV.
size_dt <- data.table::rbindlist(list(
    data.table(method = "specificity",
                cell_type = names(marker_sets),
                n_markers = vapply(marker_sets, length, integer(1L))),
    data.table(method = "legacy_ntpm",
                cell_type = names(marker_sets_legacy),
                n_markers = vapply(marker_sets_legacy, length, integer(1L)))
), use.names = TRUE, fill = TRUE)
markers_csv <- glue("{COHORT}_celltype_markers_summary.csv")
fwrite(size_dt, file.path(OUT_DIR, markers_csv))
cat(sprintf("  Saved %s (%d rows)\n", markers_csv, nrow(size_dt)))

# ---------------------------------------------------------------------------
# 2. Analyses to run
# ---------------------------------------------------------------------------
analyses <- list(
    list(name = "combined",  label = sprintf("Combined %s vs HC", DISEASE),
         file = glue("{COHORT}_combined_vs_hc.csv"), fdr_thr = 0.05),
    list(name = "pre",       label = "Pre-onset vs HC",
         file = glue("{COHORT}_pre_vs_hc.csv"),      fdr_thr = 0.20),
    list(name = "post",      label = "Post-onset vs HC",
         file = glue("{COHORT}_post_vs_hc.csv"),     fdr_thr = 0.05),
    list(name = "preVsPost", label = "Pre vs post-onset",
         file = glue("{COHORT}_pre_vs_post.csv"),    fdr_thr = 0.05)
)

# ---------------------------------------------------------------------------
# 3. Run preranked fgsea + cameraPR + legacy Fisher for each analysis
# ---------------------------------------------------------------------------
gsea_rows     <- list()
camera_rows   <- list()
legacy_rows   <- list()

for (an in analyses) {
    diff_path <- file.path(DIFF_DIR, an$file)
    if (!file.exists(diff_path)) {
        cat(sprintf("  Skipping %s: file not found (%s)\n", an$name, an$file))
        next
    }
    diff <- fread(diff_path)
    if (!all(c("protein", "t", "logFC", "P.Value") %in% names(diff))) {
        cat(sprintf("  Skipping %s: required columns missing\n", an$name))
        next
    }
    diff <- diff[!is.na(t) & !is.na(protein)]

    cat(sprintf("\n=== %s  (n_universe = %d) ===\n", an$label, nrow(diff)))

    # --- GSEA (primary) --------------------------------------------------
    rank_stat <- build_rank_stat(diff, gene_col = "protein", stat_col = "t")
    g <- run_celltype_gsea(
        stats       = rank_stat,
        marker_sets = marker_sets,
        min_size    = 5,
        max_size    = 500,
        seed        = 42
    )
    if (!is.null(g) && nrow(g) > 0L) {
        g[, analysis        := an$name]
        g[, analysis_label  := an$label]
        gsea_rows[[length(gsea_rows) + 1L]] <- g
    }

    # --- cameraPR (confirmatory) ----------------------------------------
    cm <- run_celltype_camerapr(
        stats       = rank_stat,
        marker_sets = marker_sets,
        inter_gene_cor = 0.01,
        min_size = 5, max_size = 500
    )
    if (!is.null(cm) && nrow(cm) > 0L) {
        cm[, analysis       := an$name]
        cm[, analysis_label := an$label]
        camera_rows[[length(camera_rows) + 1L]] <- cm
    }

    # --- Legacy Fisher on old nTPM>=10 markers ---------------------------
    sig <- diff[!is.na(adj.P.Val) & adj.P.Val < an$fdr_thr]
    up   <- toupper(sig[logFC > 0, protein])
    dn   <- toupper(sig[logFC < 0, protein])
    universe <- toupper(diff$protein)
    for (direction in c("up", "down")) {
        g_dep <- if (direction == "up") up else dn
        if (length(g_dep) < 3L) next
        r <- run_ora(module_genes = g_dep,
                      bg_genes     = universe,
                      gs_list      = marker_sets_legacy)
        if (is.null(r) || nrow(r) == 0L) next
        r[, analysis        := an$name]
        r[, analysis_label  := an$label]
        r[, direction       := direction]
        r[, fdr_threshold   := an$fdr_thr]
        legacy_rows[[length(legacy_rows) + 1L]] <- r
    }
}

# ---------------------------------------------------------------------------
# 4. Consolidate + write outputs
# ---------------------------------------------------------------------------
if (length(gsea_rows) == 0L)
    stop("No GSEA results produced — check that differential CSVs exist in ", DIFF_DIR)

gsea_dt <- rbindlist(gsea_rows, use.names = TRUE, fill = TRUE)
# Re-apply BH FDR within each analysis (fgsea already BH-adjusts per call; this
# is the same within-analysis correction made explicit and recomputed to be
# safe across cell types).
gsea_dt[, padj := p.adjust(pval, method = "BH"), by = analysis]
setcolorder(gsea_dt, c("analysis", "analysis_label", "cell_type",
                        "size", "NES", "ES", "pval", "padj", "leading_edge"))
gsea_csv <- glue("{COHORT}_celltype_gsea.csv")
fwrite(gsea_dt, file.path(OUT_DIR, gsea_csv))
cat(sprintf("\nSaved %s (%d rows)\n", gsea_csv, nrow(gsea_dt)))

if (length(camera_rows) > 0L) {
    camera_dt <- rbindlist(camera_rows, use.names = TRUE, fill = TRUE)
    camera_dt[, padj := p.adjust(pval, method = "BH"), by = analysis]
    camera_csv <- glue("{COHORT}_celltype_camerapr.csv")
    fwrite(camera_dt, file.path(OUT_DIR, camera_csv))
    cat(sprintf("Saved %s (%d rows)\n", camera_csv, nrow(camera_dt)))
}

if (length(legacy_rows) > 0L) {
    legacy_dt <- rbindlist(legacy_rows, use.names = TRUE, fill = TRUE)
    legacy_dt[, p_adj := p.adjust(p_twosided, method = "BH"),
              by = .(analysis, direction)]
    legacy_dt[, log2_OR := log2(pmax(odds_ratio, 0.01))]
    legacy_csv <- glue("{COHORT}_celltype_overrep.csv")
    fwrite(legacy_dt, file.path(OUT_DIR, legacy_csv))
    cat(sprintf("Saved %s (LEGACY Fisher, %d rows)\n",
                legacy_csv, nrow(legacy_dt)))
}

# ---------------------------------------------------------------------------
# 5. Summary print
# ---------------------------------------------------------------------------
cat("\nGSEA significant cell types (FDR < 0.05):\n")
print(gsea_dt[padj < 0.05, .(analysis, cell_type, size, NES = round(NES, 2),
                               padj = signif(padj, 2))])

# ---------------------------------------------------------------------------
# 6. Summary dot plot — all analyses side-by-side
# ---------------------------------------------------------------------------
cat("\nBuilding GSEA summary plot...\n")
an_order  <- c("combined", "pre", "post", "preVsPost")
an_labels <- c(combined = "Combined", pre = "Pre-onset",
               post = "Post-onset", preVsPost = "Pre vs post")

p <- make_celltype_gsea_plot(
    gsea_dt,
    contrast_col    = "analysis",
    contrast_order  = an_order,
    contrast_labels = an_labels,
    title_str       = sprintf("%s cell-type enrichment vs HC (preranked fgsea)", DISEASE),
    subtitle_str    = sprintf(
        "Preranked fgsea on limma moderated-t | NES>0: markers higher in %s vs HC | NES<0: lower | * FDR<0.05",
        DISEASE
    )
)
if (!is.null(p)) {
    plot_pdf <- glue("{COHORT}_celltype_gsea_plot.pdf")
    ggsave(file.path(OUT_DIR, plot_pdf),
           p, width = 8.5, height = 6.5, device = cairo_pdf)
    cat(sprintf("  Saved %s\n", plot_pdf))
}

cat("\nDone.\n")
