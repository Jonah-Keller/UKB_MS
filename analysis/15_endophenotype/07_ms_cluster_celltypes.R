#!/usr/bin/env Rscript
# 07_ms_cluster_celltypes.R — cell-type enrichment per MS cluster DEP set
#
# Refactored (2025) to use specificity-filtered markers + preranked fgsea
# instead of nTPM>=10 Fisher. See helpers/celltype_overrep_plot.R for the
# rationale. Legacy Fisher output is retained in parallel for reproducibility
# of pre-refactor supplementary material.
#
# Per cluster (C0, C1, C2 vs MS-None):
#   Primary   : fgsea on moderated-t from limma (signed, rank-based)
#   Secondary : cameraPR (inter-gene correlation aware)
#   Legacy    : one-sided Fisher ORA on up/down modules (P.Value < 0.05)
#
# Outputs (results/endophenotype/cluster_celltypes/):
#   celltypes_gsea_<contrast>.csv      NEW per-contrast fgsea results
#   celltypes_camerapr_<contrast>.csv  NEW per-contrast cameraPR results
#   celltypes_ora_<contrast>_up.csv / _down.csv   LEGACY Fisher
#   cluster_celltype_summary_gsea.csv  NEW consolidated GSEA
#   cluster_celltype_summary.csv       LEGACY consolidated Fisher
#   results/figures/5/panel_r/s/t.pdf  per-cluster dot plot (NES-based)
#   results/figures/5/panel_u.pdf      cross-cluster heatmap (NES-based)

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
    library(glue)
    library(scales)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))
source(here::here("analysis", "helpers", "celltype_overrep_plot.R"))

cfg <- load_disease_config()
NONE_LABEL <- glue("{cfg$disease_short_caps}-None")

DEP_DIR     <- here::here("results", "endophenotype", "cluster_proteomics")
OUT_DIR     <- here::here("results", "endophenotype", "cluster_celltypes")
FIG_DIR     <- here::here("results", "figures", "5")
WALCHLI_CSV <- here::here("data", "reference", "walchli_avg_expr_by_celltype.csv")
HPA_FILE    <- here::here("data", "reference", "hpa_rna_single_cell_type.tsv.zip")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_GSEA <- 0.05   # primary significance threshold
FDR_ORA  <- 0.20   # legacy Fisher threshold retained for back-compat
NOM_P    <- 0.05   # legacy up/down module threshold

cat("Loading specificity-filtered cell-type markers...\n")
marker_sets <- load_celltype_markers(
    walchli_csv = WALCHLI_CSV,
    hpa_file    = HPA_FILE,
    method      = "specificity",
    verbose     = TRUE
)
cat("Loading legacy nTPM>=10 markers (for Fisher reference)...\n")
marker_sets_legacy <- load_celltype_markers(
    walchli_csv = WALCHLI_CSV,
    hpa_file    = HPA_FILE,
    method      = "legacy_ntpm",
    verbose     = FALSE
)

# ---------------------------------------------------------------------------
# Run per-cluster enrichment
# ---------------------------------------------------------------------------
cluster_contrasts <- c("C0_vs_None", "C1_vs_None", "C2_vs_None")
gsea_all    <- list()
camera_all  <- list()
ora_all     <- list()

for (cname in cluster_contrasts) {
    dep_file <- file.path(DEP_DIR, sprintf("cluster_deps_%s.csv", cname))
    if (!file.exists(dep_file)) {
        cat(sprintf("  Skipping %s: file not found\n", cname)); next
    }
    tt <- fread(dep_file)
    tt <- tt[!is.na(t) & !is.na(protein)]
    cluster_tag <- sub("_vs_None", "", cname)

    cat(sprintf("\n=== %s (n=%d) ===\n", cname, nrow(tt)))

    # Primary: preranked fgsea on moderated t
    rank_stat <- build_rank_stat(tt, gene_col = "protein", stat_col = "t")
    g <- run_celltype_gsea(rank_stat, marker_sets, min_size = 5, max_size = 500,
                             seed = 42)
    if (!is.null(g) && nrow(g) > 0L) {
        g[, contrast := cname]
        g[, cluster  := cluster_tag]
        fwrite(g, file.path(OUT_DIR, sprintf("celltypes_gsea_%s.csv", cname)))
        gsea_all[[cname]] <- g
    }

    # Secondary: cameraPR confirmatory
    cm <- run_celltype_camerapr(rank_stat, marker_sets,
                                  inter_gene_cor = 0.01,
                                  min_size = 5, max_size = 500)
    if (!is.null(cm) && nrow(cm) > 0L) {
        cm[, contrast := cname]
        cm[, cluster  := cluster_tag]
        fwrite(cm, file.path(OUT_DIR, sprintf("celltypes_camerapr_%s.csv", cname)))
        camera_all[[cname]] <- cm
    }

    # Legacy: Fisher on up/down modules with nominal P-value cutoff
    bg_genes    <- toupper(tt$protein)
    up_module   <- toupper(tt[P.Value < NOM_P & logFC > 0, protein])
    down_module <- toupper(tt[P.Value < NOM_P & logFC < 0, protein])
    cat(sprintf("  legacy module sizes: up=%d  down=%d (P<%.2f)\n",
                length(up_module), length(down_module), NOM_P))

    or_up <- run_ora(up_module, bg_genes, marker_sets_legacy)
    or_dn <- run_ora(down_module, bg_genes, marker_sets_legacy)
    if (!is.null(or_up)) {
        or_up[, `:=`(direction = "up", contrast = cname, cluster = cluster_tag)]
        fwrite(or_up, file.path(OUT_DIR, sprintf("celltypes_ora_%s_up.csv", cname)))
    }
    if (!is.null(or_dn)) {
        or_dn[, `:=`(direction = "down", contrast = cname, cluster = cluster_tag)]
        fwrite(or_dn, file.path(OUT_DIR, sprintf("celltypes_ora_%s_down.csv", cname)))
    }
    ora_all[[cname]] <- rbind(if (!is.null(or_up)) or_up,
                                if (!is.null(or_dn)) or_dn)
}

if (length(gsea_all) == 0L)
    stop("No GSEA results produced - check cluster DEP files in ", DEP_DIR)

gsea_summary <- rbindlist(gsea_all, use.names = TRUE, fill = TRUE)
fwrite(gsea_summary[order(padj)], file.path(OUT_DIR, "cluster_celltype_summary_gsea.csv"))
cat(sprintf("\nGSEA summary: %d rows  |  FDR<%.2f: %d  |  FDR<0.05: %d\n",
            nrow(gsea_summary), FDR_GSEA,
            sum(gsea_summary$padj < FDR_GSEA, na.rm = TRUE),
            sum(gsea_summary$padj < 0.05,     na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Supplementary table sink (etbl_cluster_celltypes in figures.yaml)
# ---------------------------------------------------------------------------
supp_table_path <- here::here(
    "manuscript", "submission", "_build_intermediates",
    "Supplementary Table_Cluster_CellTypes.csv"
)
dir.create(dirname(supp_table_path), recursive = TRUE, showWarnings = FALSE)
supp_cols <- intersect(
    c("cluster", "cell_type", "size", "NES", "ES", "pval", "padj",
      "leading_edge", "contrast"),
    names(gsea_summary)
)
fwrite(gsea_summary[order(cluster, padj), ..supp_cols], supp_table_path)
cat(sprintf("  Supplementary table: %s\n", basename(supp_table_path)))

if (length(ora_all) > 0L) {
    ora_summary <- rbindlist(ora_all, use.names = TRUE, fill = TRUE)
    fwrite(ora_summary[order(p_adj)], file.path(OUT_DIR, "cluster_celltype_summary.csv"))
}

# ---------------------------------------------------------------------------
# Per-cluster dot panels (figure 5 r/s/t) — signed NES, one col per cluster
# ---------------------------------------------------------------------------
cat("\nBuilding cell-type panels (r, s, t)...\n")
panel_specs <- list(
    list(cname = "C0_vs_None", cluster = "C0", label = "r"),
    list(cname = "C1_vs_None", cluster = "C1", label = "s"),
    list(cname = "C2_vs_None", cluster = "C2", label = "t")
)
for (spec in panel_specs) {
    dt <- gsea_all[[spec$cname]]
    if (is.null(dt) || nrow(dt) == 0L) next
    p <- make_celltype_gsea_plot(
        dt,
        contrast_col    = "cluster",
        contrast_order  = spec$cluster,
        contrast_labels = setNames(spec$cluster, spec$cluster),
        title_str       = sprintf("%s  %s \u2014 cell-type enrichment (fgsea)",
                                    spec$label, spec$cluster),
        subtitle_str    = glue(
            "Preranked fgsea on limma moderated-t (cluster vs {NONE_LABEL}) | ",
            "NES>0: markers higher in cluster | NES<0: lower | * FDR<0.05"
        )
    )
    if (!is.null(p)) {
        ggsave(file.path(FIG_DIR, sprintf("panel_%s.pdf", spec$label)),
               p, width = 4.5, height = 6.0, device = cairo_pdf)
        cat(sprintf("  panel_%s.pdf\n", spec$label))
    }
}

# ---------------------------------------------------------------------------
# Cross-cluster heatmap (figure 5 panel u) — NES grid
# ---------------------------------------------------------------------------
cat("\nBuilding cross-cluster cell-type heatmap (panel u)...\n")
p_heat <- make_celltype_heatmap(
    gsea_summary,
    fdr_thr      = FDR_GSEA,
    cluster_col  = "cluster",
    value_limits = c(-2.5, 2.5),
    title_str    = "u  Cross-cluster cell-type enrichment",
    subtitle_str = glue(
        "fgsea NES (cluster vs {NONE_LABEL}) | specificity markers | ",
        "NES>0: higher in cluster | NES<0: lower | * FDR<0.05"
    )
)
if (!is.null(p_heat)) {
    ggsave(file.path(FIG_DIR, "panel_u.pdf"), p_heat,
           width = 4.5, height = 5.5, device = cairo_pdf)
    cat("  panel_u.pdf\n")
}

cat("\n07_ms_cluster_celltypes.R complete.\n")
