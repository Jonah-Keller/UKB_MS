#!/usr/bin/env Rscript
# 10_postms_cluster_celltypes.R — cell-type enrichment per post-onset cluster
#
# Post-onset mirror of 07_ms_cluster_celltypes.R. Refactored to use
# specificity-filtered markers + preranked fgsea (see helpers/celltype_overrep_plot.R).
# Legacy Fisher retained for back-compat.
#
# Outputs (results/endophenotype/postms_cluster_celltypes/):
#   celltypes_gsea_<contrast>.csv        NEW per-contrast fgsea
#   celltypes_camerapr_<contrast>.csv    NEW per-contrast cameraPR
#   celltypes_ora_<contrast>_up.csv/_down.csv   LEGACY Fisher
#   postms_cluster_celltype_summary_gsea.csv   consolidated GSEA
#   postms_cluster_celltype_summary.csv         consolidated Fisher (legacy)
#   results/figures/5S/panel_h/i/j.pdf  per-cluster NES dot plot
#   results/figures/5S/panel_k.pdf       cross-cluster NES heatmap

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

DEP_DIR     <- here::here("results", "endophenotype", "postms_cluster_proteomics")
OUT_DIR     <- here::here("results", "endophenotype", "postms_cluster_celltypes")
FIG_DIR     <- here::here("results", "figures", "5S")
WALCHLI_CSV <- here::here("data", "reference", "walchli_avg_expr_by_celltype.csv")
HPA_FILE    <- here::here("data", "reference", "hpa_rna_single_cell_type.tsv.zip")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_GSEA <- 0.05
FDR_ORA  <- 0.20
NOM_P    <- 0.05

cat("Loading specificity-filtered cell-type markers...\n")
marker_sets <- load_celltype_markers(
    walchli_csv = WALCHLI_CSV, hpa_file = HPA_FILE,
    method = "specificity", verbose = TRUE
)
cat("Loading legacy nTPM>=10 markers (for Fisher reference)...\n")
marker_sets_legacy <- load_celltype_markers(
    walchli_csv = WALCHLI_CSV, hpa_file = HPA_FILE,
    method = "legacy_ntpm", verbose = FALSE
)

cluster_contrasts <- c("C0_vs_None", "C1_vs_None", "C2_vs_None")
gsea_all   <- list()
camera_all <- list()
ora_all    <- list()

for (cname in cluster_contrasts) {
    dep_file <- file.path(DEP_DIR, sprintf("postms_deps_%s.csv", cname))
    if (!file.exists(dep_file)) {
        cat(sprintf("  Skipping %s: file not found\n", cname)); next
    }
    tt <- fread(dep_file)
    tt <- tt[!is.na(t) & !is.na(protein)]
    cluster_tag <- sub("_vs_None", "", cname)

    cat(sprintf("\n=== %s (n=%d) ===\n", cname, nrow(tt)))

    rank_stat <- build_rank_stat(tt, gene_col = "protein", stat_col = "t")

    g <- run_celltype_gsea(rank_stat, marker_sets, min_size = 5, max_size = 500,
                             seed = 42)
    if (!is.null(g) && nrow(g) > 0L) {
        g[, contrast := cname]
        g[, cluster  := cluster_tag]
        fwrite(g, file.path(OUT_DIR, sprintf("celltypes_gsea_%s.csv", cname)))
        gsea_all[[cname]] <- g
    }

    cm <- run_celltype_camerapr(rank_stat, marker_sets,
                                  inter_gene_cor = 0.01,
                                  min_size = 5, max_size = 500)
    if (!is.null(cm) && nrow(cm) > 0L) {
        cm[, contrast := cname]
        cm[, cluster  := cluster_tag]
        fwrite(cm, file.path(OUT_DIR, sprintf("celltypes_camerapr_%s.csv", cname)))
        camera_all[[cname]] <- cm
    }

    # Legacy Fisher on up/down modules
    bg_genes    <- toupper(tt$protein)
    up_module   <- toupper(tt[P.Value < NOM_P & logFC > 0, protein])
    down_module <- toupper(tt[P.Value < NOM_P & logFC < 0, protein])
    cat(sprintf("  legacy module sizes: up=%d  down=%d\n",
                length(up_module), length(down_module)))
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
fwrite(gsea_summary[order(padj)],
       file.path(OUT_DIR, "postms_cluster_celltype_summary_gsea.csv"))
cat(sprintf("\nGSEA summary: %d rows  |  FDR<%.2f: %d\n",
            nrow(gsea_summary), FDR_GSEA,
            sum(gsea_summary$padj < FDR_GSEA, na.rm = TRUE)))

if (length(ora_all) > 0L) {
    ora_summary <- rbindlist(ora_all, use.names = TRUE, fill = TRUE)
    fwrite(ora_summary[order(p_adj)],
           file.path(OUT_DIR, "postms_cluster_celltype_summary.csv"))
}

# ---------------------------------------------------------------------------
# Per-cluster dot panels (5S panels h/i/j) — signed NES
# ---------------------------------------------------------------------------
cat("\nBuilding cell-type panels (h, i, j)...\n")
panel_specs <- list(
    list(cname = "C0_vs_None", cluster = "C0", label = "h"),
    list(cname = "C1_vs_None", cluster = "C1", label = "i"),
    list(cname = "C2_vs_None", cluster = "C2", label = "j")
)
for (spec in panel_specs) {
    dt <- gsea_all[[spec$cname]]
    if (is.null(dt) || nrow(dt) == 0L) next
    p <- make_celltype_gsea_plot(
        dt,
        contrast_col    = "cluster",
        contrast_order  = spec$cluster,
        contrast_labels = setNames(sprintf("%s (post-onset)", spec$cluster),
                                    spec$cluster),
        title_str       = sprintf("%s  %s \u2014 cell-type enrichment (post-onset)",
                                    spec$label, spec$cluster),
        subtitle_str    = glue(
            "fgsea on limma moderated-t (cluster vs {NONE_LABEL}, post-onset) | ",
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
# Cross-cluster heatmap (5S panel k)
# ---------------------------------------------------------------------------
cat("\nBuilding cross-cluster cell-type heatmap (panel k)...\n")
p_heat <- make_celltype_heatmap(
    gsea_summary, fdr_thr = FDR_GSEA,
    cluster_col  = "cluster",
    value_limits = c(-2.5, 2.5),
    title_str    = "k  Cross-cluster cell-type enrichment (post-onset)",
    subtitle_str = glue(
        "fgsea NES (cluster vs {NONE_LABEL}, post-onset) | specificity markers | ",
        "NES>0: higher in cluster | NES<0: lower | * FDR<0.05"
    )
)
if (!is.null(p_heat)) {
    ggsave(file.path(FIG_DIR, "panel_k.pdf"), p_heat,
           width = 4.5, height = 5.5, device = cairo_pdf)
    cat("  panel_k.pdf\n")
}

cat("\n10_postms_cluster_celltypes.R complete.\n")
