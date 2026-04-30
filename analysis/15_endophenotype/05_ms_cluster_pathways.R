#!/usr/bin/env Rscript
# 05_ms_cluster_pathways.R — GO:BP enrichment for cluster-specific DEP modules
#
# Approach: identical to 04_prepost_go.R (figure 1 panel e).
#   enrichGO ORA (clusterProfiler) on nominal P<0.05 up/down DEP modules.
#   Visualisation via make_go_dotplot() helper — same style as panel e.
#
# Module threshold: P.Value < 0.05 (nominal; PSM-matched cluster vs MS-None).
# Background: all proteins with valid P.Value in that contrast (~2866).
#
# Outputs (results/endophenotype/cluster_pathways/):
#   cluster_go_<contrast>.csv   — enrichGO results per contrast
#   panel_j/k/l.pdf             — per-cluster GO dotplot (figure 1 panel-e style)
#   panel_m.pdf                 — cross-cluster GO dotplot

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "go_dotplot.R"))

DEP_DIR <- file.path(PROJ_DIR, "results", "endophenotype", "cluster_proteomics")
OUT_DIR <- file.path(PROJ_DIR, "results", "endophenotype", "cluster_pathways")
FIG_DIR <- file.path(PROJ_DIR, "results", "figures", "5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
NOM_P <- 0.05

# ── Run enrichGO per cluster contrast ────────────────────────────────────────
cat("\nRunning GO:BP enrichment per cluster...\n")
cluster_contrasts <- c("C0_vs_None", "C1_vs_None", "C2_vs_None")
all_go <- list()

for (cname in cluster_contrasts) {
    dep_file <- file.path(DEP_DIR, sprintf("cluster_deps_%s.csv", cname))
    if (!file.exists(dep_file)) { cat(sprintf("  Skipping %s (file not found)\n", cname)); next }

    tt <- fread(dep_file)
    tt <- tt[!is.na(P.Value)]

    up_syms   <- toupper(tt[P.Value < NOM_P & logFC > 0, protein])
    down_syms <- toupper(tt[P.Value < NOM_P & logFC < 0, protein])

    cat(sprintf("  %s: up=%d  down=%d  (P<%.2f)\n",
                cname, length(up_syms), length(down_syms), NOM_P))

    ego_up   <- run_go(up_syms,   sprintf("%s up",   cname))
    ego_down <- run_go(down_syms, sprintf("%s down", cname))

    go_dt <- rbindlist(Filter(Negate(is.null), list(
        collect_go(ego_up,   "up_DEPs"),
        collect_go(ego_down, "down_DEPs")
    )))
    if (nrow(go_dt) > 0L) {
        go_dt[, contrast := cname]
        go_dt[, cluster  := sub("_vs_None", "", cname)]
        fwrite(go_dt, file.path(OUT_DIR, sprintf("cluster_go_%s.csv", cname)))
    }
    all_go[[cname]] <- go_dt
}

# ── Per-cluster panels j / k / l (figure 1 panel-e style) ───────────────────
cat("\nBuilding GO dotplots (panels j, k, l)...\n")

panel_specs <- list(
    list(cname = "C0_vs_None", cluster = "C0", label = "j",
         up_label = "Up in C0 (spine/connective)",
         dn_label = "Down in C0"),
    list(cname = "C1_vs_None", cluster = "C1", label = "k",
         up_label = "Up in C1 (cranial/demyel.)",
         dn_label = "Down in C1"),
    list(cname = "C2_vs_None", cluster = "C2", label = "l",
         up_label = "Up in C2 (neurological/EBV)",
         dn_label = "Down in C2")
)

for (spec in panel_specs) {
    go_dt <- all_go[[spec$cname]]
    if (is.null(go_dt) || nrow(go_dt) == 0L) {
        cat(sprintf("  %s: no GO terms — skipping\n", spec$cname)); next
    }
    col <- CLUST_COLS[spec$cluster]

    p <- make_go_dotplot(
        go_dt        = go_dt,   # pass with original "up_DEPs"/"down_DEPs" labels
        title_str    = sprintf("%s  %s \u2014 GO:BP enrichment", spec$label, spec$cluster),
        subtitle_str = sprintf("enrichGO ORA | nominal P<%.2f module vs MS-None",  NOM_P),
        up_col     = col,
        down_col   = COL_DOWN,
        up_label   = spec$up_label,
        down_label = spec$dn_label
    )
    if (!is.null(p)) {
        ggsave(file.path(FIG_DIR, sprintf("panel_%s.pdf", spec$label)),
               p, width = 7.0, height = 6.5, device = cairo_pdf)
        cat(sprintf("  panel_%s.pdf\n", spec$label))
    }
}

# ── Panel m: cross-cluster GO dotplot ────────────────────────────────────────
cat("\nBuilding cross-cluster GO dotplot (panel m)...\n")

go_all <- rbindlist(all_go, fill = TRUE)

if (nrow(go_all) > 0L) {
    p_m <- make_cross_cluster_go_dotplot(
        go_all,
        title_str    = "m  Cross-cluster GO:BP enrichment",
        subtitle_str = sprintf("enrichGO ORA | nominal P<%.2f module | top 8 terms per cluster", NOM_P)
    )
    if (!is.null(p_m)) {
        ggsave(file.path(FIG_DIR, "panel_m.pdf"), p_m,
               width = 8.5, height = 7.0, device = cairo_pdf)
        cat("  panel_m.pdf\n")
    }
}

cat("\n05_ms_cluster_pathways.R complete.\n")
cat("Next: 06_ms_cluster_ppi.R (STRING PPI network)\n")
cat("      07_ms_cluster_celltypes.R (cell type mapping)\n")
