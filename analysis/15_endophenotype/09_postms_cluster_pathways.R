#!/usr/bin/env Rscript
# 09_postms_cluster_pathways.R — GO:BP enrichment for post-onset cluster DEP modules
#
# Mirror of 05_ms_cluster_pathways.R using post-onset DEP results from
# 08_postms_cluster_proteomics.R. Output goes to Figure 5S.
#
# Outputs (results/endophenotype/postms_cluster_pathways/):
#   postms_go_<contrast>.csv     — enrichGO results per contrast
#   panel_d/e/f.pdf              — per-cluster GO dotplot
#   panel_g.pdf                  — cross-cluster GO dotplot

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
    library(glue)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))
source(here::here("analysis", "helpers", "go_dotplot.R"))

cfg <- load_disease_config()
NONE_LABEL <- glue("{cfg$disease_short_caps}-None")

DEP_DIR <- here::here("results", "endophenotype", "postms_cluster_proteomics")
OUT_DIR <- here::here("results", "endophenotype", "postms_cluster_pathways")
FIG_DIR <- here::here("results", "figures", "5S")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
NOM_P <- 0.05

# ── Run enrichGO per cluster contrast ────────────────────────────────────────
cat("\nRunning GO:BP enrichment per cluster (post-onset)...\n")
cluster_contrasts <- c("C0_vs_None", "C1_vs_None", "C2_vs_None")
all_go <- list()

for (cname in cluster_contrasts) {
    dep_file <- file.path(DEP_DIR, sprintf("postms_deps_%s.csv", cname))
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
        fwrite(go_dt, file.path(OUT_DIR, sprintf("postms_go_%s.csv", cname)))
    }
    all_go[[cname]] <- go_dt
}

# ── Per-cluster panels d / e / f ─────────────────────────────────────────────
cat("\nBuilding GO dotplots (panels d, e, f)...\n")

panel_specs <- list(
    list(cname = "C0_vs_None", cluster = "C0", label = "d",
         up_label = "Up in C0 (spine/connective)  · post-onset",
         dn_label = "Down in C0  · post-onset"),
    list(cname = "C1_vs_None", cluster = "C1", label = "e",
         up_label = "Up in C1 (cranial/demyel.)  · post-onset",
         dn_label = "Down in C1  · post-onset"),
    list(cname = "C2_vs_None", cluster = "C2", label = "f",
         up_label = "Up in C2 (neurological/EBV)  · post-onset",
         dn_label = "Down in C2  · post-onset")
)

for (spec in panel_specs) {
    go_dt <- all_go[[spec$cname]]
    if (is.null(go_dt) || nrow(go_dt) == 0L) {
        cat(sprintf("  %s: no GO terms — skipping\n", spec$cname)); next
    }
    col <- CLUST_COLS[spec$cluster]

    p <- make_go_dotplot(
        go_dt        = go_dt,
        title_str    = sprintf("%s  %s \u2014 GO:BP enrichment (post-onset)", spec$label, spec$cluster),
        subtitle_str = glue("enrichGO ORA | nominal P<{format(NOM_P, nsmall=2)} module vs {NONE_LABEL} | post-onset only"),
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

# ── Panel g: cross-cluster GO dotplot ────────────────────────────────────────
cat("\nBuilding cross-cluster GO dotplot (panel g)...\n")

go_all <- rbindlist(all_go, fill = TRUE)

if (nrow(go_all) > 0L) {
    p_g <- make_cross_cluster_go_dotplot(
        go_all,
        title_str    = "g  Cross-cluster GO:BP enrichment (post-onset)",
        subtitle_str = sprintf("enrichGO ORA | nominal P<%.2f module | top 8 terms per cluster | post-onset only", NOM_P)
    )
    if (!is.null(p_g)) {
        ggsave(file.path(FIG_DIR, "panel_g.pdf"), p_g,
               width = 8.5, height = 7.0, device = cairo_pdf)
        cat("  panel_g.pdf\n")
    }
}

cat("\n09_postms_cluster_pathways.R complete.\n")
cat("Next: 10_postms_cluster_celltypes.R (cell type mapping)\n")
