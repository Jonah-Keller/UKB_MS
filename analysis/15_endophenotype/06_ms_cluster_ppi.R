#!/usr/bin/env Rscript
# 06_ms_cluster_ppi.R — STRING PPI networks for cluster-specific DEPs
#
# Uses the same circular chord diagram as Figure 1 panel F (see
# analysis/helpers/chord_diagram.R).  Communities = clusters (C0/C1/C2)
# plus a "Shared" community for proteins strongly connected across clusters.
#
# Per-cluster panels: seed = top 30 proteins by |rank_stat| from that
# cluster's DEP contrast.  Cross-cluster panel: union of seeds, community =
# cluster with strongest association per protein.
#
# String threshold: 0.40 (medium confidence), consistent with existing pipeline.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
    library(glue)
    library(igraph)
    library(STRINGdb)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))
source(here::here("analysis", "helpers", "chord_diagram.R"))

cfg <- load_disease_config()

DEP_DIR <- here::here("results", "endophenotype",
                      glue("{cfg$cohort_short}_cluster_proteomics"))
OUT_DIR <- here::here("results", "endophenotype",
                      glue("{cfg$cohort_short}_cluster_ppi"))
FIG_DIR <- here::here("results", "figures", "5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "string_cache"), showWarnings = FALSE, recursive = TRUE)

set.seed(42)
STRING_SCORE <- 400L
N_SEED       <- 30L

# Cluster colour palette — matches figure1.R CPAL style and CLUST_COLS in figure5.R
# Community IDs: 1=C0, 2=C1, 3=C2
CPAL <- c(
    "1" = "#E6A817",   # C0 — amber
    "2" = "#CC0066",   # C1 — magenta
    "3" = "#2B4C7E",   # C2 — navy
    "0" = "#bdbdbd"    # orphan
)
CPAL_LIGHT <- c(
    "1" = "#fdf3d9",
    "2" = "#fce4ef",
    "3" = "#dfe6f0",
    "0" = "#f5f5f5"
)

# ── 1. Initialise STRINGdb ────────────────────────────────────────────────────
cat("Initialising STRINGdb (v12, human)...\n")
string_db <- STRINGdb$new(version = "12", species = 9606,
                           score_threshold = STRING_SCORE,
                           network_type = "full",
                           input_directory = file.path(OUT_DIR, "string_cache"))

# ── 2. Build per-cluster seed sets ───────────────────────────────────────────
cat("Building seed protein sets...\n")
cluster_contrasts <- c("C0_vs_None", "C1_vs_None", "C2_vs_None")
seeds <- list()

for (cname in cluster_contrasts) {
    dep_file <- file.path(DEP_DIR, sprintf("cluster_deps_%s.csv", cname))
    if (!file.exists(dep_file)) next
    tt <- fread(dep_file)
    tt[, rank_stat := sign(logFC) * -log10(P.Value)]
    seeds[[cname]] <- tt[order(-abs(rank_stat))][1:min(N_SEED, .N)]
}

# ── 3. Helper: build igraph from STRING ───────────────────────────────────────
build_string_graph <- function(protein_genes, cname) {
    prots  <- data.frame(gene = toupper(protein_genes), stringsAsFactors = FALSE)
    mapped <- string_db$map(prots, "gene", removeUnmappedRows = TRUE)
    if (nrow(mapped) < 3L) return(NULL)

    edges_raw <- string_db$get_interactions(mapped$STRING_id)
    if (is.null(edges_raw) || nrow(edges_raw) == 0L) return(NULL)

    id2gene <- setNames(mapped$gene, mapped$STRING_id)
    el <- edges_raw[, c("from","to","combined_score")]
    el$from <- id2gene[el$from]; el$to <- id2gene[el$to]
    el <- el[!is.na(el$from) & !is.na(el$to), ]

    fwrite(as.data.table(edges_raw),
           file.path(OUT_DIR, sprintf("ppi_edges_%s.csv", cname)))

    g <- graph_from_data_frame(el, directed = FALSE,
                                vertices = data.frame(name = unique(c(el$from, el$to))))
    g <- simplify(g, edge.attr.comb = list(combined_score = "max"))
    E(g)$combined_score <- E(g)$combined_score / 1000
    cat(sprintf("  %s: %d proteins → %d STRING edges\n",
                cname, vcount(g), ecount(g)))
    g
}

# ── 4. Per-cluster chord diagrams ─────────────────────────────────────────────
cat("\nBuilding per-cluster PPI chord diagrams...\n")

panel_specs <- list(
    list(cname = "C0_vs_None", cl_id = 1L, label = "n", cl_name = "C0"),
    list(cname = "C1_vs_None", cl_id = 2L, label = "o", cl_name = "C1"),
    list(cname = "C2_vs_None", cl_id = 3L, label = "p", cl_name = "C2")
)

for (spec in panel_specs) {
    tt <- seeds[[spec$cname]]
    if (is.null(tt)) next

    g <- build_string_graph(tt$protein, spec$cname)
    if (is.null(g)) {
        cat(sprintf("  %s: no graph, skipping\n", spec$cname)); next
    }

    # Direction-based community assignment:
    # up-regulated → community 1, down-regulated → community 2, NS → orphan 0
    dep_dir <- setNames(
        fifelse(toupper(tt$protein[match(V(g)$name, toupper(tt$protein))]) %in%
                    toupper(tt$protein[tt$logFC > 0]), 1L, 2L),
        V(g)$name
    )
    dep_dir[is.na(dep_dir)] <- 0L

    dir_label <- setNames(toupper(tt$protein), toupper(tt$protein))
    dir_col <- c(
        "1" = unname(CPAL[as.character(spec$cl_id)]),
        "2" = "#56B4E9",   # sky blue — down-regulated (matches COL_PRE_DOWN)
        "0" = "#bdbdbd"
    )
    dir_col_light <- c(
        "1" = unname(CPAL_LIGHT[as.character(spec$cl_id)]),
        "2" = "#e8f5fd",
        "0" = "#f5f5f5"
    )
    dir_labels <- c("1" = sprintf("%s up", spec$cl_name),
                    "2" = sprintf("%s down", spec$cl_name))

    p <- make_chord_ppi(
        g              = g,
        comm_vec       = dep_dir,
        comm_cols      = dir_col,
        comm_cols_light = dir_col_light,
        comm_labels    = dir_labels,
        max_nodes      = 50L,
        title_str      = sprintf("%s  %s protein interaction network", spec$label, spec$cl_name),
        subtitle_str   = sprintf("STRING v12 ≥%d | top %d seeds | up = %s, down = sky blue",
                                 STRING_SCORE, N_SEED, spec$cl_name)
    )

    out <- file.path(FIG_DIR, sprintf("panel_%s.pdf", spec$label))
    ggsave(out, p, width = 7.5, height = 7.0, device = cairo_pdf)
    cat(sprintf("  panel_%s.pdf\n", spec$label))
}

# ── 5. Cross-cluster chord diagram ────────────────────────────────────────────
cat("\nBuilding cross-cluster PPI chord diagram...\n")

all_seeds <- rbindlist(lapply(cluster_contrasts, function(cn) {
    dt <- seeds[[cn]]
    if (is.null(dt)) return(NULL)
    dt[, cluster_id := match(sub("_vs_None","",cn), c("C0","C1","C2"))]
    dt[, .(gene = toupper(protein), logFC, P.Value, cluster_id)]
}))
all_seeds <- all_seeds[, .SD[which.max(abs(logFC))], by = gene]

g_cross <- build_string_graph(all_seeds$gene, "cross_cluster")

if (!is.null(g_cross)) {
    comm_cross <- setNames(
        all_seeds$cluster_id[match(V(g_cross)$name, all_seeds$gene)],
        V(g_cross)$name
    )
    comm_cross[is.na(comm_cross)] <- 0L

    cross_labels <- c("1" = "C0  spine/connective",
                      "2" = "C1  cranial/demyel.",
                      "3" = "C2  neurol/EBV")

    p_cross <- make_chord_ppi(
        g               = g_cross,
        comm_vec        = comm_cross,
        comm_cols       = CPAL,
        comm_cols_light = CPAL_LIGHT,
        comm_labels     = cross_labels,
        max_nodes       = 60L,
        title_str       = "q  Cross-cluster PPI network",
        subtitle_str    = sprintf("STRING v12 ≥%d | colour = cluster with strongest signal | n=%d proteins",
                                  STRING_SCORE, vcount(g_cross))
    )
    ggsave(file.path(FIG_DIR, "panel_q.pdf"), p_cross,
           width = 8.5, height = 8.0, device = cairo_pdf)
    cat("  panel_q.pdf\n")
}

cat("\n06_ms_cluster_ppi.R complete.\n")
