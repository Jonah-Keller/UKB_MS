#!/usr/bin/env Rscript
# figure1.R — Figure 1: The pre-onset MS plasma proteomic signature
#
# Framework: 4 analyses × 5 visualization types.
#
# Four analyses:
#   (1) combined   — All MS vs HC (combined pre + post onset)
#   (2) pre        — Pre-onset MS vs HC  [PROMOTED TO MAIN FIGURE]
#   (3) post       — Post-onset MS vs HC
#   (4) preVsPost  — Direct pre vs post-onset (MS-only limma)
#
# Five visualization types (run for each analysis):
#   V1. Volcano plot
#   V2. Ranked barplot (top 20 by nominal p-value — pre-onset has few FDR<0.05)
#   V3. GO:BP enrichment dotplot (up + down split)
#   V4. HPA scRNA-seq cell-type expression heatmap (pheatmap)
#   V5. Cell-type over-representation test (Fisher's exact, dot plot)
#
# Main figure panels (pre-onset analysis):
#   a  Study design schematic (placeholder)
#   b  UMAP cohort overview
#   c  Pre-onset volcano
#   d  Pre-onset ranked barplot (top 20 nominal p)
#   e  Pre-onset GO:BP enrichment (up + down)
#   f  Pre-onset PPI network (STRING ≥ 0.7, Louvain communities)
#   g  Pre-onset HPA cell-type expression heatmap
#   h  Pre-onset cell-type over-representation test
#   i  Pre vs post GO pathway comparison heatmap (both analyses)
#
# Supplementary (all 4 analyses, same 5 viz types) → figure1_supp.R
#
# Output: results/figures/1/panel_{a..i}.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(limma)
    library(igraph)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "limma_utils.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "celltype_overrep_plot.R"))

FIG_DIR  <- file.path(PROJ_DIR, "results", "figures", "1")
DATA_DIR <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed")
DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
NET_DIR  <- file.path(PROJ_DIR, "results", "networks")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# All colour constants come from ukb_theme.R (sourced above).
# Stage identity (UMAP, trajectories): COL_PRE / COL_POST / COL_COMB
# Directional DEP colours:             COL_PRE_UP/DOWN, COL_POST_UP/DOWN, COL_COMB_UP/DOWN

# ---------------------------------------------------------------------------
# Helper: save panel
# ---------------------------------------------------------------------------
save_panel <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(FIG_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: panel_", name, ".pdf\n", sep = "")
}

# ---------------------------------------------------------------------------
# Helper: reusable volcano builder
# ---------------------------------------------------------------------------
# Labeling policy (applied to every volcano):
#   1. Must-label set (never dropped by ggrepel) = union of
#        (a) proteins passing Bonferroni correction (P < 0.05/n),
#        (b) explicit `highlight` vector of symbols named in manuscript prose.
#   2. Remaining label slots (up to n_label) are filled by P-value ranking
#      among FDR < fdr_line proteins.
#   Must-label layer uses max.overlaps=Inf so names we commit to in text
#   are always visible; fill layer uses a normal max.overlaps.
make_volcano <- function(df, title_str, n_label = 15, fdr_line = 0.05,
                         up_col = COL_PRE_UP, down_col = COL_PRE_DOWN,
                         up_label = "Up", down_label = "Down",
                         highlight = character()) {
    df <- copy(df)
    if ("adj.P.Val" %in% names(df) && !"fdr" %in% names(df))
        df[, fdr := adj.P.Val]
    df[, direction := fcase(
        fdr < fdr_line & logFC > 0, "Up",
        fdr < fdr_line & logFC < 0, "Down",
        default = "NS"
    )]
    df[, direction := factor(direction, levels = c("Up", "Down", "NS"))]
    df[, nlp := -log10(P.Value)]
    df[, protein_upper := toupper(protein)]
    bonf          <- -log10(0.05 / nrow(df))
    highlight_up  <- toupper(highlight)

    must_label <- df[nlp >= bonf | protein_upper %in% highlight_up]
    fill_budget <- max(0L, n_label - nrow(must_label))
    fill <- df[fdr < fdr_line &
               !(protein_upper %in% must_label$protein_upper)][
        order(P.Value)][seq_len(min(fill_budget, .N))]

    n_up   <- sum(df$direction == "Up",   na.rm = TRUE)
    n_down <- sum(df$direction == "Down", na.rm = TRUE)
    ggplot(df, aes(x = logFC, y = nlp, colour = direction)) +
        geom_point(size = 0.55, alpha = 0.55) +
        geom_hline(yintercept = bonf, linetype = "dashed",
                   linewidth = 0.3, colour = "grey45") +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
        # Must-label: Bonferroni-passing + manuscript-named proteins, never dropped
        geom_text_repel(
            data = must_label, aes(label = protein_upper),
            size = 2.0, fontface = "bold", max.overlaps = Inf,
            segment.size = 0.22, segment.colour = "grey40",
            box.padding = 0.3, min.segment.length = 0,
            show.legend = FALSE
        ) +
        # Fill: remaining FDR-sig, ranked by P-value
        geom_text_repel(
            data = fill, aes(label = protein_upper),
            size = 2.0, fontface = "bold", max.overlaps = 15,
            segment.size = 0.22, segment.colour = "grey55",
            box.padding = 0.28, show.legend = FALSE
        ) +
        scale_colour_manual(
            values = c("Up" = up_col, "Down" = down_col, "NS" = "grey65"),
            labels = c("Up"   = sprintf("%s (n=%d)", up_label,   n_up),
                       "Down" = sprintf("%s (n=%d)", down_label, n_down),
                       "NS"   = "NS"),
            name = NULL
        ) +
        annotate("text",
                 x = max(df$logFC, na.rm = TRUE) * 0.92,
                 y = bonf + 0.4, label = "Bonferroni",
                 size = 2.1, hjust = 1, colour = "grey45") +
        labs(title = title_str,
             x = expression(log[2]~"fold change"),
             y = expression(-log[10]~italic(P))) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 7))
}

# ---------------------------------------------------------------------------
# Helper: ranked barplot by nominal p-value (top N per direction)
# ---------------------------------------------------------------------------
make_ranked_barplot <- function(df, title_str, n_each = 20,
                                up_col = COL_PRE_UP, down_col = COL_PRE_DOWN,
                                up_label = "Up in MS", down_label = "Down in MS") {
    df <- copy(df)
    if ("adj.P.Val" %in% names(df) && !"fdr" %in% names(df))
        df[, fdr := adj.P.Val]
    # Rank by nominal p-value (most useful for pre-onset where FDR<0.05 rare)
    df[, direction := fifelse(logFC > 0, up_label, down_label)]
    df[, protein_upper := toupper(protein)]

    top_up   <- df[logFC > 0][order(P.Value)][seq_len(min(n_each, .N))]
    top_down <- df[logFC < 0][order(P.Value)][seq_len(min(n_each, .N))]
    bar <- rbind(top_up, top_down)
    bar[, protein_f := factor(protein_upper, levels = protein_upper[order(logFC)])]

    # Significance annotation
    bar[, sig_str := fcase(
        fdr < 0.001, "***",
        fdr < 0.01,  "**",
        fdr < 0.05,  "*",
        P.Value < 0.05, ".",
        default = ""
    )]

    ggplot(bar, aes(x = logFC, y = protein_f, fill = direction)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = sig_str,
                      x = ifelse(logFC > 0, logFC + 0.01, logFC - 0.01),
                      hjust = ifelse(logFC > 0, 0, 1)),
                  size = 2.5, colour = "grey30") +
        geom_vline(xintercept = 0, linewidth = 0.3) +
        scale_fill_manual(
            values = c(setNames(up_col, up_label), setNames(down_col, down_label)),
            name = NULL
        ) +
        labs(title = title_str,
             subtitle = "Top 20 per direction by nominal p | * FDR<0.05  ** <0.01  *** <0.001",
             x = expression(log[2]~"fold change"), y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y  = element_text(size = 6),
              legend.position = "bottom",
              plot.subtitle = element_text(size = 6.5, colour = "grey40"))
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat("Loading data...\n")
qc        <- fread(file.path(DATA_DIR, "ms_olink_qc.csv"))
comb      <- fread(file.path(DIFF_DIR, "ms_combined_vs_hc.csv"))
pre_diff  <- fread(file.path(DIFF_DIR, "ms_pre_vs_hc.csv"))
post_diff <- fread(file.path(DIFF_DIR, "ms_post_vs_hc.csv"))

# Pre/post GO (from 04_prepost_go.R)
PRE_GO_FILE  <- file.path(NET_DIR, "ms_pre_go_results.csv")
POST_GO_FILE <- file.path(NET_DIR, "ms_post_go_results.csv")
go_pre  <- if (file.exists(PRE_GO_FILE))  fread(PRE_GO_FILE)  else data.table()
go_post <- if (file.exists(POST_GO_FILE)) fread(POST_GO_FILE) else data.table()

# STRING edges and cell-type connectivity
EDGE_FILE <- file.path(NET_DIR, "ms_string_edges.csv")
edges     <- if (file.exists(EDGE_FILE)) fread(EDGE_FILE) else data.table()

# Walchli+HPA integrated cell-type matrix for pre-onset (primary panel g)
WALCHLI_HPA_PRE_FILE <- file.path(NET_DIR, "ms_pre_walchli_hpa_matrix.csv")
# Fallback to HPA-only if integrated not yet generated
HPA_PRE_FILE <- file.path(NET_DIR, "ms_pre_hpa_celltype_matrix.csv")
HPA_MAT_FILE <- if (file.exists(HPA_PRE_FILE)) HPA_PRE_FILE else
                    file.path(NET_DIR, "ms_hpa_celltype_matrix.csv")

# Cell-type enrichment (new preranked-GSEA pipeline — specificity-filtered markers)
# Primary source: ms_celltype_gsea.csv (NES, padj). Falls back to the legacy
# Fisher CSV for backwards compatibility with pre-refactor runs.
GSEA_FILE    <- file.path(NET_DIR, "ms_celltype_gsea.csv")
OVERREP_FILE <- file.path(NET_DIR, "ms_celltype_overrep.csv")
gsea    <- if (file.exists(GSEA_FILE))    fread(GSEA_FILE)    else data.table()
overrep <- if (file.exists(OVERREP_FILE)) fread(OVERREP_FILE) else data.table()

n_pre  <- sum(qc$ms_status == "pre_onset",  na.rm = TRUE)
n_post <- sum(qc$ms_status == "post_onset", na.rm = TRUE)
n_hc   <- sum(qc$ms_status == "control",    na.rm = TRUE)
cat(sprintf("  QC: %d HC, %d pre-onset, %d post-onset\n", n_hc, n_pre, n_post))
cat(sprintf("  Pre-onset GO terms: %d | Post-onset GO terms: %d\n",
            nrow(go_pre), nrow(go_post)))

# ---------------------------------------------------------------------------
# Panel a — Study design schematic (placeholder)
# ---------------------------------------------------------------------------
cat("Building panel a...\n")
pA <- ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "grey96", colour = "grey70", linewidth = 0.5) +
    annotate("text", x = 0.5, y = 0.56,
             label = "Study design schematic",
             size = 5, colour = "grey40", fontface = "bold") +
    annotate("text", x = 0.5, y = 0.44,
             label = "Replace with illustration",
             size = 3.5, colour = "grey60") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = "a  Study design") +
    theme_void() +
    theme(plot.title = element_text(size = 9),
          plot.margin = margin(4, 4, 4, 4))
save_panel(pA, "a_study_design", 5, 3.5)

# ---------------------------------------------------------------------------
# Panel b — Control-derived PC1/PC2: two-facet (full cohort + MS temporal split)
# These are the same PCs used as covariates in the differential expression
# models, making panel b a direct visualization of the technical-variance
# subspace that is regressed out.
# ---------------------------------------------------------------------------
cat("Building panel b...\n")
META_COLS_B <- c("eid","ms_status","age_at_sampling","age_at_diagnosis",
                 "years_to_diagnosis","sex","olink_instance","qc_outlier",
                 "UMAP1","UMAP2","mean_npx")
prot_cols_b <- setdiff(names(qc), META_COLS_B)
qc_nonout   <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]
pc_dt_b     <- compute_control_pcs(qc_nonout, prot_cols_b)
pc_df <- merge(qc_nonout[, .(eid, ms_status)], pc_dt_b, by = "eid")
pc_df[, ms_status := factor(ms_status, levels = c("control", "post_onset", "pre_onset"))]
pc_df[, group_left := fifelse(ms_status == "control", "HC", "MS")]
pc_df[, group_left := factor(group_left, levels = c("HC", "MS"))]

pB_left <- ggplot() +
    geom_point(data = pc_df[group_left == "HC"],
               aes(PC1, PC2, colour = "HC"),
               size = 0.35, alpha = 0.25, stroke = 0) +
    geom_point(data = pc_df[group_left == "MS"],
               aes(PC1, PC2, colour = "All MS"),
               size = 0.55, alpha = 0.80, stroke = 0) +
    scale_colour_manual(
        values = c("HC" = "grey72", "All MS" = COL_COMB),
        labels = c("HC" = sprintf("HC  (n=%s)", format(sum(qc$ms_status == "control"), big.mark = ",")),
                   "All MS" = sprintf("MS  (n=%s)", format(n_pre + n_post, big.mark = ","))),
        name = NULL
    ) +
    labs(title = "All samples", x = "PC1", y = "PC2") +
    theme_ukb(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "bottom", legend.text = element_text(size = 7)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 2))

pB_right <- ggplot() +
    geom_point(data = pc_df[ms_status == "control"],
               aes(PC1, PC2),
               colour = "grey72", size = 0.35, alpha = 0.20, stroke = 0) +
    geom_point(data = pc_df[ms_status != "control"],
               aes(PC1, PC2, colour = ms_status),
               size = 0.65, alpha = 0.85, stroke = 0) +
    scale_colour_manual(
        values = c("pre_onset" = COL_PRE, "post_onset" = COL_POST),
        labels = c("pre_onset"  = sprintf("Pre-onset (n=%d)",  n_pre),
                   "post_onset" = sprintf("Post-onset (n=%d)", n_post)),
        name = NULL
    ) +
    labs(title = "MS cases by disease stage", x = "PC1", y = "PC2") +
    theme_ukb(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 2))

pB <- (pB_left | pB_right) +
    plot_annotation(title = "b  Plasma proteomic landscape (control-derived PC1/PC2)")
save_panel(pB, "b_umap_cohort", 7, 3.2)

# ---------------------------------------------------------------------------
# Panel c — Pre-onset volcano (MAIN: promoted from panel j)
# ---------------------------------------------------------------------------
cat("Building panel c...\n")
# Proteins named in the Results prose — guaranteed to be labeled regardless
# of |logFC| rank or ggrepel crowding. Keep in sync with manuscript.Rmd.
PANEL_C_HIGHLIGHT <- c("NEFL", "ERBB2",
                        "OSM", "TGFA", "MMP9",
                        "PTPRC", "BGN", "OMG")
pC <- make_volcano(
    pre_diff,
    "c  Pre-onset MS vs HC",
    fdr_line  = 0.20,
    up_label  = "Up in MS", down_label = "Down in MS",
    highlight = PANEL_C_HIGHLIGHT
)
save_panel(pC, "c_volcano_pre", 3.5, 4)

# ---------------------------------------------------------------------------
# Panel d — Pre-onset ranked barplot (top 20 by nominal p per direction)
# ---------------------------------------------------------------------------
cat("Building panel d...\n")
pD <- make_ranked_barplot(
    pre_diff,
    "d  Pre-onset DEPs ranked by nominal p"
)
save_panel(pD, "d_barplot_pre", 5.5, 5.5)

# ---------------------------------------------------------------------------
# Panel e — GO:BP enrichment: pre-onset DEPs (up + down split)
# Source: ms_pre_go_results.csv (04_prepost_go.R)
# ---------------------------------------------------------------------------
cat("Building panel e...\n")
tryCatch({
    if (nrow(go_pre) == 0) stop("ms_pre_go_results.csv empty — run 04_prepost_go.R")

    pre_up   <- go_pre[gene_set == "up_DEPs"   & p.adjust < 0.05]
    pre_down <- go_pre[gene_set == "down_DEPs" & p.adjust < 0.05]
    if (nrow(pre_up)   > 0) pre_up[,   direction := "Up in pre-onset MS"]
    if (nrow(pre_down) > 0) pre_down[, direction := "Down in pre-onset MS"]

    top_e <- rbindlist(Filter(Negate(is.null), list(
        if (nrow(pre_up)   > 0) pre_up[order(p.adjust)][seq_len(min(10, .N))]   else NULL,
        if (nrow(pre_down) > 0) pre_down[order(p.adjust)][seq_len(min(10, .N))] else NULL
    )))
    if (nrow(top_e) == 0) stop("No GO terms to plot")

    top_e[, direction := factor(direction,
                                levels = c("Up in pre-onset MS", "Down in pre-onset MS"))]
    top_e[, Description := factor(Description,
                                  levels = top_e[order(direction, FoldEnrichment), Description])]

    n_pre_up   <- sum(pre_diff$adj.P.Val < 0.20 & pre_diff$logFC > 0, na.rm = TRUE)
    n_pre_down <- sum(pre_diff$adj.P.Val < 0.20 & pre_diff$logFC < 0, na.rm = TRUE)

    pE <- ggplot(top_e, aes(x = FoldEnrichment, y = Description,
                             size = Count, colour = direction)) +
        geom_point() +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
        scale_colour_manual(
            values = c("Up in pre-onset MS" = COL_PRE_UP, "Down in pre-onset MS" = COL_PRE_DOWN),
            guide  = guide_legend(override.aes = list(size = 3.5))
        ) +
        scale_size_continuous(name = "Gene count", range = c(2, 6),
                              guide = guide_legend(override.aes = list(colour = "grey40"))) +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        labs(title    = "e  GO:BP enrichment — pre-onset MS DEPs",
             subtitle = sprintf("Pre-onset DEPs: %d up, %d down",
                                n_pre_up, n_pre_down),
             x = "Fold enrichment", y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y   = element_text(size = 7),
              strip.text    = element_text(size = 8, face = "bold"),
              legend.position = "right",
              plot.subtitle = element_text(size = 7, colour = "grey40"))
    save_panel(pE, "e_go_pre", 5.5, 6)
}, error = function(e) cat("  Panel e skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Panel f — PPI network (STRING ≥ 0.7, Louvain communities) — pre-onset DEPs
# NOTE: requires ms_string_edges.csv generated by 01_ppi_network.R.
# Currently using combined DEP network; re-run with pre-onset DEPs for final figure.
# ---------------------------------------------------------------------------
cat("Building panel f...\n")
tryCatch({
    if (nrow(edges) == 0) stop("ms_string_edges.csv not found — run 01_ppi_network.R")

    suppressPackageStartupMessages(library(enrichR))

    pre_sig <- pre_diff[adj.P.Val < 0.20]
    pre_sig[, direction_f := fifelse(logFC > 0, "Up", "Down")]
    dir_map_pre <- setNames(pre_sig$direction_f, toupper(pre_sig$protein))

    edge_filt <- edges[combined_score >= 0.7]
    edge_filt[, gene_a := toupper(gene_a)]
    edge_filt[, gene_b := toupper(gene_b)]

    # Keep only edges where both nodes are pre-onset DEPs (if enough exist)
    pre_genes <- toupper(pre_sig$protein)
    edge_pre  <- edge_filt[gene_a %in% pre_genes & gene_b %in% pre_genes]

    # Fall back to full edge set filtered to pre-onset nodes (may include some non-DEP neighbours)
    if (nrow(edge_pre) < 10) {
        cat("  Few pre-onset edges; using all edges with pre-onset nodes\n")
        edge_pre <- edge_filt[gene_a %in% pre_genes | gene_b %in% pre_genes]
    }

    g   <- graph_from_data_frame(edge_pre[, .(gene_a, gene_b, combined_score)],
                                  directed = FALSE)
    g   <- simplify(g, edge.attr.comb = list(combined_score = "max"))
    deg <- degree(g)
    if (vcount(g) > 50) {
        g   <- induced_subgraph(g, names(sort(deg, decreasing = TRUE))[1:50])
        deg <- degree(g)
    }
    cat(sprintf("  Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

    # Community detection — only show communities with ≥ 3 members to avoid slivers
    set.seed(42)
    lou <- cluster_louvain(g)
    mem <- membership(lou)
    csz <- sort(table(mem), decreasing = TRUE)
    csz <- csz[csz >= 3]                        # drop tiny communities
    n_c <- min(5, length(csz))
    top_c <- as.integer(names(csz)[seq_len(n_c)])
    cid   <- setNames(rep(0L, vcount(g)), V(g)$name)
    for (i in seq_len(n_c)) cid[names(mem)[mem == top_c[i]]] <- i

    # Direction labels
    dir_map_all <- dir_map_pre[V(g)$name]
    dir_map_all[is.na(dir_map_all)] <- "NS"

    # EnrichR community labels
    setEnrichrSite("Enrichr")
    shorten_path <- function(s, nw = 3) {
        s <- gsub(" \\(.*\\)$", "", s)
        w <- strsplit(s, " ")[[1]]
        paste(w[seq_len(min(nw, length(w)))], collapse = " ")
    }
    clabels <- setNames(paste0("M", seq_len(n_c)), as.character(seq_len(n_c)))
    for (ci in seq_len(n_c)) {
        gns <- names(cid)[cid == ci]
        if (length(gns) < 3) next
        lbl <- tryCatch({
            er <- enrichr(gns, "KEGG_2021_Human")[[1]]
            if (is.data.frame(er) && nrow(er) > 0 && er$Adjusted.P.value[1] < 0.2)
                shorten_path(er$Term[1])
            else paste0("M", ci)
        }, error = function(e) paste0("M", ci))
        clabels[as.character(ci)] <- lbl
    }

    CPAL <- c("1" = "#e6ab02", "2" = "#fd8d3c", "3" = "#74c476",
              "4" = "#66c2a5", "5" = "#9e9ac8", "0" = "#bdbdbd")
    # Node fill: tinted white using the community's background colour so nodes
    # read cleanly on the sector polygon without carrying a direction colour.
    # Each node gets a very light version of its community colour as fill.
    CPAL_LIGHT <- c("1" = "#fdf6d3", "2" = "#fde8d5", "3" = "#e5f7e5",
                    "4" = "#e0f5f0", "5" = "#ede9f7", "0" = "#f5f5f5")

    R_RING <- 3.0; GAP <- 0.012   # ring radius; tighter gaps → complete circle
    node_dt <- data.table(
        name      = V(g)$name,
        comm      = cid[V(g)$name],
        degree    = deg,
        direction = dir_map_all
    )
    setorder(node_dt, comm, -degree)

    # N_top = nodes in top communities; used for arc-fraction so the ring fills
    # exactly 360° with no empty gap from comm-0 (orphan) nodes.
    # pos is sized for all nodes; comm-0 entries remain NULL → dropped via NA filter.
    N_top <- sum(node_dt$comm > 0)
    avail <- 1.0 - GAP * n_c
    pos   <- vector("list", nrow(node_dt)); names(pos) <- node_dt$name
    cang  <- vector("list", n_c)
    cur   <- 0.0

    for (ci in seq_len(n_c)) {
        ii <- which(node_dt$comm == ci); n_ci <- length(ii); if (n_ci == 0) next
        frac <- avail * n_ci / N_top
        sa <- cur * 2 * pi - pi / 2; ea <- (cur + frac) * 2 * pi - pi / 2
        for (j in seq_len(n_ci)) {
            a <- (cur + frac * (j - 0.5) / n_ci) * 2 * pi - pi / 2
            pos[[node_dt$name[ii[j]]]] <- c(R_RING * cos(a), R_RING * sin(a))
        }
        cang[[ci]] <- list(s = sa, e = ea, m = (sa + ea) / 2)
        cur <- cur + frac + GAP
    }

    node_dt[, x := vapply(pos[name], function(p) if (is.null(p)) NA_real_ else p[1], numeric(1))]
    node_dt[, y := vapply(pos[name], function(p) if (is.null(p)) NA_real_ else p[2], numeric(1))]
    node_dt <- node_dt[!is.na(x)]
    node_dt[, comm_col       := CPAL[as.character(comm)]]
    node_dt[, comm_col_light := CPAL_LIGHT[as.character(comm)]]

    sector_list <- lapply(seq_len(n_c), function(ci) {
        if (is.null(cang[[ci]])) return(NULL)
        th <- seq(cang[[ci]]$s, cang[[ci]]$e, length.out = 55)
        lo <- R_RING - 0.55; hi <- R_RING + 0.55
        data.table(x = c(hi * cos(th), lo * cos(rev(th))),
                   y = c(hi * sin(th), lo * sin(rev(th))), ci = ci)
    })

    el  <- as_edgelist(g)
    edt <- data.table(u = el[, 1], v = el[, 2], s = E(g)$combined_score)
    edt[, cu := cid[u]]; edt[, cv := cid[v]]
    edt[, within := (cu == cv & cu > 0)]
    slo <- min(edt$s); shi <- max(edt$s); sR <- max(shi - slo, 1e-3)

    intra_list <- lapply(seq_len(n_c), function(ci) {
        sub <- edt[within == TRUE & cu == ci]; if (nrow(sub) == 0) return(NULL)
        rbindlist(lapply(seq_len(nrow(sub)), function(i) {
            pu <- pos[[sub$u[i]]]; pv <- pos[[sub$v[i]]]
            data.table(x = pu[1], y = pu[2], xend = pv[1], yend = pv[2],
                       lw = 0.5 + 1.8 * (sub$s[i] - slo) / sR, ci = ci)
        }))
    })

    bzier <- function(p0, p1, p2, p3, n = 80) {
        t <- seq(0, 1, length.out = n)
        data.frame(
            x = (1-t)^3*p0[1] + 3*(1-t)^2*t*p1[1] + 3*(1-t)*t^2*p2[1] + t^3*p3[1],
            y = (1-t)^3*p0[2] + 3*(1-t)^2*t*p1[2] + 3*(1-t)*t^2*p2[2] + t^3*p3[2]
        )
    }

    inter_pairs <- edt[within == FALSE & cu > 0 & cv > 0,
                       .(tot = sum(s)), by = .(c1 = pmin(cu, cv), c2 = pmax(cu, cv))]
    mb <- max(inter_pairs$tot, 1); Ra <- R_RING - 0.35

    chord_list <- lapply(seq_len(nrow(inter_pairs)), function(k) {
        c1 <- inter_pairs$c1[k]; c2 <- inter_pairs$c2[k]
        if (is.null(cang[[c1]]) || is.null(cang[[c2]])) return(NULL)
        m1 <- cang[[c1]]$m; m2 <- cang[[c2]]$m
        x1 <- Ra*cos(m1); y1 <- Ra*sin(m1); x2 <- Ra*cos(m2); y2 <- Ra*sin(m2)
        cx <- (x1+x2)*0.08; cy <- (y1+y2)*0.08; dx <- x2-x1; dy <- y2-y1
        bz <- bzier(c(x1,y1), c(cx-dy*0.04, cy+dx*0.04),
                    c(cx+dy*0.04, cy-dx*0.04), c(x2,y2))
        bz$nb <- inter_pairs$tot[k] / mb; bz$c1 <- c1; bz$pid <- k; bz
    })

    R_LBL <- R_RING + 1.20
    clab_dt <- rbindlist(lapply(seq_len(n_c), function(ci) {
        if (is.null(cang[[ci]])) return(NULL)
        m <- cang[[ci]]$m; lx <- R_LBL*cos(m); ly <- R_LBL*sin(m)
        data.table(x = lx, y = ly, label = clabels[as.character(ci)],
                   hjust = if (cos(m) > 0.15) 0 else if (cos(m) < -0.15) 1 else 0.5,
                   col = CPAL[as.character(ci)])
    }))

    pF <- ggplot() + coord_equal(clip = "off") + theme_void() +
        theme(plot.title = element_text(size = 9),
              plot.subtitle = element_text(size = 7, colour = "grey40"),
              legend.position = "none",
              plot.margin = margin(40, 90, 40, 60, "pt"))

    th_f   <- seq(0, 2*pi, length.out = 300)
    ring_d <- data.frame(x = R_RING*cos(th_f), y = R_RING*sin(th_f))

    for (ci in seq_len(n_c)) {
        sec <- sector_list[[ci]]; if (is.null(sec)) next
        pF <- pF + geom_polygon(data = sec, aes(x=x, y=y),
                                fill = CPAL[as.character(ci)], colour = NA, alpha = 0.30)
    }
    pF <- pF + geom_path(data = ring_d, aes(x=x, y=y),
                         colour = "#d8d8d8", linewidth = 0.35, linetype = "dashed")
    for (k in seq_along(chord_list)) {
        bz <- chord_list[[k]]; if (is.null(bz)) next
        pF <- pF + geom_path(data = bz, aes(x=x, y=y),
                             colour = CPAL[as.character(bz$c1[1])],
                             linewidth = 1.0 + 4.0*bz$nb[1],
                             alpha = 0.20 + 0.22*bz$nb[1])
    }
    for (ci in seq_len(n_c)) {
        seg <- intra_list[[ci]]; if (is.null(seg)) next
        pF <- pF + geom_segment(data = seg, aes(x=x, y=y, xend=xend, yend=yend),
                                colour = CPAL[as.character(ci)], linewidth = 0.85, alpha = 0.55)
    }
    pF <- pF +
        geom_point(data = node_dt, aes(x=x, y=y, fill=comm_col_light),
                   shape = 21, size = 6.5, stroke = 0.65,
                   colour = "grey20", alpha = 0.95) +
        scale_fill_identity() +
        geom_text(data = node_dt, aes(x=x, y=y, label=name),
                  size = 1.55, colour = "grey10", fontface = "bold")
    for (i in seq_len(nrow(clab_dt))) {
        pF <- pF + annotate("text", x=clab_dt$x[i], y=clab_dt$y[i],
                            label=clab_dt$label[i], colour=clab_dt$col[i],
                            hjust=clab_dt$hjust[i], size=2.4, fontface="bold", lineheight=0.85)
    }
    pF <- pF + labs(
        title    = "f  Pre-onset MS DEP protein interaction network",
        subtitle = sprintf("%d proteins | %d interactions | STRING ≥ 0.7 | Louvain communities",
                           vcount(g), ecount(g))
    )
    save_panel(pF, "f_network_pre", 9.5, 9.0)
}, error = function(e) cat("  Panel f skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Panel g — Walchli brain + HPA immune integrated cell-type heatmap (pre-onset DEPs)
# Rows = 16 cell types (Vascular → CNS → Immune).
# Columns = top 40 pre-onset DEPs, FIRST split into up-regulated and down-regulated
# blocks (primary grouping on the proteomic signal), then within each block ordered
# by peak cell-type compartment (Vascular → CNS → Immune) and then by peak cell
# type (in display order), with |logFC| as final tiebreaker. A validation CSV
# recording the peak cell type for each protein is written alongside the PDF so
# the Vascular/stromal-up vs. Immune-down pattern can be inspected directly.
# Run 07_walchli_hpa_integrated_heatmap.R to regenerate the input matrix.
# TODO (Minerva): If Walchli CSV needs regeneration from full Seurat object
#   (~14 GB), submit analysis/hpc/extract_walchli_avgexpr.lsf to Minerva.
# ---------------------------------------------------------------------------
cat("Building panel g...\n")
tryCatch({
    suppressPackageStartupMessages({
        library(pheatmap)
        library(RColorBrewer)
    })

    g_mat_file <- if (file.exists(WALCHLI_HPA_PRE_FILE)) WALCHLI_HPA_PRE_FILE else HPA_MAT_FILE
    if (!file.exists(g_mat_file))
        stop("Run 07_walchli_hpa_integrated_heatmap.R (or 05_hpa_celltype_heatmap.R) first")

    use_integrated <- file.exists(WALCHLI_HPA_PRE_FILE)
    cat(sprintf("  Using %s matrix\n",
                if (use_integrated) "Walchli+HPA integrated" else "HPA-only"))

    wh_mat  <- fread(g_mat_file)
    ct_cols <- setdiff(names(wh_mat), c("gene", "logFC", "direction", "abs_logFC"))

    top40_genes <- wh_mat[!is.na(logFC)][order(-abs(logFC))][seq_len(min(40, .N)), gene]
    mat_g   <- wh_mat[gene %in% top40_genes, ..ct_cols]
    genes_g <- wh_mat[gene %in% top40_genes, gene]
    mat_g   <- mat_g[match(top40_genes, genes_g)]

    expr_mat_g <- as.matrix(mat_g)
    rownames(expr_mat_g) <- top40_genes
    expr_mat_g[is.na(expr_mat_g)] <- 0

    dep_dir_g  <- wh_mat[gene %in% top40_genes, .(gene, direction, logFC)]
    dep_dir_g  <- dep_dir_g[match(top40_genes, dep_dir_g$gene)]

    unique_dirs <- unique(dep_dir_g$direction[!is.na(dep_dir_g$direction)])
    dir_col_vals <- c()
    for (d in unique_dirs) {
        dir_col_vals[d] <- if (grepl("Up|Higher pre", d, ignore.case = TRUE)) COL_PRE_UP else COL_PRE_DOWN
    }
    up_dirs   <- grep("Up|Higher pre",   unique_dirs, value = TRUE)
    down_dirs <- grep("Down|Higher post", unique_dirs, value = TRUE)

    # Row (cell type) compartment scaffolding
    CT_ORDER_DISPLAY <- c(
        "Endothelial", "Pericytes", "Smooth muscle", "Fibroblasts",
        "Oligodendrocytes", "Microglia", "Astrocytes", "Neurons",
        "T cells", "B cells", "NK cells", "Monocytes",
        "Dendritic cells", "Macrophages", "Plasma cells", "Granulocytes"
    )
    VASC_CT   <- c("Endothelial", "Pericytes", "Smooth muscle", "Fibroblasts")
    CNS_CT    <- c("Oligodendrocytes", "Microglia", "Astrocytes", "Neurons")
    IMMUNE_CT <- setdiff(CT_ORDER_DISPLAY, c(VASC_CT, CNS_CT))

    mat_display <- t(expr_mat_g)
    row_order   <- intersect(CT_ORDER_DISPLAY, rownames(mat_display))
    if (length(row_order) < 2) row_order <- rownames(mat_display)  # fallback for HPA-only
    mat_display <- mat_display[row_order, , drop = FALSE]

    # ── Per-protein peak cell-type validation ────────────────────────────────
    # Peak = argmax over the per-gene min-max-normalised expression vector.
    # (Matrix is already masked to NA → 0 for low-expression entries.)
    gene_mat <- expr_mat_g[, row_order, drop = FALSE]
    peak_idx <- apply(gene_mat, 1, function(x) {
        if (all(x == 0)) NA_integer_ else which.max(x)
    })
    peak_ct  <- ifelse(is.na(peak_idx), NA_character_, row_order[peak_idx])
    peak_val <- apply(gene_mat, 1, function(x) if (all(x == 0)) NA_real_ else max(x))

    peak_compartment <- data.table::fcase(
        peak_ct %in% VASC_CT,   "Vascular",
        peak_ct %in% CNS_CT,    "CNS",
        peak_ct %in% IMMUNE_CT, "Immune",
        default = "None"
    )

    gene_info <- data.table(
        gene             = top40_genes,
        direction        = dep_dir_g$direction,
        logFC            = dep_dir_g$logFC,
        abs_logFC        = abs(dep_dir_g$logFC),
        peak_ct          = peak_ct,
        peak_val         = peak_val,
        peak_compartment = peak_compartment
    )

    # Factor ordering: up-first, then down; compartments Vasc → CNS → Immune → None;
    # peak_ct by display order.
    dir_levels <- c(up_dirs, down_dirs)
    gene_info[, direction        := factor(direction, levels = dir_levels)]
    gene_info[, peak_compartment := factor(peak_compartment,
                                           levels = c("Vascular", "CNS", "Immune", "None"))]
    gene_info[, peak_ct          := factor(peak_ct, levels = row_order)]
    setorder(gene_info, direction, peak_compartment, peak_ct, -abs_logFC, na.last = TRUE)

    ordered_genes <- as.character(gene_info$gene)
    mat_display   <- mat_display[, ordered_genes, drop = FALSE]

    # Gap between the up and down blocks (primary split)
    n_up        <- sum(gene_info$direction %in% up_dirs)
    gaps_col    <- if (n_up > 0 && n_up < nrow(gene_info)) n_up else NULL

    # ── Column annotations: Direction (up/down) + Peak compartment ───────────
    ann_col_g <- data.frame(
        Direction           = as.character(gene_info$direction),
        `Peak compartment`  = as.character(gene_info$peak_compartment),
        row.names           = ordered_genes,
        check.names         = FALSE,
        stringsAsFactors    = FALSE
    )

    # ── Row annotations: cell-type compartment bands ─────────────────────────
    ct_group <- ifelse(row_order %in% VASC_CT, "Vascular",
                ifelse(row_order %in% CNS_CT,  "CNS", "Immune"))
    ann_row_g   <- data.frame(Compartment = ct_group, row.names = row_order)
    comp_cols_g <- c("Vascular" = "#E67E22", "CNS" = "#27AE60",
                     "Immune"   = "#2980B9", "None" = "grey85")
    gaps_row    <- c(sum(ct_group == "Vascular"),
                     sum(ct_group == "Vascular") + sum(ct_group == "CNS"))
    gaps_row    <- gaps_row[gaps_row > 0 & gaps_row < nrow(mat_display)]

    ann_cols_g <- list(
        Direction           = dir_col_vals,
        Compartment         = comp_cols_g[c("Vascular", "CNS", "Immune")],
        `Peak compartment`  = comp_cols_g
    )

    # ── Validation: where do up- vs down-regulated proteins actually peak? ──
    cat("\n  === Peak cell-type validation (pre-onset DEPs, top 40) ===\n")
    comp_dist <- gene_info[, .N, by = .(direction, peak_compartment)][
                    order(direction, peak_compartment)]
    cat("  Compartment distribution by direction:\n")
    print(comp_dist)
    ct_dist <- gene_info[!is.na(peak_ct),
                         .N, by = .(direction, peak_ct)][order(direction, peak_ct)]
    cat("\n  Cell-type peak counts:\n")
    print(ct_dist)
    val_csv <- file.path(FIG_DIR, "panel_g_peak_celltype_validation.csv")
    fwrite(gene_info, val_csv)
    cat(sprintf("  Saved per-protein peak assignments → %s\n", basename(val_csv)))

    out_g <- file.path(FIG_DIR, "panel_g_celltype_heatmap_pre.pdf")
    pdf(out_g, width = 11, height = 6, onefile = FALSE)
    pheatmap(mat_display,
        color                = colorRampPalette(c("white", COL_PRE_UP))(100),
        breaks               = seq(0, 1, length.out = 101),
        cluster_rows         = FALSE,   # enforce Vascular → CNS → Immune order
        cluster_cols         = FALSE,   # enforce up/down block + peak-CT order
        show_rownames        = TRUE,
        show_colnames        = TRUE,
        fontsize_row         = 9,
        fontsize_col         = 6,
        annotation_col       = ann_col_g,
        annotation_row       = ann_row_g,
        annotation_colors    = ann_cols_g,
        annotation_names_col = FALSE,
        annotation_names_row = FALSE,
        border_color         = NA,
        gaps_row             = if (length(gaps_row) > 0) gaps_row else NULL,
        gaps_col             = if (!is.null(gaps_col))   gaps_col else NULL,
        main = "g  Pre-onset MS DEPs: Walchli brain + HPA immune cell-type expression"
    )
    dev.off()
    cat("  Saved: panel_g_celltype_heatmap_pre.pdf\n")
}, error = function(e) cat("  Panel g skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Panel h — Cell-type enrichment (pre-onset vs HC)  [REFACTORED 2025]
# Source: ms_celltype_gsea.csv (06_celltype_overrepresentation.R)
# Method: preranked fgsea on limma moderated-t; specificity-filtered markers
#         (tau + per-cell-type z on log1p(nTPM)) via helpers/celltype_overrep_plot.R.
# Display: single column (pre-onset NES); dot colour = NES (red = up, blue = down),
#          size = -log10(FDR); asterisks mark FDR < {0.05, 0.01, 0.001}.
# Rationale: the legacy Fisher-on-up/down-lists test was biased by marker
# set size (HPA nTPM>=10 covered 35-50% of Olink universe); see helper header.
# ---------------------------------------------------------------------------
cat("Building panel h (GSEA pre-onset)...\n")
tryCatch({
    if (nrow(gsea) == 0)
        stop("ms_celltype_gsea.csv not found — run 06_celltype_overrepresentation.R first")

    pre_gsea <- gsea[analysis == "pre"]
    if (nrow(pre_gsea) == 0) stop("No pre-onset GSEA results")

    pH <- make_celltype_gsea_plot(
        pre_gsea,
        contrast_col    = "analysis",
        contrast_order  = "pre",
        contrast_labels = c(pre = "Pre-onset MS"),
        title_str       = "h  Cell-type enrichment \u2014 pre-onset MS vs HC",
        subtitle_str    = paste0(
            "Preranked fgsea on limma moderated-t | ",
            "NES>0: markers higher in pre-MS | NES<0: lower | * FDR<0.05"
        ),
        high_col        = COL_PRE,
        low_col         = COL_PRE_DOWN
    )
    if (is.null(pH)) stop("make_celltype_gsea_plot returned NULL")
    save_panel(pH, "h_celltype_gsea_pre", 4.5, 6.5)
}, error = function(e) cat("  Panel h skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Panel i — Pre vs post GO pathway comparison heatmap
# Rows: union of top GO terms from pre + post; columns: stage;
# fill: -log10(p.adj); grey = not significant
# ---------------------------------------------------------------------------
cat("Building panel i...\n")
tryCatch({
    if (nrow(go_pre) == 0 || nrow(go_post) == 0)
        stop("Missing pre or post GO results — run 04_prepost_go.R first")

    shorten_go <- function(s, n = 45) ifelse(nchar(s) > n, paste0(substr(s, 1, n), "…"), s)

    go_pre_i  <- copy(go_pre);  go_pre_i[,  stage := "Pre-onset"]
    go_post_i <- copy(go_post); go_post_i[, stage := "Post-onset"]
    both <- rbind(go_pre_i, go_post_i)

    # 5 terms per direction (up/down) per stage → ≤20 unique terms total
    N_TERMS <- 5
    top_terms_i <- both[p.adjust < 0.05][
        , .SD[order(p.adjust)][seq_len(min(N_TERMS, .N))],
        by = .(stage, gene_set)
    ]
    keep_ids <- unique(top_terms_i$ID)
    cat(sprintf("  %d unique GO terms across pre/post top lists\n", length(keep_ids)))
    if (length(keep_ids) == 0) stop("No GO terms to plot for panel i")

    heat_i <- both[ID %in% keep_ids, .(ID, Description, gene_set, stage, p.adjust, FoldEnrichment)]
    heat_i  <- heat_i[, .SD[which.min(p.adjust)], by = .(ID, stage)]

    heat_wide <- dcast(heat_i, ID + Description ~ stage, value.var = "p.adjust", fill = NA_real_)
    setnames(heat_wide, c("Pre-onset", "Post-onset"), c("pval_pre", "pval_post"), skip_absent = TRUE)
    if (!"pval_pre"  %in% names(heat_wide)) heat_wide[, pval_pre  := NA_real_]
    if (!"pval_post" %in% names(heat_wide)) heat_wide[, pval_post := NA_real_]

    heat_wide[, nlp_pre  := ifelse(is.na(pval_pre),  0, pmin(-log10(pval_pre),  8))]
    heat_wide[, nlp_post := ifelse(is.na(pval_post), 0, pmin(-log10(pval_post), 8))]

    heat_wide[, category := fcase(
        !is.na(pval_pre) & pval_pre  < 0.05 & !is.na(pval_post) & pval_post < 0.05, "Both",
        !is.na(pval_pre) & pval_pre  < 0.05,                                          "Pre-onset only",
        !is.na(pval_post) & pval_post < 0.05,                                         "Post-onset only",
        default = "Neither"
    )]
    heat_wide <- heat_wide[category != "Neither"]

    heat_wide[, cat_order := fcase(
        category == "Both",           1L,
        category == "Pre-onset only", 2L,
        default = 3L
    )]
    setorder(heat_wide, cat_order, -nlp_pre)
    heat_wide[, desc_short := shorten_go(Description)]
    heat_wide[, desc_short := make.unique(desc_short, sep = " ")]
    heat_wide[, desc_f := factor(desc_short, levels = rev(unique(desc_short)))]

    long_i <- melt(heat_wide, id.vars = c("ID", "desc_f", "category"),
                   measure.vars = c("nlp_pre", "nlp_post"),
                   variable.name = "stage_var", value.name = "nlp")
    long_i[, stage_label := fifelse(stage_var == "nlp_pre", "Pre-onset", "Post-onset")]
    long_i[, stage_label := factor(stage_label, levels = c("Pre-onset", "Post-onset"))]
    long_i[nlp == 0, nlp := NA_real_]

    cat_col <- c("Both" = "#6B3FA0", "Pre-onset only" = COL_PRE_UP, "Post-onset only" = COL_POST)
    heat_wide[, strip_colour := cat_col[category]]

    pI <- ggplot(long_i, aes(x = stage_label, y = desc_f, fill = nlp)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = ifelse(!is.na(nlp) & nlp >= 1, sprintf("%.1f", nlp), "")),
                  size = 2.1, colour = "white") +
        scale_fill_gradient(
            low = "#FDE8F2", high = COL_PRE_UP, na.value = "grey92",
            name = expression(-log[10]~italic(p)[adj]),
            limits = c(0, 8), oob = scales::squish,
            guide = guide_colourbar(barwidth = 0.5, barheight = 3.5)
        ) +
        annotate("rect",
                 xmin = 2.55, xmax = 2.75,
                 ymin = as.numeric(long_i[stage_label == "Pre-onset"][
                     order(as.integer(desc_f)), .(desc_f)]$desc_f) - 0.5,
                 ymax = as.numeric(long_i[stage_label == "Pre-onset"][
                     order(as.integer(desc_f)), .(desc_f)]$desc_f) + 0.5,
                 fill = heat_wide[order(-as.integer(desc_f)), strip_colour]) +
        scale_x_discrete(expand = expansion(add = c(0.5, 1.2))) +
        labs(title    = "i  Pre vs post-onset GO pathway comparison",
             subtitle = "Top GO:BP terms | Colour strip: stage specificity",
             x = NULL, y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y    = element_text(size = 7),
              axis.text.x    = element_text(size = 8.5, face = "bold"),
              panel.grid     = element_blank(),
              legend.position = "right",
              plot.subtitle  = element_text(size = 7, colour = "grey40"))
    h_i <- max(5, nrow(heat_wide) * 0.28 + 2)
    save_panel(pI, "i_go_comparison", 6, h_i)
}, error = function(e) cat("  Panel i skipped:", conditionMessage(e), "\n"))

cat("\nFigure 1 complete. Panels saved to:", FIG_DIR, "\n")
