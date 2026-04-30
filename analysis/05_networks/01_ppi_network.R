#!/usr/bin/env Rscript
# 01_ppi_network.R
# STRING PPI network + GO enrichment for MS combined DEPs (Abdelhak replication)
#
# Approach:
#   1. Load 173 MS combined DEPs
#   2. Fetch STRING interactions via REST API (species 9606, score >= 400)
#   3. Build igraph; Louvain community detection
#   4. ggraph network plot: nodes coloured by direction, sized by -log10(FDR)
#   5. GO BP enrichment (clusterProfiler) + dotplot; separate for up/down
#
# Outputs: results/networks/
#   ms_ppi_network.pdf
#   ms_ppi_network_hires.pdf
#   ms_go_enrichment_up.pdf
#   ms_go_enrichment_down.pdf
#   ms_go_enrichment_combined.pdf
#   ms_string_edges.csv
#   ms_go_results.csv

suppressPackageStartupMessages({
    library(data.table)
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(ggplot2)
    library(ggrepel)
    library(httr)
    library(jsonlite)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR  <- file.path(PROJ_DIR, "results", "networks")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STRING_SCORE_THR <- 400   # medium confidence (0–1000)
FDR_THR          <- 0.05
SPECIES          <- 9606  # Homo sapiens

# ── 1. Load DEPs ──────────────────────────────────────────────────────────────
cat("Loading MS combined DEPs...\n")
deps <- fread(file.path(DIFF_DIR, "ms_combined_vs_hc.csv"))
deps_sig <- deps[adj.P.Val < FDR_THR]
deps_sig[, direction := fifelse(logFC > 0, "Up in MS", "Down in MS")]
cat(sprintf("  %d DEPs: %d up, %d down\n",
            nrow(deps_sig), sum(deps_sig$direction == "Up in MS"),
            sum(deps_sig$direction == "Down in MS")))

gene_list <- toupper(deps_sig$protein)

# ── 2. STRING REST API ────────────────────────────────────────────────────────
cat("Querying STRING API...\n")
string_url <- "https://string-db.org/api/json/network"

resp <- tryCatch(
    POST(string_url,
         body = list(
             identifiers     = paste(gene_list, collapse = "\n"),
             species         = as.character(SPECIES),
             required_score  = as.character(STRING_SCORE_THR),
             caller_identity = "ukb_ms_proteomics"
         ),
         encode = "form",
         timeout(60)),
    error = function(e) {
        message("STRING API call failed: ", e$message)
        NULL
    }
)

if (is.null(resp) || resp$status_code != 200) {
    stop("STRING API unavailable or returned error. Check internet connection.")
}

edges_raw <- fromJSON(rawToChar(resp$content))
edges <- as.data.table(edges_raw)
cat(sprintf("  STRING returned %d edges before score filter\n", nrow(edges)))

# score field from REST API is already 0–1
edges <- edges[, .(
    gene_a         = preferredName_A,
    gene_b         = preferredName_B,
    combined_score = score            # already 0–1
)]
edges <- edges[combined_score >= STRING_SCORE_THR / 1000]
edges <- edges[gene_a %in% gene_list & gene_b %in% gene_list]

cat(sprintf("  After filtering to DEPs only: %d edges\n", nrow(edges)))
fwrite(edges, file.path(OUT_DIR, "ms_string_edges.csv"))

if (nrow(edges) == 0) stop("No STRING edges found for DEP set — check gene names.")

# ── 3. Build igraph ───────────────────────────────────────────────────────────
cat("Building PPI network...\n")
g <- graph_from_data_frame(
    d        = edges[, .(gene_a, gene_b, combined_score)],
    directed = FALSE,
    vertices = deps_sig[, .(
        name       = toupper(protein),
        direction  = direction,
        logFC      = logFC,
        neg_log10p = -log10(adj.P.Val)
    )]
)

# Remove isolated nodes (not connected in STRING at this threshold)
g <- delete_vertices(g, which(degree(g) == 0))
cat(sprintf("  Network: %d nodes, %d edges (after removing isolates)\n",
            vcount(g), ecount(g)))

# Louvain community detection
set.seed(42)
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)
cat(sprintf("  Louvain: %d communities\n", max(membership(communities))))

# ── 4. Network plot ────────────────────────────────────────────────────────────
cat("Plotting PPI network...\n")
tg <- as_tbl_graph(g)

# Colour scheme matching UKB theme
DIR_COLS <- c("Up in MS" = "#CC0066", "Down in MS" = "#56B4E9")

# Node label: show only high-confidence hubs
hub_threshold <- quantile(degree(g), 0.75)
V(g)$hub_label <- ifelse(degree(g) >= hub_threshold, V(g)$name, NA_character_)
tg <- tg %>% activate(nodes) %>%
    mutate(hub_label = ifelse(centrality_degree() >= hub_threshold, name, NA_character_))

p_net <- ggraph(tg, layout = "stress") +
    geom_edge_link(aes(alpha = combined_score), colour = "grey70", width = 0.4,
                   show.legend = FALSE) +
    scale_edge_alpha(range = c(0.15, 0.7)) +
    geom_node_point(aes(colour = direction, size = neg_log10p), alpha = 0.85) +
    geom_node_text(aes(label = hub_label), repel = TRUE, size = 2.2,
                   max.overlaps = 20, colour = "black",
                   bg.colour = "white", bg.r = 0.1) +
    scale_colour_manual(values = DIR_COLS, name = NULL) +
    scale_size_continuous(name = expression(-log[10](FDR)),
                          range = c(1.5, 6), breaks = c(5, 10, 20, 30)) +
    labs(title = sprintf("MS DEP PPI network (STRING score \u2265 %.1f)",
                         STRING_SCORE_THR / 1000),
         subtitle = sprintf("%d nodes, %d edges | %d Louvain communities",
                            vcount(g), ecount(g), max(membership(communities)))) +
    theme_graph(base_size = 10, base_family = "") +
    theme(
        legend.position   = "bottom",
        plot.title        = element_text(size = 10, face = "plain"),
        plot.subtitle     = element_text(size = 8, colour = "grey40"),
        legend.key.size   = unit(0.8, "lines"),
        panel.background  = element_rect(fill = "white", colour = NA),
        plot.background   = element_rect(fill = "white", colour = NA)
    )

ggsave(file.path(OUT_DIR, "ms_ppi_network.pdf"),
       p_net, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "ms_ppi_network_hires.pdf"),
       p_net, width = 12, height = 10, device = cairo_pdf)
cat("  Saved: ms_ppi_network.pdf\n")

# ── 5. GO BP enrichment ────────────────────────────────────────────────────────
cat("Running GO BP enrichment...\n")

run_enrichment <- function(gene_symbols, label) {
    if (length(gene_symbols) < 5) {
        message(sprintf("  Skipping %s enrichment: only %d genes", label, length(gene_symbols)))
        return(NULL)
    }
    # Map symbols to Entrez IDs
    entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
    if (nrow(entrez) == 0) { message("  No Entrez mappings found for ", label); return(NULL) }

    ego <- enrichGO(
        gene          = entrez$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message(sprintf("  No significant GO terms for %s", label))
        return(NULL)
    }
    cat(sprintf("  %s: %d significant GO BP terms\n", label, nrow(as.data.frame(ego))))
    ego
}

# Background: all proteins tested (all columns in the DEP file)
all_tested <- toupper(deps$protein)
up_genes   <- toupper(deps_sig[direction == "Up in MS",   protein])
down_genes <- toupper(deps_sig[direction == "Down in MS", protein])
all_genes  <- toupper(deps_sig$protein)

ego_up   <- run_enrichment(up_genes,   "Up-regulated")
ego_down <- run_enrichment(down_genes, "Down-regulated")
ego_all  <- run_enrichment(all_genes,  "Combined DEPs")

make_dotplot <- function(ego, title_str, filename, n_terms = 20) {
    if (is.null(ego)) return(invisible(NULL))
    p <- dotplot(ego, showCategory = n_terms,
                 color = "p.adjust", label_format = 40) +
        labs(title = title_str) +
        scale_colour_gradient(low = "#CC0066", high = "#56B4E9",
                              name = "Adj. p-value") +
        theme_ukb(base_size = 9) +
        theme(
            axis.text.y   = element_text(size = 7),
            legend.position = "right",
            plot.title    = element_text(size = 9)
        )
    ggsave(file.path(OUT_DIR, filename), p,
           width = 8, height = 6, device = cairo_pdf)
    cat(sprintf("  Saved: %s\n", filename))
}

make_dotplot(ego_up,   "GO BP: Up-regulated MS DEPs",   "ms_go_enrichment_up.pdf")
make_dotplot(ego_down, "GO BP: Down-regulated MS DEPs", "ms_go_enrichment_down.pdf")
make_dotplot(ego_all,  "GO BP: All MS DEPs (combined)", "ms_go_enrichment_combined.pdf")

# Save results table
if (!is.null(ego_all)) {
    go_dt <- as.data.table(ego_all@result)[p.adjust < 0.05]
    go_dt[, gene_set := "all_DEPs"]
    if (!is.null(ego_up)) {
        go_up_dt <- as.data.table(ego_up@result)[p.adjust < 0.05]
        go_up_dt[, gene_set := "up_DEPs"]
        go_dt <- rbind(go_dt, go_up_dt, fill = TRUE)
    }
    if (!is.null(ego_down)) {
        go_down_dt <- as.data.table(ego_down@result)[p.adjust < 0.05]
        go_down_dt[, gene_set := "down_DEPs"]
        go_dt <- rbind(go_dt, go_down_dt, fill = TRUE)
    }
    fwrite(go_dt, file.path(OUT_DIR, "ms_go_results.csv"))
    cat(sprintf("  GO results written: %d rows\n", nrow(go_dt)))
}

# ── 6. KEGG pathway enrichment ────────────────────────────────────────────────
cat("Running KEGG pathway enrichment...\n")

run_kegg <- function(gene_symbols, label) {
    if (length(gene_symbols) < 5) {
        message(sprintf("  Skipping %s KEGG: only %d genes", label, length(gene_symbols)))
        return(NULL)
    }
    entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
    if (nrow(entrez) == 0) return(NULL)
    ekegg <- tryCatch(
        enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
        error = function(e) { message("  KEGG failed: ", e$message); NULL }
    )
    if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
        message(sprintf("  No significant KEGG pathways for %s", label))
        return(NULL)
    }
    cat(sprintf("  %s: %d significant KEGG pathways\n", label, nrow(as.data.frame(ekegg))))
    ekegg
}

kegg_up   <- run_kegg(up_genes,   "Up-regulated")
kegg_down <- run_kegg(down_genes, "Down-regulated")
kegg_all  <- run_kegg(all_genes,  "Combined DEPs")

make_kegg_dotplot <- function(ekegg, title_str, filename, n_terms = 20) {
    if (is.null(ekegg)) return(invisible(NULL))
    # setReadable for gene symbols in KEGG result
    ekegg2 <- tryCatch(setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                       error = function(e) ekegg)
    p <- dotplot(ekegg2, showCategory = n_terms, color = "p.adjust", label_format = 40) +
        labs(title = title_str) +
        scale_colour_gradient(low = "#CC0066", high = "#56B4E9", name = "Adj. p-value") +
        theme_ukb(base_size = 9) +
        theme(axis.text.y = element_text(size = 7), legend.position = "right",
              plot.title = element_text(size = 9))
    ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 6, device = cairo_pdf)
    cat(sprintf("  Saved: %s\n", filename))
}

make_kegg_dotplot(kegg_up,   "KEGG: Up-regulated MS DEPs",   "ms_kegg_up.pdf")
make_kegg_dotplot(kegg_down, "KEGG: Down-regulated MS DEPs", "ms_kegg_down.pdf")
make_kegg_dotplot(kegg_all,  "KEGG: All MS DEPs (combined)", "ms_kegg_combined.pdf")

# ── 7. WikiPathways enrichment ─────────────────────────────────────────────────
cat("Running WikiPathways enrichment...\n")

run_wp <- function(gene_symbols, label) {
    if (length(gene_symbols) < 5) return(NULL)
    entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
    if (nrow(entrez) == 0) return(NULL)
    ewp <- tryCatch(
        enrichWP(gene = entrez$ENTREZID, organism = "Homo sapiens",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
        error = function(e) { message("  WikiPathways failed: ", e$message); NULL }
    )
    if (is.null(ewp) || nrow(as.data.frame(ewp)) == 0) {
        message(sprintf("  No significant WikiPathways for %s", label))
        return(NULL)
    }
    cat(sprintf("  %s: %d significant WikiPathways\n", label, nrow(as.data.frame(ewp))))
    ewp
}

wp_up   <- run_wp(up_genes,   "Up-regulated")
wp_down <- run_wp(down_genes, "Down-regulated")
wp_all  <- run_wp(all_genes,  "Combined DEPs")

make_wp_dotplot <- function(ewp, title_str, filename, n_terms = 20) {
    if (is.null(ewp)) return(invisible(NULL))
    ewp2 <- tryCatch(setReadable(ewp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                     error = function(e) ewp)
    p <- dotplot(ewp2, showCategory = n_terms, color = "p.adjust", label_format = 40) +
        labs(title = title_str) +
        scale_colour_gradient(low = "#CC0066", high = "#56B4E9", name = "Adj. p-value") +
        theme_ukb(base_size = 9) +
        theme(axis.text.y = element_text(size = 7), legend.position = "right",
              plot.title = element_text(size = 9))
    ggsave(file.path(OUT_DIR, filename), p, width = 8, height = 6, device = cairo_pdf)
    cat(sprintf("  Saved: %s\n", filename))
}

make_wp_dotplot(wp_up,   "WikiPathways: Up-regulated MS DEPs",   "ms_wp_up.pdf")
make_wp_dotplot(wp_down, "WikiPathways: Down-regulated MS DEPs", "ms_wp_down.pdf")
make_wp_dotplot(wp_all,  "WikiPathways: All MS DEPs (combined)", "ms_wp_combined.pdf")

# ── 8. Append KEGG + WP to combined results CSV ───────────────────────────────
path_results_csv <- file.path(OUT_DIR, "ms_go_results.csv")
pathway_rows <- list()

for (res_obj in list(
        list(obj=kegg_up,   gs="up_DEPs",   db="KEGG"),
        list(obj=kegg_down, gs="down_DEPs", db="KEGG"),
        list(obj=kegg_all,  gs="all_DEPs",  db="KEGG"),
        list(obj=wp_up,     gs="up_DEPs",   db="WikiPathways"),
        list(obj=wp_down,   gs="down_DEPs", db="WikiPathways"),
        list(obj=wp_all,    gs="all_DEPs",  db="WikiPathways")
)) {
    if (!is.null(res_obj$obj)) {
        dt_tmp <- as.data.table(res_obj$obj@result)[p.adjust < 0.05]
        dt_tmp[, gene_set := res_obj$gs]
        dt_tmp[, database := res_obj$db]
        pathway_rows[[length(pathway_rows)+1]] <- dt_tmp
    }
}

if (length(pathway_rows) > 0) {
    pathway_dt <- rbindlist(pathway_rows, fill = TRUE)
    # Append to existing GO results if present
    if (file.exists(path_results_csv)) {
        existing <- fread(path_results_csv)
        existing[, database := "GO_BP"]
        combined_dt <- rbind(existing, pathway_dt, fill = TRUE)
    } else {
        combined_dt <- pathway_dt
    }
    fwrite(combined_dt, file.path(OUT_DIR, "ms_pathway_results.csv"))
    cat(sprintf("  Pathway results written: %d rows (KEGG + WikiPathways)\n", nrow(pathway_dt)))
}

cat("\n01_ppi_network.R complete.\n")
