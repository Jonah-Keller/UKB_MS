#!/usr/bin/env Rscript
# 02_als_ppi_network.R
# STRING PPI network + GO enrichment for ALS combined DEPs (Chia 2025 replication)
#
# Outputs: results/networks/
#   als_ppi_network.pdf
#   als_string_edges.csv
#   als_go_results.csv
#   als_go_enrichment_up.pdf
#   als_go_enrichment_down.pdf

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

args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR  <- file.path(PROJ_DIR, "results", "networks")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STRING_SCORE_THR <- 400
FDR_THR          <- 0.05
SPECIES          <- 9606

# ── 1. Load DEPs ──────────────────────────────────────────────────────────────
cat("Loading ALS combined DEPs...\n")
deps <- fread(file.path(DIFF_DIR, "als_combined_vs_hc.csv"))
deps_sig <- deps[adj.P.Val < FDR_THR]
deps_sig[, direction := fifelse(logFC > 0, "Up in ALS", "Down in ALS")]
cat(sprintf("  %d DEPs: %d up, %d down\n",
            nrow(deps_sig), sum(deps_sig$direction=="Up in ALS"),
            sum(deps_sig$direction=="Down in ALS")))

gene_list <- toupper(deps_sig$protein)

# ── 2. STRING API ─────────────────────────────────────────────────────────────
cat("Querying STRING API...\n")
resp <- tryCatch(
    POST("https://string-db.org/api/json/network",
         body = list(identifiers=paste(gene_list, collapse="\n"),
                     species=as.character(SPECIES),
                     required_score=as.character(STRING_SCORE_THR),
                     caller_identity="ukb_als_proteomics"),
         encode="form", timeout(60)),
    error = function(e) { message("STRING API failed: ", e$message); NULL }
)
if (is.null(resp) || resp$status_code != 200)
    stop("STRING API unavailable.")

edges_raw <- fromJSON(rawToChar(resp$content))
edges <- as.data.table(edges_raw)
cat(sprintf("  STRING returned %d edges\n", nrow(edges)))

edges <- edges[, .(gene_a=preferredName_A, gene_b=preferredName_B, combined_score=score)]
edges <- edges[combined_score >= STRING_SCORE_THR/1000 & gene_a %in% gene_list & gene_b %in% gene_list]
cat(sprintf("  After DEP filter: %d edges\n", nrow(edges)))
fwrite(edges, file.path(OUT_DIR, "als_string_edges.csv"))

if (nrow(edges) == 0) {
    cat("No STRING edges for ALS DEPs at this threshold — skipping network plot.\n")
    q(status=0)
}

# ── 3. Build igraph ───────────────────────────────────────────────────────────
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
g <- delete_vertices(g, which(degree(g) == 0))
set.seed(42)
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)
cat(sprintf("  Network: %d nodes, %d edges, %d communities\n",
            vcount(g), ecount(g), max(membership(communities))))

# ── 4. Network plot ────────────────────────────────────────────────────────────
DIR_COLS <- c("Up in ALS" = "#CC0066", "Down in ALS" = "#56B4E9")
tg <- as_tbl_graph(g) %>% activate(nodes) %>%
    mutate(hub_label = ifelse(centrality_degree() >= quantile(degree(g), 0.6), name, NA_character_))

p_net <- ggraph(tg, layout = "stress") +
    geom_edge_link(aes(alpha = combined_score), colour = "grey70", width = 0.5,
                   show.legend = FALSE) +
    scale_edge_alpha(range = c(0.2, 0.8)) +
    geom_node_point(aes(colour = direction, size = neg_log10p), alpha = 0.9) +
    geom_node_text(aes(label = hub_label), repel = TRUE, size = 2.4,
                   max.overlaps = 30, colour = "black",
                   bg.colour = "white", bg.r = 0.1) +
    scale_colour_manual(values = DIR_COLS, name = NULL) +
    scale_size_continuous(name = expression(-log[10](FDR)),
                          range = c(2, 8), breaks = c(5, 10, 20, 40, 80)) +
    labs(title = sprintf("ALS DEP PPI network (STRING score >= %.1f)",
                         STRING_SCORE_THR/1000),
         subtitle = sprintf("%d nodes, %d edges | %d Louvain communities",
                            vcount(g), ecount(g), max(membership(communities)))) +
    theme_graph(base_size=10, base_family="") +
    theme(legend.position="bottom",
          plot.title=element_text(size=10, face="plain"),
          plot.subtitle=element_text(size=8, colour="grey40"),
          panel.background=element_rect(fill="white",colour=NA),
          plot.background=element_rect(fill="white",colour=NA))

ggsave(file.path(OUT_DIR, "als_ppi_network.pdf"),
       p_net, width=8, height=7, device=cairo_pdf)
cat("  Saved: als_ppi_network.pdf\n")

# ── 5. GO enrichment ──────────────────────────────────────────────────────────
run_enrichment <- function(gene_symbols, label) {
    if (length(gene_symbols) < 5) return(NULL)
    entrez <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    if (nrow(entrez) == 0) return(NULL)
    ego <- enrichGO(gene=entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)
    if (is.null(ego) || nrow(as.data.frame(ego))==0) return(NULL)
    cat(sprintf("  %s: %d GO BP terms\n", label, nrow(as.data.frame(ego))))
    ego
}

make_dotplot <- function(ego, title_str, filename) {
    if (is.null(ego)) return(invisible(NULL))
    p <- dotplot(ego, showCategory=15, color="p.adjust", label_format=40) +
        labs(title=title_str) +
        scale_colour_gradient(low="#CC0066", high="#56B4E9", name="Adj. p-value") +
        theme_ukb(base_size=9) +
        theme(axis.text.y=element_text(size=7), legend.position="right")
    ggsave(file.path(OUT_DIR, filename), p, width=7.5, height=5.5, device=cairo_pdf)
    cat(sprintf("  Saved: %s\n", filename))
}

cat("Running GO BP enrichment...\n")
ego_up   <- run_enrichment(toupper(deps_sig[direction=="Up in ALS",   protein]), "Up-regulated")
ego_down <- run_enrichment(toupper(deps_sig[direction=="Down in ALS", protein]), "Down-regulated")
ego_all  <- run_enrichment(toupper(deps_sig$protein), "All ALS DEPs")

make_dotplot(ego_up,   "GO BP: Up-regulated ALS DEPs",   "als_go_enrichment_up.pdf")
make_dotplot(ego_down, "GO BP: Down-regulated ALS DEPs", "als_go_enrichment_down.pdf")
make_dotplot(ego_all,  "GO BP: All ALS DEPs",            "als_go_enrichment_combined.pdf")

go_list <- Filter(Negate(is.null), list(
    if (!is.null(ego_all))  { d <- as.data.table(ego_all@result)[p.adjust<0.05];  d[,gene_set:="all"]; d },
    if (!is.null(ego_up))   { d <- as.data.table(ego_up@result)[p.adjust<0.05];   d[,gene_set:="up"];  d },
    if (!is.null(ego_down)) { d <- as.data.table(ego_down@result)[p.adjust<0.05]; d[,gene_set:="down"];d }
))
if (length(go_list) > 0) {
    fwrite(rbindlist(go_list, fill=TRUE), file.path(OUT_DIR, "als_go_results.csv"))
    cat(sprintf("  GO results written\n"))
}

cat("\n02_als_ppi_network.R complete.\n")
