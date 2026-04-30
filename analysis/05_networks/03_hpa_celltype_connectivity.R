#!/usr/bin/env Rscript
# 03_hpa_celltype_connectivity.R
# Cell-type network connectivity analysis — Abdelhak 2026 replication
#
# Approach:
#   For each key immune/CNS cell type:
#   1. Identify which disease DEPs are expressed in that cell type
#      (curated from published scRNA-seq marker lists)
#   2. Compute observed STRING connectivity density within that subset
#   3. Permutation test (n=1000): randomly draw same-size subsets from all proteins
#   4. Enrichment Z-score and p-value
#
# Cell types tested (as in Abdelhak 2026):
#   B cells, NK cells, T cells, plasma cells, oligodendrocytes,
#   microglia, monocytes (added for completeness)
#
# Cohort prefix is read from configs/disease.yaml.
#
# Output: results/networks/
#   {cohort_short}_celltype_connectivity.pdf
#   {cohort_short}_celltype_connectivity.csv

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(igraph)
    library(glue)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config(file.path(PROJ_DIR, "configs", "disease.yaml"))
COHORT  <- cfg$cohort_short
DISEASE <- cfg$disease_short_caps

NET_DIR  <- file.path(PROJ_DIR, "results", "networks")
DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR  <- NET_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERM       <- 10000   # paper uses 10,000
STRING_THR   <- 0.7     # paper uses 0.7 (high confidence); our edges file has >=0.4
set.seed(42)

# ── 1. Load STRING edges (pre-computed by 01_ppi_network.R) ──────────────────
cat("Loading STRING edges...\n")
edges_basename <- glue("{COHORT}_string_edges.csv")
edges_file     <- file.path(NET_DIR, edges_basename)
if (!file.exists(edges_file))
    stop("Run 01_ppi_network.R first to generate ", edges_basename)

edges <- fread(edges_file)
cat(sprintf("  %d edges loaded (score >= 0.4)\n", nrow(edges)))

# Apply high-confidence threshold to match paper (0.7)
edges <- edges[combined_score >= STRING_THR]
cat(sprintf("  %d edges at score >= %.1f (high confidence)\n", nrow(edges), STRING_THR))

if (nrow(edges) == 0) stop("No edges remain at score >= ", STRING_THR)

# Build the full graph for background
all_nodes <- unique(c(edges$gene_a, edges$gene_b))
g_full    <- graph_from_data_frame(edges[, .(gene_a, gene_b, combined_score)], directed = FALSE)
cat(sprintf("  Full network: %d nodes, %d edges\n", vcount(g_full), ecount(g_full)))

# ── 2. Load DEPs ──────────────────────────────────────────────────────────────
cat(sprintf("Loading %s DEPs...\n", DISEASE))
deps_sig <- fread(file.path(DIFF_DIR, glue("{COHORT}_combined_vs_hc.csv")))[adj.P.Val < 0.05]
dep_genes <- toupper(deps_sig$protein)
cat(sprintf("  %d DEPs, %d in STRING network\n", length(dep_genes),
            sum(dep_genes %in% all_nodes)))

# ── 3. Curated cell-type marker lists ─────────────────────────────────────────
# Sources:
#   B cells:        Human Cell Atlas / Seurat PBMC markers, Tirosh 2016 Cell
#   NK cells:       Crinier 2018 Immunity, HPA single-cell
#   T cells:        Zheng 2017 Nature, HPA single-cell
#   Plasma cells:   Lam 2020 Immunity, HPA single-cell
#   Oligodendrocytes: Jäkel 2019 Nature, Mathys 2019 Nature
#   Microglia:      Masuda 2019 Nat Neurosci, HPA single-cell
#   Monocytes:      Monaco 2019 Cell Reports

cell_type_markers <- list(
    B_cells = c(
        "CD19","MS4A1","CD79A","CD79B","PAX5","FCRL1","FCRL2","FCRL3","FCRL4","FCRL5",
        "CD22","BANK1","IGHM","IGHD","BLK","FCER2","CR2","CD40","BAFF","TNFSF13B",
        "TNFSF13","CXCR5","BCL11A","EBF1","SELL","PTPRC","CD37","CD38","SLAMF7",
        "HLA-DRA","HLA-DRB1","CXCL13","MZB1","JCHAIN","PRDM1","XBP1"
    ),
    NK_cells = c(
        "NCAM1","KLRB1","KLRD1","KLRK1","KLRF1","NKG7","GZMA","GZMB","GZMH","GZMM",
        "PRF1","GNLY","CD56","CD16","FCGR3A","FCGR3B","FCER1G","TYROBP","ITGA1",
        "CX3CR1","S1PR5","EOMES","TBX21","IFNG","TNF","NCR1","NCR3","KLRC1","KLRC2"
    ),
    T_cells = c(
        "CD3D","CD3E","CD3G","CD4","CD8A","CD8B","TCF7","LEF1","CCR7","SELL",
        "IL7R","CD28","ICOS","CTLA4","PDCD1","LAG3","HAVCR2","IL2RA","FOXP3",
        "RORC","TBX21","GATA3","RUNX3","CD44","LCK","ZAP70","ITK","TRAT1",
        "CD2","CD5","CD6","CD7","CD26","CD27","CXCR3","CCR5","CCR6","CXCR4"
    ),
    Plasma_cells = c(
        "PRDM1","XBP1","IRF4","MZB1","SDC1","CD38","IGHG1","IGHG2","IGHG3","IGHG4",
        "IGHA1","IGHA2","IGLC1","IGLC2","IGKC","JCHAIN","SEC61B","DERL3","HSP90B1",
        "FKBP11","FCRL5","CXCR3","CXCR4","VLA4","ITGA4","SLP65","BLNK","CD79A"
    ),
    Oligodendrocytes = c(
        "MBP","MOG","MAG","OLIG1","OLIG2","CNP","CLDN11","PLP1","MOBP","NKX2-2",
        "SOX10","GPR37","ENPP2","ASPA","HAPLN2","FA2H","ERMN","SIRT2","MYRF",
        "GPR17","UGT8","TPPP","QKI","BCAS1","HHIP","LINC01595","MOB3B",
        "NKX6-2","SH3TC2","LPAR5"
    ),
    Microglia = c(
        "CX3CR1","TMEM119","P2RY12","P2RY13","HEXB","SIGLEC1","ITGAM","TREM2",
        "OLFML3","GPR34","SLC2A5","FCRLS","SALL1","PROS1","LPCAT2","CSF1R",
        "CD68","AIF1","ADGRE1","LGMN","CTSS","CTSD","CTSL","MRC1","MSR1",
        "C1QA","C1QB","C1QC","TYROBP","FCER1G","FCGR2B","FCGR3A"
    ),
    Monocytes = c(
        "CD14","FCGR3A","LYZ","S100A8","S100A9","S100A12","VCAN","FCN1","CD16",
        "CSF1R","CX3CR1","CCR2","ITGAM","ITGAX","CD163","MRC1","MARCO","VSIG4",
        "MS4A7","LILRB2","CLEC4A","CLEC7A","CLEC10A","SIGLEC7","SIGLEC9","IFI30",
        "HLA-DRA","HLA-DRB1","CD74","MMP9","CXCL8","IL1B","NLRP3"
    )
)

# ── 4. Connectivity metric ────────────────────────────────────────────────────
# Metric: total edge weight within the subset of DEPs ∩ cell-type markers
compute_connectivity <- function(gene_subset, g) {
    nodes_in_net <- intersect(gene_subset, V(g)$name)
    if (length(nodes_in_net) < 2) return(0)
    sg <- induced_subgraph(g, nodes_in_net)
    sum(E(sg)$combined_score)  # total weighted degree
}

# ── 5. Permutation test ────────────────────────────────────────────────────────
cat("Running permutation tests...\n")
network_genes <- V(g_full)$name

results_list <- list()

for (ct in names(cell_type_markers)) {
    markers    <- toupper(cell_type_markers[[ct]])
    dep_in_ct  <- intersect(dep_genes, markers)
    n_overlap  <- length(dep_in_ct)

    if (n_overlap < 2) {
        cat(sprintf("  %-20s — %d DEPs in markers (skip)\n", ct, n_overlap))
        next
    }

    obs_conn   <- compute_connectivity(dep_in_ct, g_full)

    # Permutation: draw n_overlap random genes from all network genes
    n_net      <- length(network_genes)
    perm_conns <- vapply(seq_len(N_PERM), function(i) {
        rand_set <- sample(network_genes, n_overlap)
        compute_connectivity(rand_set, g_full)
    }, numeric(1))

    z_score  <- (obs_conn - mean(perm_conns)) / (sd(perm_conns) + 1e-9)
    perm_p   <- mean(perm_conns >= obs_conn)

    cat(sprintf("  %-20s — %d DEPs, obs=%.2f, z=%.2f, p=%.4f\n",
                ct, n_overlap, obs_conn, z_score, perm_p))

    results_list[[ct]] <- data.table(
        cell_type   = ct,
        n_dep_markers = n_overlap,
        dep_markers = paste(dep_in_ct, collapse = ";"),
        obs_connectivity = obs_conn,
        perm_mean   = mean(perm_conns),
        perm_sd     = sd(perm_conns),
        z_score     = z_score,
        perm_p      = perm_p
    )
}

if (length(results_list) == 0) stop("No cell types had >= 2 DEPs in markers")

results_dt <- rbindlist(results_list)
results_dt[, fdr := p.adjust(perm_p, method = "BH")]
results_dt <- results_dt[order(perm_p)]
fwrite(results_dt, file.path(OUT_DIR, glue("{COHORT}_celltype_connectivity.csv")))
cat(sprintf("\nResults written: %d cell types tested\n", nrow(results_dt)))
print(results_dt[, .(cell_type, n_dep_markers, z_score, perm_p, fdr)])

# ── 6. Plot ────────────────────────────────────────────────────────────────────
results_dt[, cell_f := factor(cell_type, levels = rev(results_dt$cell_type))]
results_dt[, sig_label := fcase(
    fdr < 0.01,  "**",
    fdr < 0.05,  "*",
    perm_p < 0.05, "+",
    default = ""
)]

p_conn <- ggplot(results_dt, aes(x = z_score, y = cell_f)) +
    geom_col(aes(fill = z_score > 0), colour = "#1A1A1A", width = 0.65, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_text(aes(x = z_score + sign(z_score) * 0.15, label = sig_label),
              size = 4, vjust = 0.75) +
    scale_fill_manual(values = c("TRUE" = "#CC0066", "FALSE" = "#56B4E9"),
                      labels = c("TRUE" = "Enriched", "FALSE" = "Depleted"),
                      name = NULL) +
    labs(x = "Permutation Z-score (network connectivity)",
         y = NULL,
         title = sprintf("%s DEP cell-type network connectivity", DISEASE),
         subtitle = sprintf("STRING edges within DEPs ∩ cell-type markers (n=%d permutations)",
                            N_PERM),
         caption = "* FDR<0.05; ** FDR<0.01; + perm. p<0.05 (uncorrected)") +
    theme_ukb() +
    theme(legend.position = "bottom")

conn_pdf <- glue("{COHORT}_celltype_connectivity.pdf")
ggsave(file.path(OUT_DIR, conn_pdf),
       p_conn, width = 5.5, height = 4.5, device = cairo_pdf)
cat(sprintf("  Saved: %s\n", conn_pdf))

# Panel B: dot plot showing which DEPs are in each cell type
UP_LBL   <- glue("Up in {DISEASE}")
DOWN_LBL <- glue("Down in {DISEASE}")
ct_dep_long <- rbindlist(lapply(names(cell_type_markers), function(ct) {
    markers   <- toupper(cell_type_markers[[ct]])
    dep_in_ct <- intersect(dep_genes, markers)
    if (length(dep_in_ct) == 0) return(NULL)
    # Join with logFC and direction
    dt_sub <- deps_sig[toupper(protein) %in% dep_in_ct,
                       .(protein = toupper(protein), logFC, adj.P.Val,
                         direction = fifelse(logFC > 0, UP_LBL, DOWN_LBL))]
    dt_sub[, cell_type := ct]
    dt_sub
}))

if (nrow(ct_dep_long) > 0) {
    ct_dep_long[, protein_f := reorder(protein, -logFC)]
    ct_dep_long[, cell_f := factor(cell_type, levels = rev(results_dt$cell_type))]

    DIR_COLS <- setNames(c("#CC0066", "#56B4E9"), c(UP_LBL, DOWN_LBL))

    p_dot <- ggplot(ct_dep_long, aes(x = cell_f, y = protein_f,
                                      size = -log10(adj.P.Val), colour = direction)) +
        geom_point(alpha = 0.8) +
        scale_colour_manual(values = DIR_COLS, name = NULL) +
        scale_size_continuous(name = expression(-log[10](FDR)), range = c(1, 5)) +
        labs(x = NULL, y = NULL,
             title = "DEPs expressed in each cell type") +
        coord_flip() +
        theme_ukb() +
        theme(axis.text.y = element_text(size = 6),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

    dot_pdf <- glue("{COHORT}_celltype_dep_dotplot.pdf")
    ggsave(file.path(OUT_DIR, dot_pdf),
           p_dot, width = 7, height = max(3.5, length(unique(ct_dep_long$protein)) * 0.22),
           device = cairo_pdf)
    cat(sprintf("  Saved: %s\n", dot_pdf))
}

cat("\n03_hpa_celltype_connectivity.R complete.\n")
