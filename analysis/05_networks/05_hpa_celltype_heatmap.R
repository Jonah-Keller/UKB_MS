#!/usr/bin/env Rscript
# 05_hpa_celltype_heatmap.R
# Cell-type expression mapping for MS DEPs using HPA single-cell RNA data.
# Replicates the CADASIL "Walchli atlas" approach but with the publicly
# available HPA v23 single-cell RNA consensus dataset (81 cell types).
#
# Runs for all 4 analyses (combined, pre-onset, post-onset, pre vs post)
# and saves separate CSV matrices + heatmap PDFs for each.
#
# Approach:
#   1. Load DEPs for each analysis (with appropriate FDR threshold)
#   2. Filter HPA data to 14 MS-relevant CNS + immune cell types
#   3. Pivot to gene × cell-type nTPM matrix
#   4. Min-max normalize per gene (0 = min expression, 1 = max)
#   5. Mask low-expression entries (nTPM < 1 → NA/grey)
#   6. Hierarchical clustering on both axes
#   7. Save matrix as CSV + PDF heatmap (full + top-40 versions)
#
# Output (per analysis):
#   results/networks/ms_{analysis}_hpa_celltype_matrix.csv
#   results/networks/ms_{analysis}_hpa_celltype_heatmap.pdf
#   results/networks/ms_{analysis}_hpa_celltype_top40.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR  <- file.path(PROJ_DIR, "results", "networks")
HPA_FILE <- file.path(PROJ_DIR, "data", "reference", "hpa_rna_single_cell_type.tsv.zip")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(HPA_FILE))
    stop("HPA file missing. Download: https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip")

# ---------------------------------------------------------------------------
# MS-relevant cell types (CNS glia + neurons + peripheral immune)
# ---------------------------------------------------------------------------
MS_CELL_TYPES <- c(
    "Microglial cells",
    "Astrocytes",
    "Oligodendrocytes",
    "Oligodendrocyte precursor cells",
    "Excitatory neurons",
    "Inhibitory neurons",
    "T-cells",
    "B-cells",
    "NK-cells",
    "monocytes",
    "dendritic cells",
    "Macrophages",
    "Plasma cells",
    "granulocytes"
)

# ---------------------------------------------------------------------------
# Load HPA data once (shared across all analyses)
# ---------------------------------------------------------------------------
cat("Loading HPA single-cell RNA data...\n")
hpa <- fread(cmd = paste("unzip -p", HPA_FILE))
setnames(hpa, c("Gene", "Gene name", "Cell type", "nTPM"),
               c("ensembl", "gene_name", "cell_type", "nTPM"))
cat(sprintf("  HPA: %d rows, %d cell types, %d genes\n",
            nrow(hpa), length(unique(hpa$cell_type)),
            length(unique(hpa$gene_name))))

missing_ct <- setdiff(MS_CELL_TYPES, unique(hpa$cell_type))
if (length(missing_ct) > 0)
    cat("  Note: not found in HPA:", paste(missing_ct, collapse = ", "), "\n")
MS_CELL_TYPES <- intersect(MS_CELL_TYPES, unique(hpa$cell_type))
cat(sprintf("  Using %d cell types\n", length(MS_CELL_TYPES)))

# ---------------------------------------------------------------------------
# 4 Analyses to run
# ---------------------------------------------------------------------------
analyses <- list(
    list(
        name    = "combined",
        label   = "Combined MS vs HC",
        file    = "ms_combined_vs_hc.csv",
        fdr_thr = 0.05,
        up_col  = "Up in MS",
        dn_col  = "Down in MS"
    ),
    list(
        name    = "pre",
        label   = "Pre-onset vs HC",
        file    = "ms_pre_vs_hc.csv",
        fdr_thr = 0.20,   # relaxed threshold for low-power pre-onset window
        up_col  = "Up in MS",
        dn_col  = "Down in MS"
    ),
    list(
        name    = "post",
        label   = "Post-onset vs HC",
        file    = "ms_post_vs_hc.csv",
        fdr_thr = 0.05,
        up_col  = "Up in MS",
        dn_col  = "Down in MS"
    ),
    list(
        name    = "preVsPost",
        label   = "Pre vs post-onset",
        file    = "ms_pre_vs_post.csv",
        fdr_thr = 0.05,
        up_col  = "Higher pre-onset",
        dn_col  = "Higher post-onset"
    )
)

# ---------------------------------------------------------------------------
# Helper: min-max normalize per row
# ---------------------------------------------------------------------------
minmax <- function(x) {
    r <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (r == 0) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) / r
}

# ---------------------------------------------------------------------------
# Run heatmap pipeline for each analysis
# ---------------------------------------------------------------------------
for (an in analyses) {
    cat(sprintf("\n=== %s (FDR < %.2f) ===\n", an$label, an$fdr_thr))

    diff_path <- file.path(DIFF_DIR, an$file)
    if (!file.exists(diff_path)) {
        cat(sprintf("  Skipping: %s not found\n", an$file))
        next
    }

    deps     <- fread(diff_path)
    deps_sig <- deps[adj.P.Val < an$fdr_thr]
    deps_sig[, direction := fifelse(logFC > 0, an$up_col, an$dn_col)]
    deps_sig[, abs_logFC := abs(logFC)]
    setorder(deps_sig, -abs_logFC)

    dep_genes <- toupper(deps_sig$protein)
    cat(sprintf("  %d DEPs: %d up, %d down\n",
                nrow(deps_sig),
                sum(deps_sig$logFC > 0),
                sum(deps_sig$logFC < 0)))

    if (length(dep_genes) < 5) {
        cat("  Too few DEPs — skipping\n")
        next
    }

    # Filter HPA to this gene set
    hpa_filt <- hpa[cell_type %in% MS_CELL_TYPES & toupper(gene_name) %in% dep_genes]
    cat(sprintf("  HPA filtered: %d gene-cell_type pairs\n", nrow(hpa_filt)))

    if (nrow(hpa_filt) == 0) {
        cat("  No HPA matches — skipping\n")
        next
    }

    # Pivot to gene × cell-type matrix
    expr_wide <- dcast(hpa_filt, gene_name ~ cell_type, value.var = "nTPM", fill = 0)
    gene_col  <- expr_wide$gene_name
    expr_wide[, gene_name := NULL]
    expr_mat  <- as.matrix(expr_wide)
    rownames(expr_mat) <- toupper(gene_col)

    # Ensure all cell types present
    missing_cols <- setdiff(MS_CELL_TYPES, colnames(expr_mat))
    for (mc in missing_cols)
        expr_mat <- cbind(expr_mat, setNames(rep(0, nrow(expr_mat)), mc))
    expr_mat <- expr_mat[, MS_CELL_TYPES]

    cat(sprintf("  Expression matrix: %d genes × %d cell types\n",
                nrow(expr_mat), ncol(expr_mat)))

    # Min-max normalize per gene (row-wise)
    expr_norm <- t(apply(expr_mat, 1, minmax))
    colnames(expr_norm) <- colnames(expr_mat)

    # Remove zero-expression genes
    zero_rows  <- rowSums(expr_norm) == 0
    expr_norm  <- expr_norm[!zero_rows, ]
    cat(sprintf("  After removing zero-expression: %d genes\n", nrow(expr_norm)))

    # Mask nTPM < 1 → NA (shown as grey in heatmap)
    expr_masked <- expr_norm
    for (g in rownames(expr_masked)) {
        if (g %in% rownames(expr_mat)) {
            low <- expr_mat[g, ] < 1
            expr_masked[g, low] <- NA
        }
    }

    # ---------------------------------------------------------------------------
    # Save matrix CSV (includes logFC + direction for figure use)
    # ---------------------------------------------------------------------------
    out_dt  <- as.data.table(expr_masked, keep.rownames = "gene")
    dir_map <- deps_sig[, .(gene = toupper(protein), logFC, direction)]
    out_dt  <- merge(out_dt, dir_map, by = "gene", all.x = TRUE)
    out_dt  <- out_dt[order(-logFC)]

    csv_out <- file.path(OUT_DIR, sprintf("ms_%s_hpa_celltype_matrix.csv", an$name))
    fwrite(out_dt, csv_out)
    cat(sprintf("  Saved %s_hpa_celltype_matrix.csv (%d genes)\n", an$name, nrow(out_dt)))

    # ---------------------------------------------------------------------------
    # Heatmap colours: white → project pink (#CC0066)
    # ---------------------------------------------------------------------------
    heatmap_cols <- colorRampPalette(c("white", "#CC0066"))(100)
    ann_cols <- list(Direction = c(
        setNames("#CC0066", an$up_col),
        setNames("#56B4E9", an$dn_col)
    ))

    # ---------------------------------------------------------------------------
    # Full heatmap: drop genes with >50% NA cell types
    # ---------------------------------------------------------------------------
    na_frac   <- rowMeans(is.na(expr_masked))
    expr_plot <- expr_masked[na_frac <= 0.5, , drop = FALSE]
    cat(sprintf("  Plotting %d genes × %d cell types\n", nrow(expr_plot), ncol(expr_plot)))

    if (nrow(expr_plot) >= 3) {
        # Replace NA with 0 for clustering (cannot cluster with NAs)
        expr_for_clust <- expr_plot
        expr_for_clust[is.na(expr_for_clust)] <- 0

        # Annotation column (gene direction)
        dep_dir <- deps_sig[toupper(protein) %in% rownames(expr_plot),
                            .(gene = toupper(protein), direction)]
        dep_dir <- dep_dir[match(rownames(expr_plot), dep_dir$gene)]
        ann_col_df <- data.frame(
            Direction = dep_dir$direction,
            row.names = rownames(expr_plot)
        )

        # Transpose: rows = cell types, columns = genes (CADASIL layout)
        mat_t  <- t(expr_for_clust)
        show_labels <- ncol(mat_t) <= 60

        pdf_out <- file.path(OUT_DIR, sprintf("ms_%s_hpa_celltype_heatmap.pdf", an$name))
        pdf(pdf_out,
            width  = max(8, ncol(mat_t) * 0.12 + 3),
            height = 5.5, onefile = FALSE)
        pheatmap(mat_t,
            color                = heatmap_cols,
            breaks               = seq(0, 1, length.out = 101),
            cluster_rows         = TRUE,
            cluster_cols         = TRUE,
            show_rownames        = TRUE,
            show_colnames        = show_labels,
            fontsize_row         = 9,
            fontsize_col         = 5.5,
            annotation_col       = ann_col_df,
            annotation_colors    = ann_cols,
            annotation_names_col = FALSE,
            border_color         = NA,
            main = sprintf("%s: DEP cell-type expression (%d DEPs × %d cell types)",
                           an$label, ncol(mat_t), nrow(mat_t))
        )
        dev.off()
        cat(sprintf("  Saved ms_%s_hpa_celltype_heatmap.pdf\n", an$name))
    }

    # ---------------------------------------------------------------------------
    # Top-40 heatmap with gene labels (for figure panels)
    # ---------------------------------------------------------------------------
    top40 <- rownames(expr_plot)[rownames(expr_plot) %in%
        deps_sig[order(-abs_logFC)][seq_len(min(40, .N)), toupper(protein)]]

    if (length(top40) >= 5) {
        expr_top40 <- expr_plot[top40, , drop = FALSE]
        expr_top40[is.na(expr_top40)] <- 0
        mat_top40  <- t(expr_top40)

        ann_top40 <- data.frame(
            Direction = deps_sig[match(top40, toupper(protein)), direction],
            row.names = top40
        )

        pdf_top <- file.path(OUT_DIR, sprintf("ms_%s_hpa_celltype_top40.pdf", an$name))
        pdf(pdf_top, width = 10, height = 5.5, onefile = FALSE)
        pheatmap(mat_top40,
            color                = heatmap_cols,
            breaks               = seq(0, 1, length.out = 101),
            cluster_rows         = TRUE,
            cluster_cols         = TRUE,
            show_rownames        = TRUE,
            show_colnames        = TRUE,
            fontsize_row         = 9,
            fontsize_col         = 6,
            annotation_col       = ann_top40,
            annotation_colors    = ann_cols,
            annotation_names_col = FALSE,
            border_color         = NA,
            main = sprintf("%s: Top %d DEPs by |logFC| — HPA scRNA cell-type expression",
                           an$label, length(top40))
        )
        dev.off()
        cat(sprintf("  Saved ms_%s_hpa_celltype_top40.pdf\n", an$name))
    }
}

# Also keep the legacy combined filename for backward compatibility with figure1.R
legacy_src <- file.path(OUT_DIR, "ms_combined_hpa_celltype_matrix.csv")
legacy_dst <- file.path(OUT_DIR, "ms_hpa_celltype_matrix.csv")
if (file.exists(legacy_src) && !file.exists(legacy_dst))
    file.copy(legacy_src, legacy_dst)

cat("\nDone.\n")
