#!/usr/bin/env Rscript
# 07_walchli_hpa_integrated_heatmap.R
# Integrated cell-type expression heatmap for disease DEPs combining:
#   - Walchli normal temporal-lobe brain atlas (vascular + CNS cell types)
#     Source: pre-computed AggregateExpression CSV from CADASIL project
#             (originally extracted from walchli_object.h5Seurat on Minerva HPC)
#   - HPA v23 scRNA-seq consensus (peripheral immune cell types only)
#
# Cell types used:
#   Walchli (brain/vascular):
#     Endothelial cells, Pericytes, Smooth muscle cells, Fibroblasts,
#     Oligodendrocytes, Microglia & Macrophages, Astrocytes, Neurons
#   HPA (peripheral immune):
#     T cells, B cells, NK cells, Monocytes, Dendritic cells,
#     Macrophages, Plasma cells, Granulocytes
#
# Normalization: min-max per gene across ALL cell types combined (0-1).
# Walchli values: log1p(avg_aggregate_expression) then min-max.
# HPA values: nTPM then min-max. Both normalised to same 0-1 scale per gene.
# Masking: Walchli log1p < 1 OR HPA nTPM < 1 → NA (grey in heatmap).
#
# Runs for all 4 analyses (combined, pre, post, preVsPost).
#
# TODO (Minerva): If the Walchli CSV needs to be regenerated from the full
#   Seurat object (e.g., for updated cell-type annotations or a different atlas),
#   submit analysis/hpc/extract_walchli_avgexpr.lsf to Minerva (kellej10@minerva).
#   The object is at: /sc/arion/projects/.../walchli_object.h5Seurat (~14 GB RAM req).
#   Output will be a small CSV to replace data/reference/walchli_avg_expr_by_celltype.csv.
#
# Cohort prefix is read from configs/disease.yaml.
#
# Outputs (per analysis):
#   results/networks/{cohort_short}_{analysis}_walchli_hpa_matrix.csv
#   results/networks/{cohort_short}_{analysis}_walchli_hpa_heatmap.pdf
#   results/networks/{cohort_short}_{analysis}_walchli_hpa_top40.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(pheatmap)
    library(RColorBrewer)
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
UP_LBL   <- glue("Up in {DISEASE}")
DOWN_LBL <- glue("Down in {DISEASE}")

DIFF_DIR    <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR     <- file.path(PROJ_DIR, "results", "networks")
WALCHLI_CSV <- file.path(PROJ_DIR, "data", "reference", "walchli_avg_expr_by_celltype.csv")
HPA_FILE    <- file.path(PROJ_DIR, "data", "reference", "hpa_rna_single_cell_type.tsv.zip")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(WALCHLI_CSV))
    stop("Walchli CSV missing: ", WALCHLI_CSV)
if (!file.exists(HPA_FILE))
    stop("HPA file missing: ", HPA_FILE)

# ---------------------------------------------------------------------------
# Cell type selection
# ---------------------------------------------------------------------------
# Walchli brain/vascular cell types (column names use dots for spaces)
WALCHLI_CT_COLS <- c(
    "Endothelial.cells",
    "Pericytes",
    "Smooth.muscle.cells",
    "Fibroblasts",
    "Oligodendrocytes",
    "Microglia.and.Macrophages",
    "Astrocytes",
    "Neurons"
)
WALCHLI_CT_LABELS <- c(
    "Endothelial.cells"        = "Endothelial",
    "Pericytes"                = "Pericytes",
    "Smooth.muscle.cells"      = "Smooth muscle",
    "Fibroblasts"              = "Fibroblasts",
    "Oligodendrocytes"         = "Oligodendrocytes",
    "Microglia.and.Macrophages" = "Microglia",
    "Astrocytes"               = "Astrocytes",
    "Neurons"                  = "Neurons"
)

# HPA peripheral immune cell types (brain cell types covered by Walchli)
HPA_IMMUNE_CT <- c(
    "T-cells", "B-cells", "NK-cells",
    "monocytes", "dendritic cells", "Macrophages",
    "Plasma cells", "granulocytes"
)
HPA_CT_LABELS <- c(
    "T-cells"        = "T cells",
    "B-cells"        = "B cells",
    "NK-cells"       = "NK cells",
    "monocytes"      = "Monocytes",
    "dendritic cells" = "Dendritic cells",
    "Macrophages"    = "Macrophages",
    "Plasma cells"   = "Plasma cells",
    "granulocytes"   = "Granulocytes"
)

# Display order: vascular → CNS → immune
CT_ORDER_DISPLAY <- c(
    "Endothelial", "Pericytes", "Smooth muscle", "Fibroblasts",
    "Oligodendrocytes", "Microglia", "Astrocytes", "Neurons",
    "T cells", "B cells", "NK cells", "Monocytes",
    "Dendritic cells", "Macrophages", "Plasma cells", "Granulocytes"
)

# ---------------------------------------------------------------------------
# Load data sources
# ---------------------------------------------------------------------------
cat("Loading Walchli average expression...\n")
walchli_raw <- fread(WALCHLI_CSV)
setnames(walchli_raw, "Gene", "gene_name")
walchli_raw[, gene_name := toupper(gene_name)]
missing_cols <- setdiff(WALCHLI_CT_COLS, names(walchli_raw))
if (length(missing_cols) > 0)
    cat("  Warning: missing Walchli columns:", paste(missing_cols, collapse = ", "), "\n")
WALCHLI_CT_COLS <- intersect(WALCHLI_CT_COLS, names(walchli_raw))
cat(sprintf("  %d genes, %d cell types\n",
            nrow(walchli_raw), length(WALCHLI_CT_COLS)))

cat("Loading HPA single-cell RNA data...\n")
hpa <- fread(cmd = paste("unzip -p", HPA_FILE))
setnames(hpa, c("Gene", "Gene name", "Cell type", "nTPM"),
               c("ensembl", "gene_name", "cell_type", "nTPM"))
hpa[, gene_name := toupper(gene_name)]
missing_hpa <- setdiff(HPA_IMMUNE_CT, unique(hpa$cell_type))
if (length(missing_hpa) > 0)
    cat("  Warning: missing HPA immune types:", paste(missing_hpa, collapse = ", "), "\n")
HPA_IMMUNE_CT <- intersect(HPA_IMMUNE_CT, unique(hpa$cell_type))

# ---------------------------------------------------------------------------
# Min-max helper
# ---------------------------------------------------------------------------
minmax <- function(x) {
    r <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (r == 0) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) / r
}

# ---------------------------------------------------------------------------
# 4 Analyses
# ---------------------------------------------------------------------------
analyses <- list(
    list(name = "combined", label = sprintf("Combined %s vs HC", DISEASE),
         file = glue("{COHORT}_combined_vs_hc.csv"), fdr_thr = 0.05,
         up_col = UP_LBL, dn_col = DOWN_LBL),
    list(name = "pre",      label = "Pre-onset vs HC",
         file = glue("{COHORT}_pre_vs_hc.csv"),      fdr_thr = 0.20,
         up_col = UP_LBL, dn_col = DOWN_LBL),
    list(name = "post",     label = "Post-onset vs HC",
         file = glue("{COHORT}_post_vs_hc.csv"),     fdr_thr = 0.05,
         up_col = UP_LBL, dn_col = DOWN_LBL),
    list(name = "preVsPost", label = "Pre vs post-onset",
         file = glue("{COHORT}_pre_vs_post.csv"),    fdr_thr = 0.05,
         up_col = "Higher pre-onset", dn_col = "Higher post-onset")
)

for (an in analyses) {
    cat(sprintf("\n=== %s (FDR < %.2f) ===\n", an$label, an$fdr_thr))

    diff_path <- file.path(DIFF_DIR, an$file)
    if (!file.exists(diff_path)) { cat("  Skipping: file not found\n"); next }

    deps     <- fread(diff_path)
    deps_sig <- deps[adj.P.Val < an$fdr_thr]
    if (nrow(deps_sig) < 5) { cat("  Too few DEPs — skipping\n"); next }

    deps_sig[, direction := fifelse(logFC > 0, an$up_col, an$dn_col)]
    deps_sig[, abs_logFC := abs(logFC)]
    setorder(deps_sig, -abs_logFC)
    dep_genes <- toupper(deps_sig$protein)
    cat(sprintf("  %d DEPs (%d up, %d down)\n",
                nrow(deps_sig), sum(deps_sig$logFC > 0), sum(deps_sig$logFC < 0)))

    # ── Walchli: log1p-transform, filter to DEPs ──────────────────────────────
    walchli_deps <- walchli_raw[gene_name %in% dep_genes]
    if (nrow(walchli_deps) == 0) { cat("  No Walchli matches — skipping\n"); next }

    w_mat <- as.matrix(walchli_deps[, ..WALCHLI_CT_COLS])
    rownames(w_mat) <- walchli_deps$gene_name
    w_mat_log <- log1p(w_mat)
    WALCHLI_MASK_THR <- 1.0   # log1p ≥ 1 ≈ avg_expr ≥ 1.7

    # ── HPA: pivot to gene × immune-cell matrix ────────────────────────────────
    hpa_filt <- hpa[cell_type %in% HPA_IMMUNE_CT & gene_name %in% dep_genes]
    hpa_wide  <- dcast(hpa_filt, gene_name ~ cell_type, value.var = "nTPM", fill = 0)
    h_genes   <- hpa_wide$gene_name
    hpa_wide[, gene_name := NULL]
    # Ensure all immune types present
    for (ct in setdiff(HPA_IMMUNE_CT, names(hpa_wide)))
        hpa_wide[, (ct) := 0]
    h_mat  <- as.matrix(hpa_wide[, ..HPA_IMMUNE_CT])
    rownames(h_mat) <- h_genes
    HPA_MASK_THR <- 1.0   # nTPM ≥ 1

    # ── Find gene universe (present in both or either; use union) ─────────────
    all_genes <- union(rownames(w_mat), rownames(h_mat))
    # Pad missing genes with zeros
    pad_rows <- function(mat, target_genes) {
        missing <- setdiff(target_genes, rownames(mat))
        if (length(missing) > 0) {
            pad <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                          dimnames = list(missing, colnames(mat)))
            mat <- rbind(mat, pad)
        }
        mat[target_genes, , drop = FALSE]
    }
    w_mat_log_full <- pad_rows(w_mat_log, all_genes)
    h_mat_full     <- pad_rows(h_mat,     all_genes)

    # ── Combined min-max per gene across all cell types ───────────────────────
    combined_raw <- cbind(w_mat_log_full, h_mat_full)
    combined_norm <- t(apply(combined_raw, 1, minmax))
    colnames(combined_norm) <- c(WALCHLI_CT_COLS, HPA_IMMUNE_CT)

    # ── Apply expression masks (set low-expression entries → NA) ─────────────
    combined_masked <- combined_norm
    for (g in rownames(combined_masked)) {
        # Walchli columns
        for (ct in WALCHLI_CT_COLS) {
            if (g %in% rownames(w_mat_log)) {
                if (w_mat_log[g, ct] < WALCHLI_MASK_THR)
                    combined_masked[g, ct] <- NA
            } else {
                combined_masked[g, ct] <- NA
            }
        }
        # HPA columns
        for (ct in HPA_IMMUNE_CT) {
            if (g %in% rownames(h_mat)) {
                if (h_mat[g, ct] < HPA_MASK_THR)
                    combined_masked[g, ct] <- NA
            } else {
                combined_masked[g, ct] <- NA
            }
        }
    }

    # Drop genes with zero expression across all cell types
    zero_rows <- rowSums(combined_norm, na.rm = TRUE) == 0
    combined_masked <- combined_masked[!zero_rows, ]
    cat(sprintf("  Integrated matrix: %d genes × %d cell types\n",
                nrow(combined_masked), ncol(combined_masked)))

    # ── Rename cell-type columns to display labels ────────────────────────────
    all_labels <- c(WALCHLI_CT_LABELS[WALCHLI_CT_COLS], HPA_CT_LABELS[HPA_IMMUNE_CT])
    colnames(combined_masked) <- all_labels[colnames(combined_masked)]

    # ── Save CSV ──────────────────────────────────────────────────────────────
    out_dt  <- as.data.table(combined_masked, keep.rownames = "gene")
    dir_map <- deps_sig[, .(gene = toupper(protein), logFC, direction)]
    out_dt  <- merge(out_dt, dir_map, by = "gene", all.x = TRUE)
    out_dt  <- out_dt[order(-logFC)]
    csv_name <- glue("{COHORT}_{an$name}_walchli_hpa_matrix.csv")
    fwrite(out_dt, file.path(OUT_DIR, csv_name))
    cat(sprintf("  Saved %s\n", csv_name))

    # ── Heatmap helpers ───────────────────────────────────────────────────────
    heatmap_cols <- colorRampPalette(c("white", "#CC0066"))(100)
    ann_cols_up <- "#CC0066"; ann_cols_dn <- "#56B4E9"

    # Filter: ≤50% NA across cell types
    na_frac   <- rowMeans(is.na(combined_masked))
    expr_plot <- combined_masked[na_frac <= 0.5, , drop = FALSE]
    cat(sprintf("  Plotting %d genes × %d cell types\n", nrow(expr_plot), ncol(expr_plot)))

    if (nrow(expr_plot) < 3) { cat("  Too few genes to heatmap — skipping\n"); next }

    expr_for_clust <- expr_plot; expr_for_clust[is.na(expr_for_clust)] <- 0

    # Direction annotation
    dep_dir <- deps_sig[toupper(protein) %in% rownames(expr_plot),
                        .(gene = toupper(protein), direction)]
    dep_dir <- dep_dir[match(rownames(expr_plot), dep_dir$gene)]
    ann_col_df <- data.frame(Direction = dep_dir$direction,
                              row.names = rownames(expr_plot))
    dir_col_vals <- c(setNames(ann_cols_up, an$up_col),
                      setNames(ann_cols_dn, an$dn_col))

    # Add vascular group annotation for rows (cell-type rows after transpose)
    ct_group <- ifelse(colnames(expr_for_clust) %in%
                           c("Endothelial", "Pericytes", "Smooth muscle", "Fibroblasts"),
                       "Vascular",
                       ifelse(colnames(expr_for_clust) %in%
                                  c("Oligodendrocytes", "Microglia", "Astrocytes", "Neurons"),
                              "CNS",
                              "Immune"))
    ann_row_df <- data.frame(Compartment = ct_group,
                              row.names = colnames(expr_for_clust))
    comp_cols <- c("Vascular" = "#E67E22", "CNS" = "#27AE60", "Immune" = "#2980B9")

    mat_t <- t(expr_for_clust)
    # Enforce display row order (vascular → CNS → immune)
    row_order <- intersect(CT_ORDER_DISPLAY, rownames(mat_t))
    mat_t <- mat_t[row_order, , drop = FALSE]
    ann_row_df <- ann_row_df[row_order, , drop = FALSE]

    show_gene_labels <- ncol(mat_t) <= 60

    # Full heatmap
    pdf_name <- glue("{COHORT}_{an$name}_walchli_hpa_heatmap.pdf")
    pdf_out  <- file.path(OUT_DIR, pdf_name)
    pdf(pdf_out,
        width  = max(9, ncol(mat_t) * 0.12 + 4),
        height = 6.5, onefile = FALSE)
    pheatmap(mat_t,
        color                = heatmap_cols,
        breaks               = seq(0, 1, length.out = 101),
        cluster_rows         = FALSE,   # keep vascular → CNS → immune ordering
        cluster_cols         = TRUE,
        show_rownames        = TRUE,
        show_colnames        = show_gene_labels,
        fontsize_row         = 9,
        fontsize_col         = 5.5,
        annotation_col       = ann_col_df,
        annotation_row       = ann_row_df,
        annotation_colors    = list(Direction    = dir_col_vals,
                                    Compartment  = comp_cols),
        annotation_names_col = FALSE,
        annotation_names_row = FALSE,
        border_color         = NA,
        gaps_row             = c(4, 8),   # gap after vascular, after CNS
        main = sprintf("%s: DEP expression (Walchli brain + HPA immune, %d genes)",
                       an$label, ncol(mat_t))
    )
    dev.off()
    cat(sprintf("  Saved %s\n", pdf_name))

    # Top-40 heatmap (for main figure panel)
    top40 <- rownames(expr_plot)[rownames(expr_plot) %in%
        deps_sig[order(-abs_logFC)][seq_len(min(40, .N)), toupper(protein)]]

    if (length(top40) >= 5) {
        expr_top40 <- expr_plot[top40, , drop = FALSE]
        expr_top40[is.na(expr_top40)] <- 0
        mat_top <- t(expr_top40)
        mat_top <- mat_top[intersect(CT_ORDER_DISPLAY, rownames(mat_top)), , drop = FALSE]
        ann_top_col <- data.frame(
            Direction = deps_sig[match(top40, toupper(protein)), direction],
            row.names = top40
        )
        ann_top_row <- ann_row_df[rownames(mat_top), , drop = FALSE]

        pdf_top_name <- glue("{COHORT}_{an$name}_walchli_hpa_top40.pdf")
        pdf_top      <- file.path(OUT_DIR, pdf_top_name)
        pdf(pdf_top, width = 10, height = 6, onefile = FALSE)
        pheatmap(mat_top,
            color                = heatmap_cols,
            breaks               = seq(0, 1, length.out = 101),
            cluster_rows         = FALSE,
            cluster_cols         = TRUE,
            show_rownames        = TRUE,
            show_colnames        = TRUE,
            fontsize_row         = 9,
            fontsize_col         = 6,
            annotation_col       = ann_top_col,
            annotation_row       = ann_top_row,
            annotation_colors    = list(Direction   = dir_col_vals,
                                        Compartment = comp_cols),
            annotation_names_col = FALSE,
            annotation_names_row = FALSE,
            border_color         = NA,
            gaps_row             = c(4, 8),
            main = sprintf("%s: Top %d DEPs - Walchli brain + HPA immune",
                           an$label, length(top40))
        )
        dev.off()
        cat(sprintf("  Saved %s\n", pdf_top_name))
    }
}

cat("\nDone.\n")
