#!/usr/bin/env Rscript
# =============================================================================
# 01_olink_qc.R  —  Olink NPX QC for UKB MS/ALS proteomics replication
# =============================================================================
# Steps (following Chia 2025 / template notebooks):
#   1. Load UKB Olink instance-0 NPX (wide format, 2923 proteins)
#   2. Protein-level QC: flag proteins with >20% missing
#   3. Sample-level QC: IQR outlier detection (per-sample Q1/Q3 across proteins)
#   4. Merge with MS and ALS cohort labels + age/sex
#   5. UMAP on protein data (Chia 2025: UMAP1/2 used as limma covariates)
#   6. Write per-cohort QC'd files to data/ukb/olink/processed/
#
# Outputs:
#   results/qc/ms_npx_distribution.pdf
#   results/qc/als_npx_distribution.pdf
#   results/qc/protein_missingness.pdf
#   data/ukb/olink/processed/ms_olink_qc.csv   (samples × proteins + meta + UMAP)
#   data/ukb/olink/processed/als_olink_qc.csv
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(umap)
})

# --- Paths ---
.args       <- commandArgs(trailingOnly = FALSE)
.file_arg   <- .args[grepl("^--file=", .args)]
SCRIPT_DIR  <- if (length(.file_arg) > 0) dirname(normalizePath(sub("^--file=", "", .file_arg))) else getwd()
REPO_ROOT   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
CADASIL_DIR <- file.path(dirname(REPO_ROOT), "CADASIL_Proteome_ML_Keller_2024_Rebuttal")
DATA_DIR    <- file.path(CADASIL_DIR, "data", "ukb")
OUT_QC      <- file.path(REPO_ROOT, "results", "qc")
OUT_OLINK   <- file.path(REPO_ROOT, "data", "ukb", "olink", "processed")
dir.create(OUT_QC,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_OLINK, showWarnings = FALSE, recursive = TRUE)

OLINK_PATH    <- file.path(DATA_DIR, "olink", "i0", "olink_instance_0_extracted_data.csv")
MS_COHORT_PATH  <- file.path(REPO_ROOT, "data", "ukb", "cohort", "ms_cohort.csv")
ALS_COHORT_PATH <- file.path(REPO_ROOT, "data", "ukb", "cohort", "als_cohort.csv")

# QC thresholds
MISSING_PROTEIN_THRESH <- 0.20  # exclude proteins with >20% missing
IQR_OUTLIER_FACTOR     <- 1.5   # flag samples outside Q1/Q3 ± 1.5×IQR
UMAP_N_NEIGHBORS       <- 15
UMAP_N_PCS             <- 10    # PCs of protein matrix fed into UMAP

set.seed(42)

source(file.path(REPO_ROOT, "analysis", "helpers", "ukb_theme.R"))

# =============================================================================
# 1. Load Olink NPX (wide format)
# =============================================================================
message("Loading Olink NPX (instance 0)...")
if (!file.exists(OLINK_PATH)) stop("Olink file not found: ", OLINK_PATH)
olink_raw <- fread(OLINK_PATH, showProgress = FALSE)
setnames(olink_raw, names(olink_raw),
         str_replace(names(olink_raw), "^olink_instance_0\\.", ""))
setnames(olink_raw, "eid", "eid")   # already named eid after strip

protein_cols <- setdiff(names(olink_raw), "eid")
message(sprintf("  %d participants, %d proteins", nrow(olink_raw), length(protein_cols)))

# =============================================================================
# 2. Protein-level QC: missingness
# =============================================================================
message("Protein-level QC: missingness...")
prot_mat <- as.matrix(olink_raw[, ..protein_cols])
prot_missing_frac <- colMeans(is.na(prot_mat))

n_exclude <- sum(prot_missing_frac > MISSING_PROTEIN_THRESH)
message(sprintf("  Proteins with >%.0f%% missing: %d (excluded)",
                MISSING_PROTEIN_THRESH * 100, n_exclude))
keep_proteins <- names(prot_missing_frac)[prot_missing_frac <= MISSING_PROTEIN_THRESH]

# Plot protein missingness distribution
miss_df <- data.frame(pct_missing = prot_missing_frac * 100)
p_miss <- ggplot(miss_df, aes(x = pct_missing)) +
    geom_histogram(bins = 50, fill = UKB_PALETTE[3], colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = MISSING_PROTEIN_THRESH * 100,
               colour = UKB_PALETTE[4], linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = MISSING_PROTEIN_THRESH * 100 + 0.5,
             y = Inf, vjust = 1.5, hjust = 0, size = 3, colour = UKB_PALETTE[4],
             label = sprintf("Threshold = %.0f%%", MISSING_PROTEIN_THRESH * 100)) +
    labs(x = "% Missing per protein", y = "# proteins",
         title = sprintf("Protein missingness (n=%d proteins)", length(protein_cols))) +
    theme_ukb()
ggsave(file.path(OUT_QC, "protein_missingness.pdf"), p_miss,
       width = 7, height = 4, device = cairo_pdf)

# Restrict to QC-passing proteins
olink_qc <- olink_raw[, c("eid", keep_proteins), with = FALSE]
message(sprintf("  Retaining %d proteins after missingness filter", length(keep_proteins)))

# =============================================================================
# 3. Sample-level QC: IQR outlier detection (per Chia 2025 / template)
#    For each sample: compute Q1, Q3 of its NPX distribution across proteins.
#    Flag samples where Q1 or Q3 is outside the population-level IQR bounds.
# =============================================================================
message("Sample-level QC: IQR outlier detection...")

prot_mat_qc <- as.matrix(olink_qc[, ..keep_proteins])

# Per-sample distribution statistics
sample_stats <- data.frame(
    eid         = olink_qc$eid,
    mean_npx    = rowMeans(prot_mat_qc, na.rm = TRUE),
    sd_npx      = apply(prot_mat_qc, 1, sd, na.rm = TRUE),
    lower_iqr   = apply(prot_mat_qc, 1, quantile, 0.25, na.rm = TRUE),
    upper_iqr   = apply(prot_mat_qc, 1, quantile, 0.75, na.rm = TRUE),
    iqr_npx     = apply(prot_mat_qc, 1, IQR, na.rm = TRUE)
)

# Population-level thresholds for lower_iqr and upper_iqr
Q1_lower <- quantile(sample_stats$lower_iqr, 0.25, na.rm = TRUE)
Q3_lower <- quantile(sample_stats$lower_iqr, 0.75, na.rm = TRUE)
Q1_upper <- quantile(sample_stats$upper_iqr, 0.25, na.rm = TRUE)
Q3_upper <- quantile(sample_stats$upper_iqr, 0.75, na.rm = TRUE)

min_thresh <- Q1_lower - IQR_OUTLIER_FACTOR * (Q3_lower - Q1_lower)
max_thresh <- Q3_upper + IQR_OUTLIER_FACTOR * (Q3_upper - Q1_upper)

sample_stats$qc_outlier <- (sample_stats$lower_iqr < min_thresh) |
                           (sample_stats$upper_iqr > max_thresh)

n_outliers <- sum(sample_stats$qc_outlier, na.rm = TRUE)
message(sprintf("  %d sample-level outliers flagged (IQR × %.1f)", n_outliers, IQR_OUTLIER_FACTOR))

# Merge QC flags + mean_npx back
olink_qc <- merge(olink_qc, sample_stats[, c("eid", "qc_outlier", "mean_npx")], by = "eid")

# =============================================================================
# 4. UMAP on protein data (all non-outlier samples)
#    Used as covariates in limma to capture non-specific variation (Chia 2025)
# =============================================================================
message("Computing UMAP on protein data...")
non_outlier_idx <- which(!olink_qc$qc_outlier)
prot_mat_for_umap <- as.matrix(olink_qc[non_outlier_idx, ..keep_proteins])

# Impute missing with column medians (required for PCA/UMAP)
col_medians <- apply(prot_mat_for_umap, 2, median, na.rm = TRUE)
for (j in seq_along(col_medians)) {
    nas <- is.na(prot_mat_for_umap[, j])
    if (any(nas)) prot_mat_for_umap[nas, j] <- col_medians[j]
}

# PCA → top PCs → UMAP
pca_res <- prcomp(prot_mat_for_umap, center = TRUE, scale. = TRUE,
                  rank. = UMAP_N_PCS)
umap_res <- umap(pca_res$x, n_neighbors = UMAP_N_NEIGHBORS,
                 random_state = 42)

umap_df <- data.frame(
    eid   = olink_qc$eid[non_outlier_idx],
    UMAP1 = umap_res$layout[, 1],
    UMAP2 = umap_res$layout[, 2]
)

# Assign NA UMAP for outlier samples
all_umap <- data.frame(eid = olink_qc$eid)
all_umap <- merge(all_umap, umap_df, by = "eid", all.x = TRUE)

olink_final <- merge(olink_qc, all_umap, by = "eid")
message(sprintf("  UMAP computed on %d non-outlier samples", nrow(umap_df)))

# =============================================================================
# 5. Per-cohort merge and output
# =============================================================================
process_cohort <- function(disease_name, cohort_path, status_col) {
    message(sprintf("\nProcessing %s cohort...", disease_name))

    if (!file.exists(cohort_path)) {
        message("  Cohort file not found: ", cohort_path)
        return(invisible(NULL))
    }
    cohort <- fread(cohort_path)

    # Merge Olink with cohort (inner join: only participants with both)
    dt <- merge(cohort, olink_final, by = "eid", all = FALSE)
    message(sprintf("  %d participants after inner join with Olink+cohort", nrow(dt)))

    if (nrow(dt) == 0) {
        message("  No overlap — skipping.")
        return(invisible(NULL))
    }

    # Summarise cases vs controls
    status_counts <- table(dt[[status_col]], useNA = "ifany")
    for (s in names(status_counts)) {
        message(sprintf("    %-15s: %d", s, status_counts[s]))
    }
    n_outliers <- sum(dt$qc_outlier, na.rm = TRUE)
    message(sprintf("  IQR outliers included (flagged, not removed): %d", n_outliers))

    # Plot global NPX median by status
    dt_plot <- dt[, c(status_col, "mean_npx", "qc_outlier"), with = FALSE]
    setnames(dt_plot, status_col, "status")
    p <- ggplot(dt_plot, aes(x = status, y = mean_npx, fill = status)) +
        geom_violin(alpha = 0.6, linewidth = 0.3) +
        geom_boxplot(width = 0.2, outlier.size = 0.5, linewidth = 0.35) +
        scale_fill_manual(values = setNames(
            UKB_PALETTE[seq_len(length(unique(dt_plot$status)))],
            sort(unique(dt_plot$status)))) +
        labs(title = sprintf("%s — Global NPX distribution by status", disease_name),
             y = "Mean NPX across proteins", x = NULL) +
        theme_ukb() +
        theme(legend.position = "none")
    plot_path <- file.path(OUT_QC, sprintf("%s_npx_distribution.pdf",
                                           tolower(disease_name)))
    ggsave(plot_path, p, width = 7, height = 5, device = cairo_pdf)
    message(sprintf("  QC plot saved: %s", basename(plot_path)))

    # Write output (keep qc_outlier flag; downstream scripts can filter on it)
    out_path <- file.path(OUT_OLINK,
                          sprintf("%s_olink_qc.csv", tolower(disease_name)))
    fwrite(dt, out_path)
    message(sprintf("  Written: %s  (%d rows × %d cols)",
                    basename(out_path), nrow(dt), ncol(dt)))
    return(invisible(dt))
}

process_cohort("MS",  MS_COHORT_PATH,  "ms_status")
process_cohort("ALS", ALS_COHORT_PATH, "als_status")

message("\nOlink QC complete.")
