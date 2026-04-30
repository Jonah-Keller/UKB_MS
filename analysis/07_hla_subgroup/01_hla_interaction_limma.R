#!/usr/bin/env Rscript
# =============================================================================
# 01_hla_interaction_limma.R  —  HLA×MS interaction: limma 2×2 model
# =============================================================================
# Tests whether DRB1*15:01 carrier status amplifies or attenuates the MS
# protein signal.  Uses a 4-group no-intercept parameterisation so the
# interaction term can be extracted cleanly via makeContrasts.
#
# Three cohorts:
#   all  — (pre_onset + post_onset) vs control
#   pre  — pre_onset vs control
#   post — post_onset vs control
#
# Design (per group4 cell):
#   ~ 0 + group4 + ytd + age_at_sampling + sex_f + PC1 + PC2
# PC1/PC2 are control-derived (prcomp on HC, all samples projected in).
#
# Interaction contrast (HLA×MS):
#   (MS_HLA1 - HC_HLA1) - (MS_HLA0 - HC_HLA0)
#
# Input:
#   data/ukb/olink/processed/ms_olink_qc.csv
#   data/ukb/genetics/hla_imputed.csv
#
# Output:
#   results/endophenotype/ms_hla_interaction_all.csv
#   results/endophenotype/ms_hla_interaction_pre.csv
#   results/endophenotype/ms_hla_interaction_post.csv
# =============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(data.table)
})

# --- Paths ---
.args      <- commandArgs(trailingOnly = FALSE)
.file_arg  <- .args[grepl("^--file=", .args)]
SCRIPT_DIR <- if (length(.file_arg) > 0) dirname(normalizePath(sub("^--file=", "", .file_arg))) else getwd()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

QC_FILE  <- file.path(REPO_ROOT, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
HLA_FILE <- file.path(REPO_ROOT, "data", "ukb", "genetics", "hla_imputed.csv")
OUT_DIR  <- file.path(REPO_ROOT, "results", "endophenotype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THRESH <- 0.05

# =============================================================================
# 1. Load data
# =============================================================================
message("Loading QC'd Olink data...")
if (!file.exists(QC_FILE)) stop("QC file not found: ", QC_FILE)
if (!file.exists(HLA_FILE)) stop("HLA file not found: ", HLA_FILE)

dt  <- fread(QC_FILE,  showProgress = FALSE)
hla <- fread(HLA_FILE, showProgress = FALSE)

message(sprintf("  QC: %d rows", nrow(dt)))
message(sprintf("  HLA: %d rows", nrow(hla)))

# Merge on eid
hla[, eid := as.integer(eid)]
dt[,  eid := as.integer(eid)]
dt <- merge(dt, hla[, .(eid, drb1_1501_carrier)], by = "eid", all.x = TRUE)

# Meta columns
META_COLS <- c("eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance",
               "qc_outlier", "UMAP1", "UMAP2", "mean_npx",
               "bmi", "ever_smoker", "diabetes", "alcohol_freq", "sex_num",
               "drb1_1501_carrier", "PC1", "PC2")
protein_cols <- setdiff(names(dt), META_COLS)
message(sprintf("  %d proteins", length(protein_cols)))

# =============================================================================
# 2. Filter: outliers + required covariates + HLA carrier status
# =============================================================================
dt_filt <- dt[
    qc_outlier == FALSE &
    !is.na(age_at_sampling) &
    !is.na(sex) &
    !is.na(drb1_1501_carrier)
]
message(sprintf("  After filter: %d participants", nrow(dt_filt)))
message(sprintf("  Status: %s",
    paste(names(table(dt_filt$ms_status)), table(dt_filt$ms_status),
          sep = "=", collapse = "  ")))

# =============================================================================
# 3. Protein matrix: proteins (rows) x samples (cols), mean-impute per protein
# =============================================================================
prot_mat <- t(as.matrix(dt_filt[, ..protein_cols]))
colnames(prot_mat) <- as.character(dt_filt$eid)

row_means <- rowMeans(prot_mat, na.rm = TRUE)
for (i in seq_len(nrow(prot_mat))) {
    nas <- is.na(prot_mat[i, ])
    if (any(nas)) prot_mat[i, nas] <- row_means[i]
}

message("Computing control-only PC1/PC2 (truncated SVD via irlba)...")
ctrl_eids_hla <- as.character(dt_filt[ms_status == "control", eid])
ctrl_x_hla    <- t(prot_mat[, ctrl_eids_hla])         # samples × proteins
ctrl_ctr_hla  <- colMeans(ctrl_x_hla)
ctrl_xc_hla   <- sweep(ctrl_x_hla, 2L, ctrl_ctr_hla, "-")
sv_hla        <- irlba::irlba(ctrl_xc_hla, nv = 2L)
rot_hla       <- sv_hla$v
for (k in 1:2) {
    i_max <- which.max(abs(rot_hla[, k]))
    if (rot_hla[i_max, k] < 0) rot_hla[, k] <- -rot_hla[, k]
}
total_var_hla <- sum(apply(ctrl_xc_hla, 2L, var))
pc_var_hla    <- 100 * (sv_hla$d ^ 2 / max(1L, nrow(ctrl_xc_hla) - 1L)) / total_var_hla
message(sprintf("  Control PC1: %.1f%%  PC2: %.1f%% variance explained", pc_var_hla[1], pc_var_hla[2]))
all_pcs_hla   <- sweep(t(prot_mat), 2L, ctrl_ctr_hla, "-") %*% rot_hla
for (col in c("PC1", "PC2")) if (col %in% names(dt_filt)) set(dt_filt, j = col, value = NULL)
pc_dt         <- data.table(eid = as.integer(colnames(prot_mat)),
                             PC1 = all_pcs_hla[, 1],
                             PC2 = all_pcs_hla[, 2])
dt_filt       <- merge(dt_filt, pc_dt, by = "eid", all.x = TRUE)

# =============================================================================
# Helper: run_hla_interaction()
# =============================================================================
run_hla_interaction <- function(prot_mat_sub, meta_sub,
                                ms_levels, hc_level = "control",
                                ytd_vec, label) {
    message(sprintf("\nRunning HLA interaction limma: %s  (n=%d)",
                    label, ncol(prot_mat_sub)))

    meta_sub <- copy(meta_sub)
    meta_sub[, sex_f       := factor(sex)]
    meta_sub[, hla_carrier := factor(drb1_1501_carrier)]   # 0/1
    meta_sub[, ms_binary   := as.integer(ms_status %in% ms_levels)]
    meta_sub[, ytd         := ytd_vec]

    # 4-group cell label
    meta_sub[, group4 := paste0(
        fifelse(ms_binary == 1, "MS", "HC"),
        "_HLA",
        as.character(drb1_1501_carrier)
    )]
    meta_sub[, group4 := factor(group4)]
    g_levels <- make.names(levels(meta_sub$group4))

    # Build design: no-intercept group4 + additive covariates
    meta_df <- as.data.frame(meta_sub)
    design <- model.matrix(
        ~ 0 + group4 + ytd + age_at_sampling + sex_f + PC1 + PC2,
        data = meta_df
    )
    colnames(design) <- make.names(colnames(design))

    # Map group4 column names after make.names
    hc0  <- make.names(paste0("group4HC_HLA0"))
    hc1  <- make.names(paste0("group4HC_HLA1"))
    ms0  <- make.names(paste0("group4MS_HLA0"))
    ms1  <- make.names(paste0("group4MS_HLA1"))

    # Verify all four groups exist in design columns
    needed <- c(hc0, hc1, ms0, ms1)
    missing_cols <- setdiff(needed, colnames(design))
    if (length(missing_cols) > 0) {
        message(sprintf("  WARNING: group columns missing from design: %s",
                        paste(missing_cols, collapse = ", ")))
        message(sprintf("  Available design cols: %s",
                        paste(colnames(design), collapse = ", ")))
        return(NULL)
    }

    # Interaction contrast: (MS_HLA1 - HC_HLA1) - (MS_HLA0 - HC_HLA0)
    contr_str <- sprintf("(%s - %s) - (%s - %s)", ms1, hc1, ms0, hc0)
    message(sprintf("  Contrast: %s", contr_str))

    fit   <- lmFit(prot_mat_sub, design)
    contr <- makeContrasts(contrasts = contr_str, levels = design)
    fit2  <- contrasts.fit(fit, contr)
    fit2  <- eBayes(fit2)

    res <- as.data.frame(topTable(fit2, sort.by = "P", n = Inf))
    res$protein    <- rownames(res)
    res$comparison <- label
    res$fdr        <- p.adjust(res$P.Value, method = "BH")

    n_sig <- sum(res$fdr < FDR_THRESH)
    message(sprintf("  Interaction FDR<%.2f: %d proteins", FDR_THRESH, n_sig))

    res <- res[, c("protein", "logFC", "AveExpr", "t", "P.Value", "fdr", "comparison")]
    return(res)
}

# =============================================================================
# 4. Cohort definitions and runs
# =============================================================================

# --- All: (pre + post) vs control ---
message("\n--- Cohort: all (pre+post vs control) ---")
dt_all  <- dt_filt[ms_status %in% c("pre_onset", "post_onset", "control")]
eids_all <- as.character(dt_all$eid)
prot_all <- prot_mat[, eids_all]
ytd_all  <- dt_all$years_to_diagnosis
ytd_all[is.na(ytd_all)] <- 0

res_all <- run_hla_interaction(
    prot_mat_sub = prot_all,
    meta_sub     = dt_all,
    ms_levels    = c("pre_onset", "post_onset"),
    ytd_vec      = ytd_all,
    label        = "hla_interaction_all"
)

# --- Pre: pre_onset vs control ---
message("\n--- Cohort: pre (pre_onset vs control) ---")
dt_pre  <- dt_filt[ms_status %in% c("pre_onset", "control")]
eids_pre <- as.character(dt_pre$eid)
prot_pre <- prot_mat[, eids_pre]
ytd_pre  <- dt_pre$years_to_diagnosis
ytd_pre[is.na(ytd_pre)] <- 0

res_pre <- run_hla_interaction(
    prot_mat_sub = prot_pre,
    meta_sub     = dt_pre,
    ms_levels    = "pre_onset",
    ytd_vec      = ytd_pre,
    label        = "hla_interaction_pre"
)

# --- Post: post_onset vs control ---
message("\n--- Cohort: post (post_onset vs control) ---")
dt_post  <- dt_filt[ms_status %in% c("post_onset", "control")]
eids_post <- as.character(dt_post$eid)
prot_post <- prot_mat[, eids_post]
ytd_post  <- abs(dt_post$years_to_diagnosis)
ytd_post[is.na(ytd_post)] <- 0

res_post <- run_hla_interaction(
    prot_mat_sub = prot_post,
    meta_sub     = dt_post,
    ms_levels    = "post_onset",
    ytd_vec      = ytd_post,
    label        = "hla_interaction_post"
)

# =============================================================================
# 5. Write outputs
# =============================================================================
if (!is.null(res_all)) {
    fwrite(res_all,  file.path(OUT_DIR, "ms_hla_interaction_all.csv"))
    message("\nWrote ms_hla_interaction_all.csv  — ",
            sum(res_all$fdr < FDR_THRESH), " interaction FDR<0.05")
}
if (!is.null(res_pre)) {
    fwrite(res_pre,  file.path(OUT_DIR, "ms_hla_interaction_pre.csv"))
    message("Wrote ms_hla_interaction_pre.csv  — ",
            sum(res_pre$fdr < FDR_THRESH), " interaction FDR<0.05")
}
if (!is.null(res_post)) {
    fwrite(res_post, file.path(OUT_DIR, "ms_hla_interaction_post.csv"))
    message("Wrote ms_hla_interaction_post.csv — ",
            sum(res_post$fdr < FDR_THRESH), " interaction FDR<0.05")
}

message("\nHLA interaction limma complete.")
