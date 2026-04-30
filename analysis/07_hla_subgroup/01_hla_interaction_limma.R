#!/usr/bin/env Rscript
# =============================================================================
# 01_hla_interaction_limma.R  —  HLA x disease interaction: limma 2x2 model
# =============================================================================
# Tests whether the HLA risk-allele carrier status amplifies or attenuates the
# disease protein signal.  Uses a 4-group no-intercept parameterisation so the
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
# Interaction contrast (HLA x case):
#   (CASE_HLA1 - HC_HLA1) - (CASE_HLA0 - HC_HLA0)
#
# Disease-, HLA-, and status-column constants come from configs/disease.yaml.
# =============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(data.table)
    library(glue)
    library(here)
})

source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

HLA_CARRIER_COL <- cfg$hla_carrier_col
STATUS_COL      <- cfg$cohort_status_col
PRE_ONSET       <- cfg$status_values$pre_onset
POST_ONSET      <- cfg$status_values$post_onset
CONTROL         <- cfg$status_values$control

# --- Paths ---
.args      <- commandArgs(trailingOnly = FALSE)
.file_arg  <- .args[grepl("^--file=", .args)]
SCRIPT_DIR <- if (length(.file_arg) > 0) dirname(normalizePath(sub("^--file=", "", .file_arg))) else getwd()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

QC_FILE  <- file.path(REPO_ROOT, "data", "ukb", "olink", "processed",
                      glue("{cfg$cohort_short}_olink_qc.csv"))
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
dt <- merge(dt, hla[, .SD, .SDcols = c("eid", HLA_CARRIER_COL)],
            by = "eid", all.x = TRUE)

# Meta columns
META_COLS <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance",
               "qc_outlier", "UMAP1", "UMAP2", "mean_npx",
               "bmi", "ever_smoker", "diabetes", "alcohol_freq", "sex_num",
               HLA_CARRIER_COL, "PC1", "PC2")
protein_cols <- setdiff(names(dt), META_COLS)
message(sprintf("  %d proteins", length(protein_cols)))

# =============================================================================
# 2. Filter: outliers + required covariates + HLA carrier status
# =============================================================================
dt_filt <- dt[
    qc_outlier == FALSE &
    !is.na(age_at_sampling) &
    !is.na(sex) &
    !is.na(get(HLA_CARRIER_COL))
]
message(sprintf("  After filter: %d participants", nrow(dt_filt)))
status_tab <- table(dt_filt[[STATUS_COL]])
message(sprintf("  Status: %s",
    paste(names(status_tab), status_tab, sep = "=", collapse = "  ")))

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
ctrl_eids_hla <- as.character(dt_filt[get(STATUS_COL) == CONTROL, eid])
ctrl_x_hla    <- t(prot_mat[, ctrl_eids_hla])         # samples x proteins
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
                                case_levels, hc_level = CONTROL,
                                ytd_vec, label) {
    message(sprintf("\nRunning HLA interaction limma: %s  (n=%d)",
                    label, ncol(prot_mat_sub)))

    meta_sub <- copy(meta_sub)
    meta_sub[, sex_f       := factor(sex)]
    meta_sub[, hla_carrier := factor(get(HLA_CARRIER_COL))]   # 0/1
    meta_sub[, case_binary := as.integer(get(STATUS_COL) %in% case_levels)]
    meta_sub[, ytd         := ytd_vec]

    # 4-group cell label
    meta_sub[, group4 := paste0(
        fifelse(case_binary == 1, "CASE", "HC"),
        "_HLA",
        as.character(get(HLA_CARRIER_COL))
    )]
    meta_sub[, group4 := factor(group4)]

    # Build design: no-intercept group4 + additive covariates
    meta_df <- as.data.frame(meta_sub)
    design <- model.matrix(
        ~ 0 + group4 + ytd + age_at_sampling + sex_f + PC1 + PC2,
        data = meta_df
    )
    colnames(design) <- make.names(colnames(design))

    # Map group4 column names after make.names
    hc0 <- make.names("group4HC_HLA0")
    hc1 <- make.names("group4HC_HLA1")
    ms0 <- make.names("group4CASE_HLA0")
    ms1 <- make.names("group4CASE_HLA1")

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

    # Interaction contrast: (CASE_HLA1 - HC_HLA1) - (CASE_HLA0 - HC_HLA0)
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

run_cohort <- function(label, status_levels, case_levels, abs_ytd = FALSE) {
    message(sprintf("\n--- Cohort: %s ---", label))
    sub <- dt_filt[get(STATUS_COL) %in% status_levels]
    eids <- as.character(sub$eid)
    ytd <- sub$years_to_diagnosis
    if (abs_ytd) ytd <- abs(ytd)
    ytd[is.na(ytd)] <- 0
    run_hla_interaction(
        prot_mat_sub = prot_mat[, eids],
        meta_sub     = sub,
        case_levels  = case_levels,
        ytd_vec      = ytd,
        label        = sprintf("hla_interaction_%s", label)
    )
}

res_all  <- run_cohort("all",  c(PRE_ONSET, POST_ONSET, CONTROL),
                       c(PRE_ONSET, POST_ONSET))
res_pre  <- run_cohort("pre",  c(PRE_ONSET, CONTROL),  PRE_ONSET)
res_post <- run_cohort("post", c(POST_ONSET, CONTROL), POST_ONSET, abs_ytd = TRUE)

# =============================================================================
# 5. Write outputs
# =============================================================================
write_result <- function(res, suffix) {
    if (is.null(res)) return(invisible(NULL))
    fname <- glue("{cfg$cohort_short}_hla_interaction_{suffix}.csv")
    fwrite(res, file.path(OUT_DIR, fname))
    message(sprintf("Wrote %s — %d interaction FDR<%.2f",
                    fname, sum(res$fdr < FDR_THRESH), FDR_THRESH))
}

write_result(res_all,  "all")
write_result(res_pre,  "pre")
write_result(res_post, "post")

message("\nHLA interaction limma complete.")
