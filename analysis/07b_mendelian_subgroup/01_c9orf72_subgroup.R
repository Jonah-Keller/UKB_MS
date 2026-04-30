#!/usr/bin/env Rscript
# =============================================================================
# 01_c9orf72_subgroup.R — Mendelian carrier × ALS interaction (limma 2x2)
# =============================================================================
# Tests whether the configured Mendelian carrier status (C9orf72 hexanucleotide
# expansion >= 30 repeats by default) modifies the ALS proteomic signature.
# Same 2x2 no-intercept design as analysis/07_hla_subgroup/01_hla_interaction_limma.R
# but stratifying on cfg$mendelian_carrier_col instead of cfg$hla_carrier_col.
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
# Interaction contrast (carrier x case):
#   (CASE_C9pos - HC_C9pos) - (CASE_C9neg - HC_C9neg)
#
# Sex × C9orf72 sub-model (per cfg$sex_handling$interaction_terms_to_test):
#   ~ 0 + group4 * sex_f + ytd + age_at_sampling + PC1 + PC2
#
# Bails out gracefully if the carrier column is not present in the data
# (i.e. the C9orf72 extraction has not yet been run on UKB-RAP).
# =============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(data.table)
    library(glue)
    library(here)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "limma_utils.R"))
cfg <- load_disease_config()

CARRIER_COL    <- cfg$mendelian_carrier_col
STRICT_CARRIER <- paste0(CARRIER_COL, "_strict")
STATUS_COL     <- cfg$cohort_status_col
PRE_ONSET      <- cfg$status_values$pre_onset
POST_ONSET     <- cfg$status_values$post_onset
CONTROL        <- cfg$status_values$control

# --- Paths ---
QC_FILE  <- here("data", "ukb", "olink", "processed",
                 glue("{cfg$cohort_short}_olink_qc.csv"))
MND_FILE <- here("data", "ukb", "genetics", "c9orf72_expansion.csv")
OUT_DIR  <- here("results", "mendelian")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THRESH <- 0.05

# =============================================================================
# 1. Load data — bail out gracefully if extracts missing
# =============================================================================
message("Loading QC'd Olink data...")
if (!file.exists(QC_FILE)) stop("QC file not found: ", QC_FILE)

dt <- fread(QC_FILE, showProgress = FALSE)
message(sprintf("  QC: %d rows", nrow(dt)))

if (!file.exists(MND_FILE)) {
    message(sprintf("\nMendelian-carrier stage skipped: %s does not exist.",
                    basename(MND_FILE)))
    message("  See EXTRACT_TODO.md section 1 for the C9orf72 extraction recipe.")
    quit(save = "no", status = 0)
}

mnd <- fread(MND_FILE, showProgress = FALSE)
message(sprintf("  Mendelian carrier file: %d rows", nrow(mnd)))

if (!CARRIER_COL %in% names(mnd)) {
    message(sprintf("\nMendelian-carrier stage skipped: column '%s' not in %s",
                    CARRIER_COL, basename(MND_FILE)))
    message(sprintf("  Available columns: %s",
                    paste(names(mnd), collapse = ", ")))
    quit(save = "no", status = 0)
}

dt[, eid := as.integer(eid)]
mnd[, eid := as.integer(eid)]
keep_cols <- intersect(c("eid", CARRIER_COL, STRICT_CARRIER), names(mnd))
dt <- merge(dt, mnd[, ..keep_cols], by = "eid", all.x = TRUE)

# =============================================================================
# 2. Filter: outliers + required covariates + carrier status
# =============================================================================
dt_filt <- dt[
    qc_outlier == FALSE &
    !is.na(age_at_sampling) &
    !is.na(sex) &
    !is.na(get(CARRIER_COL))
]
message(sprintf("  After filter: %d participants", nrow(dt_filt)))

status_tab <- table(dt_filt[[STATUS_COL]], dt_filt[[CARRIER_COL]],
                    useNA = "ifany",
                    dnn = c(STATUS_COL, CARRIER_COL))
print(status_tab)

n_case_carrier <- sum(dt_filt[[STATUS_COL]] %in% c(PRE_ONSET, POST_ONSET) &
                       dt_filt[[CARRIER_COL]] == 1, na.rm = TRUE)
n_case_total   <- sum(dt_filt[[STATUS_COL]] %in% c(PRE_ONSET, POST_ONSET),
                       na.rm = TRUE)
message(sprintf("  Case carrier rate: %d / %d = %.1f%%",
                n_case_carrier, n_case_total,
                100 * n_case_carrier / max(1L, n_case_total)))

if (n_case_carrier < 5L) {
    message(sprintf("\nMendelian-carrier stage halted: only %d case carriers.\n",
                    n_case_carrier),
            "  Need >=5 to run a meaningful 2x2 limma interaction test.\n",
            "  Inspect the call quality of c9orf72_expansion.csv first.")
    quit(save = "no", status = 0)
}

# =============================================================================
# 3. Protein matrix + control-derived PC1/PC2
# =============================================================================
META_COLS <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance",
               "qc_outlier", "UMAP1", "UMAP2", "mean_npx",
               "bmi", "ever_smoker", "diabetes", "alcohol_freq", "sex_num",
               CARRIER_COL, STRICT_CARRIER, "PC1", "PC2")
protein_cols <- setdiff(names(dt_filt), META_COLS)
message(sprintf("  %d proteins", length(protein_cols)))

message("Computing control-derived PC1/PC2...")
for (col in c("PC1", "PC2"))
    if (col %in% names(dt_filt)) set(dt_filt, j = col, value = NULL)
pc_dt <- compute_control_pcs(dt_filt, protein_cols,
                             control_value = CONTROL,
                             status_col    = STATUS_COL)
dt_filt <- merge(dt_filt, pc_dt, by = "eid", all.x = TRUE)
dt_filt <- dt_filt[!is.na(PC1) & !is.na(PC2)]

# =============================================================================
# 4. 2x2 limma per cohort (all / pre / post)
# =============================================================================
prot_mat <- t(as.matrix(dt_filt[, ..protein_cols]))
colnames(prot_mat) <- as.character(dt_filt$eid)
row_means <- rowMeans(prot_mat, na.rm = TRUE)
for (i in seq_len(nrow(prot_mat))) {
    nas <- is.na(prot_mat[i, ])
    if (any(nas)) prot_mat[i, nas] <- row_means[i]
}

run_cohort <- function(case_status_filter, label) {
    keep_eids <- as.character(dt_filt[
        get(STATUS_COL) %in% c(CONTROL, case_status_filter), eid])
    if (length(keep_eids) < 30L) {
        message(sprintf("  [%s] insufficient sample (n=%d), skipping",
                        label, length(keep_eids)))
        return(NULL)
    }
    sub <- dt_filt[as.character(eid) %in% keep_eids]
    sub[, group4 := factor(fcase(
        get(STATUS_COL) == CONTROL & get(CARRIER_COL) == 0, "HC_NEG",
        get(STATUS_COL) == CONTROL & get(CARRIER_COL) == 1, "HC_POS",
        get(STATUS_COL) %in% case_status_filter & get(CARRIER_COL) == 0, "CASE_NEG",
        get(STATUS_COL) %in% case_status_filter & get(CARRIER_COL) == 1, "CASE_POS"
    ), levels = c("HC_NEG", "HC_POS", "CASE_NEG", "CASE_POS"))]
    sub[, sex_f := factor(sex)]
    sub[, ytd   := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]

    table_g <- table(sub$group4)
    message(sprintf("  [%s] group4 cells: %s", label,
                    paste(names(table_g), table_g, sep = "=", collapse = " | ")))
    if (any(table_g < 3L)) {
        message(sprintf("  [%s] one or more group4 cells <3 — skipping interaction",
                        label))
        return(NULL)
    }

    keep_idx <- as.character(sub$eid)
    expr <- prot_mat[, keep_idx]

    design <- model.matrix(~ 0 + group4 + ytd + age_at_sampling + sex_f + PC1 + PC2,
                           data = sub)
    fit <- lmFit(expr, design)
    contrast_int <- makeContrasts(
        interaction = (group4CASE_POS - group4HC_POS) - (group4CASE_NEG - group4HC_NEG),
        levels = design)
    fit2 <- contrasts.fit(fit, contrast_int)
    fit2 <- eBayes(fit2)
    res <- topTable_safe(fit2, coef = "interaction", n = Inf)
    res$cohort <- label
    res
}

write_result <- function(res, label) {
    if (is.null(res) || nrow(res) == 0) return(invisible(NULL))
    out <- file.path(OUT_DIR,
                     glue("{cfg$cohort_short}_{tolower(CARRIER_COL)}_interaction_{label}.csv"))
    fwrite(res, out)
    n_sig <- sum(res$adj.P.Val < FDR_THRESH, na.rm = TRUE)
    message(sprintf("  [%s] %d FDR-sig interaction DEPs -> %s",
                    label, n_sig, basename(out)))
}

message("\n--- Running 2x2 limma per cohort ---")
res_all  <- run_cohort(c(PRE_ONSET, POST_ONSET), "all")
res_pre  <- run_cohort(PRE_ONSET,                "pre")
res_post <- run_cohort(POST_ONSET,               "post")
write_result(res_all,  "all")
write_result(res_pre,  "pre")
write_result(res_post, "post")

# =============================================================================
# 5. Sex × carrier sub-model (per cfg$sex_handling)
# =============================================================================
# Generic match: any cfg interaction term of the form sex_x_<col> where
# <col> matches CARRIER_COL or the literal "mendelian_carrier" alias.
.iact_terms <- cfg$sex_handling$interaction_terms_to_test %||% character(0)
.iact_targets <- sub("^sex_x_", "", grep("^sex_x_", .iact_terms, value = TRUE))
do_sex_x_carrier <- any(.iact_targets %in% c(CARRIER_COL, "mendelian_carrier"))

if (do_sex_x_carrier && !is.null(res_all)) {
    message("\n--- Sex x carrier sub-model (combined cases) ---")
    sub <- dt_filt[get(STATUS_COL) %in% c(CONTROL, PRE_ONSET, POST_ONSET)]
    sub[, is_case  := as.integer(get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET))]
    sub[, carrier  := as.integer(get(CARRIER_COL))]
    sub[, sex_f    := factor(sex)]
    sub[, ytd      := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
    keep_idx <- as.character(sub$eid)
    expr <- prot_mat[, keep_idx]
    design <- model.matrix(~ is_case * carrier * sex_f +
                              ytd + age_at_sampling + PC1 + PC2,
                           data = sub)
    fit <- lmFit(expr, design)
    fit <- eBayes(fit)
    coef_name <- "is_case:carrier:sex_f1"
    if (coef_name %in% colnames(fit$coefficients)) {
        res_sx <- topTable_safe(fit, coef = coef_name, n = Inf)
        out <- file.path(OUT_DIR, glue(
            "{cfg$cohort_short}_{tolower(CARRIER_COL)}_sex_3way.csv"))
        fwrite(res_sx, out)
        n_sig <- sum(res_sx$adj.P.Val < FDR_THRESH, na.rm = TRUE)
        message(sprintf("  3-way sex×case×carrier: %d FDR-sig -> %s",
                        n_sig, basename(out)))
    } else {
        message("  is_case:carrier:sex_f1 not in design (probably collinear)")
    }
}

message("\n01_c9orf72_subgroup.R complete.")
