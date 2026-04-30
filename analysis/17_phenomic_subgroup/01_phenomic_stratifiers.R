#!/usr/bin/env Rscript
# =============================================================================
# 01_phenomic_stratifiers.R — phenomic / exposure interaction with disease
# =============================================================================
# For each entry in cfg$phenomic_stratifiers with use_as == "stratifier",
# fit a 2x2 limma interaction (carrier x case) on the QC'd Olink + covariate
# data.  Same no-intercept four-cell design as 07_hla_subgroup and 07b but
# with cohort-level phenomic stratifiers (smoking, head trauma, etc.) instead
# of genetic carriers.
#
# Sex is included as a covariate per cfg$sex_handling$primary_use, plus a
# dedicated sex × stratifier interaction sub-model when listed in
# cfg$sex_handling$interaction_terms_to_test.
#
# Bails out gracefully:
#   - if cfg$phenomic_stratifiers is null or empty: skip whole stage
#   - per stratifier: if the column is missing from QC + covariates, skip
#     that one and continue to the next
#   - per stratifier: if <5 carriers among ALS cases, skip (under-powered)
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

STATUS_COL <- cfg$cohort_status_col
PRE_ONSET  <- cfg$status_values$pre_onset
POST_ONSET <- cfg$status_values$post_onset
CONTROL    <- cfg$status_values$control

QC_FILE  <- here("data", "ukb", "olink", "processed",
                 glue("{cfg$cohort_short}_olink_qc.csv"))
COV_FILE <- here("data", "ukb", "covariates", "covariate_table.csv")
OUT_DIR  <- here("results", "phenomic")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THRESH <- 0.05

# =============================================================================
# Skip stage entirely if no stratifiers configured
# =============================================================================
strats <- cfg$phenomic_stratifiers
if (is.null(strats) || length(strats) == 0L) {
    message("Phenomic-stratifier stage skipped: cfg$phenomic_stratifiers is empty.")
    quit(save = "no", status = 0)
}

stratifier_keys <- names(strats)[
    vapply(strats, function(s) {
        if (!is.list(s)) return(FALSE)   # documentation strings (e.g. "NOT_CAPTURED")
        identical(s$use_as, "stratifier")
    }, logical(1))
]
if (length(stratifier_keys) == 0L) {
    message("Phenomic-stratifier stage skipped: no entries with use_as: stratifier.")
    quit(save = "no", status = 0)
}
message(sprintf("Configured stratifiers: %s",
                paste(stratifier_keys, collapse = ", ")))

# =============================================================================
# Load data
# =============================================================================
message("Loading QC'd Olink data...")
if (!file.exists(QC_FILE)) stop("QC file not found: ", QC_FILE)
dt <- fread(QC_FILE, showProgress = FALSE)

if (file.exists(COV_FILE)) {
    cov <- fread(COV_FILE, showProgress = FALSE)
    message(sprintf("  Covariate table: %d rows, %d cols",
                    nrow(cov), ncol(cov)))
    keep_cov <- intersect(c("eid", stratifier_keys, "bmi", "head_trauma_history",
                            "smoking_status", "ever_smoker"),
                          names(cov))
    cov_sub <- cov[, ..keep_cov]
    dt[, eid := as.integer(eid)]
    cov_sub[, eid := as.integer(eid)]
    new_cols <- setdiff(names(cov_sub), names(dt))
    if (length(new_cols) > 0L) {
        dt <- merge(dt, cov_sub[, c("eid", new_cols), with = FALSE],
                    by = "eid", all.x = TRUE)
    }
} else {
    message("  Covariate file not found; only QC columns available.")
}

META_COLS <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance", "qc_outlier",
               "UMAP1", "UMAP2", "mean_npx", "bmi", "ever_smoker",
               "smoking_status", "head_trauma_history", "diabetes",
               "alcohol_freq", "sex_num",
               cfg$hla_carrier_col, cfg$hla_dosage_col,
               cfg$mendelian_carrier_col, "PC1", "PC2")
META_COLS <- intersect(META_COLS, names(dt))
protein_cols <- setdiff(names(dt), META_COLS)
message(sprintf("  %d proteins, %d meta columns", length(protein_cols),
                length(META_COLS)))

# =============================================================================
# Filter base
# =============================================================================
dt_filt <- dt[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]
message(sprintf("  After base filter: %d participants", nrow(dt_filt)))

# Compute control-derived PCs once (shared across stratifiers)
message("Computing control-derived PC1/PC2...")
for (col in c("PC1", "PC2"))
    if (col %in% names(dt_filt)) set(dt_filt, j = col, value = NULL)
pc_dt <- compute_control_pcs(dt_filt, protein_cols,
                             control_value = CONTROL,
                             status_col    = STATUS_COL)
dt_filt <- merge(dt_filt, pc_dt, by = "eid", all.x = TRUE)
dt_filt <- dt_filt[!is.na(PC1) & !is.na(PC2)]

prot_mat <- t(as.matrix(dt_filt[, ..protein_cols]))
colnames(prot_mat) <- as.character(dt_filt$eid)
row_means <- rowMeans(prot_mat, na.rm = TRUE)
for (i in seq_len(nrow(prot_mat))) {
    nas <- is.na(prot_mat[i, ])
    if (any(nas)) prot_mat[i, nas] <- row_means[i]
}

# =============================================================================
# Per-stratifier 2x2 limma
# =============================================================================
# Generic match: each cfg interaction term of the form "sex_x_<col>" enables
# the 3-way sex × case × <col> sub-model when that stratifier is run.
.iact_terms <- cfg$sex_handling$interaction_terms_to_test %||% character(0)
.iact_targets <- sub("^sex_x_", "", grep("^sex_x_", .iact_terms, value = TRUE))

run_strat <- function(strat_key) {
    cat(sprintf("\n=== Stratifier: %s ===\n", strat_key))

    strat_col <- strat_key
    if (!strat_col %in% names(dt_filt)) {
        cat(sprintf("  Skipped: column '%s' not in QC + covariates.\n", strat_col))
        cat("  Run the corresponding extraction (see EXTRACT_TODO.md).\n")
        return(NULL)
    }

    sub <- dt_filt[!is.na(get(strat_col))]
    if (is.character(sub[[strat_col]])) {
        # smoking_status: levels never/former/current → binary "ever vs never"
        # (more powerful at n=369 than 3-way; the cfg.levels list documents
        # the source values).
        if (identical(strat_key, "smoking_status")) {
            sub[, .carrier := as.integer(get(strat_col) %in% c("former", "current",
                                                                "previous", "Previous",
                                                                "Current"))]
        } else {
            cat(sprintf("  Skipped: cannot binarise character stratifier '%s'.\n",
                        strat_key))
            return(NULL)
        }
    } else {
        sub[, .carrier := as.integer(get(strat_col))]
    }

    n_case_carrier <- sum(sub[[STATUS_COL]] %in% c(PRE_ONSET, POST_ONSET) &
                           sub$.carrier == 1, na.rm = TRUE)
    if (n_case_carrier < 5L) {
        cat(sprintf("  Skipped: only %d case carriers (need >=5).\n",
                    n_case_carrier))
        return(NULL)
    }

    keep_eids <- as.character(sub[get(STATUS_COL) %in%
                                    c(CONTROL, PRE_ONSET, POST_ONSET), eid])
    sub <- sub[as.character(eid) %in% keep_eids]
    sub[, group4 := factor(fcase(
        get(STATUS_COL) == CONTROL & .carrier == 0, "HC_NEG",
        get(STATUS_COL) == CONTROL & .carrier == 1, "HC_POS",
        get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET) & .carrier == 0, "CASE_NEG",
        get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET) & .carrier == 1, "CASE_POS"
    ), levels = c("HC_NEG", "HC_POS", "CASE_NEG", "CASE_POS"))]
    sub[, sex_f := factor(sex)]
    sub[, ytd   := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]

    cells <- table(sub$group4)
    cat(sprintf("  group4 cells: %s\n",
                paste(names(cells), cells, sep = "=", collapse = " | ")))
    if (any(cells < 3L)) {
        cat("  Skipped: one or more group4 cells <3 — interaction underpowered.\n")
        return(NULL)
    }

    keep_idx <- as.character(sub$eid)
    expr <- prot_mat[, keep_idx]

    # Add BMI as covariate when available (per cfg phenomic_stratifiers.bmi)
    bmi_term <- if ("bmi" %in% names(sub) && any(!is.na(sub$bmi))) " + bmi" else ""
    rhs <- glue("~ 0 + group4 + ytd + age_at_sampling + sex_f + PC1 + PC2{bmi_term}")
    design <- model.matrix(as.formula(rhs), data = sub)
    fit <- lmFit(expr, design)
    contrast_int <- makeContrasts(
        interaction = (group4CASE_POS - group4HC_POS) - (group4CASE_NEG - group4HC_NEG),
        levels = design)
    fit2 <- contrasts.fit(fit, contrast_int)
    fit2 <- eBayes(fit2)
    res <- topTable_safe(fit2, coef = "interaction", n = Inf)
    res$stratifier <- strat_key
    n_sig <- sum(res$adj.P.Val < FDR_THRESH, na.rm = TRUE)
    out_csv <- file.path(OUT_DIR, glue(
        "{cfg$cohort_short}_{strat_key}_interaction.csv"))
    fwrite(res, out_csv)
    cat(sprintf("  %d FDR-sig interaction DEPs -> %s\n",
                n_sig, basename(out_csv)))

    # Optional sex × stratifier 3-way submodel: enabled when cfg lists
    # sex_x_<strat_key> in sex_handling$interaction_terms_to_test.
    if (strat_key %in% .iact_targets) {
        cat("  Running sex × smoking 3-way interaction sub-model...\n")
        sub_x <- copy(sub)
        sub_x[, is_case := as.integer(get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET))]
        design_x <- model.matrix(~ is_case * .carrier * sex_f +
                                    ytd + age_at_sampling + PC1 + PC2,
                                 data = sub_x)
        fit_x <- lmFit(expr, design_x)
        fit_x <- eBayes(fit_x)
        coef_3way <- "is_case:.carrier:sex_f1"
        if (coef_3way %in% colnames(fit_x$coefficients)) {
            res_x <- topTable_safe(fit_x, coef = coef_3way, n = Inf)
            out_csv_x <- file.path(OUT_DIR, glue(
                "{cfg$cohort_short}_{strat_key}_sex_3way.csv"))
            fwrite(res_x, out_csv_x)
            cat(sprintf("  3-way sex×case×%s: %d FDR-sig -> %s\n",
                        strat_key,
                        sum(res_x$adj.P.Val < FDR_THRESH, na.rm = TRUE),
                        basename(out_csv_x)))
        } else {
            cat("  3-way coef not in design (probably collinear); skipping.\n")
        }
    }
    res
}

results <- lapply(stratifier_keys, run_strat)
names(results) <- stratifier_keys

cat("\n01_phenomic_stratifiers.R complete.\n")
ran <- vapply(results, function(r) !is.null(r), logical(1))
cat(sprintf("  Stratifiers run:     %d / %d\n",
            sum(ran), length(stratifier_keys)))
cat(sprintf("  Stratifiers skipped: %s\n",
            paste(names(results)[!ran], collapse = ", ")))
