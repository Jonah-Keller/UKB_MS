#!/usr/bin/env Rscript
# 02a_prems_phewas.R — Pre-diagnostic PheWAS: ICD enrichment before disease onset
#
# Outcome is defined per participant:
#   Cases: present = 1 iff first-occurrence age < (age_at_diagnosis - 1 yr)
#   HC:    present = 1 iff first-occurrence age < age_at_sampling
#
# This tests "was this code enriched in cases before their diagnosis compared
# to same-aged HC at the same time point?"
# No post-hoc timing filter needed — timing is baked into the outcome.
#
# Model: present ~ is_case + age_at_sampling + sex  (logistic, Wald CI)
# Correction: BH FDR across all tested codes.
# Both directions retained (OR>1 enriched, OR<1 depleted before Dx).
#
# Disease-, status-, and exclusion-code constants come from configs/disease.yaml.
#
# Output: results/comorbidity/{cohort_short}_prems_phewas.csv

suppressPackageStartupMessages({
    library(data.table)
    library(here)
    library(glue)
})

source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

STATUS_COL  <- cfg$cohort_status_col
PRE_ONSET   <- cfg$status_values$pre_onset
POST_ONSET  <- cfg$status_values$post_onset
COHORT      <- cfg$cohort_short
DISEASE_LBL <- cfg$disease_short_caps
EXCL_CODES  <- cfg$control_exclusion_codes$demyelinating

ICD_DIR <- file.path(dirname(here::here()),
                     "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                     "data", "ukb", "diagnoses", "icd_codes")
QC_FILE <- here::here("data", "ukb", "olink", "processed",
                      glue::glue("{COHORT}_olink_qc.csv"))
OUT_DIR <- here::here("results", "comorbidity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_N_CASE <- 5L
MIN_N_HC   <- 5L

# ── 1. Build analysis cohort ──────────────────────────────────────────────────
cat("Loading cohort...\n")
qc <- fread(QC_FILE, showProgress = FALSE)
qc <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]
qc[, is_case := as.integer(get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET))]

# Cases with known diagnosis age — cutoff is 1yr before Dx
cases <- qc[is_case == 1L & !is.na(age_at_diagnosis),
            .(eid, is_case, age_at_sampling, sex,
              cutoff_age = age_at_diagnosis - 1)]
# HC — cutoff is age at sampling (all available record)
hc <- qc[is_case == 0L,
         .(eid, is_case, age_at_sampling, sex, cutoff_age = age_at_sampling)]

cohort <- rbind(cases, hc)
n_case <- nrow(cases); n_hc <- nrow(hc)
cat(sprintf("  %s (with Dx age): %d  |  HC: %d\n", DISEASE_LBL, n_case, n_hc))

# ── 2. Discover code files ────────────────────────────────────────────────────
cat("Scanning ICD code files...\n")
all_files <- list.files(ICD_DIR, pattern = "[A-Z][0-9]+[.]csv",
                        recursive = TRUE, full.names = TRUE)
all_files <- all_files[!tools::file_path_sans_ext(basename(all_files)) %in% EXCL_CODES]
cat(sprintf("  %d code files (excluding %s)\n", length(all_files),
            paste(EXCL_CODES, collapse = ", ")))

# ── 3. Test one code ──────────────────────────────────────────────────────────
test_one_code <- function(fp) {
    tryCatch({
        code <- tools::file_path_sans_ext(basename(fp))
        ch   <- substr(code, 1, 1)
        blk  <- substr(code, 1, 2)

        raw <- fread(fp, showProgress = FALSE, select = 1:3, fill = TRUE)
        if (nrow(raw) == 0L) return(NULL)
        names(raw) <- c("eid", "ever_present", "age_at_event")
        raw[, eid          := as.integer(eid)]
        raw[, ever_present := as.integer(ever_present)]
        raw[, age_at_event := as.numeric(age_at_event)]

        dt <- merge(cohort, raw[, .(eid, ever_present, age_at_event)],
                    by = "eid", all.x = TRUE)
        dt[is.na(ever_present), ever_present := 0L]
        dt[is.na(age_at_event), age_at_event  := Inf]

        # Present = had code before participant-specific cutoff
        dt[, present := as.integer(ever_present == 1L & age_at_event < cutoff_age)]

        n_case_pre <- sum(dt$present[dt$is_case == 1L])
        n_hc_pre   <- sum(dt$present[dt$is_case == 0L])
        if (n_case_pre < MIN_N_CASE || n_hc_pre < MIN_N_HC) return(NULL)

        fit <- suppressWarnings(
            tryCatch(
                glm(present ~ is_case + age_at_sampling + sex,
                    data = dt, family = binomial),
                error = function(e) NULL
            )
        )
        if (is.null(fit)) return(NULL)
        cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
        if (is.null(cf) || !"is_case" %in% rownames(cf)) return(NULL)

        beta <- cf["is_case", "Estimate"]
        se   <- cf["is_case", "Std. Error"]
        pval <- cf["is_case", "Pr(>|z|)"]

        data.table(
            icd10      = code,
            block      = blk,
            chapter    = ch,
            OR         = exp(beta),
            OR_lo      = exp(beta - 1.96 * se),
            OR_hi      = exp(beta + 1.96 * se),
            beta       = beta,
            se         = se,
            pval       = pval,
            n_case_pre = n_case_pre,
            n_hc_pre   = n_hc_pre,
            pct_case   = 100 * n_case_pre / n_case,
            pct_hc     = 100 * n_hc_pre / n_hc
        )
    }, error = function(e) NULL)
}

# ── 4. Run PheWAS ─────────────────────────────────────────────────────────────
cat(sprintf("Running pre-%s PheWAS...\n", DISEASE_LBL))
N_FILES <- length(all_files)
res_list <- lapply(seq_len(N_FILES), function(i) {
    r <- test_one_code(all_files[[i]])
    if (i %% 200L == 0L) cat(sprintf("  %d / %d\n", i, N_FILES))
    r
})
phewas <- rbindlist(Filter(Negate(is.null), res_list))
cat(sprintf("\nConditions tested: %d\n", nrow(phewas)))

# ── 5. FDR correction ─────────────────────────────────────────────────────────
phewas[, fdr := p.adjust(pval, method = "BH")]
phewas <- phewas[order(pval)]

n_sig <- sum(phewas$fdr < 0.05, na.rm = TRUE)
cat(sprintf("FDR < 0.05: %d  (OR>1: %d, OR<1: %d)\n",
            n_sig,
            sum(phewas$fdr < 0.05 & phewas$OR > 1, na.rm = TRUE),
            sum(phewas$fdr < 0.05 & phewas$OR < 1, na.rm = TRUE)))

cat("\nTop 25:\n")
print(phewas[seq_len(min(25L, .N)),
             .(icd10, block, OR = round(OR, 2),
               pct_case = round(pct_case, 1), pct_hc = round(pct_hc, 1),
               pval = signif(pval, 3), fdr = round(fdr, 3))])

# ── 6. Save ───────────────────────────────────────────────────────────────────
out_path <- file.path(OUT_DIR, glue::glue("{COHORT}_prems_phewas.csv"))
fwrite(phewas, out_path)
cat(sprintf("\nSaved: %s  (%d rows)\n", basename(out_path), nrow(phewas)))
cat("02a_prems_phewas.R complete.\n")
