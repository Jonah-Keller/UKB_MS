#!/usr/bin/env Rscript
# 02a_ms_prems_phewas.R — Pre-diagnostic PheWAS: ICD enrichment before MS onset
#
# Outcome is defined per participant:
#   MS cases:  present = 1 iff first-occurrence age < (age_at_diagnosis - 1 yr)
#   HC:        present = 1 iff first-occurrence age < age_at_sampling
#
# This tests "was this code enriched in MS patients before their diagnosis
# compared to same-aged HC at the same time point?"
# No post-hoc timing filter needed — timing is baked into the outcome.
#
# Model: present ~ ms_status + age_at_sampling + sex  (logistic, Wald CI)
# Correction: BH FDR across all tested codes.
# Both directions retained (OR>1 enriched, OR<1 depleted before Dx).
#
# Output: results/comorbidity/ms_prems_phewas.csv

suppressPackageStartupMessages(library(data.table))

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

ICD_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                     "data", "ukb", "diagnoses", "icd_codes")
QC_FILE <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
OUT_DIR <- file.path(PROJ_DIR, "results", "comorbidity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_N_MS <- 5L
MIN_N_HC <- 5L
MS_EXCL  <- c("G35", "G36", "G37")

# ── 1. Build analysis cohort ──────────────────────────────────────────────────
cat("Loading cohort...\n")
qc <- fread(QC_FILE, showProgress = FALSE)
qc <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]
qc[, is_ms := as.integer(ms_status %in% c("pre_onset", "post_onset"))]

# MS cases with known diagnosis age — cutoff is 1yr before Dx
ms <- qc[is_ms == 1L & !is.na(age_at_diagnosis),
         .(eid, is_ms, age_at_sampling, sex, cutoff_age = age_at_diagnosis - 1)]
# HC — cutoff is age at sampling (all available record)
hc <- qc[is_ms == 0L,
         .(eid, is_ms, age_at_sampling, sex, cutoff_age = age_at_sampling)]

cohort <- rbind(ms, hc)
n_ms <- nrow(ms); n_hc <- nrow(hc)
cat(sprintf("  MS (with Dx age): %d  |  HC: %d\n", n_ms, n_hc))

# ── 2. Discover code files ────────────────────────────────────────────────────
cat("Scanning ICD code files...\n")
all_files <- list.files(ICD_DIR, pattern = "[A-Z][0-9]+[.]csv",
                        recursive = TRUE, full.names = TRUE)
all_files <- all_files[!tools::file_path_sans_ext(basename(all_files)) %in% MS_EXCL]
cat(sprintf("  %d code files (excluding %s)\n", length(all_files),
            paste(MS_EXCL, collapse = ", ")))

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

        n_ms_pre <- sum(dt$present[dt$is_ms == 1L])
        n_hc_pre <- sum(dt$present[dt$is_ms == 0L])
        if (n_ms_pre < MIN_N_MS || n_hc_pre < MIN_N_HC) return(NULL)

        fit <- suppressWarnings(
            tryCatch(
                glm(present ~ is_ms + age_at_sampling + sex,
                    data = dt, family = binomial),
                error = function(e) NULL
            )
        )
        if (is.null(fit)) return(NULL)
        cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
        if (is.null(cf) || !"is_ms" %in% rownames(cf)) return(NULL)

        beta <- cf["is_ms", "Estimate"]
        se   <- cf["is_ms", "Std. Error"]
        pval <- cf["is_ms", "Pr(>|z|)"]

        data.table(
            icd10    = code,
            block    = blk,
            chapter  = ch,
            OR       = exp(beta),
            OR_lo    = exp(beta - 1.96 * se),
            OR_hi    = exp(beta + 1.96 * se),
            beta     = beta,
            se       = se,
            pval     = pval,
            n_ms_pre = n_ms_pre,
            n_hc_pre = n_hc_pre,
            pct_ms   = 100 * n_ms_pre / n_ms,
            pct_hc   = 100 * n_hc_pre / n_hc
        )
    }, error = function(e) NULL)
}

# ── 4. Run PheWAS ─────────────────────────────────────────────────────────────
cat("Running pre-MS PheWAS...\n")
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
               pct_ms = round(pct_ms, 1), pct_hc = round(pct_hc, 1),
               pval = signif(pval, 3), fdr = round(fdr, 3))])

# ── 6. Save ───────────────────────────────────────────────────────────────────
out_path <- file.path(OUT_DIR, "ms_prems_phewas.csv")
fwrite(phewas, out_path)
cat(sprintf("\nSaved: %s  (%d rows)\n", basename(out_path), nrow(phewas)))
cat("02a_ms_prems_phewas.R complete.\n")
