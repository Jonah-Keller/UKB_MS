#!/usr/bin/env Rscript
# 02_ms_full_phewas.R — Comprehensive ICD-10 PheWAS: MS vs HC
#
# Scans all available ICD-10 first-occurrence CSVs in the UKB diagnoses directory.
# For each code with ≥5 MS cases and ≥5 HC cases, tests:
#   present ~ ms_status + age_at_sampling + sex   (logistic regression)
# Also computes timing for MS cases (% with code before MS diagnosis).
#
# Chapters available: A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, Y, Z
#
# Output:
#   results/comorbidity/ms_full_phewas.csv      — all tested conditions + OR, p, FDR, timing

suppressPackageStartupMessages({
    library(data.table)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

ICD_DIR   <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                        "data", "ukb", "diagnoses", "icd_codes")
QC_FILE   <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
OUT_DIR   <- file.path(PROJ_DIR, "results", "comorbidity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
MIN_CASES_MS <- 5L   # minimum MS cases with condition
MIN_CASES_HC <- 5L   # minimum HC cases with condition

# ICD-10 codes that define or are synonymous with the case-defining disease
# (G35 = MS; G36/G37 = other demyelinating diseases used in MS cohort curation)
MS_EXCLUDE_CODES <- c("G35", "G36", "G37")

# ── 1. Load cohort ─────────────────────────────────────────────────────────────
cat("Loading cohort...\n")
qc     <- fread(QC_FILE, showProgress = FALSE)
cohort <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]
cohort[, is_ms := as.integer(ms_status %in% c("pre_onset", "post_onset"))]
n_ms <- sum(cohort$is_ms); n_hc <- sum(!cohort$is_ms)
cat(sprintf("  MS: %d  |  HC: %d\n", n_ms, n_hc))

# MS cases with known diagnosis age (for timing)
ms_dx <- cohort[is_ms == 1L & !is.na(age_at_diagnosis), .(eid, age_at_diagnosis)]

# ── 2. Discover all ICD code files ─────────────────────────────────────────────
cat("Scanning ICD code files...\n")
all_files <- list.files(ICD_DIR, pattern = "^[A-Z][0-9]+\\.csv$",
                         recursive = TRUE, full.names = TRUE)
# Remove MS-defining codes
all_files <- all_files[!tools::file_path_sans_ext(basename(all_files)) %in% MS_EXCLUDE_CODES]
cat(sprintf("  Found %d ICD-10 code files (after excluding %s)\n",
            length(all_files), paste(MS_EXCLUDE_CODES, collapse=", ")))

# ── 3. Helper: load one code file and run logistic regression ─────────────────
COHORT_EIDS <- cohort$eid

test_one_code <- function(fp) {
    tryCatch({
        code <- tools::file_path_sans_ext(basename(fp))
        ch   <- substr(code, 1, 1)

        # Load file: eid (col 1), present flag (col 2), age_at_event (col 3)
        ncols <- tryCatch(ncol(fread(fp, showProgress = FALSE, nrows = 1)),
                          error = function(e) 0L)
        if (ncols < 2L) return(NULL)

        raw <- fread(fp, showProgress = FALSE,
                     select = seq_len(min(ncols, 3L)), fill = TRUE)
        if (nrow(raw) == 0L) return(NULL)

        # Standardize column names
        names(raw)[1:2] <- c("eid", "present_flag")
        if (ncol(raw) >= 3) names(raw)[3] <- "age_at_event" else raw[, age_at_event := NA_real_]
        raw[, eid          := as.integer(eid)]
        raw[, present_flag := as.integer(present_flag)]
        raw[, age_at_event := as.numeric(age_at_event)]

        raw_present <- raw[present_flag == 1L, .(eid, age_at_event)]

        # n present in each group
        n_ms_with <- sum(raw_present$eid %in% cohort[is_ms == 1L, eid])
        n_hc_with <- sum(raw_present$eid %in% cohort[is_ms == 0L, eid])
        if (n_ms_with < MIN_CASES_MS || n_hc_with < MIN_CASES_HC) return(NULL)

        # Build analysis dataset
        dt <- copy(cohort[, .(eid, is_ms, age_at_sampling, sex, age_at_diagnosis)])
        dt[, present := as.integer(eid %in% raw_present$eid)]

        # Logistic regression — Wald CI to avoid profiling slowdown
        fit <- suppressWarnings(
            tryCatch(
                glm(present ~ is_ms + age_at_sampling + sex, data = dt, family = binomial),
                error = function(e) NULL
            )
        )
        if (is.null(fit)) return(NULL)
        cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
        if (is.null(cf) || !"is_ms" %in% rownames(cf)) return(NULL)

        beta <- cf["is_ms", "Estimate"]
        se   <- cf["is_ms", "Std. Error"]
        pval <- cf["is_ms", "Pr(>|z|)"]

        or    <- exp(beta)
        or_lo <- exp(beta - 1.96 * se)
        or_hi <- exp(beta + 1.96 * se)

        # Timing: % of MS cases with this code recorded > 1yr before MS diagnosis
        timing_pct_pre <- NA_real_
        n_with_timing  <- 0L
        if (nrow(ms_dx) > 0L && any(!is.na(raw_present$age_at_event))) {
            ms_p <- raw_present[eid %in% ms_dx$eid]
            if (nrow(ms_p) > 0L) {
                tdt <- merge(ms_p, ms_dx, by = "eid", all.x = TRUE)
                tdt[, timing := age_at_event - age_at_diagnosis]
                n_with_timing  <- nrow(tdt)
                timing_pct_pre <- 100 * mean(tdt$timing < -1, na.rm = TRUE)
            }
        }

        data.table(
            icd10         = code,
            chapter       = ch,
            OR            = or,
            OR_lo         = or_lo,
            OR_hi         = or_hi,
            beta          = beta,
            se            = se,
            pval          = pval,
            n_ms_with     = n_ms_with,
            n_hc_with     = n_hc_with,
            pct_ms        = 100 * n_ms_with / n_ms,
            pct_hc        = 100 * n_hc_with / n_hc,
            n_timing_ms   = n_with_timing,
            pct_pre_onset = timing_pct_pre
        )
    }, error = function(e) NULL)
}

# ── 4. Run PheWAS across all codes ─────────────────────────────────────────────
cat("Running PheWAS...\n")
N_FILES <- length(all_files)
results_list <- lapply(seq_len(N_FILES), function(i) {
    res <- test_one_code(all_files[[i]])
    if (i %% 200L == 0L)
        cat(sprintf("  %d / %d files\n", i, N_FILES))
    res
})

phewas_full <- rbindlist(Filter(Negate(is.null), results_list))
cat(sprintf("\nTotal conditions tested: %d\n", nrow(phewas_full)))

# ── 5. Multiple testing correction ─────────────────────────────────────────────
phewas_full[, fdr := p.adjust(pval, method = "BH")]
phewas_full <- phewas_full[order(pval)]

n_fdr <- sum(phewas_full$fdr < 0.05, na.rm = TRUE)
cat(sprintf("FDR < 0.05: %d conditions\n", n_fdr))
cat(sprintf("Nominal p < 0.05: %d conditions\n", sum(phewas_full$pval < 0.05, na.rm = TRUE)))

# Top hits
cat("\nTop 20 associations:\n")
print(phewas_full[1:min(20L, .N), .(icd10, chapter, OR = round(OR, 2),
                                     pct_ms = round(pct_ms, 1),
                                     pval = signif(pval, 3),
                                     fdr = round(fdr, 3),
                                     pct_pre_onset = round(pct_pre_onset, 0))])

# ── 6. Save ────────────────────────────────────────────────────────────────────
out_path <- file.path(OUT_DIR, "ms_full_phewas.csv")
fwrite(phewas_full, out_path)
cat(sprintf("\nSaved: %s  (%d rows)\n", basename(out_path), nrow(phewas_full)))
cat("02_ms_full_phewas.R complete.\n")
