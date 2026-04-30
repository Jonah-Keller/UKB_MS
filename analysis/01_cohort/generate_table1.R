#!/usr/bin/env Rscript
# generate_table1.R — Cohort characteristics table for manuscript Table 1
# (config-driven; disease/HLA/PRS labels read from configs/disease.yaml)
#
# Output: manuscript/submission/_build_intermediates/Table_1_Cohort.csv
#
# Columns: Healthy Controls | Pre-onset {disease} | Post-onset {disease}
# Rows:    N, age at blood draw, sex, BMI, smoking, diabetes,
#          HLA carrier, disease PRS, timing variables

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(glue)
})

source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

STATUS_COL     <- cfg$cohort_status_col
STATUS_CONTROL <- cfg$status_values$control
STATUS_PRE     <- cfg$status_values$pre_onset
STATUS_POST    <- cfg$status_values$post_onset

QC_FILE  <- here::here("data", "ukb", "olink", "processed",
                       glue::glue("{cfg$cohort_short}_olink_qc.csv"))
HLA_FILE <- here::here("data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE <- here::here("data", "ukb", "genetics",
                       glue::glue("{cfg$cohort_short}_prs_scores.csv"))
OUT_FILE <- here::here("manuscript", "submission", "_build_intermediates",
                       "Table_1_Cohort.csv")

# Sibling repo (CADASIL) provides covariate tables; resolve relative to project root.
CADASIL_DIR <- normalizePath(file.path(here::here(), "..",
                                       "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                                       "data", "ukb"),
                             mustWork = FALSE)

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

# ── 1. Load QC data ──────────────────────────────────────────────────────────
cat("Loading QC data...\n")
keep_cols <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance", "qc_outlier")
qc <- fread(QC_FILE, select = intersect(keep_cols, names(fread(QC_FILE, nrows = 0))))
qc <- qc[qc_outlier == FALSE]

# ── 2. Merge covariates (BMI, smoking, diabetes) ─────────────────────────────
cat("Loading comorbidity covariates...\n")
misc <- file.path(CADASIL_DIR, "misc")

# BMI (field p21001_i0)
bmi <- fread(file.path(misc, "ukb_bmi_age.csv"), showProgress = FALSE)
setnames(bmi, c("participant.eid", "participant.p21001_i0"), c("eid", "bmi"),
         skip_absent = TRUE)
bmi <- bmi[, .(eid, bmi)]

# Smoking (field p20116_i0: 0=never, 1=previous, 2=current)
smoke <- fread(file.path(misc, "ukb_smoking_status.csv"), showProgress = FALSE)
setnames(smoke, c("participant.eid", "participant.p20116_i0"),
         c("eid", "smoking_code"), skip_absent = TRUE)
smoke[, ever_smoker := as.integer(smoking_code %in% c(1L, 2L))]
smoke <- smoke[smoking_code != -3L, .(eid, ever_smoker)]

# Diabetes (field p2443_i0)
wmh <- fread(file.path(misc, "wmh_htn_dm.csv"), showProgress = FALSE)
wmh[, diabetes := as.integer(p2443_i0 == "Yes")]
wmh <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes)]

qc <- merge(qc, bmi,   by = "eid", all.x = TRUE)
qc <- merge(qc, smoke, by = "eid", all.x = TRUE)
qc <- merge(qc, wmh,   by = "eid", all.x = TRUE)
cat(sprintf("  BMI: %d/%d  Smoking: %d/%d  Diabetes: %d/%d\n",
            sum(!is.na(qc$bmi)), nrow(qc),
            sum(!is.na(qc$ever_smoker)), nrow(qc),
            sum(!is.na(qc$diabetes)), nrow(qc)))

# ── 3. Merge HLA carrier ─────────────────────────────────────────────────────
hla_col <- cfg$hla_carrier_col
if (file.exists(HLA_FILE)) {
    cat(sprintf("Loading HLA-%s...\n", cfg$hla_allele))
    hla <- fread(HLA_FILE, showProgress = FALSE)
    qc <- merge(qc, hla[, .SD, .SDcols = c("eid", hla_col)], by = "eid", all.x = TRUE)
    cat(sprintf("  HLA: %d/%d\n", sum(!is.na(qc[[hla_col]])), nrow(qc)))
} else {
    qc[, (hla_col) := NA_integer_]
    cat("  HLA file not found — skipping\n")
}

# ── 4. Merge disease PRS ─────────────────────────────────────────────────────
prs_col <- cfg$prs_combined_col
if (file.exists(PRS_FILE)) {
    cat(sprintf("Loading %s...\n", cfg$prs_label))
    prs <- fread(PRS_FILE, showProgress = FALSE)
    qc <- merge(qc, prs[, .SD, .SDcols = c("eid", prs_col)], by = "eid", all.x = TRUE)
    cat(sprintf("  PRS: %d/%d\n", sum(!is.na(qc[[prs_col]])), nrow(qc)))
} else {
    qc[, (prs_col) := NA_real_]
    cat("  PRS file not found — skipping\n")
}

# sex: 0 = female, 1 = male (UKB field 31 convention)
groups <- list(
  hc   = qc[get(STATUS_COL) == STATUS_CONTROL],
  pre  = qc[get(STATUS_COL) == STATUS_PRE],
  post = qc[get(STATUS_COL) == STATUS_POST]
)

# ── Formatting helpers ────────────────────────────────────────────────────────
miqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  q <- quantile(x, c(0.25, 0.5, 0.75))
  sprintf("%.1f (%.1f\u2013%.1f)", q[2], q[1], q[3])
}

npct <- function(x, cond) {
  n <- sum(cond, na.rm = TRUE)
  total <- sum(!is.na(x))
  if (total == 0) return("\u2014")
  p <- round(100 * n / total, 1)
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), p)
}

# ── Build table rows ─────────────────────────────────────────────────────────
rows <- list(
  c("N",
    format(nrow(groups$hc),  big.mark = ","),
    format(nrow(groups$pre), big.mark = ","),
    format(nrow(groups$post),big.mark = ",")),

  c("Age at blood draw, median (IQR), years",
    miqr(groups$hc$age_at_sampling),
    miqr(groups$pre$age_at_sampling),
    miqr(groups$post$age_at_sampling)),

  c("Female sex, n (%)",
    npct(groups$hc$sex,   groups$hc$sex   == 0),
    npct(groups$pre$sex,  groups$pre$sex  == 0),
    npct(groups$post$sex, groups$post$sex == 0)),

  c("BMI, median (IQR), kg/m\u00B2",
    miqr(groups$hc$bmi),
    miqr(groups$pre$bmi),
    miqr(groups$post$bmi)),

  c("Ever smoker, n (%)",
    npct(groups$hc$ever_smoker,   groups$hc$ever_smoker   == 1),
    npct(groups$pre$ever_smoker,  groups$pre$ever_smoker  == 1),
    npct(groups$post$ever_smoker, groups$post$ever_smoker == 1)),

  c("Diabetes, n (%)",
    npct(groups$hc$diabetes,   groups$hc$diabetes   == 1),
    npct(groups$pre$diabetes,  groups$pre$diabetes  == 1),
    npct(groups$post$diabetes, groups$post$diabetes == 1)),

  c(glue::glue("HLA-{cfg$hla_allele} carrier, n (%)"),
    npct(groups$hc[[hla_col]],   groups$hc[[hla_col]]   == 1),
    npct(groups$pre[[hla_col]],  groups$pre[[hla_col]]  == 1),
    npct(groups$post[[hla_col]], groups$post[[hla_col]] == 1)),

  c(glue::glue("{cfg$disease_short_caps} polygenic risk score, median (IQR)"),
    miqr(groups$hc[[prs_col]]),
    miqr(groups$pre[[prs_col]]),
    miqr(groups$post[[prs_col]])),

  c(glue::glue("Age at {cfg$disease_short_caps} diagnosis, median (IQR), years"),
    "\u2014",
    miqr(groups$pre$age_at_diagnosis),
    miqr(groups$post$age_at_diagnosis)),

  c(glue::glue("Years before {cfg$disease_short_caps} diagnosis (pre-onset), median (IQR)"),
    "\u2014",
    miqr(groups$pre$years_to_diagnosis),
    "\u2014"),

  c(glue::glue("Years since {cfg$disease_short_caps} diagnosis (post-onset), median (IQR)"),
    "\u2014",
    "\u2014",
    miqr(abs(groups$post$years_to_diagnosis)))
)

col_names <- c(
  "Characteristic",
  sprintf("Healthy Controls (n=%s)", format(nrow(groups$hc), big.mark = ",")),
  sprintf("Pre-onset %s (n=%s)",  cfg$disease_short_caps, format(nrow(groups$pre), big.mark = ",")),
  sprintf("Post-onset %s (n=%s)", cfg$disease_short_caps, format(nrow(groups$post), big.mark = ","))
)

tbl <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tbl) <- col_names

write.csv(tbl, OUT_FILE, row.names = FALSE)
cat("Table 1 written to:", OUT_FILE, "\n")
print(tbl)
