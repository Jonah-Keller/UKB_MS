#!/usr/bin/env Rscript
# generate_table1.R — Cohort characteristics table for manuscript Table 1
#
# Output: manuscript/submission/_build_intermediates/Table_1_Cohort.csv
#
# Columns: Healthy Controls | Pre-onset MS | Post-onset MS
# Rows:    N, age at blood draw, sex, BMI, smoking, diabetes,
#          HLA-DRB1*15:01 carrier, MS PRS, timing variables

suppressPackageStartupMessages(library(data.table))

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
CADASIL_DIR <- file.path(dirname(PROJ_DIR),
                         "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                         "data", "ukb")

QC_FILE   <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
HLA_FILE  <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE  <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
OUT_FILE  <- file.path(PROJ_DIR, "manuscript", "submission", "_build_intermediates",
                       "Table_1_Cohort.csv")

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

# ── 1. Load QC data ──────────────────────────────────────────────────────────
cat("Loading QC data...\n")
keep_cols <- c("eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
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

# ── 3. Merge HLA-DRB1*15:01 ──────────────────────────────────────────────────
if (file.exists(HLA_FILE)) {
    cat("Loading HLA-DRB1*15:01...\n")
    hla <- fread(HLA_FILE, showProgress = FALSE)
    qc <- merge(qc, hla[, .(eid, drb1_1501_carrier)], by = "eid", all.x = TRUE)
    cat(sprintf("  HLA: %d/%d\n", sum(!is.na(qc$drb1_1501_carrier)), nrow(qc)))
} else {
    qc[, drb1_1501_carrier := NA_integer_]
    cat("  HLA file not found — skipping\n")
}

# ── 4. Merge MS PRS ──────────────────────────────────────────────────────────
if (file.exists(PRS_FILE)) {
    cat("Loading MS PRS...\n")
    prs <- fread(PRS_FILE, showProgress = FALSE)
    qc <- merge(qc, prs[, .(eid, prs_score)], by = "eid", all.x = TRUE)
    cat(sprintf("  PRS: %d/%d\n", sum(!is.na(qc$prs_score)), nrow(qc)))
} else {
    qc[, prs_score := NA_real_]
    cat("  PRS file not found — skipping\n")
}

# sex: 0 = female, 1 = male (UKB field 31 convention)
groups <- list(
  hc   = qc[ms_status == "control"],
  pre  = qc[ms_status == "pre_onset"],
  post = qc[ms_status == "post_onset"]
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

  c("HLA-DRB1*15:01 carrier, n (%)",
    npct(groups$hc$drb1_1501_carrier,   groups$hc$drb1_1501_carrier   == 1),
    npct(groups$pre$drb1_1501_carrier,  groups$pre$drb1_1501_carrier  == 1),
    npct(groups$post$drb1_1501_carrier, groups$post$drb1_1501_carrier == 1)),

  c("MS polygenic risk score, median (IQR)",
    miqr(groups$hc$prs_score),
    miqr(groups$pre$prs_score),
    miqr(groups$post$prs_score)),

  c("Age at MS diagnosis, median (IQR), years",
    "\u2014",
    miqr(groups$pre$age_at_diagnosis),
    miqr(groups$post$age_at_diagnosis)),

  c("Years before MS diagnosis (pre-onset), median (IQR)",
    "\u2014",
    miqr(groups$pre$years_to_diagnosis),
    "\u2014"),

  c("Years since MS diagnosis (post-onset), median (IQR)",
    "\u2014",
    "\u2014",
    miqr(abs(groups$post$years_to_diagnosis)))
)

col_names <- c(
  "Characteristic",
  sprintf("Healthy Controls (n=%s)", format(nrow(groups$hc), big.mark = ",")),
  sprintf("Pre-onset MS (n=%s)",      format(nrow(groups$pre), big.mark = ",")),
  sprintf("Post-onset MS (n=%s)",     format(nrow(groups$post), big.mark = ","))
)

tbl <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tbl) <- col_names

write.csv(tbl, OUT_FILE, row.names = FALSE)
cat("Table 1 written to:", OUT_FILE, "\n")
print(tbl)
