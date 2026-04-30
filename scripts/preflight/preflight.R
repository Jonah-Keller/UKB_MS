#!/usr/bin/env Rscript
# ============================================================================
# preflight.R — validate disease.yaml + cohort + extracts before running
# ============================================================================
# Run this FIRST after editing configs/disease.yaml for a new disease.
# It tells you (and the agent) exactly which stages will run, which will
# skip with a clear reason, and which will crash.  No analysis is started.
#
# Usage:
#   Rscript scripts/preflight/preflight.R
#
# Exit code: 0 if every required stage will run; 1 if any required stage
# is blocked (cohort too small, missing schema field).  Stages with merely-
# missing optional data (HLA dosages, PRS profiles, C9-style mendelian
# carriers) are reported as "skip" but do not fail preflight.
#
# Output: a Markdown table summarising each stage's expected disposition,
# plus the cohort sample-size table.
# ============================================================================

suppressPackageStartupMessages({
    library(here)
    library(yaml)
    library(data.table)
})

source(here("analysis", "helpers", "disease_config.R"))

# ── 1. disease.yaml validation ───────────────────────────────────────────────
.required <- c(
    "cohort_short", "disease_long", "disease_short_caps",
    "icd_codes", "first_occurrence_field", "cohort_status_col", "status_values",
    "control_exclusion_codes",
    "top_n_volcano", "top_n_heatmap", "top_n_trajectory", "min_sig_for_fdr",
    "sex_handling"
)

cat("══ UKB proteomics replication — preflight ══\n\n")

cfg <- tryCatch(load_disease_config(), error = function(e) {
    cat("FAIL: configs/disease.yaml does not load.\n")
    cat("  Error: ", conditionMessage(e), "\n", sep = "")
    quit(save = "no", status = 1)
})

cat(sprintf("Loaded configs/disease.yaml for cohort_short = '%s'\n",
            cfg$cohort_short %||% "<missing>"))

missing <- setdiff(.required, names(cfg))
if (length(missing) > 0L) {
    cat(sprintf("FAIL: %d required disease.yaml field(s) missing: %s\n",
                length(missing), paste(missing, collapse = ", ")))
    quit(save = "no", status = 1)
}
cat(sprintf("disease.yaml: all %d required fields present\n", length(.required)))

# ── 2. cohort sample-size check ──────────────────────────────────────────────
QC_FILE <- here("data", "ukb", "olink", "processed",
                paste0(cfg$cohort_short, "_olink_qc.csv"))
cat("\n## Cohort\n\n")

if (!file.exists(QC_FILE)) {
    cat(sprintf("FAIL: %s not found.\n", basename(QC_FILE)))
    cat("  Run analysis/01_cohort/01_build_ms_cohort.py (or the disease-specific\n")
    cat("  cohort builder), then analysis/02_olink_qc/01_olink_qc.R first.\n")
    quit(save = "no", status = 1)
}

qc <- fread(QC_FILE, showProgress = FALSE)
status_col <- cfg$cohort_status_col
if (!status_col %in% names(qc)) {
    cat(sprintf("FAIL: cfg$cohort_status_col='%s' not in %s\n",
                status_col, basename(QC_FILE)))
    quit(save = "no", status = 1)
}

qc_pass <- qc[qc_outlier == FALSE]
sv <- cfg$status_values
n_pre  <- sum(qc_pass[[status_col]] == sv$pre_onset,  na.rm = TRUE)
n_post <- sum(qc_pass[[status_col]] == sv$post_onset, na.rm = TRUE)
n_ctrl <- sum(qc_pass[[status_col]] == sv$control,    na.rm = TRUE)

cat(sprintf("Total (post-QC):       %d\n", nrow(qc_pass)))
cat(sprintf("Pre-onset cases:       %d\n", n_pre))
cat(sprintf("Post-onset cases:      %d\n", n_post))
cat(sprintf("Controls:              %d\n", n_ctrl))

if (n_pre  < 30L) cat("WARN: pre-onset N < 30 — limma will be underpowered\n")
if (n_post < 20L) cat("WARN: post-onset N < 20 — post-onset contrast may not run\n")
if (n_ctrl < 500L) cat("WARN: controls N < 500 — control PC1/PC2 may be unstable\n")

# ── 3. stage-by-stage disposition ────────────────────────────────────────────
cat("\n## Stage disposition\n\n")
cat("| stage | run? | reason |\n|---|---|---|\n")

check_file <- function(path) file.exists(here(path))

stages <- list(
    list("01_cohort",                  TRUE, "always"),
    list("02_olink_qc",                TRUE, "always"),
    list("04_differential limma",      n_pre  >= 5L || n_post >= 5L,
         if (n_pre + n_post >= 5L) "ok" else "no cases pass QC"),
    list("04_differential CNS traj",   n_pre  >= 30L,
         if (n_pre >= 30L) "ok" else "pre-onset N too small for trajectory"),
    list("05_networks PPI/GO/cell",    TRUE, "depends on 04 outputs"),
    list("06_ml classifiers",          n_pre  >= 50L,
         if (n_pre >= 50L) "ok" else "pre-onset N too small for ML"),
    list("07_hla_subgroup",
         check_file(file.path("data", "ukb", "genetics", "hla_imputed.csv")) &&
             !is.null(cfg$hla_carrier_col) &&
             cfg$hla_allele %||% "NONE" != "NONE",
         if (cfg$hla_allele %||% "NONE" == "NONE")
             "skip — cfg$hla_allele set to NONE"
         else if (!check_file(file.path("data", "ukb", "genetics", "hla_imputed.csv")))
             "skip — hla_imputed.csv missing"
         else "ok"),
    list("07b_mendelian_subgroup",
         !is.null(cfg$mendelian_carrier_col) &&
             check_file(file.path("data", "ukb", "genetics", "c9orf72_expansion.csv")),
         if (is.null(cfg$mendelian_carrier_col))
             "skip — cfg$mendelian_carrier_col not set"
         else if (!check_file(file.path("data", "ukb", "genetics", "c9orf72_expansion.csv")))
             "skip — mendelian carrier file missing (see EXTRACT_TODO.md)"
         else "ok"),
    list("08_mri",                     TRUE, "depends on UKB MRI fields"),
    list("09_sex_stratified",
         identical(cfg$sex_handling$primary_use, "stratifier") &&
             !isTRUE(cfg$sex_handling$do_not_stratify_on_sex_alone),
         if (!identical(cfg$sex_handling$primary_use, "stratifier"))
             "skip — cfg$sex_handling$primary_use is 'covariate'"
         else if (isTRUE(cfg$sex_handling$do_not_stratify_on_sex_alone))
             "skip — do_not_stratify_on_sex_alone is TRUE"
         else "ok"),
    list("10_prs",
         length(cfg$prs_pgs_ids %||% character(0)) > 0L &&
             check_file(file.path("data", "ukb", "genetics",
                                  paste0(cfg$cohort_short, "_prs_scores.csv"))),
         if (length(cfg$prs_pgs_ids %||% character(0)) == 0L)
             "skip — cfg$prs_pgs_ids is empty"
         else "ok or skip — depends on prs_scores.csv"),
    list("11_cross_disease",           !is.null(cfg$comparison_cohort_short) &&
             nzchar(cfg$comparison_cohort_short),
         if (!is.null(cfg$comparison_cohort_short) &&
             nzchar(cfg$comparison_cohort_short))
             paste0("ok if results/differential/",
                    cfg$comparison_cohort_short, "_*.csv exist")
         else "skip — cfg$comparison_cohort_short not set"),
    list("13_survival",                n_pre >= 5L,
         if (n_pre >= 5L) "ok" else "pre-onset N too small for Cox PH"),
    list("14_comorbidity (PheWAS)",    n_pre >= 30L,
         if (n_pre >= 30L) "ok" else "pre-onset N too small for PheWAS"),
    list("15_endophenotype clustering", n_pre >= 50L,
         if (n_pre >= 50L) "ok" else "pre-onset N too small for k-means"),
    list("16_cluster_gwas",
         check_file(file.path("results", "endophenotype",
                              paste0(cfg$cohort_short, "_prems_cluster_assignments.csv"))),
         "phenotype build runs locally; REGENIE requires UKB-RAP"),
    list("17_phenomic_subgroup",
         length(Filter(function(s) is.list(s) &&
                                    identical(s$use_as, "stratifier"),
                       cfg$phenomic_stratifiers %||% list())) > 0L,
         if (length(cfg$phenomic_stratifiers %||% list()) == 0L)
             "skip — cfg$phenomic_stratifiers is empty"
         else "ok or skip per stratifier — depends on covariate columns")
)

required_stages <- c("01_cohort", "02_olink_qc", "04_differential limma")
hard_fail <- FALSE
for (s in stages) {
    sym <- if (isTRUE(s[[2]])) "✓" else "·"
    cat(sprintf("| %s | %s | %s |\n", s[[1]], sym, s[[3]]))
    if (s[[1]] %in% required_stages && !isTRUE(s[[2]])) hard_fail <- TRUE
}

# ── 4. summary ───────────────────────────────────────────────────────────────
cat("\n## Summary\n\n")
n_run  <- sum(vapply(stages, function(s) isTRUE(s[[2]]), logical(1)))
n_skip <- length(stages) - n_run
cat(sprintf("Stages that will run:  %d\n", n_run))
cat(sprintf("Stages that will skip: %d\n", n_skip))

if (hard_fail) {
    cat("\nFAIL: at least one required stage is blocked.\n")
    quit(save = "no", status = 1)
}

cat("\nPreflight passed.  Run the pipeline in the order documented in\n")
cat("AGENT_HANDOFF_disease_replication.md §5.\n")
