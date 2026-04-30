#!/usr/bin/env Rscript
# 01_hla_ms_subgroup.R
# HLA risk-allele subgroup analysis — UKB-unique extension
#
# Approach:
#   - Stratify disease DEPs by HLA risk-allele carrier status
#   - Test: are DEP signatures stronger in carriers?
#   - Interaction model: DEP ~ status * hla_carrier + covariates
#
# Disease-, HLA-, and status-column constants come from configs/disease.yaml.
#
# DATA REQUIREMENTS:
#   UKB field 22182: Imputed HLA alleles (HIBAG imputation)
#   Expected: file with column cfg$hla_carrier_col (0/1) or cfg$hla_dosage_col
#   OR raw HIBAG p22182_a* columns (parsed using cfg$hla_search_patterns)
#
# Output: results/hla/
#   <cohort>_hla_subgroup_volcano.pdf
#   <cohort>_hla_subgroup_results.csv

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(limma)
    library(ggplot2)
    library(ggrepel)
    library(MatchIt)
    library(patchwork)
    library(glue)
    library(here)
})

source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

HLA_CARRIER_COL <- cfg$hla_carrier_col
HLA_DOSAGE_COL  <- cfg$hla_dosage_col
HLA_ALLELE      <- cfg$hla_allele
STATUS_COL      <- cfg$cohort_status_col
PRE_ONSET       <- cfg$status_values$pre_onset
POST_ONSET      <- cfg$status_values$post_onset
SEARCH_PATTERNS <- cfg$hla_search_patterns

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed",
                         glue("{cfg$cohort_short}_olink_qc.csv"))
DIFF_FILE   <- file.path(PROJ_DIR, "results", "differential",
                         glue("{cfg$cohort_short}_combined_vs_hc.csv"))
OUT_DIR     <- file.path(PROJ_DIR, "results", "hla")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Helper: does any pattern in `patterns` match string `x`?
match_any_pattern <- function(x, patterns) {
    if (is.null(x)) return(FALSE)
    Reduce(`|`, lapply(patterns, function(p) grepl(p, x, fixed = TRUE)))
}

# ── Data file candidates ──────────────────────────────────────────────────────
hla_candidates <- c(
    file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv"),
    file.path(CADASIL_DIR, "genotype", "ukb_hla_imputed.csv"),
    file.path(CADASIL_DIR, "genotype", "ukb_hla_field22182.csv"),
    file.path(CADASIL_DIR, "genotype", "hla_drb1.csv")
)

hla_file <- Filter(file.exists, hla_candidates)[1]

if (is.null(hla_file) || is.na(hla_file)) {
    message(strrep("=", 65))
    message(sprintf("DATA REQUIRED: HLA imputed alleles for %s not found.", HLA_ALLELE))
    message("Expected one of:")
    for (f in hla_candidates) message("  ", f)
    message("")
    message("To extract from UKB:")
    message("  Field 22182 = Imputed HLA alleles (HIBAG)")
    message("  https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22182")
    message("")
    message(sprintf("Expected column: %s (0/1) or %s (dosage)",
                    HLA_CARRIER_COL, HLA_DOSAGE_COL))
    message(strrep("=", 65))
    quit(save = "no", status = 0)
}

cat("Loading HLA data from:", hla_file, "\n")
hla <- fread(hla_file, showProgress = FALSE)

# ── Parse risk-allele carrier status ─────────────────────────────────────────
# Handle three formats:
#   1. Pre-computed carrier column (cfg$hla_carrier_col, already 0/1)
#   2. Dosage column (cfg$hla_dosage_col), threshold 0.5
#   3. Allele strings (HLA_DRB1_1 / HLA_DRB1_2 or raw p22182_a*) matched
#      against cfg$hla_search_patterns
if (HLA_CARRIER_COL %in% names(hla)) {
    cat(sprintf("  Using pre-computed %s column\n", HLA_CARRIER_COL))
    hla[, (HLA_CARRIER_COL) := as.integer(get(HLA_CARRIER_COL))]
} else if (HLA_DOSAGE_COL %in% names(hla)) {
    hla[, (HLA_CARRIER_COL) := as.integer(get(HLA_DOSAGE_COL) >= 0.5)]
} else if ("HLA_DRB1_1" %in% names(hla)) {
    hla[, (HLA_CARRIER_COL) := as.integer(
        match_any_pattern(HLA_DRB1_1, SEARCH_PATTERNS) |
        match_any_pattern(HLA_DRB1_2, SEARCH_PATTERNS)
    )]
} else if (any(grepl("^p22182_a", names(hla)))) {
    a_cols <- grep("^p22182_a", names(hla), value = TRUE)
    hla[, (HLA_CARRIER_COL) := as.integer(
        apply(as.data.frame(hla[, ..a_cols]), 1, function(row) {
            any(vapply(row, function(v)
                match_any_pattern(as.character(v), SEARCH_PATTERNS),
                logical(1)), na.rm = TRUE)
        })
    )]
} else {
    stop(sprintf(
        "HLA file format not recognised. Expected columns: %s, %s, HLA_DRB1_1/2, or p22182_a*",
        HLA_CARRIER_COL, HLA_DOSAGE_COL))
}

carrier_vec <- hla[[HLA_CARRIER_COL]]
cat(sprintf("  HLA-%s carriers: %d / %d (%.1f%%)\n",
            HLA_ALLELE,
            sum(carrier_vec, na.rm = TRUE),
            sum(!is.na(carrier_vec)),
            100 * mean(carrier_vec, na.rm = TRUE)))

# ── Load QC data ─────────────────────────────────────────────────────────────
cat("Loading QC data...\n")
qc <- fread(QC_FILE, showProgress = FALSE)
qc <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
         !is.na(UMAP1) & !is.na(UMAP2)]

dt <- merge(qc, hla[, .SD, .SDcols = c("eid", HLA_CARRIER_COL)],
            by = "eid", all.x = TRUE)
dt <- dt[!is.na(get(HLA_CARRIER_COL))]

cases <- dt[get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET)]
carrier_n    <- sum(cases[[HLA_CARRIER_COL]] == 1)
noncarrier_n <- sum(cases[[HLA_CARRIER_COL]] == 0)
cat(sprintf("  %s cases with HLA: %d  (carrier=%d, non-carrier=%d)\n",
            cfg$disease_short_caps, nrow(cases), carrier_n, noncarrier_n))

# ── Load DEPs ────────────────────────────────────────────────────────────────
deps        <- fread(DIFF_FILE)
dep_prots   <- tolower(deps[adj.P.Val < 0.05, protein])
avail_prots <- intersect(dep_prots, names(dt))

# ── Carrier vs non-carrier comparison within cases ───────────────────────────
cat("\nCarrier vs non-carrier comparison (cases)...\n")
if (carrier_n < 20 || noncarrier_n < 20) {
    cat("  WARNING: Insufficient carriers or non-carriers for comparison\n")
} else {
    cmp_list <- lapply(avail_prots, function(prot) {
        sub <- cases[!is.na(get(prot))]
        carrier <- sub[[HLA_CARRIER_COL]]
        wt <- wilcox.test(sub[[prot]] ~ carrier, exact = FALSE)
        data.table(
            protein         = prot,
            carrier_mean    = mean(sub[carrier == 1, get(prot)], na.rm = TRUE),
            noncarrier_mean = mean(sub[carrier == 0, get(prot)], na.rm = TRUE),
            pval            = wt$p.value
        )
    })
    cmp_dt <- rbindlist(cmp_list)
    cmp_dt[, fdr   := p.adjust(pval, method = "BH")]
    cmp_dt[, delta := carrier_mean - noncarrier_mean]
    cmp_dt <- cmp_dt[order(pval)]
    cmp_dt[, logFC := delta]  # approximate
    cmp_dt[, sig_label := fcase(fdr < 0.05, "*", pval < 0.05, "+", default = "")]
    fwrite(cmp_dt, file.path(OUT_DIR,
                             glue("{cfg$cohort_short}_hla_subgroup_results.csv")))
    cat(sprintf("  FDR<0.05: %d DEPs differ by HLA carrier status\n",
                sum(cmp_dt$fdr < 0.05)))

    # Volcano
    cmp_dt[, direction := fifelse(delta > 0, "Higher in carriers", "Lower in carriers")]
    cmp_dt[, highlight := fdr < 0.05 |
                          (pval < 0.05 & abs(delta) > quantile(abs(delta), 0.75))]
    cmp_dt[, label_prot := fifelse(highlight, toupper(protein), NA_character_)]

    DIR_COLS <- c("Higher in carriers" = "#CC0066", "Lower in carriers" = "#56B4E9")

    p_vol <- ggplot(cmp_dt, aes(x = delta, y = -log10(pval), colour = direction)) +
        geom_point(aes(alpha = highlight), size = 1.5) +
        scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3), guide = "none") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                   colour = "grey50", linewidth = 0.3) +
        geom_text_repel(aes(label = label_prot), size = 2.2, na.rm = TRUE,
                        max.overlaps = 15) +
        scale_colour_manual(values = DIR_COLS, name = NULL) +
        labs(x = glue("Δ NPX (HLA-{HLA_ALLELE} carriers - non-carriers)"),
             y = expression(-log[10](p)),
             title = glue("{cfg$disease_short_caps} DEPs: HLA-{HLA_ALLELE} carrier vs non-carrier"),
             subtitle = sprintf("Cases only: carrier n=%d, non-carrier n=%d",
                                carrier_n, noncarrier_n)) +
        theme_ukb()

    ggsave(file.path(OUT_DIR,
                     glue("{cfg$cohort_short}_hla_subgroup_volcano.pdf")),
           p_vol, width = 6, height = 5, device = cairo_pdf)
    cat(sprintf("  Saved: %s_hla_subgroup_volcano.pdf\n", cfg$cohort_short))
}

cat("\n01_hla_ms_subgroup.R complete.\n")
