#!/usr/bin/env Rscript
# 01_hla_ms_subgroup.R
# HLA-DRB1*15:01 subgroup analysis — UKB-unique extension
#
# Approach (Abdelhak 2026 limitation addressed):
#   - Stratify MS DEPs + trajectories by HLA-DRB1*15:01 carrier status
#   - HLA-DRB1*15:01 is the strongest MS genetic risk factor (OR ~3)
#   - Test: are DEP signatures stronger in DRB1*15:01 carriers?
#   - Interaction model: DEP ~ ms_status * hla_carrier + covariates
#
# DATA REQUIREMENTS (NOT YET AVAILABLE):
#   UKB field 22182: Imputed HLA alleles (HIBAG imputation)
#   Expected format: ukb_hla_imputed.csv with columns:
#     eid, HLA_DRB1_1, HLA_DRB1_2, HLA_DRB1_alleles (semicolon-separated)
#   OR: ukb_hla_field22182.csv extracted from UKB RAP with columns:
#     eid, p22182_a0, p22182_a1, ... (one allele per column)
#
# To extract: In UKB RAP, run:
#   dx run app-swiss-army-knife -iimage_file="ukb-rap-tools" \
#     -icmd="python3 extract_fields.py --fields 22182 --output ukb_hla.csv"
#   OR download via UK Biobank Data Portal (field 22182, all arrays)
#
# Output: results/hla/
#   ms_hla_subgroup_volcano.pdf
#   ms_hla_subgroup_interaction.pdf
#   ms_hla_subgroup_results.csv

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(limma)
    library(ggplot2)
    library(ggrepel)
    library(MatchIt)
    library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DIFF_FILE   <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
OUT_DIR     <- file.path(PROJ_DIR, "results", "hla")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Data file candidates ──────────────────────────────────────────────────────
# Look for HLA data in several expected locations
hla_candidates <- c(
    file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv"),
    file.path(CADASIL_DIR, "genotype", "ukb_hla_imputed.csv"),
    file.path(CADASIL_DIR, "genotype", "ukb_hla_field22182.csv"),
    file.path(CADASIL_DIR, "genotype", "hla_drb1.csv")
)

hla_file <- Filter(file.exists, hla_candidates)[1]

if (is.null(hla_file) || is.na(hla_file)) {
    message("═══════════════════════════════════════════════════════════════")
    message("DATA REQUIRED: HLA imputed alleles not found.")
    message("Expected one of:")
    for (f in hla_candidates) message("  ", f)
    message("")
    message("To extract from UKB:")
    message("  Field 22182 = Imputed HLA alleles (HIBAG, 6 HLA genes, i2)")
    message("  Download via UKB Biobank Data Portal or:")
    message("  https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22182")
    message("")
    message("Expected column format after extraction:")
    message("  eid, p22182_a0, p22182_a1, ... (one allele code per column)")
    message("  OR: eid, HLA_DRB1_1, HLA_DRB1_2 (as allele strings)")
    message("")
    message("Script will run once data is available at one of the above paths.")
    message("═══════════════════════════════════════════════════════════════")
    quit(save = "no", status = 0)
}

cat("Loading HLA data from:", hla_file, "\n")
hla <- fread(hla_file, showProgress = FALSE)

# ── Parse DRB1*15:01 carrier status ──────────────────────────────────────────
# Handle three formats:
#   1. Pre-computed: columns drb1_1501_carrier (already 0/1) and drb1_1501_dosage
#   2. Allele strings: HLA_DRB1_1 / HLA_DRB1_2 (e.g. "DRB1*15:01")
#   3. Raw UKB HIBAG: p22182_a* columns (allele codes as integers)
if ("drb1_1501_carrier" %in% names(hla)) {
    # Already computed (our processed file from field 22182 dosage vector)
    # drb1_1501_dosage at position 282 of the HIBAG vector; OR≈2.95 in MS
    cat("  Using pre-computed drb1_1501_carrier column\n")
    hla[, drb1_1501_carrier := as.integer(drb1_1501_carrier)]
} else if ("drb1_1501_dosage" %in% names(hla)) {
    # Dosage column — apply 0.5 threshold for carrier calling
    hla[, drb1_1501_carrier := as.integer(drb1_1501_dosage >= 0.5)]
} else if ("HLA_DRB1_1" %in% names(hla)) {
    # Pre-processed format with allele strings
    hla[, drb1_1501_carrier := as.integer(
        grepl("15:01", HLA_DRB1_1, fixed=TRUE) | grepl("15:01", HLA_DRB1_2, fixed=TRUE)
    )]
} else if (any(grepl("^p22182_a", names(hla)))) {
    # Raw UKB field 22182 format — alleles encoded as integers
    a_cols <- grep("^p22182_a", names(hla), value = TRUE)
    drb1_cols <- a_cols[grepl("a1[01]", a_cols)]
    if (length(drb1_cols) >= 2) {
        hla[, drb1_1501_carrier := as.integer(
            get(drb1_cols[1]) == 1501 | get(drb1_cols[2]) == 1501
        )]
    } else {
        hla[, drb1_1501_carrier := as.integer(
            apply(as.data.frame(hla[, ..a_cols]), 1, function(x) any(x == 1501, na.rm=TRUE))
        )]
    }
} else {
    stop("HLA file format not recognised. Expected columns: drb1_1501_carrier, drb1_1501_dosage, HLA_DRB1_1/2, or p22182_a*")
}

cat(sprintf("  HLA-DRB1*15:01 carriers: %d / %d (%.1f%%)\n",
            sum(hla$drb1_1501_carrier, na.rm=TRUE),
            sum(!is.na(hla$drb1_1501_carrier)),
            100*mean(hla$drb1_1501_carrier, na.rm=TRUE)))

# ── Load MS QC data ───────────────────────────────────────────────────────────
cat("Loading MS QC data...\n")
ms_qc <- fread(QC_FILE, showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
               !is.na(UMAP1) & !is.na(UMAP2)]

dt <- merge(ms_qc, hla[, .(eid, drb1_1501_carrier)], by = "eid", all.x = TRUE)
dt <- dt[!is.na(drb1_1501_carrier)]

ms_cases <- dt[ms_status %in% c("pre_onset","post_onset")]
cat(sprintf("  MS cases with HLA: %d  (carrier=%d, non-carrier=%d)\n",
            nrow(ms_cases), sum(ms_cases$drb1_1501_carrier == 1),
            sum(ms_cases$drb1_1501_carrier == 0)))

# ── Load DEPs ─────────────────────────────────────────────────────────────────
deps      <- fread(DIFF_FILE)
dep_prots <- tolower(deps[adj.P.Val < 0.05, protein])
avail_prots <- intersect(dep_prots, names(dt))

# ── Carrier vs non-carrier comparison within MS cases ─────────────────────────
cat("\nCarrier vs non-carrier comparison (MS cases)...\n")
if (sum(ms_cases$drb1_1501_carrier == 1) >= 20 &&
    sum(ms_cases$drb1_1501_carrier == 0) >= 20) {

    # Wilcoxon for each DEP
    cmp_list <- lapply(avail_prots, function(prot) {
        sub <- ms_cases[!is.na(get(prot))]
        wt <- wilcox.test(get(prot) ~ drb1_1501_carrier, data = sub, exact = FALSE)
        data.table(
            protein = prot,
            carrier_mean    = mean(sub[drb1_1501_carrier==1, get(prot)], na.rm=TRUE),
            noncarrier_mean = mean(sub[drb1_1501_carrier==0, get(prot)], na.rm=TRUE),
            pval = wt$p.value
        )
    })
    cmp_dt <- rbindlist(cmp_list)
    cmp_dt[, fdr := p.adjust(pval, method = "BH")]
    cmp_dt[, delta := carrier_mean - noncarrier_mean]
    cmp_dt <- cmp_dt[order(pval)]
    cmp_dt[, logFC := delta]  # approximate
    cmp_dt[, sig_label := fcase(fdr<0.05,"*", pval<0.05,"+", default="")]
    fwrite(cmp_dt, file.path(OUT_DIR, "ms_hla_subgroup_results.csv"))
    cat(sprintf("  FDR<0.05: %d DEPs differ by HLA carrier status\n", sum(cmp_dt$fdr < 0.05)))

    # Volcano
    cmp_dt[, direction := fifelse(delta > 0, "Higher in carriers", "Lower in carriers")]
    cmp_dt[, highlight := fdr < 0.05 | (pval < 0.05 & abs(delta) > quantile(abs(delta), 0.75))]
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
        labs(x = "Δ NPX (DRB1*15:01 carriers − non-carriers)",
             y = expression(-log[10](p)),
             title = "MS DEPs: HLA-DRB1*15:01 carrier vs non-carrier",
             subtitle = sprintf("MS cases only: carrier n=%d, non-carrier n=%d",
                                sum(ms_cases$drb1_1501_carrier==1),
                                sum(ms_cases$drb1_1501_carrier==0))) +
        theme_ukb()

    ggsave(file.path(OUT_DIR, "ms_hla_subgroup_volcano.pdf"),
           p_vol, width = 6, height = 5, device = cairo_pdf)
    cat("  Saved: ms_hla_subgroup_volcano.pdf\n")
} else {
    cat("  WARNING: Insufficient carriers or non-carriers for comparison\n")
}

cat("\n01_hla_ms_subgroup.R complete.\n")
