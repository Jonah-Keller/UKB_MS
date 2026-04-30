#!/usr/bin/env Rscript
# =============================================================================
# 02_limma_als_vs_hc.R  —  Limma differential expression: ALS vs healthy controls
# =============================================================================
# Replicating Chia 2025 approach in UKB Olink data.
#
# Three contrasts:
#   A. Combined  — (pre_onset + post_onset) vs control
#   B. Pre-onset — pre_onset vs control           (presymptomatic)
#   C. Post-onset — post_onset vs control         (after clinical onset)
#
# Model (per Chia 2025 template):
#   ~0 + als_status + age_at_sampling + sex + UMAP1 + UMAP2
#
# Input:
#   data/ukb/olink/processed/als_olink_qc.csv
#
# Output:
#   results/differential/als_combined_vs_hc.csv
#   results/differential/als_pre_vs_hc.csv
#   results/differential/als_post_vs_hc.csv
#   results/differential/als_volcano_combined.pdf
#   results/differential/als_volcano_pre.pdf
#   results/differential/als_volcano_post.pdf
# =============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(ggrepel)
    library(stringr)
    library(MatchIt)
    library(splines)
})

# --- Paths ---
.args      <- commandArgs(trailingOnly = FALSE)
.file_arg  <- .args[grepl("^--file=", .args)]
SCRIPT_DIR <- if (length(.file_arg) > 0) dirname(normalizePath(sub("^--file=", "", .file_arg))) else getwd()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
OUT_DIR    <- file.path(REPO_ROOT, "results", "differential")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(SCRIPT_DIR, "..", "helpers", "ukb_theme.R"))

QC_FILE <- file.path(REPO_ROOT, "data", "ukb", "olink", "processed", "als_olink_qc.csv")

CADASIL_DIR <- file.path(dirname(REPO_ROOT), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")

# If TRUE, restrict ALS cases to those with als_confirmed_at_death == 1
# (MND confirmed at death; excludes G12.x with uncertain coding).
# Default FALSE keeps all G12 cases.
ALS_CONFIRMED_ONLY <- FALSE

# PSM match ratio
PSM_RATIO   <- 10   # 1:10 nearest-neighbour matching
PSM_CALIPER <- 0.2  # caliper in SD of logit propensity score

# Thresholds
FDR_THRESH  <- 0.05
LOGFC_LABEL <- 0.3   # label DEPs on volcano with |logFC| > this AND FDR < 0.05

set.seed(42)

# =============================================================================
# Helper: load_covariates() — BMI, smoking, diabetes, alcohol from CADASIL_DIR
# =============================================================================
load_covariates <- function(cadasil_dir) {
    misc <- file.path(cadasil_dir, "misc")

    # BMI (field p21001_i0)
    bmi <- fread(file.path(misc, "ukb_bmi_age.csv"), showProgress = FALSE)
    setnames(bmi, c("participant.eid", "participant.p21001_i0"), c("eid", "bmi"),
             skip_absent = TRUE)
    bmi <- bmi[, .(eid, bmi)]

    # Smoking (0=never, 1=previous, 2=current; −3=prefer not to answer)
    smoke <- fread(file.path(misc, "ukb_smoking_status.csv"), showProgress = FALSE)
    setnames(smoke, c("participant.eid", "participant.p20116_i0"),
             c("eid", "smoking_code"), skip_absent = TRUE)
    smoke[, ever_smoker := as.integer(smoking_code %in% c(1L, 2L))]
    smoke <- smoke[smoking_code != -3L, .(eid, ever_smoker)]

    # Diabetes (p2443_i0: "Yes" / "No" / "Borderline" / "Prefer not to answer")
    wmh <- fread(file.path(misc, "wmh_htn_dm.csv"), showProgress = FALSE)
    wmh[, diabetes := as.integer(p2443_i0 == "Yes")]
    wmh <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes)]

    # Alcohol intake frequency (1–6 ordinal; −3=prefer not to answer)
    alc <- fread(file.path(misc, "ukb_alcohol_intake.csv"), showProgress = FALSE)
    setnames(alc, c("participant.eid", "participant.p1558_i0"),
             c("eid", "alcohol_freq"), skip_absent = TRUE)
    alc <- alc[alcohol_freq > 0L, .(eid, alcohol_freq)]

    Reduce(function(x, y) merge(x, y, by = "eid", all.x = TRUE),
           list(bmi, smoke, wmh, alc))
}

# =============================================================================
# Helper: match_cohort() — propensity-score match cases to controls
# Returns a data.table containing only the matched participants.
# =============================================================================
match_cohort <- function(dt_sub, case_levels,
                         ratio   = PSM_RATIO,
                         caliper = PSM_CALIPER) {
    dt_sub <- copy(dt_sub)
    dt_sub[, is_case := as.integer(als_status %in% case_levels)]

    PSM_VARS <- c("age_at_sampling", "sex_num", "bmi", "ever_smoker",
                  "diabetes", "alcohol_freq")
    dt_sub[, sex_num := as.integer(factor(sex)) - 1L]  # 0/1

    complete_mask <- complete.cases(
        as.data.frame(dt_sub[, c("is_case", ..PSM_VARS)])
    )
    dt_c <- dt_sub[complete_mask]
    message(sprintf("  PSM input: %d cases, %d controls (%d excluded: missing covariates)",
                    sum(dt_c$is_case), sum(!dt_c$is_case),
                    nrow(dt_sub) - nrow(dt_c)))

    psm_df <- as.data.frame(dt_c[, c("eid", "is_case", ..PSM_VARS)])

    m <- matchit(
        is_case ~ age_at_sampling + sex_num + bmi + ever_smoker +
                  diabetes + alcohol_freq,
        data        = psm_df,
        method      = "nearest",
        ratio       = ratio,
        caliper     = caliper,
        std.caliper = TRUE
    )

    matched_eids <- as.integer(match.data(m)$eid)
    dt_matched   <- dt_c[eid %in% matched_eids]
    message(sprintf("  After PSM: n=%d  (%d cases, %d controls; ratio ~1:%g)",
                    nrow(dt_matched),
                    sum(dt_matched$is_case),
                    sum(!dt_matched$is_case),
                    sum(!dt_matched$is_case) / max(1, sum(dt_matched$is_case))))
    return(dt_matched)
}

# =============================================================================
# 1. Load QC'd data
# =============================================================================
message("Loading QC'd Olink + ALS cohort data...")
if (!file.exists(QC_FILE)) stop("QC file not found: ", QC_FILE, "\nRun 02_olink_qc/01_olink_qc.R first.")
dt <- fread(QC_FILE, showProgress = FALSE)

# Meta columns
META_COLS <- c("eid", "als_status", "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance",
               "qc_outlier", "UMAP1", "UMAP2", "mean_npx",
               "als_confirmed_at_death",
               "bmi", "ever_smoker", "diabetes", "alcohol_freq", "sex_num")
protein_cols <- setdiff(names(dt), META_COLS)
message(sprintf("  %d participants, %d proteins", nrow(dt), length(protein_cols)))
message(sprintf("  Status counts: %s", paste(capture.output(print(table(dt$als_status))), collapse = " ")))

# =============================================================================
# 1b. Optional: restrict to MND-confirmed cases only
# =============================================================================
if (ALS_CONFIRMED_ONLY) {
    dt[als_status != "control" & als_confirmed_at_death != 1, als_status := "control"]
    message("  ALS_CONFIRMED_ONLY: restricted to MND-confirmed cases")
}

# =============================================================================
# 2. Filter: remove flagged outliers, require complete covariates
# =============================================================================
dt_filt <- dt[qc_outlier == FALSE &
              !is.na(age_at_sampling) &
              !is.na(sex) &
              !is.na(UMAP1) &
              !is.na(UMAP2)]

message(sprintf("  After outlier/covariate filter: %d participants", nrow(dt_filt)))
message(sprintf("  Status: %s", paste(names(table(dt_filt$als_status)),
                                       table(dt_filt$als_status), sep = "=", collapse = "  ")))

# =============================================================================
# 2b. Load and merge comorbidity covariates for PSM
# =============================================================================
message("Loading comorbidity covariates for PSM...")
covs <- load_covariates(CADASIL_DIR)
covs[, eid := as.integer(eid)]
dt_filt <- merge(dt_filt, covs, by = "eid", all.x = TRUE)
message(sprintf("  BMI available:      %d / %d",   sum(!is.na(dt_filt$bmi)),          nrow(dt_filt)))
message(sprintf("  Smoking available:  %d / %d",   sum(!is.na(dt_filt$ever_smoker)),   nrow(dt_filt)))
message(sprintf("  Diabetes available: %d / %d",   sum(!is.na(dt_filt$diabetes)),      nrow(dt_filt)))
message(sprintf("  Alcohol available:  %d / %d",   sum(!is.na(dt_filt$alcohol_freq)),  nrow(dt_filt)))

# =============================================================================
# 3. Protein matrix: proteins (rows) × samples (cols)
#    Mean-impute any remaining missing values per protein
# =============================================================================
prot_mat <- t(as.matrix(dt_filt[, ..protein_cols]))
colnames(prot_mat) <- as.character(dt_filt$eid)  # character so subsetting uses names not indices

# Impute missing with per-protein mean
row_means <- rowMeans(prot_mat, na.rm = TRUE)
for (i in seq_len(nrow(prot_mat))) {
    nas <- is.na(prot_mat[i, ])
    if (any(nas)) prot_mat[i, nas] <- row_means[i]
}

# =============================================================================
# Helper: run_limma() — design + fit + contrasts + topTable
# =============================================================================
run_limma <- function(prot_mat_sub, meta_sub, contrast_str, label,
                      include_mean_npx = FALSE) {
    message(sprintf("\nRunning limma: %s  (n=%d)", label, ncol(prot_mat_sub)))

    meta_sub$als_status_f <- relevel(factor(meta_sub$als_status), ref = "control")
    meta_sub$sex_f         <- factor(meta_sub$sex)

    if (include_mean_npx) {
        design <- model.matrix(
            ~ 0 + als_status_f + mean_npx + age_at_sampling + sex_f + UMAP1 + UMAP2,
            data = meta_sub
        )
    } else {
        design <- model.matrix(
            ~ 0 + als_status_f + age_at_sampling + sex_f + UMAP1 + UMAP2,
            data = meta_sub
        )
    }
    colnames(design) <- make.names(colnames(design))

    fit  <- lmFit(prot_mat_sub, design)
    contr <- makeContrasts(contrasts = contrast_str, levels = design)
    fit2  <- contrasts.fit(fit, contr)
    fit2  <- eBayes(fit2)

    sig_counts <- summary(decideTests(fit2, p.value = FDR_THRESH))
    message(sprintf("  DEPs at FDR<%.2f: %d up, %d down",
                    FDR_THRESH, sig_counts["Up", 1], sig_counts["Down", 1]))

    res <- as.data.frame(topTable(fit2, sort.by = "P", n = Inf))
    res$protein    <- rownames(res)
    res$comparison <- label
    res <- res[, c("protein", "comparison", "logFC", "AveExpr", "t",
                   "P.Value", "adj.P.Val", "B")]
    return(res)
}

# =============================================================================
# 4a. Contrast A: Combined ALS (pre + post) vs control  —  PSM + mean_npx
# =============================================================================
message("\n--- Contrast A: Combined ALS vs HC ---")
dt_combined_raw <- dt_filt[als_status %in% c("pre_onset", "post_onset", "control")]
dt_combined_raw <- copy(dt_combined_raw)
dt_combined_raw[als_status != "control", als_status := "als"]

dt_combined <- match_cohort(dt_combined_raw, case_levels = "als")
prot_comb   <- prot_mat[, as.character(dt_combined$eid)]
meta_comb   <- as.data.frame(dt_combined[, .(eid, als_status, age_at_sampling,
                                              sex, UMAP1, UMAP2, mean_npx)])
rownames(meta_comb) <- as.character(meta_comb$eid)

res_combined <- run_limma(prot_comb, meta_comb,
                          contrast_str     = "als_status_fals - als_status_fcontrol",
                          label            = "combined_vs_hc",
                          include_mean_npx = TRUE)

# =============================================================================
# 4b. Contrast B: Pre-onset vs control  —  PSM
# =============================================================================
message("\n--- Contrast B: Pre-onset vs HC ---")
dt_pre_raw <- dt_filt[als_status %in% c("pre_onset", "control")]
dt_pre     <- match_cohort(dt_pre_raw, case_levels = "pre_onset")
prot_pre   <- prot_mat[, as.character(dt_pre$eid)]
meta_pre   <- as.data.frame(dt_pre[, .(eid, als_status, age_at_sampling,
                                        sex, UMAP1, UMAP2)])
rownames(meta_pre) <- as.character(meta_pre$eid)

res_pre <- run_limma(prot_pre, meta_pre,
                     contrast_str = "als_status_fpre_onset - als_status_fcontrol",
                     label        = "pre_onset_vs_hc")

# =============================================================================
# 4c. Contrast C: Post-onset vs control  —  PSM
# =============================================================================
message("\n--- Contrast C: Post-onset vs HC ---")
dt_post_raw <- dt_filt[als_status %in% c("post_onset", "control")]
dt_post     <- match_cohort(dt_post_raw, case_levels = "post_onset")
prot_post   <- prot_mat[, as.character(dt_post$eid)]
meta_post   <- as.data.frame(dt_post[, .(eid, als_status, age_at_sampling,
                                          sex, UMAP1, UMAP2)])
rownames(meta_post) <- as.character(meta_post$eid)

res_post <- run_limma(prot_post, meta_post,
                      contrast_str = "als_status_fpost_onset - als_status_fcontrol",
                      label        = "post_onset_vs_hc")

# =============================================================================
# 5. Write result CSVs
# =============================================================================
fwrite(res_combined, file.path(OUT_DIR, "als_combined_vs_hc.csv"))
fwrite(res_pre,      file.path(OUT_DIR, "als_pre_vs_hc.csv"))
fwrite(res_post,     file.path(OUT_DIR, "als_post_vs_hc.csv"))
message("\nResult files written:")
message("  als_combined_vs_hc.csv  — ", sum(res_combined$adj.P.Val < FDR_THRESH), " DEPs at FDR<0.05")
message("  als_pre_vs_hc.csv       — ", sum(res_pre$adj.P.Val < FDR_THRESH),      " DEPs at FDR<0.05")
message("  als_post_vs_hc.csv      — ", sum(res_post$adj.P.Val < FDR_THRESH),     " DEPs at FDR<0.05")

# =============================================================================
# 6. Volcano plots
# =============================================================================
make_volcano <- function(res, title, out_file) {
    highlight <- c("nefl", "gfap", "chit1", "stmn2", "nefh", "nefm",
                   "tardbp", "fus", "sod1", "tff3", "cd163", "lgals3")

    # Use shared VOLCANO_COLORS (redefined here with ALS label for axis clarity)
    VOLCANO_COLORS <- c("Up in ALS" = "#CC0066", "Down in ALS" = "#56B4E9", "NS" = "grey75")

    res$sig       <- res$adj.P.Val < FDR_THRESH
    res$label     <- ifelse(res$sig & (abs(res$logFC) > LOGFC_LABEL |
                                        res$protein %in% highlight),
                            res$protein, NA)
    res$direction <- case_when(
        res$sig & res$logFC > 0 ~ "Up in ALS",
        res$sig & res$logFC < 0 ~ "Down in ALS",
        TRUE                    ~ "NS"
    )

    n_up   <- sum(res$direction == "Up in ALS",   na.rm = TRUE)
    n_down <- sum(res$direction == "Down in ALS", na.rm = TRUE)

    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value),
                          colour = direction, label = label)) +
        geom_point(size = 0.55, alpha = 0.65) +
        geom_text_repel(size = 2.4, max.overlaps = 25, na.rm = TRUE,
                        colour = "grey20", segment.size = 0.3) +
        scale_colour_manual(
            values = VOLCANO_COLORS,
            labels = c(
                "Up in ALS"   = sprintf("Up in ALS (%d)", n_up),
                "Down in ALS" = sprintf("Down in ALS (%d)", n_down),
                "NS"          = "NS"
            )
        ) +
        geom_hline(yintercept = -log10(max(res$P.Value[res$sig], na.rm = TRUE)),
                   linetype = "dashed", colour = "grey50", linewidth = 0.35) +
        theme_classic(base_size = 11) +
        theme(
            legend.position  = "bottom",
            legend.key.size  = unit(0.9, "lines"),
            strip.background = element_blank()
        ) +
        labs(title = title, x = expression(log[2]~FC~(ALS~vs~HC)),
             y = expression(-log[10](italic(P))), colour = NULL)

    ggsave(out_file, p, width = 7, height = 6, device = cairo_pdf)
    message(sprintf("  Volcano saved: %s", basename(out_file)))
}

make_volcano(res_combined, "ALS combined vs HC",
             file.path(OUT_DIR, "als_volcano_combined.pdf"))
make_volcano(res_pre, "ALS pre-onset vs HC",
             file.path(OUT_DIR, "als_volcano_pre.pdf"))
make_volcano(res_post, "ALS post-onset vs HC",
             file.path(OUT_DIR, "als_volcano_post.pdf"))

message("\nLimma ALS vs HC complete.")
