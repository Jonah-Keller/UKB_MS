#!/usr/bin/env Rscript
# =============================================================================
# 01_limma_ms_vs_hc.R  —  Limma differential expression: cases vs healthy controls
# =============================================================================
# Replicating Abdelhak et al. 2026 (Nature Medicine) approach in UKB Olink data.
# Disease cohort, status-value labels, and output filename prefix are read from
# configs/disease.yaml — no disease-specific values are hardcoded here.
#
# Four contrasts:
#   A. Combined   — (pre_onset + post_onset) vs control  [no ytd; mean_npx adjusted]
#   B. Pre-onset  — pre_onset vs control                 [no ytd; mean_npx adjusted]
#   C. Post-onset — post_onset vs control                [abs(ytd); mean_npx adjusted]
#   D. Pre-onset temporal model (same PSM cohort as B)
#        protein ~ case_binary + ytd_case + mean_npx + age + sex + PC1 + PC2
#        ytd_case = case_binary × years_to_diagnosis (0 for HC)
#        case_binary coefficient → case–HC gap at diagnosis (ytd = 0)
#        ytd_case  coefficient   → per-year change in gap going back from Dx
#                                  positive = gap larger earlier (early/persistent)
#                                  negative = gap grows near Dx (late-emerging)
#
# Covariates: age_at_sampling, sex, mean_npx (technical), PC1/PC2 (control-derived)
# PC derivation: prcomp on healthy controls only; all samples projected in.
#   Rationale: control-derived PCs capture technical variance without disease
#   signal contamination, preventing collinearity with the outcome.
#
# Input:
#   data/ukb/olink/processed/<cohort_short>_olink_qc.csv
#
# Output (filenames templated by cfg$cohort_short):
#   results/differential/<cohort_short>_combined_vs_hc.csv
#   results/differential/<cohort_short>_pre_vs_hc.csv
#   results/differential/<cohort_short>_post_vs_hc.csv
#   results/differential/<cohort_short>_pre_temporal.csv  [Contrast D]
#   results/differential/<cohort_short>_volcano_combined.pdf
#   results/differential/<cohort_short>_volcano_pre.pdf
#   results/differential/<cohort_short>_volcano_post.pdf
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
    library(here)
    library(glue)
})

source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

STATUS_COL  <- cfg$cohort_status_col
PRE         <- cfg$status_values$pre_onset
POST        <- cfg$status_values$post_onset
CONTROL     <- cfg$status_values$control
COMBINED    <- cfg$cohort_short                              # combined-cases factor level

OUT_DIR <- here::here("results", "differential")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(here::here("analysis", "helpers", "ukb_theme.R"))

QC_FILE <- here::here("data", "ukb", "olink", "processed",
                       glue::glue("{cfg$cohort_short}_olink_qc.csv"))

CADASIL_DIR <- file.path(dirname(here::here()),
                          "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")

# PSM match ratio
PSM_RATIO   <- 10   # 1:10 nearest-neighbour matching
PSM_CALIPER <- 0.2  # caliper in SD of logit propensity score

# Thresholds
FDR_THRESH <- 0.05
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
    dt_sub[, is_case := as.integer(get(STATUS_COL) %in% case_levels)]

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
message("Loading QC'd Olink + cohort data...")
if (!file.exists(QC_FILE)) stop("QC file not found: ", QC_FILE, "\nRun 02_olink_qc/01_olink_qc.R first.")
dt <- fread(QC_FILE, showProgress = FALSE)

# Meta columns (UMAP1/UMAP2 kept here only to exclude them from protein_cols)
META_COLS <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance",
               "qc_outlier", "UMAP1", "UMAP2", "mean_npx",
               "bmi", "ever_smoker", "diabetes", "alcohol_freq", "sex_num",
               "PC1", "PC2")
protein_cols <- setdiff(names(dt), META_COLS)
message(sprintf("  %d participants, %d proteins", nrow(dt), length(protein_cols)))
message(sprintf("  Status counts: %s",
                paste(capture.output(print(table(dt[[STATUS_COL]]))), collapse = " ")))

# =============================================================================
# 2. Filter: remove flagged outliers, require complete covariates
# =============================================================================
dt_filt <- dt[qc_outlier == FALSE &
              !is.na(age_at_sampling) &
              !is.na(sex)]

message(sprintf("  After outlier/covariate filter: %d participants", nrow(dt_filt)))
message(sprintf("  Status: %s",
                paste(names(table(dt_filt[[STATUS_COL]])),
                      table(dt_filt[[STATUS_COL]]), sep = "=", collapse = "  ")))

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

# Compute PC1/PC2 from healthy controls only, then project all participants into
# that space.  Control-derived PCs capture technical/batch variance without
# disease-signal contamination, preventing collinearity with the outcome variable.
message("Computing control-only PC1/PC2 (truncated SVD via irlba)...")
ctrl_eids <- as.character(dt_filt[get(STATUS_COL) == CONTROL, eid])
ctrl_x    <- t(prot_mat[, ctrl_eids])                 # samples × proteins
ctrl_ctr  <- colMeans(ctrl_x)
ctrl_xc   <- sweep(ctrl_x, 2L, ctrl_ctr, "-")
sv        <- irlba::irlba(ctrl_xc, nv = 2L)
rotation  <- sv$v
for (k in 1:2) {
    i_max <- which.max(abs(rotation[, k]))
    if (rotation[i_max, k] < 0) rotation[, k] <- -rotation[, k]
}
total_var <- sum(apply(ctrl_xc, 2L, var))
pc_var    <- 100 * (sv$d ^ 2 / max(1L, nrow(ctrl_xc) - 1L)) / total_var
message(sprintf("  Control PC1: %.1f%%  PC2: %.1f%% variance explained", pc_var[1], pc_var[2]))
all_pcs   <- sweep(t(prot_mat), 2L, ctrl_ctr, "-") %*% rotation
for (col in c("PC1", "PC2")) if (col %in% names(dt_filt)) set(dt_filt, j = col, value = NULL)
pc_dt     <- data.table(eid = as.integer(colnames(prot_mat)),
                        PC1 = all_pcs[, 1],
                        PC2 = all_pcs[, 2])
dt_filt   <- merge(dt_filt, pc_dt, by = "eid", all.x = TRUE)

# =============================================================================
# Helper: run_limma() — design + fit + contrasts + topTable
# =============================================================================
run_limma <- function(prot_mat_sub, meta_sub, contrast_str, label,
                      include_mean_npx = FALSE, ytd = NULL) {
    message(sprintf("\nRunning limma: %s  (n=%d)", label, ncol(prot_mat_sub)))

    meta_sub$status_f <- relevel(factor(meta_sub[[STATUS_COL]]), ref = CONTROL)
    meta_sub$sex_f    <- factor(meta_sub$sex)

    if (!is.null(ytd)) {
        meta_sub$ytd <- ytd
        if (include_mean_npx) {
            design <- model.matrix(
                ~ 0 + status_f + ytd + mean_npx + age_at_sampling + sex_f + PC1 + PC2,
                data = meta_sub
            )
        } else {
            design <- model.matrix(
                ~ 0 + status_f + ytd + age_at_sampling + sex_f + PC1 + PC2,
                data = meta_sub
            )
        }
    } else if (include_mean_npx) {
        design <- model.matrix(
            ~ 0 + status_f + mean_npx + age_at_sampling + sex_f + PC1 + PC2,
            data = meta_sub
        )
    } else {
        design <- model.matrix(
            ~ 0 + status_f + age_at_sampling + sex_f + PC1 + PC2,
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
# 4a. Contrast A: Combined MS (pre + post) vs control  —  PSM + mean_npx
# =============================================================================
message("\n--- Contrast A: Combined cases vs HC ---")
dt_combined_raw <- dt_filt[get(STATUS_COL) %in% c(PRE, POST, CONTROL)]
dt_combined_raw <- copy(dt_combined_raw)
dt_combined_raw[get(STATUS_COL) != CONTROL, (STATUS_COL) := COMBINED]

dt_combined <- match_cohort(dt_combined_raw, case_levels = COMBINED)
prot_comb   <- prot_mat[, as.character(dt_combined$eid)]
meta_comb   <- as.data.frame(dt_combined[, c("eid", STATUS_COL, "age_at_sampling",
                                              "sex", "PC1", "PC2", "mean_npx"),
                                          with = FALSE])
rownames(meta_comb) <- as.character(meta_comb$eid)

res_combined <- run_limma(prot_comb, meta_comb,
                          contrast_str     = glue::glue("status_f{COMBINED} - status_f{CONTROL}"),
                          label            = "combined_vs_hc",
                          include_mean_npx = TRUE)

# =============================================================================
# 4b. Contrast B: Pre-onset vs control  —  PSM
# =============================================================================
message("\n--- Contrast B: Pre-onset vs HC ---")
dt_pre_raw <- dt_filt[get(STATUS_COL) %in% c(PRE, CONTROL)]
dt_pre     <- match_cohort(dt_pre_raw, case_levels = PRE)
prot_pre   <- prot_mat[, as.character(dt_pre$eid)]
meta_pre   <- as.data.frame(dt_pre[, c("eid", STATUS_COL, "age_at_sampling",
                                        "sex", "PC1", "PC2", "mean_npx"),
                                    with = FALSE])
rownames(meta_pre) <- as.character(meta_pre$eid)

res_pre <- run_limma(prot_pre, meta_pre,
                     contrast_str     = glue::glue("status_f{PRE} - status_f{CONTROL}"),
                     label            = "pre_onset_vs_hc",
                     include_mean_npx = TRUE)

# =============================================================================
# 4c. Contrast C: Post-onset vs control  —  PSM
# =============================================================================
message("\n--- Contrast C: Post-onset vs HC ---")
dt_post_raw <- dt_filt[get(STATUS_COL) %in% c(POST, CONTROL)]
dt_post     <- match_cohort(dt_post_raw, case_levels = POST)
prot_post   <- prot_mat[, as.character(dt_post$eid)]
meta_post   <- as.data.frame(dt_post[, c("eid", STATUS_COL, "age_at_sampling",
                                          "sex", "PC1", "PC2", "mean_npx",
                                          "years_to_diagnosis"),
                                      with = FALSE])
rownames(meta_post) <- as.character(meta_post$eid)
ytd_c <- abs(meta_post$years_to_diagnosis)
ytd_c[is.na(ytd_c)] <- 0

res_post <- run_limma(prot_post, meta_post,
                      contrast_str     = glue::glue("status_f{POST} - status_f{CONTROL}"),
                      label            = "post_onset_vs_hc",
                      ytd              = ytd_c,
                      include_mean_npx = TRUE)

# =============================================================================
# 4d. Contrast D: Pre-onset temporal model — ytd_case interaction
# Reuses the PSM-matched pre-onset + HC cohort from Contrast B.
#
# ytd_case = case_binary × years_to_diagnosis
#          = years_to_diagnosis  for pre-onset cases (positive = further from Dx)
#          = 0                   for HC
#
# Two coefficients extracted per protein:
#   case_binary → case–HC gap at time of diagnosis (ytd = 0)
#   ytd_case    → per-year change in gap going back from Dx
#                 +ve: gap is larger further from Dx  (early/persistent signal)
#                 −ve: gap grows as Dx approaches     (late-emerging signal)
# =============================================================================
message("\n--- Contrast D: Pre-onset temporal model ---")

dt_pre_t  <- copy(dt_pre)
dt_pre_t[, case_binary := as.integer(get(STATUS_COL) == PRE)]
dt_pre_t[, ytd_case    := case_binary * fifelse(is.na(years_to_diagnosis), 0,
                                                  years_to_diagnosis)]

prot_pre_t <- prot_mat[, as.character(dt_pre_t$eid)]
meta_pre_t <- as.data.frame(dt_pre_t[, .(eid, case_binary, ytd_case, mean_npx,
                                           age_at_sampling, sex, PC1, PC2)])
meta_pre_t$sex_f <- factor(meta_pre_t$sex)
rownames(meta_pre_t) <- as.character(meta_pre_t$eid)

design_t <- model.matrix(
    ~ case_binary + ytd_case + mean_npx + age_at_sampling + sex_f + PC1 + PC2,
    data = meta_pre_t
)
colnames(design_t) <- make.names(colnames(design_t))

fit_t <- lmFit(prot_pre_t, design_t)
fit_t <- eBayes(fit_t)

ms_coef          <- as.data.frame(topTable(fit_t, coef = "case_binary",
                                            sort.by = "P", n = Inf))
ms_coef$protein  <- rownames(ms_coef)
ms_coef$fdr_ms   <- p.adjust(ms_coef$P.Value, method = "BH")

ytd_coef         <- as.data.frame(topTable(fit_t, coef = "ytd_case",
                                             sort.by = "P", n = Inf))
ytd_coef$protein <- rownames(ytd_coef)
ytd_coef$fdr_ytd <- p.adjust(ytd_coef$P.Value, method = "BH")

res_temporal <- merge(
    ms_coef[,  c("protein", "logFC", "t", "P.Value", "fdr_ms")],
    ytd_coef[, c("protein", "logFC", "t", "P.Value", "fdr_ytd")],
    by = "protein", suffixes = c("_ms", "_ytd")
)
setDT(res_temporal)
# Predicted MS–HC gap 5 years before Dx (key summary for figure 2 landscape panel)
res_temporal[, gap_5yr := logFC_ms + 5 * logFC_ytd]

n_ms_sig  <- sum(res_temporal$fdr_ms  < FDR_THRESH)
n_ytd_sig <- sum(res_temporal$fdr_ytd < FDR_THRESH)
message(sprintf("  MS-effect  FDR<%.2f: %d proteins", FDR_THRESH, n_ms_sig))
message(sprintf("  Temporal slope FDR<%.2f: %d proteins", FDR_THRESH, n_ytd_sig))

# =============================================================================
# 5. Write result CSVs
# =============================================================================
fname_combined <- glue::glue("{cfg$cohort_short}_combined_vs_hc.csv")
fname_pre      <- glue::glue("{cfg$cohort_short}_pre_vs_hc.csv")
fname_post     <- glue::glue("{cfg$cohort_short}_post_vs_hc.csv")
fname_temporal <- glue::glue("{cfg$cohort_short}_pre_temporal.csv")

fwrite(res_combined, file.path(OUT_DIR, fname_combined))
fwrite(res_pre,      file.path(OUT_DIR, fname_pre))
fwrite(res_post,     file.path(OUT_DIR, fname_post))
fwrite(res_temporal, file.path(OUT_DIR, fname_temporal))
message("\nResult files written:")
message("  ", fname_combined, "  -- ", sum(res_combined$adj.P.Val < FDR_THRESH), " DEPs at FDR<0.05")
message("  ", fname_pre,      "  -- ", sum(res_pre$adj.P.Val      < FDR_THRESH), " DEPs at FDR<0.05")
message("  ", fname_post,     "  -- ", sum(res_post$adj.P.Val     < FDR_THRESH), " DEPs at FDR<0.05")
message("  ", fname_temporal, "  -- ", n_ms_sig, " case-effect | ", n_ytd_sig, " temporal-slope DEPs")

# =============================================================================
# 6. Volcano plots
# =============================================================================
make_volcano <- function(res, title, out_file) {
    highlight <- cfg$volcano_highlight_proteins
    up_label   <- glue::glue("Up in {cfg$disease_short_caps}")
    down_label <- glue::glue("Down in {cfg$disease_short_caps}")
    res$sig       <- res$adj.P.Val < FDR_THRESH
    res$label     <- ifelse(res$sig & (abs(res$logFC) > LOGFC_LABEL |
                                        res$protein %in% highlight),
                            res$protein, NA)
    res$direction <- case_when(
        res$sig & res$logFC > 0 ~ up_label,
        res$sig & res$logFC < 0 ~ down_label,
        TRUE                    ~ "NS"
    )

    n_up   <- sum(res$direction == up_label,   na.rm = TRUE)
    n_down <- sum(res$direction == down_label, na.rm = TRUE)

    colour_values <- setNames(c(COL_PRE_UP, COL_DOWN, "grey75"),
                               c(up_label, down_label, "NS"))
    colour_labels <- setNames(
        c(sprintf("%s (%d)", up_label, n_up),
          sprintf("%s (%d)", down_label, n_down),
          "NS"),
        c(up_label, down_label, "NS")
    )

    p <- ggplot(res, aes(x = logFC, y = -log10(P.Value),
                          colour = direction, label = label)) +
        geom_point(size = 0.55, alpha = 0.65) +
        geom_text_repel(size = 2.4, max.overlaps = 25, na.rm = TRUE,
                        colour = "grey20", segment.size = 0.3) +
        scale_colour_manual(values = colour_values, labels = colour_labels) +
        geom_hline(yintercept = -log10(max(res$P.Value[res$sig], na.rm = TRUE)),
                   linetype = "dashed", colour = "grey50", linewidth = 0.35) +
        theme_ukb(base_size = 11) +
        theme(
            legend.position  = "bottom",
            legend.key.size  = unit(0.9, "lines"),
            strip.background = element_blank()
        ) +
        labs(title = title,
             x = bquote(log[2]~FC~(.(cfg$disease_short_caps)~vs~HC)),
             y = expression(-log[10](italic(P))),
             colour = NULL)

    ggsave(out_file, p, width = 7, height = 6, device = cairo_pdf)
    message(sprintf("  Volcano saved: %s", basename(out_file)))
}

make_volcano(res_combined,
             glue::glue("{cfg$disease_short_caps} combined vs HC"),
             file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_volcano_combined.pdf")))
make_volcano(res_pre,
             glue::glue("{cfg$disease_short_caps} pre-onset vs HC"),
             file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_volcano_pre.pdf")))
make_volcano(res_post,
             glue::glue("{cfg$disease_short_caps} post-onset vs HC"),
             file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_volcano_post.pdf")))

message("\nLimma cases vs HC complete.")
