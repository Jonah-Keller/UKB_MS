#!/usr/bin/env Rscript
# 03_als_shap_prodromal.R
# SHAP feature importance + prodromal risk-score analysis for ALS — Chia 2025 replication
#
# Approach:
#   - Load best ALS model (glmnet, highest test AUC)
#   - Compute approximate SHAP values via fastshap on test + full pre-onset set
#   - SHAP beeswarm summary plot (top 20 features)
#   - Risk score = P(ALS | model) for all participants
#   - Prodromal trajectory: risk score ~ years_to_diagnosis (pre-onset only)
#   - Linear regression + loess smoother
#
# Input:
#   results/ml/models/als_glmnet_model.rds
#   data/ukb/olink/processed/als_olink_qc.csv
#   results/differential/als_combined_vs_hc.csv
#
# Output: results/ml/
#   als_shap_beeswarm.pdf
#   als_prodromal_risk_score.pdf
#   als_risk_score_all.csv

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(caret)
    library(ggplot2)
    library(ggrepel)
    library(MatchIt)
    library(patchwork)
    library(fastshap)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

QC_FILE    <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "als_olink_qc.csv")
DIFF_FILE  <- file.path(PROJ_DIR, "results", "differential", "als_combined_vs_hc.csv")
CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
MODEL_DIR  <- file.path(PROJ_DIR, "results", "ml", "models")
OUT_DIR    <- file.path(PROJ_DIR, "results", "ml")

FDR_THR   <- 0.05
PSM_RATIO <- 10
PSM_CAL   <- 0.2
set.seed(42)

# ── Load covariates (same as ML script) ───────────────────────────────────────
load_covariates <- function(cadasil_dir) {
    misc <- file.path(cadasil_dir, "misc")
    bmi   <- fread(file.path(misc, "ukb_bmi_age.csv"), showProgress = FALSE)
    setnames(bmi, c("participant.eid","participant.p21001_i0"), c("eid","bmi"), skip_absent=TRUE)
    bmi <- bmi[, .(eid, bmi)]
    smoke <- fread(file.path(misc, "ukb_smoking_status.csv"), showProgress = FALSE)
    setnames(smoke, c("participant.eid","participant.p20116_i0"), c("eid","smoking_code"), skip_absent=TRUE)
    smoke[, ever_smoker := as.integer(smoking_code %in% c(1L,2L))]
    smoke <- smoke[smoking_code != -3L, .(eid, ever_smoker)]
    wmh <- fread(file.path(misc, "wmh_htn_dm.csv"), showProgress = FALSE)
    wmh[, diabetes := as.integer(p2443_i0 == "Yes")]
    wmh <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes)]
    alc <- fread(file.path(misc, "ukb_alcohol_intake.csv"), showProgress = FALSE)
    setnames(alc, c("participant.eid","participant.p1558_i0"), c("eid","alcohol_freq"), skip_absent=TRUE)
    alc <- alc[alcohol_freq > 0L, .(eid, alcohol_freq)]
    Reduce(function(x,y) merge(x,y,by="eid",all.x=TRUE), list(bmi,smoke,wmh,alc))
}

# ── 1. Load model ─────────────────────────────────────────────────────────────
cat("Loading best ALS model (glmnet)...\n")
model_path <- file.path(MODEL_DIR, "als_glmnet_model.rds")
if (!file.exists(model_path)) stop("ALS glmnet model not found at: ", model_path)
model <- readRDS(model_path)
cat(sprintf("  Model loaded. Best CV ROC: %.3f\n", max(model$results$ROC, na.rm=TRUE)))

# ── 2. Rebuild feature set (must match training data) ─────────────────────────
cat("Loading DEPs and QC data...\n")
deps        <- fread(DIFF_FILE)
dep_proteins <- deps[adj.P.Val < FDR_THR, tolower(protein)]

dt   <- fread(QC_FILE, showProgress = FALSE)
covs <- load_covariates(CADASIL_DIR)
dt   <- merge(dt, covs, by = "eid", all.x = TRUE)
dt[, sex_num := as.integer(factor(sex)) - 1L]

PSM_VARS <- c("age_at_sampling","sex_num","bmi","ever_smoker","diabetes","alcohol_freq")
dt_filt  <- dt[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
               !is.na(UMAP1) & !is.na(UMAP2)]
dt_filt[, is_case := as.integer(als_status %in% c("pre_onset","post_onset"))]
dt_c <- dt_filt[complete.cases(as.data.frame(dt_filt[, c("is_case", ..PSM_VARS)]))]

psm_df <- as.data.frame(dt_c[, c("eid","is_case", ..PSM_VARS)])
m <- matchit(is_case ~ age_at_sampling + sex_num + bmi + ever_smoker + diabetes + alcohol_freq,
             data = psm_df, method = "nearest", ratio = PSM_RATIO,
             caliper = PSM_CAL, std.caliper = TRUE)
matched_eids <- as.integer(match.data(m)$eid)
dt_ml <- dt_c[eid %in% matched_eids]

avail_prots <- intersect(dep_proteins, names(dt_ml))
dt_ml[, outcome := factor(fifelse(is_case == 1L, "ALS", "Control"), levels = c("ALS","Control"))]
feat_cols  <- c("eid", "outcome", avail_prots)
ml_dt      <- dt_ml[, ..feat_cols]
ml_dt_full <- ml_dt[complete.cases(ml_dt[, -"eid"])]  # track eids
ml_eids    <- ml_dt_full$eid
ml_dt      <- ml_dt_full[, -"eid"]  # drop eid for modeling

set.seed(42)
ml_df     <- as.data.frame(ml_dt)
train_idx <- createDataPartition(ml_df$outcome, p = 0.80, list = FALSE)
train_df  <- ml_df[ train_idx, ]
test_df   <- ml_df[-train_idx, ]

feature_cols <- avail_prots
cat(sprintf("  %d features, train=%d, test=%d\n",
            length(feature_cols), nrow(train_df), nrow(test_df)))

# ── 3. SHAP values (fastshap) ─────────────────────────────────────────────────
cat("Computing SHAP values on test set...\n")
pfun <- function(object, newdata) {
    predict(object, newdata = as.data.frame(newdata), type = "prob")[, "ALS"]
}

# Use train_df as background, test_df as explain set (for speed)
X_train <- as.matrix(train_df[, feature_cols])
X_test  <- as.matrix(test_df[, feature_cols])

set.seed(2024)
shap_vals <- tryCatch(
    fastshap::explain(
        object      = model,
        X           = X_train,
        pred_wrapper = pfun,
        newdata     = X_test,
        nsim        = 50,
        adjust      = TRUE
    ),
    error = function(e) { message("SHAP failed: ", e$message); NULL }
)

if (!is.null(shap_vals)) {
    shap_dt <- as.data.table(shap_vals)
    shap_long <- melt(shap_dt, measure.vars = feature_cols,
                      variable.name = "protein", value.name = "shap")
    # Feature values for colour
    feat_dt   <- as.data.table(X_test)
    feat_long <- melt(feat_dt, measure.vars = feature_cols,
                      variable.name = "protein", value.name = "feature_val")
    shap_long[, feature_val := feat_long$feature_val]

    # Mean |SHAP| rank
    shap_importance <- shap_long[, .(mean_abs_shap = mean(abs(shap))), by = protein]
    shap_importance <- shap_importance[order(-mean_abs_shap)]
    top_prots <- shap_importance[1:min(20, .N), protein]

    shap_plot <- shap_long[protein %in% top_prots]
    shap_plot[, protein_f := factor(protein,
                                    levels = rev(shap_importance[protein %in% top_prots, protein]))]

    p_shap <- ggplot(shap_plot, aes(x = shap, y = protein_f, colour = feature_val)) +
        geom_jitter(height = 0.2, size = 0.7, alpha = 0.6) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
        scale_colour_gradient2(low = "#56B4E9", mid = "grey90", high = "#CC0066",
                               midpoint = 0, name = "Feature\nvalue") +
        labs(x = "SHAP value (impact on ALS probability)",
             y = NULL, title = "ALS classifier: SHAP feature importance (test set)") +
        theme_ukb() +
        theme(axis.text.y = element_text(size = 7))

    ggsave(file.path(OUT_DIR, "als_shap_beeswarm.pdf"),
           p_shap, width = 6, height = 5.5, device = cairo_pdf)
    cat("  Saved: als_shap_beeswarm.pdf\n")
} else {
    cat("  Skipping SHAP beeswarm (fastshap failed)\n")
}

# ── 4. Risk score for all participants ────────────────────────────────────────
cat("\nComputing risk scores for all matched participants...\n")
X_all   <- as.matrix(ml_df[, feature_cols])
risk_all <- predict(model, newdata = as.data.frame(X_all), type = "prob")[, "ALS"]

risk_dt <- data.table(
    eid        = ml_eids,
    outcome    = ml_df$outcome,
    risk_score = risk_all
)
# Attach years_to_diagnosis
ytd <- dt_ml[, .(eid, years_to_diagnosis, als_status)]
risk_dt <- merge(risk_dt, ytd, by = "eid", all.x = TRUE)

fwrite(risk_dt, file.path(OUT_DIR, "als_risk_score_all.csv"))
cat("  Saved: als_risk_score_all.csv\n")
cat(sprintf("  ALS median risk: %.3f  Control median risk: %.3f\n",
            median(risk_dt[outcome == "ALS", risk_score]),
            median(risk_dt[outcome == "Control", risk_score])))

# ── 5. Prodromal trajectory ────────────────────────────────────────────────────
cat("\nBuilding prodromal trajectory...\n")
pre_dt <- risk_dt[als_status == "pre_onset" & !is.na(years_to_diagnosis)]
cat(sprintf("  Pre-onset cases with years_to_diagnosis: %d\n", nrow(pre_dt)))

if (nrow(pre_dt) >= 10) {
    # Linear regression
    lm_fit   <- lm(risk_score ~ years_to_diagnosis, data = pre_dt)
    lm_coef  <- coef(summary(lm_fit))
    lm_r2    <- summary(lm_fit)$r.squared
    slope_p  <- lm_coef["years_to_diagnosis", "Pr(>|t|)"]
    cat(sprintf("  Linear regression: beta=%.4f, p=%.3g, R2=%.3f\n",
                lm_coef["years_to_diagnosis","Estimate"], slope_p, lm_r2))

    # X-axis: years relative to diagnosis (negative = before)
    # Reference line: control median
    ctrl_med <- median(risk_dt[outcome == "Control", risk_score])

    p_prod <- ggplot(pre_dt, aes(x = years_to_diagnosis, y = risk_score)) +
        geom_hline(yintercept = ctrl_med, linetype = "dashed",
                   colour = "grey50", linewidth = 0.4) +
        annotate("text", x = min(pre_dt$years_to_diagnosis, na.rm=TRUE) + 0.5,
                 y = ctrl_med + 0.01, label = "Ctrl median",
                 size = 2.5, colour = "grey50", hjust = 0) +
        geom_point(colour = "#CC0066", alpha = 0.65, size = 1.8) +
        geom_smooth(method = "loess", span = 0.9,
                    colour = "#2B4C7E", fill = "#2B4C7E", alpha = 0.15,
                    linewidth = 0.8, se = TRUE) +
        geom_smooth(method = "lm", se = FALSE,
                    colour = "#E6A817", linewidth = 0.6, linetype = "dashed") +
        annotate("text",
                 x = max(pre_dt$years_to_diagnosis, na.rm=TRUE) - 0.5,
                 y = max(pre_dt$risk_score, na.rm=TRUE) - 0.02,
                 label = sprintf("β=%.3f, p=%s", lm_coef["years_to_diagnosis","Estimate"],
                                 ifelse(slope_p < 0.001, "<0.001", sprintf("%.3f", slope_p))),
                 size = 2.8, hjust = 1) +
        labs(x = "Years to ALS diagnosis (negative = pre-onset)",
             y = "ALS risk score (P[ALS | model])",
             title = "ALS prodromal risk trajectory",
             subtitle = sprintf("n=%d pre-onset cases", nrow(pre_dt))) +
        coord_cartesian(ylim = c(0, 1)) +
        theme_ukb()

    # Panel B: boxplot comparing pre/post vs control
    compare_dt <- risk_dt[!is.na(als_status)]
    compare_dt[, group := fcase(
        als_status == "pre_onset",  "Pre-onset ALS",
        als_status == "post_onset", "Post-onset ALS",
        als_status == "control",    "Control"
    )]
    compare_dt <- compare_dt[!is.na(group)]
    compare_dt[, group_f := factor(group, levels = c("Control","Pre-onset ALS","Post-onset ALS"))]

    group_cols <- c("Control" = "grey60", "Pre-onset ALS" = "#E6A817", "Post-onset ALS" = "#CC0066")

    ctrl_sample <- compare_dt[group == "Control"][sample(.N, min(500L, .N))]
    plot_dt     <- rbind(compare_dt[group != "Control"], ctrl_sample)

    p_box <- ggplot(plot_dt, aes(x = group_f, y = risk_score, fill = group)) +
        geom_boxplot(outlier.size = 0.5, width = 0.5, colour = "#1A1A1A") +
        scale_fill_manual(values = group_cols, guide = "none") +
        labs(x = NULL, y = "Risk score", title = "ALS risk score by group") +
        coord_cartesian(ylim = c(0, 1)) +
        theme_ukb()

    p_final <- p_prod + p_box + plot_layout(widths = c(2.5, 1))
    ggsave(file.path(OUT_DIR, "als_prodromal_risk_score.pdf"),
           p_final, width = 9, height = 4, device = cairo_pdf)
    cat("  Saved: als_prodromal_risk_score.pdf\n")
} else {
    cat("  Too few pre-onset cases for trajectory plot\n")
}

cat("\n03_als_shap_prodromal.R complete.\n")
