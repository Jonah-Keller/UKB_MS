#!/usr/bin/env Rscript
# 04_ms_preonset_ml.R
# Presymptomatic disease ML classifier — Abdelhak 2026 replication
#
# Approach (per Abdelhak Methods):
#   - Features: top 50 combined DEPs (FDR < 0.05)
#   - Restricted to pre-onset cases + PSM-matched controls (1:5)
#   - Stratified 80/20 split
#   - 5-fold × 10 repeats, ROC metric
#   - Models: glmnet, gbm, rf, nnet
#   - Paper AUC: ~0.81 (random forest)
#
# Cohort, status column, and filename prefixes are read from configs/disease.yaml.
#
# Input:
#   data/ukb/olink/processed/{cohort_short}_olink_qc.csv
#   results/differential/{cohort_short}_combined_vs_hc.csv
#
# Output: results/ml/
#   {cohort_short}_preonset_ml_roc_curves.pdf
#   {cohort_short}_preonset_ml_model_comparison.pdf
#   {cohort_short}_preonset_ml_feature_importance.pdf
#   {cohort_short}_preonset_ml_results.csv
#   models/  ({cohort_short}_pre_*.rds)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(caret)
    library(pROC)
    library(ggplot2)
    library(MatchIt)
    library(patchwork)
    library(Boruta)
    library(glue)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config(file.path(PROJ_DIR, "configs", "disease.yaml"))
COHORT     <- cfg$cohort_short
DISEASE    <- cfg$disease_short_caps
STATUS_COL <- cfg$cohort_status_col
PRE_VAL    <- cfg$status_values$pre_onset
CTRL_VAL   <- cfg$status_values$control
POS_LBL    <- DISEASE
NEG_LBL    <- "Control"

QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed",
                         glue("{COHORT}_olink_qc.csv"))
DIFF_FILE   <- file.path(PROJ_DIR, "results", "differential",
                         glue("{COHORT}_combined_vs_hc.csv"))
CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
OUT_DIR     <- file.path(PROJ_DIR, "results", "ml")
MODEL_DIR   <- file.path(OUT_DIR, "models")
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR    <- 0.05
PSM_RATIO  <- 5          # 1:5 matching (per Abdelhak Methods)
PSM_CAL    <- 0.2
TRAIN_PROP <- 0.80
set.seed(42)

# Demographics always forced into the model regardless of Boruta protein selection
# Only age + sex: guaranteed non-missing for all participants
# BMI/smoking kept only for PSM matching (not as model features)
DEMO_FEATS <- c("age_at_sampling", "sex_num")

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
    wmh[, sbp      := as.numeric(p4080_i0_a0)]   # systolic blood pressure
    wmh <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes, sbp)]
    alc <- fread(file.path(misc, "ukb_alcohol_intake.csv"), showProgress = FALSE)
    setnames(alc, c("participant.eid","participant.p1558_i0"), c("eid","alcohol_freq"), skip_absent=TRUE)
    alc <- alc[alcohol_freq > 0L, .(eid, alcohol_freq)]
    Reduce(function(x,y) merge(x,y,by="eid",all.x=TRUE), list(bmi,smoke,wmh,alc))
}

# ── 1. DEPs ────────────────────────────────────────────────────────────────────
cat("Loading DEPs...\n")
deps         <- fread(DIFF_FILE)
dep_proteins <- deps[adj.P.Val < FDR_THR][order(adj.P.Val)][1:min(50,.N), tolower(protein)]
cat(sprintf("  %d top DEPs\n", length(dep_proteins)))

# ── 2. QC data ─────────────────────────────────────────────────────────────────
cat("Loading QC data...\n")
dt   <- fread(QC_FILE, showProgress = FALSE)
covs <- load_covariates(CADASIL_DIR)
dt   <- merge(dt, covs, by = "eid", all.x = TRUE)
dt[, sex_num := as.integer(factor(sex)) - 1L]

dt_filt <- dt[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
              !is.na(UMAP1) & !is.na(UMAP2)]

# ── 3. RESTRICT to pre-onset cases + controls ──────────────────────────────────
cat("\nRestricting to pre-onset cases...\n")
dt_pre <- dt_filt[get(STATUS_COL) %in% c(PRE_VAL, CTRL_VAL)]
cat(sprintf("  Pre-onset %s: %d, Controls: %d\n", DISEASE,
            sum(dt_pre[[STATUS_COL]] == PRE_VAL),
            sum(dt_pre[[STATUS_COL]] == CTRL_VAL)))

PSM_VARS <- c("age_at_sampling","sex_num","bmi","ever_smoker","diabetes","alcohol_freq")
# SBP excluded: not in paper's PSM spec, and its missingness costs ~11 pre-onset cases

dt_pre[, is_case := as.integer(get(STATUS_COL) == PRE_VAL)]

# Keep PSM vars that are available; drop any that are all-NA
psm_avail <- PSM_VARS[PSM_VARS %in% names(dt_pre) &
                       sapply(PSM_VARS[PSM_VARS %in% names(dt_pre)],
                              function(v) mean(is.na(dt_pre[[v]])) < 0.5)]
dt_c <- dt_pre[complete.cases(as.data.frame(dt_pre[, c("is_case", ..psm_avail)]))]
cat(sprintf("  After covariate filtering: %d cases, %d controls\n",
            sum(dt_c$is_case), sum(!dt_c$is_case)))
cat(sprintf("  PSM variables (%d): %s\n", length(psm_avail), paste(psm_avail, collapse=", ")))

psm_formula <- as.formula(paste("is_case ~", paste(psm_avail, collapse = " + ")))
psm_df <- as.data.frame(dt_c[, c("eid","is_case", ..psm_avail)])
m <- matchit(psm_formula, data = psm_df, method = "nearest", ratio = PSM_RATIO,
             caliper = PSM_CAL, std.caliper = TRUE)
matched_eids <- as.integer(match.data(m)$eid)
dt_ml <- dt_c[eid %in% matched_eids]
cat(sprintf("  After PSM 1:%d: %d cases, %d controls\n",
            PSM_RATIO, sum(dt_ml$is_case), sum(!dt_ml$is_case)))

# ── 4. Feature matrix ──────────────────────────────────────────────────────────
avail_prots <- intersect(dep_proteins, names(dt_ml))
avail_demo  <- intersect(DEMO_FEATS,   names(dt_ml))   # demographics always in model
cat(sprintf("  Proteins: %d DEPs | Demographics: %s\n",
            length(avail_prots), paste(avail_demo, collapse=", ")))

dt_ml[, outcome := factor(fifelse(is_case == 1L, POS_LBL, NEG_LBL),
                          levels = c(POS_LBL, NEG_LBL))]
feat_cols <- c("outcome", avail_prots, avail_demo)
ml_dt     <- dt_ml[, ..feat_cols]

# Median-impute protein NAs to avoid dropping MS cases due to protein missingness.
# Demographics (age_at_sampling, sex_num) are always complete; only proteins need imputation.
for (col in avail_prots) {
    if (col %in% names(ml_dt) && any(is.na(ml_dt[[col]]))) {
        med_val <- median(ml_dt[[col]], na.rm = TRUE)
        ml_dt[is.na(get(col)), (col) := med_val]
    }
}
ml_dt <- ml_dt[complete.cases(ml_dt)]   # now only removes rows with missing demographics
cat(sprintf("  Complete cases: %d  (%s=%d, %s=%d)\n",
            nrow(ml_dt), POS_LBL, sum(ml_dt$outcome == POS_LBL),
            NEG_LBL, sum(ml_dt$outcome == NEG_LBL)))

if (sum(ml_dt$outcome == POS_LBL) < 20) {
    stop("Too few pre-onset ", DISEASE, " cases for ML (n=",
         sum(ml_dt$outcome == POS_LBL), ")")
}

# ── 5. 80/20 split ─────────────────────────────────────────────────────────────
ml_df     <- as.data.frame(ml_dt)
train_idx <- createDataPartition(ml_df$outcome, p = TRAIN_PROP, list = FALSE)
train_df  <- ml_df[ train_idx, ]
test_df   <- ml_df[-train_idx, ]
cat(sprintf("  Train: %d  (%s=%d, %s=%d)\n",
            nrow(train_df), POS_LBL, sum(train_df$outcome == POS_LBL),
            NEG_LBL, sum(train_df$outcome == NEG_LBL)))
cat(sprintf("  Test : %d  (%s=%d, %s=%d)\n",
            nrow(test_df), POS_LBL, sum(test_df$outcome == POS_LBL),
            NEG_LBL, sum(test_df$outcome == NEG_LBL)))

# ── 5b. Boruta feature selection (paper: LOO 5-fold, features in ≥2 folds) ────
BORUTA_FILE <- file.path(MODEL_DIR, glue("{COHORT}_pre_boruta_features.rds"))

if (file.exists(BORUTA_FILE)) {
    cat("\nLoading cached Boruta features...\n")
    boruta_feats <- readRDS(BORUTA_FILE)
} else {
    cat("\nRunning Boruta feature selection (5-fold LOO on training set)...\n")
    set.seed(42)
    fold_idx <- createFolds(train_df$outcome, k = 5, list = TRUE)
    boruta_confirmed <- lapply(seq_along(fold_idx), function(k) {
        fold_df <- train_df[-fold_idx[[k]], ]
        X_fold  <- fold_df[, avail_prots, drop = FALSE]
        y_fold  <- fold_df$outcome
        b <- tryCatch(
            Boruta(x = X_fold, y = y_fold, doTrace = 0, maxRuns = 100),
            error = function(e) { message("  Boruta fold ", k, " failed: ", e$message); NULL }
        )
        if (is.null(b)) return(character(0))
        names(b$finalDecision[b$finalDecision == "Confirmed"])
    })
    feat_votes   <- table(unlist(boruta_confirmed))
    boruta_feats <- names(feat_votes[feat_votes >= 1])  # ≥1/5 folds threshold
    cat(sprintf("  Boruta confirmed %d features (≥1/5 folds)\n", length(boruta_feats)))
    saveRDS(boruta_feats, BORUTA_FILE)
}

if (length(boruta_feats) >= 2) {
    boruta_avail <- intersect(boruta_feats, names(train_df))
    # Always add demographics back regardless of Boruta result
    demo_avail   <- intersect(avail_demo, names(train_df))
    final_feats  <- unique(c(boruta_avail, demo_avail))
    cat(sprintf("  Boruta proteins: %d | + demographics: %s | total: %d\n",
                length(boruta_avail), paste(demo_avail, collapse=", "), length(final_feats)))
    train_df    <- train_df[, c("outcome", final_feats), drop = FALSE]
    test_df     <- test_df[,  c("outcome", final_feats), drop = FALSE]
    avail_prots <- final_feats
    fwrite(data.table(feature = boruta_avail, type = "protein"),
           file.path(OUT_DIR, glue("{COHORT}_preonset_boruta_features.csv")))
} else {
    # No Boruta: use all proteins + demographics
    demo_avail  <- intersect(avail_demo, names(train_df))
    final_feats <- unique(c(setdiff(names(train_df), "outcome"), demo_avail))
    cat(sprintf("  WARNING: Boruta confirmed %d — using all %d DEPs + demographics\n",
                length(boruta_feats), length(final_feats)))
}

# ── 6. caret ──────────────────────────────────────────────────────────────────
myControl <- trainControl(
    method = "repeatedcv", number = 5, repeats = 10,
    summaryFunction = twoClassSummary,
    classProbs = TRUE, savePredictions = "final", verboseIter = FALSE,
    sampling = "down"    # downsample controls to 1:1 per fold; works well with Boruta-reduced features
)

train_model <- function(method, label, extra_args = list()) {
    rds_path <- file.path(MODEL_DIR,
                          glue("{COHORT}_pre_boruta_{method}_model.rds"))
    if (file.exists(rds_path)) {
        cat(sprintf("  Loading cached: %s\n", label))
        model <- readRDS(rds_path)
        cat(sprintf("    Best ROC: %.3f\n", max(model$results$ROC, na.rm=TRUE)))
        return(model)
    }
    cat(sprintf("  Training: %s...\n", label))
    set.seed(2121)
    args_base <- list(form = outcome ~ ., data = train_df, method = method,
                      metric = "ROC", preProc = c("medianImpute","center","scale"), trControl = myControl)
    model <- tryCatch(do.call(train, c(args_base, extra_args)),
                      error = function(e) { message("  FAILED: ", e$message); NULL })
    if (!is.null(model)) {
        saveRDS(model, rds_path)
        cat(sprintf("    Best ROC: %.3f\n", max(model$results$ROC, na.rm=TRUE)))
    }
    model
}

glmnet_model <- train_model("glmnet", "Elastic Net")

model_list <- Filter(Negate(is.null), list(glmnet = glmnet_model))
cat(sprintf("\n%d models trained.\n", length(model_list)))

# ── 7. Model comparison ────────────────────────────────────────────────────────
resamp     <- resamples(model_list)
resamp_dt  <- as.data.table(resamp$values, keep.rownames = "Resample")
resamp_long <- melt(resamp_dt, id.vars = "Resample",
                    variable.name = "model_metric", value.name = "value")
resamp_long[, value := as.numeric(value)]
resamp_long[, c("model","metric") := tstrsplit(model_metric, "~", fixed=TRUE)]
resamp_long <- resamp_long[metric == "ROC"]

p_compare <- ggplot(resamp_long, aes(x = reorder(model, -value, median), y = value)) +
    geom_boxplot(fill = "#56B4E9", colour = "#1A1A1A", outlier.size = 0.8, width = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = 0.81, linetype = "dotted", colour = "#CC0066", linewidth = 0.5) +
    annotate("text", x = 0.6, y = 0.82, label = "Abdelhak et al. AUC 0.81",
             size = 2.5, colour = "#CC0066", hjust = 0) +
    labs(x = NULL, y = "Cross-validated AUC (ROC)",
         title = sprintf("Pre-onset %s: elastic net classifier", DISEASE),
         subtitle = "5-fold CV × 10 repeats | PSM-matched HC 1:5 | Boruta features") +
    coord_cartesian(ylim = c(0.4, 1.0)) +
    theme_ukb()

cmp_pdf <- glue("{COHORT}_preonset_ml_model_comparison.pdf")
ggsave(file.path(OUT_DIR, cmp_pdf),
       p_compare, width = 4.5, height = 4, device = cairo_pdf)
cat(sprintf("  Saved: %s\n", cmp_pdf))

# ── 8. ROC on test set ─────────────────────────────────────────────────────────
MODEL_COLS  <- c("#CC0066","#2B4C7E","#E6A817","#56B4E9")
roc_dt_list <- list()
test_results <- list()

for (i in seq_along(model_list)) {
    nm  <- names(model_list)[i]
    mod <- model_list[[i]]
    probs <- tryCatch(predict(mod, newdata = test_df, type = "prob")[, POS_LBL],
                      error = function(e) NULL)
    if (is.null(probs)) next
    roc_obj <- roc(test_df$outcome, probs, levels = c(NEG_LBL, POS_LBL),
                   direction = "<", quiet = TRUE)
    roc_dt_list[[nm]] <- data.table(model = nm, fpr = 1 - roc_obj$specificities,
                                    tpr = roc_obj$sensitivities, auc = as.numeric(auc(roc_obj)))
    test_results[[nm]] <- data.table(model = nm, AUC = round(as.numeric(auc(roc_obj)), 4))
    cat(sprintf("  %s  test AUC = %.3f\n", nm, as.numeric(auc(roc_obj))))
}

if (length(roc_dt_list) > 0) {
    roc_all <- rbindlist(roc_dt_list)
    roc_all[, model_label := sprintf("%s (AUC=%.3f)", model, auc)]
    p_roc <- ggplot(roc_all, aes(x = fpr, y = tpr, colour = model_label)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    colour = "grey60", linewidth = 0.3) +
        geom_line(linewidth = 0.8) +
        scale_colour_manual(values = setNames(MODEL_COLS[seq_along(model_list)],
                                              unique(roc_all$model_label)), name = NULL) +
        labs(x = "1 – Specificity (FPR)", y = "Sensitivity (TPR)",
             title = sprintf("Pre-onset %s: ROC curves (held-out test set)", DISEASE)) +
        coord_equal(xlim = c(0,1), ylim = c(0,1)) +
        theme_ukb()
    roc_pdf <- glue("{COHORT}_preonset_ml_roc_curves.pdf")
    ggsave(file.path(OUT_DIR, roc_pdf),
           p_roc, width = 5.5, height = 5, device = cairo_pdf)
    cat(sprintf("  Saved: %s\n", roc_pdf))
}

if (length(test_results) > 0) {
    results_dt <- rbindlist(test_results)
    fwrite(results_dt, file.path(OUT_DIR, glue("{COHORT}_preonset_ml_results.csv")))
}

# ── 9. Feature importance ──────────────────────────────────────────────────────
imp_list <- list()
for (nm in names(model_list)) {
    imp <- tryCatch(varImp(model_list[[nm]])$importance, error = function(e) NULL)
    if (is.null(imp)) next
    imp_dt  <- as.data.table(imp, keep.rownames = "protein")
    val_col <- intersect(c(POS_LBL, "Overall", names(imp_dt)[-1]), names(imp_dt))[1]
    imp_dt[, importance := get(val_col)][, model := nm]
    imp_list[[nm]] <- imp_dt[, .(protein, importance, model)]
}

if (length(imp_list) > 0) {
    imp_all <- rbindlist(imp_list)
    imp_avg <- imp_all[, .(mean_imp = mean(importance, na.rm=TRUE)), by = protein]
    imp_avg <- imp_avg[order(-mean_imp)][1:min(25, .N)]
    imp_avg[, protein_f := factor(protein, levels = rev(protein))]
    p_imp <- ggplot(imp_avg, aes(x = mean_imp, y = protein_f)) +
        geom_col(fill = "#CC0066", width = 0.7) +
        labs(x = "Mean variable importance", y = NULL,
             title = sprintf("Top features: pre-onset %s classifier", DISEASE)) +
        theme_ukb() +
        theme(axis.text.y = element_text(size = 7))
    imp_pdf <- glue("{COHORT}_preonset_ml_feature_importance.pdf")
    ggsave(file.path(OUT_DIR, imp_pdf),
           p_imp, width = 5, height = 5.5, device = cairo_pdf)
    cat(sprintf("  Saved: %s\n", imp_pdf))
    fwrite(imp_avg, file.path(OUT_DIR, glue("{COHORT}_preonset_ml_feature_importance.csv")))
}

cat("\n04_ms_preonset_ml.R complete.\n")
