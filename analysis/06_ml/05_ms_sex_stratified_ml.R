#!/usr/bin/env Rscript
# 05_ms_sex_stratified_ml.R
# Sex-stratified MS classifier: 6 glmnet models (female/male × combined/pre-onset/post-onset)
#
# Approach:
#   - Per-sex PSM (no sex_num — already stratified by sex)
#   - glmnet elastic net, 5-fold × 3 repeats CV, ROC metric
#   - Linear SHAP for combined models: mean |SHAP| = |coef| × sd(feature)
#
# Output: results/ml/
#   ms_sex_ml_roc_data.csv        (fpr/tpr/auc for all 6 models)
#   ms_sex_ml_shap.csv            (protein, female_shap, male_shap — combined models)
#   ms_sex_ml_roc_6curves.pdf
#   ms_sex_ml_shap_paired.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(caret)
    library(pROC)
    library(ggplot2)
    library(ggrepel)
    library(MatchIt)
    library(glmnet)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
SEX_DIR     <- file.path(PROJ_DIR, "results", "sex_stratified")
CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
OUT_DIR     <- file.path(PROJ_DIR, "results", "ml")
MODEL_DIR   <- file.path(OUT_DIR, "models_sex")
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

N_FEATURES <- 50
PSM_RATIO  <- 10
PSM_CAL    <- 0.2
TRAIN_PROP <- 0.80
set.seed(42)

# 6 model configurations
CONFIGS <- list(
    female_combined  = list(sex=0L, sex_lbl="Female", stage="Combined",   stage_tag="combined",
                            statuses=c("pre_onset","post_onset")),
    female_preonset  = list(sex=0L, sex_lbl="Female", stage="Pre-onset",  stage_tag="preonset",
                            statuses="pre_onset"),
    female_postonset = list(sex=0L, sex_lbl="Female", stage="Post-onset", stage_tag="postonset",
                            statuses="post_onset"),
    male_combined    = list(sex=1L, sex_lbl="Male",   stage="Combined",   stage_tag="combined",
                            statuses=c("pre_onset","post_onset")),
    male_preonset    = list(sex=1L, sex_lbl="Male",   stage="Pre-onset",  stage_tag="preonset",
                            statuses="pre_onset"),
    male_postonset   = list(sex=1L, sex_lbl="Male",   stage="Post-onset", stage_tag="postonset",
                            statuses="post_onset")
)

# Visual encoding: sex = colour family, stage = linetype
ROC_COLS <- c(
    "Female Combined"   = "#CC0066",
    "Female Post-onset" = "#E8779D",
    "Female Pre-onset"  = "#F2B8CE",
    "Male Combined"     = "#2B4C7E",
    "Male Post-onset"   = "#5B7FA6",
    "Male Pre-onset"    = "#A8BFCF"
)
ROC_LTY <- c(
    "Female Combined"   = "solid",
    "Female Post-onset" = "dashed",
    "Female Pre-onset"  = "dotted",
    "Male Combined"     = "solid",
    "Male Post-onset"   = "dashed",
    "Male Pre-onset"    = "dotted"
)

# ── 1. Load QC data ────────────────────────────────────────────────────────────
load_covariates <- function(cadasil_dir) {
    misc  <- file.path(cadasil_dir, "misc")
    bmi   <- fread(file.path(misc, "ukb_bmi_age.csv"), showProgress = FALSE)
    setnames(bmi, c("participant.eid","participant.p21001_i0"), c("eid","bmi"), skip_absent = TRUE)
    bmi   <- bmi[, .(eid, bmi)]
    smoke <- fread(file.path(misc, "ukb_smoking_status.csv"), showProgress = FALSE)
    setnames(smoke, c("participant.eid","participant.p20116_i0"), c("eid","smoking_code"), skip_absent = TRUE)
    smoke[, ever_smoker := as.integer(smoking_code %in% c(1L, 2L))]
    smoke <- smoke[smoking_code != -3L, .(eid, ever_smoker)]
    wmh   <- fread(file.path(misc, "wmh_htn_dm.csv"), showProgress = FALSE)
    wmh[, diabetes := as.integer(p2443_i0 == "Yes")]
    wmh   <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes)]
    alc   <- fread(file.path(misc, "ukb_alcohol_intake.csv"), showProgress = FALSE)
    setnames(alc, c("participant.eid","participant.p1558_i0"), c("eid","alcohol_freq"), skip_absent = TRUE)
    alc   <- alc[alcohol_freq > 0L, .(eid, alcohol_freq)]
    Reduce(function(x, y) merge(x, y, by = "eid", all.x = TRUE),
           list(bmi, smoke, wmh, alc))
}

cat("Loading data...\n")
dt   <- fread(QC_FILE, showProgress = FALSE)
covs <- load_covariates(CADASIL_DIR)
dt   <- merge(dt, covs, by = "eid", all.x = TRUE)
dt[, sex_num := as.integer(factor(sex)) - 1L]

META_COLS <- c("eid","ms_status","age_at_sampling","age_at_diagnosis","years_to_diagnosis",
               "sex","olink_instance","qc_outlier","UMAP1","UMAP2","mean_npx",
               "bmi","ever_smoker","diabetes","alcohol_freq","sex_num")
protein_cols <- setdiff(names(dt), META_COLS)

dt_filt <- dt[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
              !is.na(UMAP1) & !is.na(UMAP2)]
cat(sprintf("  After QC: %d (Female=%d, Male=%d)\n",
            nrow(dt_filt), sum(dt_filt$sex == 0), sum(dt_filt$sex == 1)))

PSM_VARS <- c("age_at_sampling", "bmi", "ever_smoker", "diabetes", "alcohol_freq")

# ── 2. Train 6 models ──────────────────────────────────────────────────────────
run_config <- function(cfg) {
    tag <- sprintf("%s_%s", tolower(cfg$sex_lbl), cfg$stage_tag)
    cat(sprintf("\n=== %s %s ===\n", cfg$sex_lbl, cfg$stage))

    dep_file <- file.path(SEX_DIR,
        fifelse(cfg$sex == 0L, "ms_female_vs_hc_all.csv", "ms_male_vs_hc_all.csv"))
    deps     <- fread(dep_file)
    features <- deps[order(adj.P.Val)][seq_len(min(N_FEATURES, .N)), tolower(protein)]

    dt_sx <- dt_filt[sex == cfg$sex & ms_status %in% c(cfg$statuses, "control")]
    dt_sx[, is_case := as.integer(ms_status %in% cfg$statuses)]
    dt_c  <- dt_sx[complete.cases(as.data.frame(dt_sx[, c("is_case", ..PSM_VARS)]))]

    n_case <- sum(dt_c$is_case)
    n_ctrl <- sum(!dt_c$is_case)
    cat(sprintf("  PSM input: %d cases, %d controls\n", n_case, n_ctrl))
    if (n_case < 20) { cat("  Too few cases — skipping\n"); return(NULL) }

    psm_df <- as.data.frame(dt_c[, c("eid", "is_case", ..PSM_VARS)])
    m <- tryCatch(
        matchit(is_case ~ age_at_sampling + bmi + ever_smoker + diabetes + alcohol_freq,
                data = psm_df, method = "nearest", ratio = PSM_RATIO,
                caliper = PSM_CAL, std.caliper = TRUE),
        error = function(e) { cat("  PSM failed:", e$message, "\n"); NULL }
    )
    if (is.null(m)) return(NULL)

    dt_ml <- dt_c[eid %in% as.integer(match.data(m)$eid)]
    cat(sprintf("  After PSM: %d cases, %d controls\n",
                sum(dt_ml$is_case), sum(!dt_ml$is_case)))

    avail <- intersect(features, names(dt_ml))
    dt_ml[, outcome := factor(fifelse(is_case == 1L, "MS", "Control"),
                               levels = c("MS", "Control"))]
    ml_df <- as.data.frame(dt_ml[, c("outcome", avail), with = FALSE])
    ml_df <- ml_df[complete.cases(ml_df), ]
    cat(sprintf("  Complete cases: %d  (MS=%d, Control=%d)\n",
                nrow(ml_df), sum(ml_df$outcome=="MS"), sum(ml_df$outcome=="Control")))

    train_idx <- createDataPartition(ml_df$outcome, p = TRAIN_PROP, list = FALSE)
    train_df  <- ml_df[ train_idx, ]
    test_df   <- ml_df[-train_idx, ]

    # Fewer CV repeats for small pre-onset cohorts
    n_cv_rep <- if (n_case < 100) 3L else 5L
    ctrl <- trainControl(
        method = "repeatedcv", number = 5, repeats = n_cv_rep,
        summaryFunction = twoClassSummary, classProbs = TRUE,
        savePredictions = "final", verboseIter = FALSE
    )

    rds <- file.path(MODEL_DIR, sprintf("glmnet_%s.rds", tag))
    if (file.exists(rds)) {
        cat("  Loading cached model\n")
        model <- readRDS(rds)
    } else {
        cat(sprintf("  Training glmnet (%d-fold × %d)...\n", 5L, n_cv_rep))
        set.seed(2121)
        model <- tryCatch(
            train(outcome ~ ., data = train_df, method = "glmnet",
                  metric = "ROC", preProc = c("center", "scale"), trControl = ctrl),
            error = function(e) { cat("  Training failed:", e$message, "\n"); NULL }
        )
        if (!is.null(model)) saveRDS(model, rds)
    }
    if (is.null(model)) return(NULL)
    cat(sprintf("  Best CV ROC: %.3f\n", max(model$results$ROC, na.rm = TRUE)))

    # Test-set ROC
    probs <- tryCatch(
        predict(model, newdata = test_df, type = "prob")[, "MS"],
        error = function(e) NULL
    )
    if (is.null(probs)) return(NULL)
    roc_obj  <- roc(test_df$outcome, probs, levels = c("Control","MS"),
                    direction = "<", quiet = TRUE)
    test_auc <- as.numeric(auc(roc_obj))
    cat(sprintf("  Test AUC: %.3f\n", test_auc))

    # ROC curve points
    roc_dt <- data.table(
        config    = tag,
        sex_lbl   = cfg$sex_lbl,
        stage     = cfg$stage,
        curve_key = paste(cfg$sex_lbl, cfg$stage),
        fpr       = 1 - roc_obj$specificities,
        tpr       = roc_obj$sensitivities,
        auc       = test_auc
    )

    # Linear SHAP for combined models: |coef| * sd(feature in training set)
    shap_dt <- NULL
    if (cfg$stage_tag == "combined") {
        coef_mat   <- coef(model$finalModel, s = model$bestTune$lambda)
        coef_vec   <- as.numeric(coef_mat[-1, 1])
        feat_names <- rownames(coef_mat)[-1]
        feat_sd    <- apply(train_df[, feat_names, drop = FALSE], 2, sd, na.rm = TRUE)
        shap_vals  <- abs(coef_vec) * feat_sd
        shap_dt    <- data.table(
            protein  = feat_names,
            shap     = shap_vals,
            sex_lbl  = cfg$sex_lbl
        )
        shap_dt <- shap_dt[shap > 0][order(-shap)]
    }

    list(tag = tag, config = cfg, model = model, roc_dt = roc_dt, shap_dt = shap_dt,
         test_auc = test_auc, train_df = train_df, test_df = test_df)
}

results <- lapply(CONFIGS, run_config)
results <- Filter(Negate(is.null), results)
cat(sprintf("\n%d / 6 models trained successfully\n", length(results)))

# ── 3. Save ROC data ───────────────────────────────────────────────────────────
roc_all <- rbindlist(lapply(results, `[[`, "roc_dt"), fill = TRUE)
fwrite(roc_all, file.path(OUT_DIR, "ms_sex_ml_roc_data.csv"))
cat("  Saved: ms_sex_ml_roc_data.csv\n")

# AUC summary
auc_summary <- unique(roc_all[, .(curve_key, auc)])[order(curve_key)]
cat("\nTest AUC summary:\n"); print(auc_summary, row.names = FALSE)

# ── 4. Save SHAP data (combined models only) ──────────────────────────────────
shap_list <- Filter(Negate(is.null), lapply(results, `[[`, "shap_dt"))
if (length(shap_list) > 0) {
    shap_long <- rbindlist(shap_list)
    shap_wide <- dcast(shap_long, protein ~ sex_lbl, value.var = "shap", fill = 0)
    setnames(shap_wide,
             intersect(c("Female","Male"), names(shap_wide)),
             intersect(c("female_shap","male_shap"),
                       paste0(tolower(intersect(c("Female","Male"), names(shap_wide))), "_shap")))
    # Handle single-sex case gracefully
    if (!"female_shap" %in% names(shap_wide)) shap_wide[, female_shap := 0]
    if (!"male_shap"   %in% names(shap_wide)) shap_wide[, male_shap   := 0]
    shap_wide[, max_shap := pmax(female_shap, male_shap)]
    shap_wide <- shap_wide[order(-max_shap)]
    fwrite(shap_wide[, .(protein, female_shap, male_shap)],
           file.path(OUT_DIR, "ms_sex_ml_shap.csv"))
    cat("  Saved: ms_sex_ml_shap.csv\n")
}

# ── 5. 6-curve ROC plot ────────────────────────────────────────────────────────
cat("\nPlotting 6-curve ROC...\n")

auc_dt <- unique(roc_all[, .(curve_key, auc)])
roc_all[, curve_label := sprintf("%s (AUC=%.3f)", curve_key,
                                  auc_dt[curve_key == .BY$curve_key, auc]),
        by = curve_key]

# Ensure consistent level ordering: Female first, then Male; within each: Combined, Post, Pre
lvl_order <- c("Female Combined","Female Post-onset","Female Pre-onset",
               "Male Combined","Male Post-onset","Male Pre-onset")
roc_all[, curve_key := factor(curve_key, levels = lvl_order)]
roc_all[, curve_label := factor(
    sprintf("%s (AUC=%.3f)", as.character(curve_key), auc),
    levels = sprintf("%s (AUC=%.3f)", lvl_order,
                     sapply(lvl_order, function(k) {
                         v <- unique(roc_all[curve_key == k, auc])
                         if (length(v) == 0) NA else v[1]
                     }))
)]
roc_all <- roc_all[!is.na(curve_key)]

# Colours and linetypes keyed to curve_label (which includes AUC)
lab_to_key <- setNames(as.character(unique(roc_all$curve_key)),
                        as.character(unique(roc_all$curve_label)))
lbl_cols <- setNames(
    ROC_COLS[match(lab_to_key[levels(roc_all$curve_label)], names(ROC_COLS))],
    levels(roc_all$curve_label)
)
lbl_lty <- setNames(
    ROC_LTY[match(lab_to_key[levels(roc_all$curve_label)], names(ROC_LTY))],
    levels(roc_all$curve_label)
)

p_roc6 <- ggplot(roc_all, aes(x = fpr, y = tpr,
                                colour = curve_label, linetype = curve_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid",
                colour = "grey80", linewidth = 0.3) +
    geom_line(linewidth = 0.75) +
    scale_colour_manual(values = lbl_cols, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = lbl_lty, name = NULL, drop = FALSE) +
    guides(colour   = guide_legend(keywidth = 1.6, keyheight = 0.75, ncol = 1),
           linetype = guide_legend(keywidth = 1.6, keyheight = 0.75, ncol = 1)) +
    labs(
        title    = "j  Sex-stratified MS classifier ROC curves",
        subtitle = "glmnet elastic net | held-out test set | solid=combined, dashed=post-onset, dotted=pre-onset",
        x = "1 \u2013 Specificity (FPR)",
        y = "Sensitivity (TPR)"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_ukb(base_size = 9) +
    theme(legend.position  = c(0.67, 0.24),
          legend.text      = element_text(size = 6.5),
          legend.key.size  = unit(0.85, "lines"),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.subtitle    = element_text(size = 6.5, colour = "grey40"))

ggsave(file.path(OUT_DIR, "ms_sex_ml_roc_6curves.pdf"),
       p_roc6, width = 5.5, height = 5.5, device = cairo_pdf)
cat("  Saved: ms_sex_ml_roc_6curves.pdf\n")

# ── 6. SHAP paired bar chart ───────────────────────────────────────────────────
cat("\nPlotting SHAP paired bar chart...\n")
shap_file <- file.path(OUT_DIR, "ms_sex_ml_shap.csv")
if (file.exists(shap_file)) {
    shap <- fread(shap_file)

    # Load dep_cat annotation
    fa_all <- fread(file.path(SEX_DIR, "ms_female_vs_hc_all.csv"))
    ma_all <- fread(file.path(SEX_DIR, "ms_male_vs_hc_all.csv"))
    wide_ann <- merge(
        fa_all[, .(protein, logFC_f=logFC, fdr_f=adj.P.Val)],
        ma_all[, .(protein, logFC_m=logFC, fdr_m=adj.P.Val)],
        by = "protein"
    )
    wide_ann[, dep_cat := fcase(
        fdr_f<0.05 & fdr_m<0.05 & sign(logFC_f)==sign(logFC_m), "Both concordant",
        fdr_f<0.05 & fdr_m<0.05 & sign(logFC_f)!=sign(logFC_m), "Both discordant",
        fdr_f<0.05 & fdr_m>=0.05,                                "Female-specific",
        fdr_f>=0.05 & fdr_m<0.05,                                "Male-specific",
        default                                                  = "NS"
    )]
    shap <- merge(shap, wide_ann[, .(protein, dep_cat)], by="protein", all.x=TRUE)
    shap[is.na(dep_cat), dep_cat := "NS"]

    shap[, max_shap := pmax(female_shap, male_shap)]
    top_shap <- shap[order(-max_shap)][seq_len(min(35, .N))]
    top_shap[, prot_f := factor(toupper(protein), levels=rev(toupper(protein)))]

    shap_long <- melt(top_shap, id.vars=c("protein","prot_f","dep_cat","max_shap"),
                      measure.vars=c("female_shap","male_shap"),
                      variable.name="sex", value.name="shap_val")
    shap_long[, sex_lbl := fifelse(sex=="female_shap","Female","Male")]
    shap_long[, sex_lbl := factor(sex_lbl, levels=c("Female","Male"))]

    SEX_FILL <- c("Female"="#CC0066","Male"="#2B4C7E")

    p_shap <- ggplot(shap_long, aes(x=shap_val, y=prot_f, fill=sex_lbl)) +
        geom_col(position=position_dodge(width=0.7), width=0.65, alpha=0.9) +
        scale_fill_manual(values=SEX_FILL, name=NULL) +
        scale_x_continuous(expand=expansion(mult=c(0,0.06))) +
        labs(
            title    = "e  Sex-stratified MS classifier: SHAP feature importance",
            subtitle = "Linear SHAP = |coef| \u00d7 sd(feature) | combined models only | top 35",
            x = "Mean |SHAP| (log-odds scale)", y = NULL
        ) +
        theme_ukb(base_size=9) +
        theme(axis.text.y     = element_text(size=6.5),
              legend.position = "top",
              plot.subtitle   = element_text(size=7, colour="grey40"))

    ggsave(file.path(OUT_DIR, "ms_sex_ml_shap_paired.pdf"),
           p_shap, width=5.5, height=8, device=cairo_pdf)
    cat("  Saved: ms_sex_ml_shap_paired.pdf\n")
}

cat("\n05_ms_sex_stratified_ml.R complete.\n")
