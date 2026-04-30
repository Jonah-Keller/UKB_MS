#!/usr/bin/env Rscript
# 02_als_ml_caret.R
# Binary ML classifier: ALS vs healthy controls — Chia 2025 replication
#
# Approach (per Chia template):
#   - Features: 31 combined DEPs (FDR < 0.05) from ALS limma
#   - PSM-matched cohort (1:10, same as limma)
#   - Stratified 80/20 split
#   - 5-fold × 10 repeats, ROC metric
#   - Models: glmnet, gbm, rf, nnet
#   - ROC curves + variable importance
#
# Input:
#   data/ukb/olink/processed/als_olink_qc.csv
#   results/differential/als_combined_vs_hc.csv
#
# Output: results/ml/
#   als_ml_roc_curves.pdf
#   als_ml_model_comparison.pdf
#   als_ml_feature_importance.pdf
#   als_ml_results.csv

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(caret)
    library(pROC)
    library(ggplot2)
    library(ggrepel)
    library(MatchIt)
    library(patchwork)
})

args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "als_olink_qc.csv")
DIFF_FILE   <- file.path(PROJ_DIR, "results", "differential", "als_combined_vs_hc.csv")
CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
OUT_DIR     <- file.path(PROJ_DIR, "results", "ml")
MODEL_DIR   <- file.path(OUT_DIR, "models")
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR    <- 0.05
PSM_RATIO  <- 10
PSM_CAL    <- 0.2
TRAIN_PROP <- 0.80
set.seed(42)

# Load PSM covariates
load_covariates <- function(cadasil_dir) {
    misc <- file.path(cadasil_dir, "misc")
    bmi   <- fread(file.path(misc, "ukb_bmi_age.csv"), showProgress=FALSE)
    setnames(bmi, c("participant.eid","participant.p21001_i0"), c("eid","bmi"), skip_absent=TRUE)
    bmi <- bmi[, .(eid, bmi)]
    smoke <- fread(file.path(misc, "ukb_smoking_status.csv"), showProgress=FALSE)
    setnames(smoke, c("participant.eid","participant.p20116_i0"), c("eid","smoking_code"), skip_absent=TRUE)
    smoke[, ever_smoker := as.integer(smoking_code %in% c(1L,2L))]
    smoke <- smoke[smoking_code != -3L, .(eid, ever_smoker)]
    wmh <- fread(file.path(misc, "wmh_htn_dm.csv"), showProgress=FALSE)
    wmh[, diabetes := as.integer(p2443_i0=="Yes")]
    wmh <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes)]
    alc <- fread(file.path(misc, "ukb_alcohol_intake.csv"), showProgress=FALSE)
    setnames(alc, c("participant.eid","participant.p1558_i0"), c("eid","alcohol_freq"), skip_absent=TRUE)
    alc <- alc[alcohol_freq > 0L, .(eid, alcohol_freq)]
    Reduce(function(x,y) merge(x,y,by="eid",all.x=TRUE), list(bmi,smoke,wmh,alc))
}

# ── 1. Load DEPs ──────────────────────────────────────────────────────────────
cat("Loading ALS DEP list...\n")
deps <- fread(DIFF_FILE)
dep_proteins <- deps[adj.P.Val < FDR_THR, tolower(protein)]
cat(sprintf("  %d DEPs as features\n", length(dep_proteins)))

# ── 2. Load QC'd Olink ───────────────────────────────────────────────────────
cat("Loading QC'd ALS Olink data...\n")
dt <- fread(QC_FILE, showProgress=FALSE)
covs <- load_covariates(CADASIL_DIR)
dt <- merge(dt, covs, by="eid", all.x=TRUE)
dt[, sex_num := as.integer(factor(sex)) - 1L]

dt_filt <- dt[qc_outlier == FALSE &
              !is.na(age_at_sampling) & !is.na(sex) &
              !is.na(UMAP1) & !is.na(UMAP2)]
cat(sprintf("  After QC filter: %d participants\n", nrow(dt_filt)))
cat(sprintf("  Status: %s\n",
            paste(names(table(dt_filt$als_status)),
                  as.integer(table(dt_filt$als_status)), sep="=", collapse="  ")))

# ── 3. PSM ────────────────────────────────────────────────────────────────────
cat("\nPSM matching (combined ALS vs HC)...\n")
PSM_VARS <- c("age_at_sampling","sex_num","bmi","ever_smoker","diabetes","alcohol_freq")
dt_filt[, is_case := as.integer(als_status %in% c("pre_onset","post_onset"))]
dt_c <- dt_filt[complete.cases(as.data.frame(dt_filt[, c("is_case", ..PSM_VARS)]))]
cat(sprintf("  PSM input: %d cases, %d controls (%d excluded)\n",
            sum(dt_c$is_case), sum(!dt_c$is_case), nrow(dt_filt)-nrow(dt_c)))

psm_df <- as.data.frame(dt_c[, c("eid","is_case", ..PSM_VARS)])
m <- matchit(is_case ~ age_at_sampling + sex_num + bmi + ever_smoker + diabetes + alcohol_freq,
             data=psm_df, method="nearest", ratio=PSM_RATIO, caliper=PSM_CAL, std.caliper=TRUE)
matched_eids <- as.integer(match.data(m)$eid)
dt_ml <- dt_c[eid %in% matched_eids]
cat(sprintf("  After PSM: %d cases, %d controls\n",
            sum(dt_ml$is_case), sum(!dt_ml$is_case)))

# ── 4. Feature matrix ─────────────────────────────────────────────────────────
cat("\nBuilding feature matrix...\n")
avail_prots <- intersect(dep_proteins, names(dt_ml))
cat(sprintf("  Using %d / %d DEPs\n", length(avail_prots), length(dep_proteins)))

dt_ml[, outcome := factor(fifelse(is_case==1L, "ALS", "Control"), levels=c("ALS","Control"))]
feat_cols <- c("outcome", avail_prots)
ml_dt <- dt_ml[, ..feat_cols]
ml_dt <- ml_dt[complete.cases(ml_dt)]
cat(sprintf("  Complete cases: %d  (ALS=%d, Control=%d)\n",
            nrow(ml_dt), sum(ml_dt$outcome=="ALS"), sum(ml_dt$outcome=="Control")))

# ── 5. 80/20 split ────────────────────────────────────────────────────────────
ml_df <- as.data.frame(ml_dt)
train_idx <- createDataPartition(ml_df$outcome, p=TRAIN_PROP, list=FALSE)
train_df  <- ml_df[ train_idx, ]
test_df   <- ml_df[-train_idx, ]
cat(sprintf("  Train: %d  (ALS=%d, Control=%d)\n",
            nrow(train_df), sum(train_df$outcome=="ALS"), sum(train_df$outcome=="Control")))
cat(sprintf("  Test : %d  (ALS=%d, Control=%d)\n",
            nrow(test_df),  sum(test_df$outcome=="ALS"),  sum(test_df$outcome=="Control")))

# ── 6. caret setup ────────────────────────────────────────────────────────────
myControl <- trainControl(
    method="repeatedcv", number=5, repeats=10,
    summaryFunction=twoClassSummary,
    classProbs=TRUE, savePredictions="final", verboseIter=FALSE
)

train_model <- function(method, label, extra_args=list()) {
    rds_path <- file.path(MODEL_DIR, paste0("als_", method, "_model.rds"))
    if (file.exists(rds_path)) {
        cat(sprintf("  Loading cached: %s\n", label))
        model <- readRDS(rds_path)
        cat(sprintf("    Best ROC: %.3f\n", max(model$results$ROC, na.rm=TRUE)))
        return(model)
    }
    cat(sprintf("  Training: %s...\n", label))
    set.seed(2121)
    args_base <- list(form=outcome~., data=train_df, method=method,
                      metric="ROC", preProc=c("center","scale"), trControl=myControl)
    model <- tryCatch(do.call(train, c(args_base, extra_args)),
                      error=function(e) { message("  FAILED: ", e$message); NULL })
    if (!is.null(model)) {
        saveRDS(model, rds_path)
        cat(sprintf("    Best ROC: %.3f\n", max(model$results$ROC, na.rm=TRUE)))
    }
    model
}

glmnet_model <- train_model("glmnet", "Elastic Net GLM")
gbm_model    <- train_model("gbm",    "Gradient Boosting")
rf_model     <- train_model("rf",     "Random Forest",
                             list(importance=TRUE, ntree=100,
                                  tuneGrid=expand.grid(mtry=c(5,10,15))))
nnet_model   <- train_model("nnet",   "Neural Network",
                             list(tuneGrid=expand.grid(size=c(5,10,15), decay=c(0.01,0.1))))

model_list <- Filter(Negate(is.null),
                     list(glmnet=glmnet_model, gbm=gbm_model, rf=rf_model, nnet=nnet_model))
cat(sprintf("\n%d models trained.\n", length(model_list)))

# ── 7. Model comparison ───────────────────────────────────────────────────────
resamp <- resamples(model_list)
resamp_dt <- as.data.table(resamp$values, keep.rownames="Resample")
resamp_long <- melt(resamp_dt, id.vars="Resample", variable.name="model_metric", value.name="value")
resamp_long[, value := as.numeric(value)]
resamp_long[, c("model","metric") := tstrsplit(model_metric, "~", fixed=TRUE)]
resamp_long <- resamp_long[metric=="ROC"]

p_compare <- ggplot(resamp_long, aes(x=reorder(model,-value,median), y=value)) +
    geom_boxplot(fill="#2B4C7E", colour="#1A1A1A", outlier.size=0.8, width=0.5) +
    geom_hline(yintercept=0.5, linetype="dashed", colour="grey50", linewidth=0.3) +
    labs(x=NULL, y="Cross-validated AUC (ROC)",
         title="ALS classifier: 5-fold CV × 10 repeats (train set)") +
    coord_cartesian(ylim=c(0.4,1.0)) +
    theme_ukb()

ggsave(file.path(OUT_DIR, "als_ml_model_comparison.pdf"),
       p_compare, width=4.5, height=4, device=cairo_pdf)
cat("  Saved: als_ml_model_comparison.pdf\n")

# ── 8. ROC on test set ────────────────────────────────────────────────────────
MODEL_COLS <- c("#CC0066","#2B4C7E","#E6A817","#56B4E9")
roc_dt_list <- list()
test_results <- list()

for (i in seq_along(model_list)) {
    nm  <- names(model_list)[i]
    mod <- model_list[[i]]
    probs <- tryCatch(predict(mod, newdata=test_df, type="prob")[,"ALS"], error=function(e) NULL)
    if (is.null(probs)) next
    roc_obj <- roc(test_df$outcome, probs, levels=c("Control","ALS"), direction="<", quiet=TRUE)
    roc_dt_list[[nm]] <- data.table(model=nm, fpr=1-roc_obj$specificities,
                                     tpr=roc_obj$sensitivities, auc=as.numeric(auc(roc_obj)))
    test_results[[nm]] <- data.table(model=nm, AUC=round(as.numeric(auc(roc_obj)),4))
    cat(sprintf("  %s  test AUC = %.3f\n", nm, as.numeric(auc(roc_obj))))
}

if (length(roc_dt_list) > 0) {
    roc_all <- rbindlist(roc_dt_list)
    roc_all[, model_label := sprintf("%s (AUC=%.3f)", model, auc)]
    p_roc <- ggplot(roc_all, aes(x=fpr, y=tpr, colour=model_label)) +
        geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60", linewidth=0.3) +
        geom_line(linewidth=0.8) +
        scale_colour_manual(values=setNames(MODEL_COLS[seq_along(model_list)],
                                            unique(roc_all$model_label)), name=NULL) +
        labs(x="1 – Specificity (FPR)", y="Sensitivity (TPR)",
             title="ALS ROC curves — held-out test set") +
        coord_equal(xlim=c(0,1), ylim=c(0,1)) +
        theme_ukb()
    ggsave(file.path(OUT_DIR, "als_ml_roc_curves.pdf"),
           p_roc, width=5.5, height=5, device=cairo_pdf)
    cat("  Saved: als_ml_roc_curves.pdf\n")
}

if (length(test_results) > 0) {
    results_dt <- rbindlist(test_results)
    fwrite(results_dt, file.path(OUT_DIR, "als_ml_results.csv"))
}

# ── 9. Feature importance ─────────────────────────────────────────────────────
imp_list <- list()
for (nm in names(model_list)) {
    imp <- tryCatch(varImp(model_list[[nm]])$importance, error=function(e) NULL)
    if (is.null(imp)) next
    imp_dt <- as.data.table(imp, keep.rownames="protein")
    val_col <- intersect(c("ALS","Overall",names(imp_dt)[-1]), names(imp_dt))[1]
    imp_dt[, importance:=get(val_col)][, model:=nm]
    imp_list[[nm]] <- imp_dt[, .(protein, importance, model)]
}
if (length(imp_list) > 0) {
    imp_all <- rbindlist(imp_list)
    imp_avg <- imp_all[, .(mean_imp=mean(importance,na.rm=TRUE)), by=protein]
    imp_avg <- imp_avg[order(-mean_imp)][1:min(30,.N)]
    imp_avg[, protein_f := factor(protein, levels=rev(protein))]
    p_imp <- ggplot(imp_avg, aes(x=mean_imp, y=protein_f)) +
        geom_col(fill="#2B4C7E", width=0.7) +
        labs(x="Mean variable importance", y=NULL,
             title="Top features: ALS classifier") +
        theme_ukb() +
        theme(axis.text.y=element_text(size=7))
    ggsave(file.path(OUT_DIR, "als_ml_feature_importance.pdf"),
           p_imp, width=5, height=5.5, device=cairo_pdf)
    cat("  Saved: als_ml_feature_importance.pdf\n")
    fwrite(imp_avg, file.path(OUT_DIR, "als_ml_feature_importance.csv"))
}

cat("\n02_als_ml_caret.R complete.\n")
