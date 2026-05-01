#!/usr/bin/env Rscript
# 11_ms_cluster_prems_ml.R — Endophenotype-stratified pre-MS prediction models
#
# Compares subtype-specific pre-MS classification performance across C0/C1/C2.
#
# Approach:
#   - Pooled "All pre-onset" baseline: loaded directly from Figure 1 results
#     ({cohort_short}_preonset_ml_results.csv). No re-fitting — ensures exact consistency.
#   - Cluster-specific models: same Boruta feature set as Figure 1
#     ({cohort_short}_pre_boruta_features.rds), glmnet elastic net.
#     PSM: 1:5 NN on age+sex, caliper 0.2 SD.
#     CV: proper LOOCV with savePredictions="final" (n<30) or 5×10 repeat CV
#         with held-out 20% test set (n>=30). Train/test split done BEFORE fitting.
#   - Power note: C0~11, C1~9, C2~21 — all use LOOCV. Treat as exploratory.
#
# Outputs (results/endophenotype/cluster_prems_ml/):
#   cluster_prems_auc_summary.csv   — per-cluster AUC + CI
#   cluster_prems_features_<c>.csv  — glmnet coefficients per cluster
#   panel_v.pdf                     — AUC comparison bar chart
#   panel_v2_roc.pdf                — ROC overlay

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
    library(glue)
    library(caret)
    library(pROC)
    library(MatchIt)
    library(glmnet)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))

cfg <- load_disease_config()

QC_FILE       <- here::here("data", "ukb", "olink", "processed",
                            glue("{cfg$cohort_short}_olink_qc.csv"))
CLUST_FILE    <- here::here("results", "endophenotype",
                            glue("{cfg$cohort_short}_prems_cluster_assignments.csv"))
FIG1_AUC_FILE <- here::here("results", "ml",
                            glue("{cfg$cohort_short}_preonset_ml_results.csv"))
BORUTA_FILE   <- here::here("results", "ml", "models",
                            glue("{cfg$cohort_short}_pre_boruta_features.rds"))
OUT_DIR       <- here::here("results", "endophenotype",
                            glue("{cfg$cohort_short}_cluster_prems_ml"))
FIG_DIR       <- here::here("results", "figures", "5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
PSM_RATIO  <- 5L
PSM_CAL    <- 0.2
MIN_CASES  <- 8L          # minimum pre-onset cases to attempt ML
DEMO_FEATS <- c("age_at_sampling", "sex_num")
CLUST_COLS["All pre-onset"] <- "grey45"

status_col   <- cfg$cohort_status_col
pre_state    <- cfg$status_values$pre_onset
ctrl_state   <- cfg$status_values$control
case_label   <- cfg$disease_short_caps

# ── 1. Load cohort ────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc    <- fread(QC_FILE,    showProgress = FALSE)
clust <- fread(CLUST_FILE)

meta_cols <- c("eid", status_col, "age_at_sampling","sex","age_at_diagnosis",
               "years_to_diagnosis","qc_outlier","plate_id","well_position",
               "hla_carrier","prs_score","comorbidity_cluster")
prot_cols <- setdiff(names(qc), meta_cols)
prot_cols <- prot_cols[vapply(prot_cols, function(x) is.numeric(qc[[x]]), logical(1))]

qc[, sex_num := as.integer(sex == "Female")]

dt <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
         get(status_col) %in% c(pre_state, ctrl_state)]
dt <- merge(dt, clust[, .(eid, cluster)], by = "eid", all.x = TRUE)
.unique_clusters <- sort(unique(dt$cluster[!is.na(dt$cluster) &
                                            dt[[status_col]] == pre_state]))
.cluster_labels  <- paste0("C", .unique_clusters)
dt[get(status_col) == pre_state, cluster_f := factor(
    ifelse(is.na(cluster), "None", paste0("C", cluster)),
    levels = c("None", .cluster_labels)
)]

n_pre <- dt[get(status_col) == pre_state, .N, by = cluster_f]
cat(glue("  Pre-onset {case_label} by cluster:\n")); print(n_pre)
cat(sprintf("  Total controls: %d\n", sum(dt[[status_col]] == ctrl_state)))

# ── 2. Load Figure 1 Boruta features (ensures consistency) ───────────────────
cat("\nLoading Figure 1 Boruta feature set...\n")
if (!file.exists(BORUTA_FILE)) stop(glue("Boruta features not found at {BORUTA_FILE} — run upstream pre-onset ML script first"))
boruta_feats <- readRDS(BORUTA_FILE)
cat(sprintf("  %d Boruta-confirmed proteins\n", length(boruta_feats)))

# Feature set: Boruta proteins + demographics, lowercased to match QC column names
feat_set <- unique(c(tolower(boruta_feats), DEMO_FEATS))
feat_set <- intersect(feat_set, names(dt))
cat(sprintf("  Feature set: %d proteins + %d demo = %d total\n",
            length(intersect(tolower(boruta_feats), names(dt))),
            length(intersect(DEMO_FEATS, names(dt))),
            length(feat_set)))

# ── 3. Pooled baseline from Figure 1 (no re-fitting) ─────────────────────────
cat(glue("\nLoading Figure 1 pooled pre-onset AUC ({basename(FIG1_AUC_FILE)})...\n"))
fig1_auc <- fread(FIG1_AUC_FILE)
cat("  Figure 1 results:\n"); print(fig1_auc)
pooled_auc    <- fig1_auc[model == "glmnet", AUC]
pooled_auc_lo <- NA_real_
pooled_auc_hi <- NA_real_
if (length(pooled_auc) == 0) {
    pooled_auc <- max(fig1_auc$AUC, na.rm = TRUE)
    cat("  glmnet not found — using best model AUC\n")
}
cat(sprintf("  Pooled baseline AUC (glmnet): %.3f\n", pooled_auc))

# ── 4. Cluster-specific models ────────────────────────────────────────────────
fit_cluster_model <- function(cluster_label, pre_eids, ctrl_dt) {
    n_cases <- length(pre_eids)
    cat(sprintf("\n  === %s: %d pre-onset cases ===\n", cluster_label, n_cases))

    if (n_cases < MIN_CASES) {
        cat(sprintf("  Skipping: %d cases < minimum %d\n", n_cases, MIN_CASES))
        return(NULL)
    }

    cases_dt <- dt[eid %in% pre_eids]
    cases_dt[, is_case := 1L]

    # PSM on age + sex (drop sex_num if zero-variance in this cluster)
    ctrl_sub <- ctrl_dt[, .(eid, is_case = 0L, age_at_sampling, sex_num)]
    sub_psm  <- rbind(
        cases_dt[, .(eid, is_case = 1L, age_at_sampling, sex_num)],
        ctrl_sub
    )
    sub_psm  <- sub_psm[complete.cases(sub_psm)]

    psm_vars <- c("age_at_sampling", "sex_num")
    # Remove zero-variance PSM vars (e.g. all-female cluster)
    psm_vars <- psm_vars[vapply(psm_vars, function(v)
        var(sub_psm[is_case == 1, get(v)], na.rm = TRUE) > 0 ||
        var(sub_psm[[v]], na.rm = TRUE) > 0,
        logical(1))]
    psm_form <- as.formula(paste("is_case ~", paste(psm_vars, collapse = " + ")))

    m <- tryCatch(
        matchit(psm_form, data = as.data.frame(sub_psm),
                method = "nearest", ratio = PSM_RATIO,
                caliper = PSM_CAL, std.caliper = TRUE),
        error = function(e) {
            cat(sprintf("  PSM (caliper) failed: %s — retrying without caliper\n",
                        conditionMessage(e)))
            tryCatch(
                matchit(psm_form, data = as.data.frame(sub_psm),
                        method = "nearest", ratio = PSM_RATIO),
                error = function(e2) { cat(sprintf("  PSM failed: %s\n", conditionMessage(e2))); NULL }
            )
        }
    )
    if (is.null(m)) return(NULL)

    matched_eids <- sub_psm$eid[m$weights > 0]
    ml_dt <- dt[eid %in% matched_eids]
    ml_dt[, is_case := as.integer(eid %in% pre_eids)]
    ml_dt[, outcome := factor(fifelse(is_case == 1L, case_label, "Control"),
                              levels = c(case_label, "Control"))]

    avail_feats <- intersect(feat_set, names(ml_dt))
    # Remove zero-variance features in this cluster's dataset
    nz_feats <- avail_feats[vapply(avail_feats, function(v) {
        vals <- ml_dt[[v]]
        !all(is.na(vals)) && var(vals, na.rm = TRUE) > 0
    }, logical(1))]
    if (length(nz_feats) < 2L) {
        cat("  Too few non-zero-variance features — skipping\n"); return(NULL)
    }

    for (col in nz_feats) {
        if (any(is.na(ml_dt[[col]]))) {
            med_val <- median(ml_dt[[col]], na.rm = TRUE)
            ml_dt[is.na(get(col)), (col) := med_val]
        }
    }
    ml_dt <- ml_dt[complete.cases(ml_dt[, c("outcome", ..nz_feats)])]
    n_cases  <- sum(ml_dt$outcome == case_label)
    cat(sprintf("  PSM dataset: %d %s, %d Control | features: %d\n",
                n_cases, case_label, sum(ml_dt$outcome == "Control"), length(nz_feats)))

    ml_df <- as.data.frame(ml_dt[, c("outcome", ..nz_feats)])

    # LOOCV for all cluster models (n always small)
    # savePredictions="final" is required to get per-fold predictions for AUC
    tc <- trainControl(
        method          = "LOOCV",
        classProbs      = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "final"
    )

    fit <- tryCatch(
        train(outcome ~ ., data = ml_df, method = "glmnet",
              trControl = tc, metric = "ROC",
              tuneGrid = expand.grid(alpha  = c(0, 0.5, 1),
                                     lambda = 10^seq(-4, 1, length.out = 20)),
              preProcess = c("center", "scale")),
        error = function(e) {
            cat(sprintf("  glmnet failed: %s\n", conditionMessage(e))); NULL
        }
    )
    if (is.null(fit)) return(NULL)

    # AUC from LOOCV held-out predictions
    loocv_preds <- fit$pred
    if (is.null(loocv_preds) || nrow(loocv_preds) == 0 || !(case_label %in% names(loocv_preds))) {
        cat("  LOOCV predictions unavailable\n"); return(NULL)
    }
    # Use best-tune predictions
    bt <- fit$bestTune
    loocv_best <- loocv_preds[loocv_preds$alpha == bt$alpha &
                               abs(loocv_preds$lambda - bt$lambda) < 1e-9, ]
    if (nrow(loocv_best) == 0) {
        warning(sprintf("%s: bestTune params not matched in fold predictions — using all folds",
                        cluster_label))
        loocv_best <- loocv_preds
    }

    roc_obj <- tryCatch(
        roc(loocv_best$obs, loocv_best[[case_label]],
            levels = c("Control", case_label), direction = "<", quiet = TRUE),
        error = function(e) NULL
    )
    if (is.null(roc_obj)) return(NULL)

    auc_val <- as.numeric(auc(roc_obj))
    ci_obj  <- tryCatch(ci.auc(roc_obj, conf.level = 0.95, method = "delong"),
                        error = function(e) NULL)
    auc_lo  <- if (!is.null(ci_obj)) ci_obj[1] else NA_real_
    auc_hi  <- if (!is.null(ci_obj)) ci_obj[3] else NA_real_
    cat(sprintf("  LOOCV AUC = %.3f (95%% CI: %.3f–%.3f)\n", auc_val, auc_lo, auc_hi))

    # Save glmnet coefficients at best lambda
    coef_mat <- coef(fit$finalModel, s = fit$bestTune$lambda)
    coef_dt  <- data.table(
        feature = rownames(coef_mat)[-1],
        coef    = as.numeric(coef_mat[-1]),
        cluster = cluster_label
    )
    coef_dt <- coef_dt[coef != 0][order(-abs(coef))]
    fwrite(coef_dt, file.path(OUT_DIR, sprintf("cluster_prems_features_%s.csv", cluster_label)))

    list(cluster = cluster_label, n_cases = n_cases,
         auc = auc_val, auc_lo = auc_lo, auc_hi = auc_hi, roc = roc_obj)
}

ctrl_dt     <- dt[get(status_col) == ctrl_state]
all_results <- list()

for (cl in .cluster_labels) {
    pre_eids_cl <- dt[get(status_col) == pre_state & cluster_f == cl, eid]
    all_results[[cl]] <- fit_cluster_model(cl, pre_eids_cl, ctrl_dt)
}

# ── 5. AUC summary (cluster-specific + Figure 1 pooled baseline) ─────────────
auc_dt <- rbindlist(lapply(.cluster_labels, function(nm) {
    r <- all_results[[nm]]
    if (is.null(r)) data.table(cluster = nm, n_cases = NA_integer_,
                               auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_)
    else data.table(cluster = nm, n_cases = r$n_cases,
                    auc = r$auc, auc_lo = r$auc_lo, auc_hi = r$auc_hi)
}))

# Prepend pooled baseline from Figure 1
pooled_row <- data.table(
    cluster  = "All pre-onset",
    n_cases  = dt[get(status_col) == pre_state, .N],
    auc      = pooled_auc,
    auc_lo   = pooled_auc_lo,
    auc_hi   = pooled_auc_hi
)
auc_dt <- rbind(pooled_row, auc_dt)
fwrite(auc_dt, file.path(OUT_DIR, "cluster_prems_auc_summary.csv"))
cat("\nAUC summary:\n"); print(auc_dt)

# ── 6. Panel v: AUC forest plot ───────────────────────────────────────────────
cat("\nBuilding AUC forest plot (panel v)...\n")

plot_dt <- auc_dt[!is.na(auc)]
plot_dt[, cluster := factor(cluster,
                             levels = rev(c("All pre-onset", .cluster_labels)))]
.forest_palette <- c(cluster_palette(length(.cluster_labels)),
                     "All pre-onset" = "#1A1A1A")

p_v <- ggplot(plot_dt, aes(y = cluster, x = auc, colour = cluster)) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               linewidth = 0.35, colour = "grey55") +
    geom_errorbar(aes(xmin = auc_lo, xmax = auc_hi),
                  orientation = "y", width = 0.22, linewidth = 0.55, na.rm = TRUE) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = .forest_palette, guide = "none") +
    scale_x_continuous(limits = c(0.3, 1.0), breaks = seq(0.3, 1.0, 0.1),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(y        = NULL,
         x        = "AUC (95% CI — DeLong)",
         title    = glue("v  Pre-{case_label} classifier AUC by subtype"),
         subtitle = "glmnet | Boruta features (Figure 1) | LOOCV | dashed = chance") +
    theme_ukb(base_size = 9) +
    theme(axis.text.y       = element_text(size = 8.5, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3))

ggsave(file.path(FIG_DIR, "panel_v.pdf"), p_v,
       width = 3.5, height = 3.0, device = cairo_pdf)
cat("  panel_v.pdf\n")

# ── 7. Panel n2: ROC overlay for clusters with valid AUC ─────────────────────
roc_results <- Filter(Negate(is.null), all_results)
if (length(roc_results) >= 1) {
    roc_dt <- rbindlist(lapply(names(roc_results), function(nm) {
        r <- roc_results[[nm]]
        data.table(cluster = nm,
                   fpr = 1 - r$roc$specificities,
                   tpr = r$roc$sensitivities)
    }))
    roc_dt[, cluster := factor(cluster, levels = .cluster_labels)]
    auc_labels <- setNames(
        vapply(names(roc_results), function(nm)
            sprintf("%s  AUC=%.2f", nm, roc_results[[nm]]$auc),
            character(1)),
        names(roc_results)
    )

    p_roc <- ggplot(roc_dt, aes(x = fpr, y = tpr, colour = cluster)) +
        geom_line(linewidth = 0.8, alpha = 0.9) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    linewidth = 0.3, colour = "grey60") +
        scale_colour_manual(values = CLUST_COLS,
                            labels = auc_labels[levels(roc_dt$cluster)]) +
        scale_x_continuous(limits = c(0,1), expand = expansion(mult = 0.02)) +
        scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0.02)) +
        labs(x        = "1 - Specificity",
             y        = "Sensitivity",
             title    = glue("n2  ROC — subtype-specific pre-{case_label} (LOOCV)"),
             subtitle = "glmnet | Boruta features (Figure 1) | exploratory",
             colour   = NULL) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom",
              legend.text     = element_text(size = 7.5))

    ggsave(file.path(FIG_DIR, "panel_v2_roc.pdf"), p_roc,
           width = 4.0, height = 4.0, device = cairo_pdf)
    cat("  panel_v2_roc.pdf\n")
}

cat("\n11_ms_cluster_prems_ml.R complete.\n")
cat(glue("Note: pooled AUC is from Figure 1 ({basename(FIG1_AUC_FILE)}).\n"))
cat("      Cluster AUCs use LOOCV — treat as exploratory given small n.\n")
