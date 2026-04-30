#!/usr/bin/env Rscript
# 03_cognitive_ms.R
# Cognitive function in disease cohort — UKB-unique extension
#
# Approach:
#   1. Case vs HC comparison: g-factor, fluid intelligence, reaction time, symbol digit
#   2. Pre-onset vs post-onset vs HC
#   3. Spearman correlation: DEPs vs cognitive scores
#   4. Linear regression: g-factor ~ NEFL + age + sex (cases)
#
# Data: cognitive_gcsvd_merged.csv from UKB cognitive_function/
#   gcsvd_score           — general cognitive factor (higher = better)
#   fluid_intelligence    — fluid intelligence score
#   reaction_time         — mean reaction time ms (higher = slower)
#   symbol_digit_i2       — symbol digit substitution (processing speed)
#   trail_making_TMT_B_i2 — trail making B duration (executive function)
#
# Output: results/cognitive/<cohort>_cognitive_*.{csv,pdf}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(glue)
    library(here)
})

source(here::here("analysis", "helpers", "ukb_theme.R"))
source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

STATUS_COL <- cfg$cohort_status_col
PRE        <- cfg$status_values$pre_onset
POST       <- cfg$status_values$post_onset
CTRL       <- cfg$status_values$control
CASE_VALS  <- c(PRE, POST)

CADASIL_DIR <- file.path(dirname(here::here()),
                         "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
COG_FILE    <- file.path(CADASIL_DIR, "cognitive_function", "cognitive_gcsvd_merged.csv")
QC_FILE     <- here::here("data", "ukb", "olink", "processed",
                          glue::glue("{cfg$cohort_short}_olink_qc.csv"))
DIFF_FILE   <- here::here("results", "differential",
                          glue::glue("{cfg$cohort_short}_combined_vs_hc.csv"))
OUT_DIR     <- here::here("results", "cognitive")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR <- 0.05

# ── 1. Load cognitive data ─────────────────────────────────────────────────────
cat("Loading cognitive data...\n")
cog_cols <- c("eid","gcsvd_score",
              "fluid_intelligence__fi_score__i0",
              "reaction_time__mean_reaction_time_ms__i0",
              "symbol_digit__n_correct__i2",
              "symbol_digit__n_attempted__i2",
              "trail_making__duration_alphanumeric_trail_s__i2",
              "trail_making__duration_numeric_trail_s__i2",
              "matrix_pattern__n_puzzles_correct__i2")

cog <- fread(COG_FILE, select = intersect(cog_cols, fread(COG_FILE, nrows=0) |> names()),
             showProgress = FALSE)
cat(sprintf("  Cognitive data: %d participants\n", nrow(cog)))

# ── 2. Load disease cohort ────────────────────────────────────────────────────
cat(sprintf("Loading %s cohort...\n", cfg$disease_short_caps))
case_qc <- fread(QC_FILE, showProgress = FALSE)
case_qc <- case_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]

dt <- merge(case_qc, cog, by = "eid", all.x = FALSE)
cat(sprintf("  %s cases with cog data: %d (pre=%d, post=%d)\n",
            cfg$disease_short_caps,
            sum(dt[[STATUS_COL]] %in% CASE_VALS),
            sum(dt[[STATUS_COL]] == PRE),
            sum(dt[[STATUS_COL]] == POST)))
cat(sprintf("  Controls with cog data: %d\n", sum(dt[[STATUS_COL]] == CTRL)))

PRE_LABEL  <- glue::glue("Pre-onset {cfg$disease_short_caps}")
POST_LABEL <- glue::glue("Post-onset {cfg$disease_short_caps}")
GROUP_LEVELS <- c("Control", PRE_LABEL, POST_LABEL)

dt[, is_case := as.integer(get(STATUS_COL) %in% CASE_VALS)]
dt[, case_group := fcase(
    get(STATUS_COL) == PRE,  PRE_LABEL,
    get(STATUS_COL) == POST, POST_LABEL,
    get(STATUS_COL) == CTRL, "Control"
)]
dt <- dt[!is.na(case_group)]
dt[, case_group_f := factor(case_group, levels = GROUP_LEVELS)]

# ── 3. Case vs HC: cognitive score comparisons ────────────────────────────────
cat(sprintf("\n%s vs HC cognitive comparisons...\n", cfg$disease_short_caps))

cog_tests <- list(
    list(col = "gcsvd_score",             label = "General cognitive factor (g)", higher = TRUE),
    list(col = "fluid_intelligence__fi_score__i0",    label = "Fluid intelligence",          higher = TRUE),
    list(col = "reaction_time__mean_reaction_time_ms__i0", label = "Reaction time (ms)",    higher = FALSE),
    list(col = "symbol_digit__n_correct__i2",         label = "Symbol digit correct (i2)",   higher = TRUE),
    list(col = "trail_making__duration_alphanumeric_trail_s__i2",
         label = "Trail-making B duration (i2)", higher = FALSE),
    list(col = "matrix_pattern__n_puzzles_correct__i2",
         label = "Matrix pattern correct (i2)", higher = TRUE)
)

# Keep only available columns
cog_tests <- Filter(function(x) x$col %in% names(dt), cog_tests)

case_vs_hc <- rbindlist(lapply(cog_tests, function(ct) {
    sub <- dt[!is.na(get(ct$col)) & get(STATUS_COL) %in% c(CASE_VALS, CTRL)]
    if (sum(sub$is_case) < 10) return(NULL)
    wt <- wilcox.test(get(ct$col) ~ is_case, data = sub, exact = FALSE)
    data.table(
        measure       = ct$label,
        col           = ct$col,
        higher_better = ct$higher,
        case_n        = sum(!is.na(sub[is_case==1, get(ct$col)])),
        hc_n          = sum(!is.na(sub[is_case==0, get(ct$col)])),
        case_mean     = mean(sub[is_case==1, get(ct$col)], na.rm=TRUE),
        hc_mean       = mean(sub[is_case==0, get(ct$col)], na.rm=TRUE),
        cohen_d       = (mean(sub[is_case==1, get(ct$col)], na.rm=TRUE) -
                         mean(sub[is_case==0, get(ct$col)], na.rm=TRUE)) /
                        sd(sub[[ct$col]], na.rm=TRUE),
        pval          = wt$p.value
    )
}))
case_vs_hc[, fdr := p.adjust(pval, method = "BH")]
case_vs_hc[, sig_label := fcase(fdr<0.01,"**", fdr<0.05,"*", pval<0.05,"+", default="")]
cat(sprintf("  FDR<0.05: %d / %d cognitive measures\n",
            sum(case_vs_hc$fdr < 0.05), nrow(case_vs_hc)))
print(case_vs_hc[, .(measure, case_n, case_mean, hc_mean, cohen_d, pval, fdr)])
fwrite(case_vs_hc, file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_cognitive_results.csv")))

# ── 4. Boxplot comparison ──────────────────────────────────────────────────────
# Sample controls for display
ctrl_ids <- dt[case_group == "Control", eid]
ctrl_samp <- sample(ctrl_ids, min(2000L, length(ctrl_ids)))
plot_dt <- dt[case_group != "Control" | eid %in% ctrl_samp]

GROUP_COLS <- setNames(c("grey65", "#E6A817", "#CC0066"), GROUP_LEVELS)

box_plots <- lapply(cog_tests[1:min(4, length(cog_tests))], function(ct) {
    sub <- plot_dt[!is.na(get(ct$col))]
    p_fdr <- case_vs_hc[col == ct$col, fdr]
    sig   <- ifelse(length(p_fdr) && p_fdr < 0.05,
                    sprintf("p[FDR]=%.3g", p_fdr), "n.s.")
    ggplot(sub, aes(x = case_group_f, y = get(ct$col), fill = case_group)) +
        geom_boxplot(outlier.size = 0.4, width = 0.5, colour = "#1A1A1A") +
        scale_fill_manual(values = GROUP_COLS, guide = "none") +
        annotate("text", x = 2, y = Inf,
                 label = sig, size = 2.5, vjust = 1.5, colour = "grey30") +
        labs(x = NULL, y = ct$label) +
        theme_ukb(base_size = 8) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
})

p_cog <- wrap_plots(box_plots, ncol = 2) +
    plot_annotation(title = glue::glue("Cognitive function: {cfg$disease_short_caps} vs healthy controls"),
                    subtitle = sprintf("%s n=%d, HC n=%d",
                                       cfg$disease_short_caps,
                                       sum(dt[[STATUS_COL]] %in% CASE_VALS),
                                       sum(dt[[STATUS_COL]] == CTRL)))
ggsave(file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_cognitive_comparison.pdf")),
       p_cog, width = 7, height = 5.5, device = cairo_pdf)
cat(sprintf("  Saved: %s_cognitive_comparison.pdf\n", cfg$cohort_short))

# ── 5. DEP–cognitive correlations ─────────────────────────────────────────────
cat(sprintf("\nDEP–cognitive correlations (%s cases)...\n", cfg$disease_short_caps))
deps        <- fread(DIFF_FILE)
dep_prots   <- tolower(deps[adj.P.Val < FDR_THR, protein])
avail_prots <- intersect(dep_prots, names(dt))
case_only   <- dt[get(STATUS_COL) %in% CASE_VALS]
cat(sprintf("  %s cases: %d, proteins: %d\n",
            cfg$disease_short_caps, nrow(case_only), length(avail_prots)))

# Focus on g-factor (most power) + symbol digit if available
cog_outcomes <- intersect(c("gcsvd_score","symbol_digit__n_correct__i2"),
                           names(case_only))

corr_list <- list()
for (cog_col in cog_outcomes) {
    for (prot in avail_prots) {
        vals <- case_only[!is.na(get(cog_col)) & !is.na(get(prot)),
                          .(cog = get(cog_col), npx = get(prot))]
        if (nrow(vals) < 15) next
        ct <- cor.test(vals$cog, vals$npx, method = "spearman", exact = FALSE)
        corr_list[[paste(cog_col, prot)]] <- data.table(
            cog_measure = cog_col,
            protein     = prot,
            rho         = ct$estimate,
            pval        = ct$p.value,
            n           = nrow(vals)
        )
    }
}

if (length(corr_list) > 0) {
    corr_dt <- rbindlist(corr_list)
    corr_dt[, fdr := p.adjust(pval, method = "BH"), by = cog_measure]
    corr_dt <- corr_dt[order(pval)]
    fwrite(corr_dt, file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_cognitive_dep_correlations.csv")))
    cat(sprintf("  %d protein–cognitive pairs; FDR<0.05: %d\n",
                nrow(corr_dt), sum(corr_dt$fdr < 0.05)))

    # Volcano-style plot for g-factor
    gcsvd_corr <- corr_dt[cog_measure == "gcsvd_score"]
    if (nrow(gcsvd_corr) > 0) {
        gcsvd_corr[, direction := fifelse(rho > 0, "Positive", "Negative")]
        gcsvd_corr[, sig_flag := fdr < 0.05]
        gcsvd_corr[, label_prot := fifelse(
            rank(pval) <= 10 | abs(rho) > 0.25, toupper(protein), NA_character_
        )]

        p_vol <- ggplot(gcsvd_corr, aes(x = rho, y = -log10(pval),
                                         colour = interaction(direction, sig_flag))) +
            geom_point(size = 1.5, alpha = 0.7) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                       colour = "grey50", linewidth = 0.3) +
            geom_text_repel(aes(label = label_prot), size = 2.2,
                            na.rm = TRUE, max.overlaps = 15) +
            scale_colour_manual(
                values = c("Positive.FALSE" = "grey70", "Negative.FALSE" = "grey70",
                           "Positive.TRUE"  = "#CC0066", "Negative.TRUE"  = "#56B4E9"),
                labels = c("Positive.FALSE" = "Pos (n.s.)", "Negative.FALSE" = "Neg (n.s.)",
                           "Positive.TRUE"  = "Pos (FDR<0.05)", "Negative.TRUE"  = "Neg (FDR<0.05)"),
                name = NULL
            ) +
            labs(x = "Spearman ρ (protein ~ g-factor)",
                 y = expression(-log[10](p)),
                 title = glue::glue("{cfg$disease_short_caps} DEPs correlated with general cognitive factor (g)"),
                 subtitle = sprintf("%s cases n≈%d",
                                    cfg$disease_short_caps, max(gcsvd_corr$n))) +
            theme_ukb()

        ggsave(file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_cognitive_dep_correlations.pdf")),
               p_vol, width = 6, height = 5, device = cairo_pdf)
        cat(sprintf("  Saved: %s_cognitive_dep_correlations.pdf\n", cfg$cohort_short))
    }
}

# ── 6. NEFL ~ cognitive regression ────────────────────────────────────────────
cat("\nNEFL ~ cognitive regression...\n")
nefl_col <- intersect(c("nefl","nfl"), names(case_only))[1]

if (!is.na(nefl_col) && "gcsvd_score" %in% names(case_only)) {
    reg_dt <- case_only[!is.na(get(nefl_col)) & !is.na(gcsvd_score) &
                        !is.na(age_at_sampling) & !is.na(sex)]
    reg_dt[, sex_num := as.integer(factor(sex)) - 1L]
    if (nrow(reg_dt) >= 20) {
        lm_adj <- lm(gcsvd_score ~ get(nefl_col) + age_at_sampling + sex_num, data = reg_dt)
        lm_coef <- coef(summary(lm_adj))
        nefl_beta <- lm_coef[2, "Estimate"]
        nefl_p    <- lm_coef[2, "Pr(>|t|)"]
        cat(sprintf("  NEFL ~ g-factor (adjusted): β=%.3f, p=%.3g\n", nefl_beta, nefl_p))

        # Scatter with jitter by status
        col_by_status   <- setNames(c("#E6A817", "#CC0066"), CASE_VALS)
        label_by_status <- setNames(c(PRE_LABEL, POST_LABEL), CASE_VALS)
        p_nefl_cog <- ggplot(reg_dt,
                             aes(x = get(nefl_col), y = gcsvd_score, colour = get(STATUS_COL))) +
            geom_point(alpha = 0.6, size = 1.8) +
            geom_smooth(method = "lm", se = TRUE, colour = "#2B4C7E",
                        linewidth = 0.8, fill = "#2B4C7E", alpha = 0.12) +
            scale_colour_manual(values = col_by_status,
                                labels = label_by_status,
                                name = NULL) +
            annotate("text", x = Inf, y = Inf,
                     label = sprintf("β=%.3f, p=%s (age/sex adj.)",
                                     nefl_beta,
                                     ifelse(nefl_p < 0.001, "<0.001",
                                            sprintf("%.3f", nefl_p))),
                     hjust = 1.05, vjust = 1.5, size = 2.8) +
            labs(x = "NfL NPX (log2)", y = "General cognitive factor (g-score)",
                 title = glue::glue("NfL vs cognitive function in {cfg$disease_short_caps}")) +
            theme_ukb()

        ggsave(file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_cognitive_nefl_regression.pdf")),
               p_nefl_cog, width = 5.5, height = 4.5, device = cairo_pdf)
        cat(sprintf("  Saved: %s_cognitive_nefl_regression.pdf\n", cfg$cohort_short))
    }
}

cat("\n03_cognitive_ms.R complete.\n")
