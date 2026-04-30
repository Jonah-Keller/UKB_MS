#!/usr/bin/env Rscript
# 03_cognitive_ms.R
# Cognitive function in MS — UKB-unique extension
#
# Approach:
#   1. MS vs HC comparison: g-factor, fluid intelligence, reaction time, symbol digit
#   2. Pre-onset vs post-onset vs HC
#   3. Spearman correlation: 173 MS DEPs vs cognitive scores
#   4. Linear regression: g-factor ~ NEFL + age + sex + years_to_dx (MS cases)
#
# Data: cognitive_gcsvd_merged.csv from UKB cognitive_function/
#   gcsvd_score           — general cognitive factor (higher = better)
#   fluid_intelligence    — fluid intelligence score
#   reaction_time         — mean reaction time ms (higher = slower)
#   symbol_digit_i2       — symbol digit substitution (processing speed, MS-sensitive)
#   trail_making_TMT_B_i2 — trail making B duration (executive function, higher = slower)
#
# Output: results/cognitive/
#   ms_cognitive_comparison.pdf
#   ms_cognitive_dep_correlations.pdf
#   ms_cognitive_nefl_regression.pdf
#   ms_cognitive_results.csv

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
COG_FILE    <- file.path(CADASIL_DIR, "cognitive_function", "cognitive_gcsvd_merged.csv")
QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DIFF_FILE   <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
OUT_DIR     <- file.path(PROJ_DIR, "results", "cognitive")
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

# ── 2. Load MS cohort ─────────────────────────────────────────────────────────
cat("Loading MS cohort...\n")
ms_qc <- fread(QC_FILE, showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]

dt <- merge(ms_qc, cog, by = "eid", all.x = FALSE)
cat(sprintf("  MS cases with cog data: %d (pre=%d, post=%d)\n",
            sum(dt$ms_status %in% c("pre_onset","post_onset")),
            sum(dt$ms_status == "pre_onset"),
            sum(dt$ms_status == "post_onset")))
cat(sprintf("  Controls with cog data: %d\n", sum(dt$ms_status == "control")))

dt[, is_ms := as.integer(ms_status %in% c("pre_onset","post_onset"))]
dt[, ms_group := fcase(
    ms_status == "pre_onset",  "Pre-onset MS",
    ms_status == "post_onset", "Post-onset MS",
    ms_status == "control",    "Control"
)]
dt <- dt[!is.na(ms_group)]
dt[, ms_group_f := factor(ms_group, levels = c("Control","Pre-onset MS","Post-onset MS"))]

# ── 3. MS vs HC: cognitive score comparisons ──────────────────────────────────
cat("\nMS vs HC cognitive comparisons...\n")

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

ms_vs_hc <- rbindlist(lapply(cog_tests, function(ct) {
    sub <- dt[!is.na(get(ct$col)) & ms_status %in% c("pre_onset","post_onset","control")]
    if (sum(sub$is_ms) < 10) return(NULL)
    wt <- wilcox.test(get(ct$col) ~ is_ms, data = sub, exact = FALSE)
    data.table(
        measure    = ct$label,
        col        = ct$col,
        higher_better = ct$higher,
        ms_n       = sum(!is.na(sub[is_ms==1, get(ct$col)])),
        hc_n       = sum(!is.na(sub[is_ms==0, get(ct$col)])),
        ms_mean    = mean(sub[is_ms==1, get(ct$col)], na.rm=TRUE),
        hc_mean    = mean(sub[is_ms==0, get(ct$col)], na.rm=TRUE),
        cohen_d    = (mean(sub[is_ms==1, get(ct$col)], na.rm=TRUE) -
                      mean(sub[is_ms==0, get(ct$col)], na.rm=TRUE)) /
                     sd(sub[[ct$col]], na.rm=TRUE),
        pval       = wt$p.value
    )
}))
ms_vs_hc[, fdr := p.adjust(pval, method = "BH")]
ms_vs_hc[, sig_label := fcase(fdr<0.01,"**", fdr<0.05,"*", pval<0.05,"+", default="")]
cat(sprintf("  FDR<0.05: %d / %d cognitive measures\n",
            sum(ms_vs_hc$fdr < 0.05), nrow(ms_vs_hc)))
print(ms_vs_hc[, .(measure, ms_n, ms_mean, hc_mean, cohen_d, pval, fdr)])
fwrite(ms_vs_hc, file.path(OUT_DIR, "ms_cognitive_results.csv"))

# ── 4. Boxplot comparison ──────────────────────────────────────────────────────
# Sample controls for display
ctrl_ids <- dt[ms_group == "Control", eid]
ctrl_samp <- sample(ctrl_ids, min(2000L, length(ctrl_ids)))
plot_dt <- dt[ms_group != "Control" | eid %in% ctrl_samp]

GROUP_COLS <- c("Control" = "grey65", "Pre-onset MS" = "#E6A817", "Post-onset MS" = "#CC0066")

box_plots <- lapply(cog_tests[1:min(4, length(cog_tests))], function(ct) {
    sub <- plot_dt[!is.na(get(ct$col))]
    p_fdr <- ms_vs_hc[col == ct$col, fdr]
    sig   <- ifelse(length(p_fdr) && p_fdr < 0.05,
                    sprintf("p[FDR]=%.3g", p_fdr), "n.s.")
    ggplot(sub, aes(x = ms_group_f, y = get(ct$col), fill = ms_group)) +
        geom_boxplot(outlier.size = 0.4, width = 0.5, colour = "#1A1A1A") +
        scale_fill_manual(values = GROUP_COLS, guide = "none") +
        annotate("text", x = 2, y = Inf,
                 label = sig, size = 2.5, vjust = 1.5, colour = "grey30") +
        labs(x = NULL, y = ct$label) +
        theme_ukb(base_size = 8) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
})

p_cog <- wrap_plots(box_plots, ncol = 2) +
    plot_annotation(title = "Cognitive function: MS vs healthy controls",
                    subtitle = sprintf("MS n=%d, HC n=%d",
                                       sum(dt$ms_status %in% c("pre_onset","post_onset")),
                                       sum(dt$ms_status == "control")))
ggsave(file.path(OUT_DIR, "ms_cognitive_comparison.pdf"),
       p_cog, width = 7, height = 5.5, device = cairo_pdf)
cat("  Saved: ms_cognitive_comparison.pdf\n")

# ── 5. DEP–cognitive correlations ─────────────────────────────────────────────
cat("\nDEP–cognitive correlations (MS cases)...\n")
deps      <- fread(DIFF_FILE)
dep_prots <- tolower(deps[adj.P.Val < FDR_THR, protein])
avail_prots <- intersect(dep_prots, names(dt))
ms_only   <- dt[ms_status %in% c("pre_onset","post_onset")]
cat(sprintf("  MS cases: %d, proteins: %d\n", nrow(ms_only), length(avail_prots)))

# Focus on g-factor (most power) + symbol digit if available
cog_outcomes <- intersect(c("gcsvd_score","symbol_digit__n_correct__i2"),
                           names(ms_only))

corr_list <- list()
for (cog_col in cog_outcomes) {
    for (prot in avail_prots) {
        vals <- ms_only[!is.na(get(cog_col)) & !is.na(get(prot)),
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
    fwrite(corr_dt, file.path(OUT_DIR, "ms_cognitive_dep_correlations.csv"))
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
                 title = "MS DEPs correlated with general cognitive factor (g)",
                 subtitle = sprintf("MS cases n≈%d", max(gcsvd_corr$n))) +
            theme_ukb()

        ggsave(file.path(OUT_DIR, "ms_cognitive_dep_correlations.pdf"),
               p_vol, width = 6, height = 5, device = cairo_pdf)
        cat("  Saved: ms_cognitive_dep_correlations.pdf\n")
    }
}

# ── 6. NEFL ~ cognitive regression ────────────────────────────────────────────
cat("\nNEFL ~ cognitive regression...\n")
nefl_col <- intersect(c("nefl","nfl"), names(ms_only))[1]

if (!is.na(nefl_col) && "gcsvd_score" %in% names(ms_only)) {
    reg_dt <- ms_only[!is.na(get(nefl_col)) & !is.na(gcsvd_score) &
                      !is.na(age_at_sampling) & !is.na(sex)]
    reg_dt[, sex_num := as.integer(factor(sex)) - 1L]
    if (nrow(reg_dt) >= 20) {
        lm_adj <- lm(gcsvd_score ~ get(nefl_col) + age_at_sampling + sex_num, data = reg_dt)
        lm_coef <- coef(summary(lm_adj))
        nefl_beta <- lm_coef[2, "Estimate"]
        nefl_p    <- lm_coef[2, "Pr(>|t|)"]
        cat(sprintf("  NEFL ~ g-factor (adjusted): β=%.3f, p=%.3g\n", nefl_beta, nefl_p))

        # Scatter with jitter by status
        col_by_status <- c(pre_onset = "#E6A817", post_onset = "#CC0066")
        p_nefl_cog <- ggplot(reg_dt,
                             aes(x = get(nefl_col), y = gcsvd_score, colour = ms_status)) +
            geom_point(alpha = 0.6, size = 1.8) +
            geom_smooth(method = "lm", se = TRUE, colour = "#2B4C7E",
                        linewidth = 0.8, fill = "#2B4C7E", alpha = 0.12) +
            scale_colour_manual(values = col_by_status,
                                labels = c(pre_onset = "Pre-onset MS",
                                           post_onset = "Post-onset MS"),
                                name = NULL) +
            annotate("text", x = Inf, y = Inf,
                     label = sprintf("β=%.3f, p=%s (age/sex adj.)",
                                     nefl_beta,
                                     ifelse(nefl_p < 0.001, "<0.001",
                                            sprintf("%.3f", nefl_p))),
                     hjust = 1.05, vjust = 1.5, size = 2.8) +
            labs(x = "NfL NPX (log2)", y = "General cognitive factor (g-score)",
                 title = "NfL vs cognitive function in MS") +
            theme_ukb()

        ggsave(file.path(OUT_DIR, "ms_cognitive_nefl_regression.pdf"),
               p_nefl_cog, width = 5.5, height = 4.5, device = cairo_pdf)
        cat("  Saved: ms_cognitive_nefl_regression.pdf\n")
    }
}

cat("\n03_cognitive_ms.R complete.\n")
