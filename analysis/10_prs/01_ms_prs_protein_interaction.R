#!/usr/bin/env Rscript
# 01_ms_prs_protein_interaction.R
# PRS × protein trajectory interaction
#
# Does higher genetic risk (PRS) accelerate the proteomic signal pre-onset?
#
# Design:
#   - All cases (pre + post) with PRS scores
#   - Split by PRS quartile (Q1=low, Q4=high)
#   - For top CNS + immune proteins: fit spline trajectory per quartile
#   - Formal interaction test: protein ~ ns(ytd) * prs_quartile + age + sex
#
# Output: results/prs/
#   {cohort}_prs_trajectories.pdf      — PRS-stratified spline trajectories
#   {cohort}_prs_interaction_tests.csv — interaction p-values per protein

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(splines)
    library(patchwork)
    library(here)
    library(glue)
})
source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))
cfg <- load_disease_config()

cohort       <- cfg$cohort_short
status_col   <- cfg$cohort_status_col
sv           <- cfg$status_values
prs_col      <- cfg$prs_combined_col
disease_caps <- cfg$disease_short_caps

QC_FILE   <- here::here("data", "ukb", "olink", "processed",
                        glue::glue("{cohort}_olink_qc.csv"))
PRS_FILE  <- here::here("data", "ukb", "genetics",
                        glue::glue("{cohort}_prs_scores.csv"))
DIFF_FILE <- here::here("results", "differential",
                        glue::glue("{cohort}_combined_vs_hc.csv"))
OUT_DIR   <- here::here("results", "prs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

YTD_MIN <- -20; YTD_MAX <- 16
set.seed(42)

# ── 1. Load & merge ───────────────────────────────────────────────────────────
cat("Loading data...\n")
qc  <- fread(QC_FILE,  showProgress = FALSE)
prs <- fread(PRS_FILE, showProgress = FALSE)
dt  <- merge(qc[qc_outlier == FALSE], prs[, c("eid", prs_col), with = FALSE],
             by = "eid")

case_states <- c(sv$pre_onset, sv$post_onset)
ms <- dt[get(status_col) %in% case_states &
         !is.na(years_to_diagnosis) & !is.na(get(prs_col)) &
         years_to_diagnosis >= YTD_MIN & years_to_diagnosis <= YTD_MAX]
cat(sprintf("  %s cases with PRS: %d (pre=%d, post=%d)\n",
            disease_caps, nrow(ms),
            sum(ms[[status_col]] == sv$pre_onset),
            sum(ms[[status_col]] == sv$post_onset)))

# ── 2. PRS quartile assignment ────────────────────────────────────────────────
ms[, prs_q := cut(get(prs_col),
                   breaks = quantile(get(prs_col), c(0, .25, .50, .75, 1), na.rm = TRUE),
                   labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
                   include.lowest = TRUE)]
cat("PRS quartile sizes:\n")
print(ms[, .N, by = prs_q][order(prs_q)])

# ── 3. Protein list: top CNS + immune DEPs ────────────────────────────────────
deps <- fread(DIFF_FILE)
deps[, fdr := p.adjust(P.Value, method = "BH")]
top_prots <- deps[fdr < 0.05][order(-abs(logFC)), tolower(protein)][seq_len(min(6, .N))]
# Always include NEFL and GFAP if available
core <- intersect(c("nefl", "gfap", top_prots), names(ms))
core <- unique(core)[seq_len(min(6, length(core)))]
cat("  Proteins:", paste(toupper(core), collapse = ", "), "\n")

# ── 4. Trajectory splines per quartile ───────────────────────────────────────
x_grid <- seq(YTD_MIN, YTD_MAX, length.out = 100)
PRS_COLS <- c("Q1 (low)" = "#56B4E9", "Q2" = "#74c476",
              "Q3" = "#E6A817",        "Q4 (high)" = "#CC0066")

plot_list <- list()

for (prot in core) {
    d <- ms[!is.na(get(prot)), .(
        npx = get(prot), ytd = years_to_diagnosis,
        prs_q, age = age_at_sampling, sex
    )]
    if (nrow(d) < 30) next

    # Fit per-quartile splines
    traj_all <- rbindlist(lapply(levels(ms$prs_q), function(q) {
        dq <- d[prs_q == q]
        if (nrow(dq) < 8) return(NULL)
        fit <- tryCatch(
            lm(npx ~ ns(ytd, df = 3) + age + sex, data = dq),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        nd  <- data.table(ytd = x_grid,
                          age = median(dq$age, na.rm=TRUE),
                          sex = as.numeric(names(sort(table(dq$sex), decreasing=TRUE))[1]))
        pr  <- predict(fit, newdata = nd, se.fit = TRUE)
        data.table(ytd = x_grid, pred = pr$fit,
                   lci = pr$fit - 1.96*pr$se.fit,
                   uci = pr$fit + 1.96*pr$se.fit,
                   prs_q = q, n_q = nrow(dq))
    }))
    if (nrow(traj_all) == 0) next
    traj_all[, prs_q := factor(prs_q, levels = levels(ms$prs_q))]

    p <- ggplot(traj_all, aes(x = ytd, colour = prs_q, fill = prs_q)) +
        geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.12, colour = NA) +
        geom_line(aes(y = pred), linewidth = 0.75) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
        annotate("text", x = 0.25, y = Inf, label = "Dx",
                 hjust = 0, vjust = 1.5, size = 2.2, colour = "grey40") +
        scale_colour_manual(values = PRS_COLS, name = "PRS quartile") +
        scale_fill_manual(values = PRS_COLS,   name = "PRS quartile") +
        scale_x_continuous(breaks = seq(-16, 16, by = 4)) +
        labs(title = toupper(prot),
             x = "Years relative to MS diagnosis",
             y = paste0(toupper(prot), " (NPX)")) +
        theme_ukb(base_size = 9) +
        guides(colour = guide_legend(override.aes = list(linewidth = 1.2)),
               fill   = "none")

    plot_list[[prot]] <- p
}

if (length(plot_list) > 0) {
    n_p   <- length(plot_list)
    n_col <- min(3, n_p)
    fig_prs <- wrap_plots(plot_list, ncol = n_col) +
        plot_annotation(
            title    = glue::glue("PRS-stratified protein trajectories across {disease_caps} disease course"),
            subtitle = "Splines per PRS quartile (Q1=lowest, Q4=highest genetic risk) | age + sex adjusted",
            theme    = theme_ukb(base_size = 9)
        )
    ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_trajectories.pdf")),
           fig_prs, width = 4 * n_col, height = 3.8 * ceiling(n_p / n_col),
           device = cairo_pdf)
    cat(sprintf("  Saved: %s_prs_trajectories.pdf\n", cohort))
}

# ── 5. Formal interaction tests ───────────────────────────────────────────────
cat("\nRunning interaction tests (protein ~ ns(ytd)*prs_quartile + age + sex)...\n")
int_results <- rbindlist(lapply(core, function(prot) {
    d <- ms[!is.na(get(prot)), .(
        npx = get(prot), ytd = years_to_diagnosis,
        prs_q = as.integer(prs_q),   # ordinal: 1-4
        age = age_at_sampling, sex
    )]
    if (nrow(d) < 30) return(NULL)
    # Additive model
    fit0 <- tryCatch(
        lm(npx ~ ns(ytd, df=3) + prs_q + age + sex, data = d),
        error = function(e) NULL
    )
    # Interaction model
    fit1 <- tryCatch(
        lm(npx ~ ns(ytd, df=3) * prs_q + age + sex, data = d),
        error = function(e) NULL
    )
    if (is.null(fit0) || is.null(fit1)) return(NULL)
    av <- tryCatch(anova(fit0, fit1), error = function(e) NULL)
    if (is.null(av)) return(NULL)
    data.table(protein = toupper(prot),
               F_stat  = round(av$F[2], 3),
               pval    = signif(av$`Pr(>F)`[2], 4),
               n       = nrow(d))
}))

if (!is.null(int_results) && nrow(int_results) > 0) {
    int_results[, fdr := p.adjust(pval, method = "BH")]
    int_results <- int_results[order(pval)]
    cat("\nInteraction test results:\n")
    print(int_results)
    fwrite(int_results, file.path(OUT_DIR,
                                  glue::glue("{cohort}_prs_interaction_tests.csv")))
    cat(sprintf("  Saved: %s_prs_interaction_tests.csv\n", cohort))
}

cat("\n01_ms_prs_protein_interaction.R complete.\n")
