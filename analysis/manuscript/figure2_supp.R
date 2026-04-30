#!/usr/bin/env Rscript
# figure2_supp.R — Supplementary Figure 2: Demoted pre-onset DEP + literature biomarkers
#
# BGN is a nominal pre-onset DEP (included here to make room in main Figure 2
# for the promoted ERBB2 trajectory). GFAP and MOG are established MS biomarkers
# that do not reach FDR<0.20 in the pre-onset limma analysis (GFAP FDR=0.964,
# MOG FDR=0.396) and are retained here as supplementary context alongside their
# post-onset trajectories.
#
# Panels:
#   a  GFAP  trajectory (not significant pre-onset; full disease course)
#   b  MOG   trajectory (not significant pre-onset; full disease course)
#   c  BGN   trajectory (nominal pre-onset DEP; demoted from main Fig 2)
#
# Output: results/figures/2_supp/

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(splines)
    library(MatchIt)
    library(glue)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))

cfg          <- load_disease_config()
COHORT       <- cfg$cohort_short
DISEASE_CAPS <- cfg$disease_short_caps
STATUS_COL   <- cfg$cohort_status_col
SV           <- cfg$status_values
PRS_COL      <- cfg$prs_combined_col
PRS_LABEL    <- cfg$prs_label
PRE_LBL      <- glue("pre-{DISEASE_CAPS}")

FIG_DIR  <- file.path(PROJ_DIR, "results", "figures", "2_supp")
DATA_DIR <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed")
DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

save_panel <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(FIG_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: panel_", name, ".pdf\n", sep = "")
}

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc <- fread(file.path(DATA_DIR, glue("{COHORT}_olink_qc.csv")))
if (STATUS_COL != "ms_status" && STATUS_COL %in% names(qc))
    setnames(qc, STATUS_COL, "ms_status")
qc <- qc[is.na(age_at_diagnosis) | (age_at_diagnosis >= 5 & age_at_diagnosis <= 90)]

TRAJ_FILE <- file.path(DIFF_DIR, "cns_trajectories.csv")
traj <- fread(TRAJ_FILE)
traj[, protein := toupper(protein)]

# PSM-matched HC reference (mirrors figure2.R)
cat("  Computing PSM-matched HC reference...\n")
ms_pool <- qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE &
               !is.na(age_at_sampling) & !is.na(sex) & !is.na(years_to_diagnosis)]
hc_pool <- qc[ms_status == SV$control & qc_outlier == FALSE &
               !is.na(age_at_sampling) & !is.na(sex)]
pool_dt <- rbind(
    ms_pool[, .(eid, is_ms = 1L, age_at_sampling, sex, years_to_diagnosis)],
    hc_pool[, .(eid, is_ms = 0L, age_at_sampling, sex, years_to_diagnosis = NA_real_)]
)
set.seed(42)
m_obj <- tryCatch(
    matchit(is_ms ~ age_at_sampling + sex, data = as.data.frame(pool_dt),
            method = "nearest", ratio = 3, replace = FALSE),
    error = function(e) NULL
)
psm_hc_ref <- NULL
if (!is.null(m_obj)) {
    md      <- as.data.table(match.data(m_obj))
    ms_sub  <- md[is_ms == 1, .(subclass, ytd_ms = years_to_diagnosis)]
    hc_md   <- merge(md[is_ms == 0, .(eid, subclass)], ms_sub, by = "subclass")
    hc_rows <- merge(hc_md[, .(eid, pseudo_ytd = ytd_ms)], hc_pool, by = "eid")
    hc_rows[, years_to_diagnosis := pseudo_ytd]
    psm_hc_ref <- hc_rows
    cat(sprintf("  PSM: %d HC observations\n", nrow(psm_hc_ref)))
}

hc_psm_spline <- function(col_nm, ytd_grid) {
    if (is.null(psm_hc_ref) || !col_nm %in% names(psm_hc_ref)) return(NULL)
    d <- psm_hc_ref[, .(ytd = years_to_diagnosis, npx = get(col_nm))][!is.na(npx) & !is.na(ytd)]
    if (nrow(d) < 10) return(NULL)
    fit <- tryCatch(lm(npx ~ ns(ytd, df = 3), data = d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pr <- predict(fit, newdata = data.frame(ytd = ytd_grid), se.fit = TRUE)
    data.table(ytd = ytd_grid,
               hc_ref     = pr$fit,
               hc_ref_lci = pr$fit - 1.96 * pr$se.fit,
               hc_ref_uci = pr$fit + 1.96 * pr$se.fit)
}

# ─────────────────────────────────────────────────────────────────────────────
# Inline BGN trajectory fit (BGN is not in cns_trajectories.csv; fit here in
# the same style as the figure2.R immune panels so efig2a is self-contained).
# ─────────────────────────────────────────────────────────────────────────────
fit_traj_inline <- function(col_nm, qc_in, ytd_grid = seq(-8, 12, by = 0.25)) {
    d <- qc_in[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE &
               !is.na(years_to_diagnosis) & !is.na(age_at_sampling) & !is.na(sex)]
    if (!col_nm %in% names(d)) return(NULL)
    d <- d[, .(npx = get(col_nm), years_to_diagnosis, age_at_sampling, sex)][!is.na(npx)]
    if (nrow(d) < 20) return(NULL)
    fit <- tryCatch(
        lm(npx ~ ns(years_to_diagnosis, df = 3) + age_at_sampling + sex, data = d),
        error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    nd <- data.frame(
        years_to_diagnosis = ytd_grid,
        age_at_sampling    = median(d$age_at_sampling, na.rm = TRUE),
        sex                = as.numeric(names(sort(table(d$sex), decreasing = TRUE))[1])
    )
    pr <- predict(fit, newdata = nd, se.fit = TRUE)
    hc_vals <- qc_in[ms_status == SV$control & qc_outlier == FALSE, get(col_nm)]
    hc_mean <- mean(hc_vals, na.rm = TRUE)
    hc_sd   <- sd(hc_vals,   na.rm = TRUE)
    data.table(
        protein     = toupper(col_nm),
        years_to_dx = ytd_grid,
        pred_npx    = pr$fit,
        pred_lci    = pr$fit - 1.96 * pr$se.fit,
        pred_uci    = pr$fit + 1.96 * pr$se.fit,
        hc_mean     = hc_mean,
        hc_lci      = hc_mean - hc_sd,
        hc_uci      = hc_mean + hc_sd
    )
}

bgn_traj <- fit_traj_inline("bgn", qc)
if (!is.null(bgn_traj)) traj <- rbind(traj, bgn_traj, fill = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Panels a, b, c — GFAP, MOG, and BGN trajectories
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panels a, b, c (GFAP, MOG, BGN trajectories)...\n")

SUPP_META <- list(
    GFAP = list(pid = "a_gfap_trajectory", col = "#E6A817",
                note = "GFAP \u2014 not significant pre-onset (FDR=0.964)"),
    MOG  = list(pid = "b_mog_trajectory",  col = COL_POST,
                note = "MOG \u2014 not significant pre-onset (FDR=0.396, nominal p=0.038)"),
    BGN  = list(pid = "c_bgn_trajectory",  col = COL_PRE_DOWN,
                note = "BGN \u2014 nominal pre-onset DEP (female-stratified FDR<0.05)")
)

for (prot in names(SUPP_META)) {
    m   <- SUPP_META[[prot]]
    t_p <- traj[protein == prot]
    if (nrow(t_p) == 0) { cat("  Skipping", prot, ": no trajectory data\n"); next }

    ytd_g   <- t_p$years_to_dx
    ref_p   <- hc_psm_spline(tolower(prot), ytd_g)
    use_ref <- !is.null(ref_p)
    col_map <- setNames(c(m$col, "grey45"), c(DISEASE_CAPS, "HC (PSM-matched)"))

    p <- ggplot(t_p, aes(x = years_to_dx)) +
        {if (use_ref) list(
            geom_ribbon(data = ref_p, aes(x = ytd, ymin = hc_ref_lci, ymax = hc_ref_uci),
                        inherit.aes = FALSE, fill = "grey70", alpha = 0.35),
            geom_line(data = ref_p, aes(x = ytd, y = hc_ref, colour = "HC (PSM-matched)"),
                      inherit.aes = FALSE, linewidth = 0.6, linetype = "dashed")
        ) else list(
            geom_hline(data = t_p[1], aes(yintercept = hc_mean, colour = "HC (PSM-matched)"),
                       linetype = "dashed", linewidth = 0.5)
        )} +
        geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), fill = m$col, alpha = 0.18) +
        geom_line(aes(y = pred_npx, colour = DISEASE_CAPS), linewidth = 0.9) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
        annotate("text", x = 0.3, y = Inf, label = "Dx",
                 hjust = 0, vjust = 1.5, size = 2.3, colour = "grey40") +
        scale_colour_manual(values = col_map, name = NULL) +
        scale_x_continuous(breaks = seq(-8, 12, by = 4)) +
        labs(
            title    = glue("{substr(m$pid, 1, 1)}  {prot} — {DISEASE_CAPS} disease course"),
            subtitle = m$note,
            x        = glue("Years relative to {DISEASE_CAPS} diagnosis"),
            y        = paste0(prot, " (NPX)")
        ) +
        theme_ukb(base_size = 9) +
        guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))

    save_panel(p, m$pid, 3.5, 3.5)
}

# ─────────────────────────────────────────────────────────────────────────────
# Panel d — ERBB2 PRS-stratified trajectory (counterpart to main Fig 2o NEFL)
#
# ERBB2 is a headline CNS marker in Figure 2 but is not tested in the existing
# `ms_prs_interaction_tests.csv` (which only contains the top-6 DEPs). We
# compute the PRS×ytd interaction test inline so this panel is self-contained.
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel d (ERBB2 PRS-stratified trajectory)...\n")

tryCatch({
    PRS_FILE_D <- file.path(PROJ_DIR, "data", "ukb", "genetics", glue("{COHORT}_prs_scores.csv"))
    if (!file.exists(PRS_FILE_D)) stop("PRS file not found: ", PRS_FILE_D)
    if (!"erbb2" %in% names(qc))  stop("ERBB2 not in QC columns")

    prs_d <- fread(PRS_FILE_D, showProgress = FALSE)
    if (PRS_COL != "prs_score" && PRS_COL %in% names(prs_d))
        setnames(prs_d, PRS_COL, "prs_score")
    ms_d  <- merge(
        qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE &
           !is.na(years_to_diagnosis) & !is.na(erbb2) &
           years_to_diagnosis >= -20 & years_to_diagnosis <= 16],
        prs_d[, .(eid, prs_score)], by = "eid"
    )
    cat(sprintf("  ERBB2 PRS cases: n=%d\n", nrow(ms_d)))

    ms_d[, prs_q := cut(prs_score,
                         breaks = quantile(prs_score, c(0, .25, .50, .75, 1), na.rm = TRUE),
                         labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
                         include.lowest = TRUE)]

    x_grid_d  <- seq(-20, 16, length.out = 100)

    traj_prs_d <- rbindlist(lapply(levels(ms_d$prs_q), function(q) {
        dq <- ms_d[prs_q == q & !is.na(erbb2) & !is.na(age_at_sampling)]
        if (nrow(dq) < 8) return(NULL)
        fit <- tryCatch(
            lm(erbb2 ~ ns(years_to_diagnosis, df = 3) + age_at_sampling + sex, data = dq),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        nd <- data.table(years_to_diagnosis = x_grid_d,
                         age_at_sampling    = median(dq$age_at_sampling, na.rm = TRUE),
                         sex                = as.numeric(names(sort(table(dq$sex),
                                                                    decreasing = TRUE))[1]))
        pr <- predict(fit, newdata = nd, se.fit = TRUE)
        data.table(ytd   = x_grid_d,
                   pred  = pr$fit,
                   lci   = pr$fit - 1.96 * pr$se.fit,
                   uci   = pr$fit + 1.96 * pr$se.fit,
                   prs_q = q)
    }))
    traj_prs_d[, prs_q := factor(prs_q, levels = levels(ms_d$prs_q))]

    # Inline PRS×time interaction test (nested ANOVA, mirrors
    # analysis/10_prs/01_ms_prs_protein_interaction.R)
    d_int <- ms_d[!is.na(erbb2),
                  .(npx = erbb2, ytd = years_to_diagnosis,
                    prs_q = as.integer(prs_q),
                    age = age_at_sampling, sex)]
    int_label_d <- ""
    if (nrow(d_int) >= 30) {
        fit0 <- tryCatch(lm(npx ~ ns(ytd, df = 3) + prs_q + age + sex,  data = d_int),
                         error = function(e) NULL)
        fit1 <- tryCatch(lm(npx ~ ns(ytd, df = 3) * prs_q + age + sex,  data = d_int),
                         error = function(e) NULL)
        if (!is.null(fit0) && !is.null(fit1)) {
            av <- tryCatch(anova(fit0, fit1), error = function(e) NULL)
            if (!is.null(av)) {
                pv <- av$`Pr(>F)`[2]
                int_label_d <- sprintf("\nPRS\u00d7time interaction p = %s",
                                       formatC(pv, format = "g", digits = 2))
                cat(sprintf("  ERBB2 PRS\u00d7time interaction p = %.4g (F=%.3f)\n",
                            pv, av$F[2]))
            }
        }
    }

    traj_ext_d <- traj_prs_d[prs_q %in% c("Q1 (low)", "Q4 (high)")]
    traj_mid_d <- traj_prs_d[prs_q %in% c("Q2", "Q3")]

    pD <- ggplot() +
        geom_line(data = traj_mid_d,
                  aes(x = ytd, y = pred, group = prs_q),
                  linewidth = 0.35, colour = "grey70", linetype = "dashed") +
        geom_ribbon(data = traj_ext_d,
                    aes(x = ytd, ymin = lci, ymax = uci, fill = prs_q),
                    alpha = 0.18, colour = NA) +
        geom_line(data = traj_ext_d,
                  aes(x = ytd, y = pred, colour = prs_q), linewidth = 0.95) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.35, colour = "grey40") +
        annotate("text", x = -0.3, y = Inf,
                 label = "\u2190 post-Dx  |  pre-Dx \u2192",
                 hjust = 0.5, vjust = 1.6, size = 1.9, colour = "grey50") +
        scale_colour_manual(values = c("Q1 (low)" = "#56B4E9",
                                        "Q4 (high)" = "#CC0066"),
                            name = "PRS\nquartile") +
        scale_fill_manual(values = c("Q1 (low)" = "#56B4E9",
                                      "Q4 (high)" = "#CC0066"),
                          guide = "none") +
        scale_x_continuous(breaks = seq(-16, 16, by = 4)) +
        labs(
            title    = paste0("d  ERBB2 trajectory by ", PRS_LABEL, int_label_d),
            subtitle = "Q1 (lowest) vs Q4 (highest) PRS shown with 95% CI | Q2/Q3 dashed | age+sex adjusted",
            x        = glue("Years relative to {DISEASE_CAPS} diagnosis"),
            y        = "ERBB2 (NPX)"
        ) +
        theme_ukb(base_size = 9) +
        guides(colour = guide_legend(override.aes = list(linewidth = 1.5)))

    save_panel(pD, "d_prs_erbb2_trajectory", 4, 3.8)

}, error = function(e) cat("  Panel d skipped:", conditionMessage(e), "\n"))

cat("\nSupplementary Figure 2 complete. Panels saved to:", FIG_DIR, "\n")
