#!/usr/bin/env Rscript
# =============================================================================
# 03_cns_trajectories.R — CNS biomarker trajectories: NEFL, MOG, GFAP, ERBB2
# Cross-sectional pseudo-longitudinal analysis replicating Abdelhak 2026 Fig. 3
# Approach: RCS regression within MS cases across years_to_diagnosis
# =============================================================================
# For each protein (NEFL, MOG, GFAP):
#   1. Fit restricted cubic spline (ns, df=3) on MS cases vs years_to_diagnosis,
#      adjusting for age_at_sampling, sex, PC1, PC2 (control-derived)
#   2. Predict over grid −12 to +6 years (0.25-yr intervals)
#   3. HC reference: age/sex-adjusted mean ± 1.96 SE
#   4. Identify earliest divergence year (predicted MS CI outside HC CI)
#   5. Plot trajectory + HC band + annotation; combine into 3-panel figure
#   6. Write per-protein PDFs, combined PDF, and summary CSV
#
# Input:
#   data/ukb/olink/processed/ms_olink_qc.csv
#
# Outputs:
#   results/differential/cns_trajectory_nefl.pdf
#   results/differential/cns_trajectory_mog.pdf
#   results/differential/cns_trajectory_gfap.pdf
#   results/differential/cns_trajectories_combined.pdf
#   results/differential/cns_trajectories.csv
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(splines)
    library(boot)
    library(patchwork)
})

# =============================================================================
# Paths
# =============================================================================
.args      <- commandArgs(trailingOnly = FALSE)
.file_arg  <- .args[grepl("^--file=", .args)]
SCRIPT_DIR <- if (length(.file_arg) > 0) dirname(normalizePath(sub("^--file=", "", .file_arg))) else getwd()
REPO_ROOT  <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
OUT_DIR    <- file.path(REPO_ROOT, "results", "differential")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

QC_FILE <- file.path(REPO_ROOT, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")

# =============================================================================
# Constants
# =============================================================================
PROTEINS   <- c("NEFL", "MOG", "GFAP", "ERBB2")   # target proteins (will resolve to actual col names)
# Sign convention (Abdelhak 2026 standard): negative = before onset, positive = after onset.
# Our years_to_diagnosis column uses the OPPOSITE sign (>0 = presymptomatic), so we negate
# inside analyse_protein() when constructing years_to_dx for the model and plot.
YR_MIN     <- -8    # 8 years before onset
YR_MAX     <- 12    # 12 years after onset
YR_STEP    <-  0.25
PRED_GRID  <- seq(YR_MIN, YR_MAX, by = YR_STEP)
CI_MULT    <- 1.96

MS_COLOR   <- "#2B4C7E"   # UKB_NAVY for MS predicted line
HC_COLOR   <- "grey60"    # grey for HC reference band

set.seed(42)

source(file.path(REPO_ROOT, "analysis", "helpers", "ukb_theme.R"))

# =============================================================================
# 1. Load QC'd data
# =============================================================================
message("Loading QC'd Olink + MS cohort data...")
if (!file.exists(QC_FILE)) stop("QC file not found: ", QC_FILE, "\nRun 02_olink_qc/01_olink_qc.R first.")
dt <- fread(QC_FILE, showProgress = FALSE)

META_COLS <- c("eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance",
               "qc_outlier", "UMAP1", "UMAP2", "mean_npx",
               "PC1", "PC2")
protein_cols <- setdiff(names(dt), META_COLS)

# Source control-derived PC helper (matches 01_limma_ms_vs_hc.R pattern)
source(file.path(REPO_ROOT, "analysis", "helpers", "limma_utils.R"))
message(sprintf("  %d participants, %d proteins", nrow(dt), length(protein_cols)))
message(sprintf("  Status counts: %s",
                paste(names(table(dt$ms_status)),
                      as.integer(table(dt$ms_status)),
                      sep = "=", collapse = "  ")))

# =============================================================================
# 2. Helper: resolve protein column name (case-insensitive fallback)
# =============================================================================
resolve_col <- function(target, available) {
    # Try exact match first (lowercase convention used in UKB Olink data)
    exact <- target
    if (exact %in% available) return(exact)
    # Try lowercase
    lo <- tolower(target)
    if (lo %in% available) return(lo)
    # Try uppercase
    up <- toupper(target)
    if (up %in% available) return(up)
    # Try title case
    ti <- paste0(toupper(substr(target, 1, 1)), tolower(substr(target, 2, nchar(target))))
    if (ti %in% available) return(ti)
    stop(sprintf("Column for protein '%s' not found. Tried: %s, %s, %s, %s",
                 target, exact, lo, up, ti))
}

# =============================================================================
# 3. Master filter: non-outliers, valid covariates, relevant statuses
# =============================================================================
message("\nApplying master QC filter...")
dt_filt <- dt[qc_outlier == FALSE &
              ms_status %in% c("pre_onset", "post_onset", "control") &
              !is.na(age_at_sampling) &
              !is.na(sex)]

# Compute control-only PC1/PC2 and merge onto dt_filt.
# Rationale: PCs derived from HC protein matrix capture technical variance
# without disease-signal contamination; used as covariates in the RCS model.
message("Computing control-only PC1/PC2...")
for (col in c("PC1", "PC2")) if (col %in% names(dt_filt)) set(dt_filt, j = col, value = NULL)
pc_dt   <- compute_control_pcs(dt_filt, protein_cols)
dt_filt <- merge(dt_filt, pc_dt, by = "eid", all.x = TRUE)
dt_filt <- dt_filt[!is.na(PC1) & !is.na(PC2)]
message(sprintf("  After filter: %d participants", nrow(dt_filt)))
message(sprintf("  Status: %s",
                paste(names(table(dt_filt$ms_status)),
                      as.integer(table(dt_filt$ms_status)),
                      sep = "=", collapse = "  ")))

# Factor encoding (sex) — used consistently below
dt_filt[, sex_f := factor(sex)]

# =============================================================================
# 4. Per-protein analysis function
# =============================================================================
analyse_protein <- function(prot_name, dt_in) {

    message(sprintf("\n--- Protein: %s ---", prot_name))

    # Resolve actual column name
    col_nm <- resolve_col(prot_name, names(dt_in))
    message(sprintf("  Column resolved: '%s'", col_nm))

    # ------------------------------------------------------------------
    # 4a. MS subset: must have years_to_diagnosis, restricted to window.
    # Negate so that years_to_dx < 0 = presymptomatic (before onset),
    # years_to_dx > 0 = post-onset (after onset) — Abdelhak 2026 convention.
    # ------------------------------------------------------------------
    ms_sub <- dt_in[!is.na(years_to_diagnosis) &
                    -years_to_diagnosis >= YR_MIN &
                    -years_to_diagnosis <= YR_MAX &
                    !is.na(get(col_nm))]
    n_ms <- nrow(ms_sub)
    message(sprintf("  MS cases in model: %d (years_to_dx range: %.1f to %.1f before/after onset)",
                    n_ms,
                    min(-ms_sub$years_to_diagnosis),
                    max(-ms_sub$years_to_diagnosis)))
    if (n_ms < 30) warning(sprintf("  Very few MS cases for %s (n=%d)", prot_name, n_ms))

    # ------------------------------------------------------------------
    # 4b. HC subset: controls with valid PCs (already filtered above)
    # ------------------------------------------------------------------
    hc_sub <- dt_in[ms_status == "control" & !is.na(get(col_nm))]
    n_hc   <- nrow(hc_sub)
    message(sprintf("  HC in reference: %d", n_hc))

    # ------------------------------------------------------------------
    # 4c. Fit RCS model on MS cases
    # ------------------------------------------------------------------
    ms_df <- as.data.frame(ms_sub[, .(
        npx             = get(col_nm),
        years_to_dx     = -years_to_diagnosis,   # negate: negative=pre-onset, positive=post-onset
        age_at_sampling = age_at_sampling,
        sex_f           = sex_f,
        PC1             = PC1,
        PC2             = PC2
    )])

    fit_ms <- lm(
        npx ~ ns(years_to_dx, df = 3) + age_at_sampling + sex_f + PC1 + PC2,
        data = ms_df
    )
    message(sprintf("  RCS model R2 = %.3f, RMSE = %.3f",
                    summary(fit_ms)$r.squared,
                    sqrt(mean(residuals(fit_ms)^2))))

    # ------------------------------------------------------------------
    # 4d. Prediction grid (hold covariates at MS-case medians)
    # ------------------------------------------------------------------
    # Majority sex level for prediction
    majority_sex <- names(which.max(table(ms_df$sex_f)))

    newdata_grid <- data.frame(
        years_to_dx     = PRED_GRID,
        age_at_sampling = median(ms_df$age_at_sampling, na.rm = TRUE),
        sex_f           = factor(majority_sex, levels = levels(ms_df$sex_f)),
        PC1             = median(ms_df$PC1, na.rm = TRUE),
        PC2             = median(ms_df$PC2, na.rm = TRUE)
    )

    pred_out <- predict(fit_ms, newdata = newdata_grid, interval = "confidence", level = 0.95)
    pred_df  <- data.frame(
        protein     = prot_name,
        years_to_dx = PRED_GRID,
        pred_npx    = pred_out[, "fit"],
        pred_lci    = pred_out[, "lwr"],
        pred_uci    = pred_out[, "upr"]
    )

    # ------------------------------------------------------------------
    # 4e. HC reference: age/sex-adjusted residuals + grand mean
    # ------------------------------------------------------------------
    hc_df <- as.data.frame(hc_sub[, .(
        npx             = get(col_nm),
        age_at_sampling = age_at_sampling,
        sex_f           = sex_f
    )])

    # Adjust HC NPX for age and sex, preserve scale via grand mean
    fit_hc_adj <- lm(npx ~ age_at_sampling + sex_f, data = hc_df)
    hc_resid   <- residuals(fit_hc_adj) + mean(hc_df$npx, na.rm = TRUE)

    hc_mean <- mean(hc_resid)
    hc_sd   <- sd(hc_resid)
    hc_se   <- hc_sd / sqrt(length(hc_resid))

    # Display band: ±1 SD (shows the population distribution, not CI of mean)
    # With n~45K, ±1.96 SE is ~0.005 wide — essentially invisible on the plot.
    hc_disp_lci <- hc_mean - hc_sd
    hc_disp_uci <- hc_mean + hc_sd

    # Divergence detection band: ±1.96 SE (statistical significance vs HC mean)
    hc_lci <- hc_mean - CI_MULT * hc_se
    hc_uci <- hc_mean + CI_MULT * hc_se

    message(sprintf("  HC reference: mean=%.3f, SD=%.3f, SE=%.4f | display ±1SD [%.3f, %.3f]",
                    hc_mean, hc_sd, hc_se, hc_disp_lci, hc_disp_uci))

    # Add HC stats to prediction df (display values for ribbon; SE-based stored separately)
    pred_df$hc_mean <- hc_mean
    pred_df$hc_lci  <- hc_disp_lci   # ±1 SD for CSV/display
    pred_df$hc_uci  <- hc_disp_uci
    pred_df$hc_lci_se <- hc_lci      # ±1.96 SE for divergence reference
    pred_df$hc_uci_se <- hc_uci

    # ------------------------------------------------------------------
    # 4f. Find earliest divergence timepoint
    #     Elevation: MS lower CI > HC upper CI
    #     Suppression: MS upper CI < HC lower CI
    # ------------------------------------------------------------------
    elevated   <- pred_df$pred_lci > pred_df$hc_uci
    suppressed <- pred_df$pred_uci < pred_df$hc_lci

    diverge_yr   <- NA_real_
    diverge_dir  <- NA_character_

    # Sign convention: years_to_dx < 0 = presymptomatic (before onset), > 0 = post-onset.
    # PRED_GRID runs from YR_MIN (most presymptomatic) to YR_MAX (most post-onset).
    # Use SE-based HC CI for divergence test (statistical significance).
    # "First divergence" = most presymptomatic point = smallest (most negative) index where
    # the condition holds → which(...)[1] gives the first match scanning from YR_MIN.

    elevated   <- pred_df$pred_lci > pred_df$hc_uci_se
    suppressed <- pred_df$pred_uci < pred_df$hc_lci_se

    recovery_yr  <- NA_real_

    if (any(elevated, na.rm = TRUE)) {
        diverge_yr  <- pred_df$years_to_dx[which(elevated)[1]]
        diverge_dir <- "elevated"
        if (diverge_yr < 0) {
            message(sprintf("  First presymptomatic elevation above HC: %.2f yrs before onset", abs(diverge_yr)))
        } else {
            message(sprintf("  Elevation above HC starts %.2f yrs after onset (not presymptomatic)", diverge_yr))
        }
    } else if (any(suppressed, na.rm = TRUE)) {
        diverge_yr  <- pred_df$years_to_dx[which(suppressed)[1]]
        diverge_dir <- "suppressed"
        message(sprintf("  First suppression below HC: %.2f yrs (negative = before onset)", diverge_yr))
    } else {
        message("  No significant divergence from HC detected in the prediction window.")
    }

    # Recovery: first post-onset grid point where the divergence condition
    # no longer holds. Only meaningful when pre-onset divergence is present.
    if (!is.na(diverge_yr) && diverge_yr < 0) {
        post_idx <- which(pred_df$years_to_dx > 0)
        cond_post <- if (diverge_dir == "elevated") elevated[post_idx] else suppressed[post_idx]
        first_end <- which(!cond_post)[1]
        if (!is.na(first_end)) {
            recovery_yr <- pred_df$years_to_dx[post_idx][first_end]
            message(sprintf("  Recovery toward HC: %.2f yrs after onset", recovery_yr))
        }
    }

    # ------------------------------------------------------------------
    # Return all results
    # ------------------------------------------------------------------
    list(
        pred_df      = pred_df,
        fit_ms       = fit_ms,
        hc_mean      = hc_mean,
        hc_se        = hc_se,
        hc_lci       = hc_lci,
        hc_uci       = hc_uci,
        diverge_yr   = diverge_yr,
        diverge_dir  = diverge_dir,
        recovery_yr  = recovery_yr,
        n_ms         = n_ms,
        n_hc         = n_hc
    )
}

# =============================================================================
# 5. Plot function
# =============================================================================
make_traj_plot <- function(res, prot_name) {

    pd  <- res$pred_df
    dyr <- res$diverge_yr
    ddir <- res$diverge_dir

    # y-axis range: encompass both MS CI and HC CI with a little padding
    y_lo <- min(pd$pred_lci, pd$hc_lci, na.rm = TRUE)
    y_hi <- max(pd$pred_uci, pd$hc_uci, na.rm = TRUE)
    y_pad <- (y_hi - y_lo) * 0.12
    y_lo <- y_lo - y_pad
    y_hi <- y_hi + y_pad

    p <- ggplot(pd, aes(x = years_to_dx)) +

        # --- HC reference band ---
        geom_ribbon(aes(ymin = hc_lci, ymax = hc_uci),
                    fill  = HC_COLOR, alpha = 0.15) +
        geom_hline(aes(yintercept = hc_mean),
                   colour   = HC_COLOR,
                   linetype = "dashed",
                   linewidth = 0.6) +

        # --- MS trajectory ---
        geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
                    fill  = MS_COLOR, alpha = 0.20) +
        geom_line(aes(y = pred_npx),
                  colour    = MS_COLOR,
                  linewidth = 0.9) +

        # --- Diagnosis line ---
        geom_vline(xintercept = 0,
                   linetype   = "dashed",
                   colour     = "grey30",
                   linewidth  = 0.5) +
        annotate("text",
                 x      = 0.15,
                 y      = y_hi - y_pad * 0.5,
                 label  = "Diagnosis",
                 hjust  = 0,
                 vjust  = 1,
                 size   = 3,
                 colour = "grey30") +

        # --- Axes and labels ---
        scale_x_continuous(
            limits = c(YR_MIN, YR_MAX),
            breaks = seq(YR_MIN, YR_MAX, by = 2)
        ) +
        coord_cartesian(ylim = c(y_lo, y_hi)) +
        labs(
            title = prot_name,
            x     = "Years relative to MS diagnosis (negative = before onset)",
            y     = "NPX (normalized)"
        ) +
        theme_ukb()

    # --- Annotation: earliest divergence ---
    if (!is.na(dyr)) {
        arrow_y_top  <- y_lo + (y_hi - y_lo) * 0.18
        arrow_y_tip  <- y_lo + (y_hi - y_lo) * 0.06
        label_text   <- if (ddir == "elevated") "First above HC" else "First below HC"

        p <- p +
            annotate("segment",
                     x    = dyr, xend = dyr,
                     y    = arrow_y_top, yend = arrow_y_tip,
                     arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
                     colour = "#CC0066",
                     linewidth = 0.6) +
            annotate("text",
                     x      = dyr,
                     y      = arrow_y_top + (y_hi - y_lo) * 0.04,
                     label  = sprintf("%s\n(yr %.1f)", label_text, dyr),
                     hjust  = 0.5,
                     vjust  = 0,
                     size   = 2.5,
                     colour = "#CC0066")
    }

    # --- Sample size footnote in subtitle ---
    p <- p + labs(
        subtitle = sprintf("MS n=%d  |  HC n=%d", res$n_ms, res$n_hc)
    )

    p
}

# =============================================================================
# 6. Run analyses for each protein
# =============================================================================
message("\n=== Running per-protein trajectory analyses ===")

results_list <- lapply(PROTEINS, function(prot) {
    tryCatch(
        analyse_protein(prot, dt_filt),
        error = function(e) {
            message(sprintf("  ERROR for %s: %s", prot, conditionMessage(e)))
            NULL
        }
    )
})
names(results_list) <- PROTEINS

# =============================================================================
# 7. Save per-protein PDFs
# =============================================================================
message("\n=== Saving per-protein PDFs ===")

plot_list <- vector("list", length(PROTEINS))
names(plot_list) <- PROTEINS

for (prot in PROTEINS) {
    res <- results_list[[prot]]
    if (is.null(res)) {
        message(sprintf("  Skipping %s (analysis failed)", prot))
        next
    }

    p   <- make_traj_plot(res, prot)
    plot_list[[prot]] <- p

    out_pdf <- file.path(OUT_DIR, sprintf("cns_trajectory_%s.pdf", tolower(prot)))
    ggsave(out_pdf, p, width = 6, height = 4.5, device = cairo_pdf)
    message(sprintf("  Saved: %s", basename(out_pdf)))
}

# =============================================================================
# 8. Combined 3-panel figure
# =============================================================================
message("\n=== Building combined 3-panel figure ===")

valid_plots <- Filter(Negate(is.null), plot_list)

if (length(valid_plots) >= 1) {
    # Remove individual x-axis labels for cleaner panel, keep on bottom row only
    # (patchwork handles this via collect_guides or axis labelling)
    combined <- patchwork::wrap_plots(valid_plots, nrow = 1) +
        patchwork::plot_annotation(
            title   = "CNS biomarker trajectories relative to MS diagnosis",
            caption = "RCS spline (df=3); blue ribbon = 95% prediction CI (MS cases). Grey band = HC mean +/- 1 SD (age/sex adjusted, n~45K). Arrow = earliest presymptomatic divergence.",
            theme   = theme(
                plot.title   = element_text(size = 12, face = "plain"),
                plot.caption = element_text(size = 8, colour = "grey40")
            )
        ) &
        theme_ukb()

    out_combined <- file.path(OUT_DIR, "cns_trajectories_combined.pdf")
    ggsave(out_combined,
           combined,
           width  = 6 * length(valid_plots),
           height = 5,
           device = cairo_pdf)
    message(sprintf("  Saved combined figure: %s", basename(out_combined)))
} else {
    message("  No valid plots — combined figure not written.")
}

# =============================================================================
# 9. Write CSV: predicted values + HC stats for all proteins
# =============================================================================
message("\n=== Writing summary CSV ===")

csv_list <- lapply(PROTEINS, function(prot) {
    res <- results_list[[prot]]
    if (is.null(res)) return(NULL)
    res$pred_df[, c("protein", "years_to_dx",
                    "pred_npx", "pred_lci", "pred_uci",
                    "hc_mean", "hc_lci", "hc_uci")]
})
csv_list <- Filter(Negate(is.null), csv_list)

if (length(csv_list) > 0) {
    out_csv <- rbindlist(lapply(csv_list, as.data.table))
    fwrite(out_csv, file.path(OUT_DIR, "cns_trajectories.csv"))
    message(sprintf("  Saved: cns_trajectories.csv  (%d rows)", nrow(out_csv)))
} else {
    message("  No data to write to CSV.")
}

# =============================================================================
# 10. Summary of divergence findings
# =============================================================================
message("\n=== Divergence summary ===")
for (prot in PROTEINS) {
    res <- results_list[[prot]]
    if (is.null(res)) {
        message(sprintf("  %s: FAILED", prot))
    } else if (is.na(res$diverge_yr)) {
        message(sprintf("  %s: no significant divergence from HC in [%.0f, %.0f] years",
                        prot, YR_MIN, YR_MAX))
    } else {
        message(sprintf("  %s: first %s vs HC at %.2f years relative to diagnosis",
                        prot, res$diverge_dir, res$diverge_yr))
    }
}

summary_rows <- lapply(PROTEINS, function(prot) {
    res <- results_list[[prot]]
    if (is.null(res)) return(NULL)
    data.table(
        protein     = prot,
        diverge_yr  = res$diverge_yr,
        diverge_dir = res$diverge_dir,
        recovery_yr = res$recovery_yr,
        n_ms        = res$n_ms,
        n_hc        = res$n_hc
    )
})
summary_rows <- Filter(Negate(is.null), summary_rows)
if (length(summary_rows) > 0) {
    fwrite(rbindlist(summary_rows),
           file.path(OUT_DIR, "cns_trajectories_divergence.csv"))
    message("  Saved: cns_trajectories_divergence.csv")
}

message("\n03_cns_trajectories.R complete.")
