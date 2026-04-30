#!/usr/bin/env Rscript
# 01_ms_protein_cox.R
# Per-protein Cox regression for MS incidence + rolling HR by time-to-diagnosis
#
# Design:
#   Population: pre-onset MS cases (n~127) + healthy controls (n~49,573)
#   Outcome:    MS diagnosis (event=1; time = years_to_diagnosis for cases,
#               16 = max follow-up for controls)
#   Model:      Surv(time, event) ~ scale(protein) + age_at_sampling + sex +
#               PC1 + PC2  (control-only PCA, all samples projected)
#
# Rolling HR:
#   Divide pre-onset cases into time-to-diagnosis bins:
#     0–4 yr, 4–8 yr, 8–12 yr, 12+ yr before diagnosis
#   For each bin + all controls, re-run Cox per protein
#   Plot HR trajectories across time windows (when does each protein "turn on"?)
#
# Output: results/survival/
#   ms_protein_cox_forest.pdf     — top 30 proteins, overall Cox HR
#   ms_protein_cox_rolling.pdf    — HR x time-window heatmap for top proteins
#   ms_protein_cox_results.csv    — full table of all proteins tested

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "limma_utils.R"))

QC_FILE   <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DIFF_FILE <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
OUT_DIR   <- file.path(PROJ_DIR, "results", "survival")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MAX_FOLLOW_UP <- 16   # years; max years_to_diagnosis in pre-onset cases
set.seed(42)

# ── 1. Prepare survival dataset ───────────────────────────────────────────────
cat("Preparing survival dataset...\n")
ms_qc <- fread(QC_FILE, showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]

# Compute control-only PC1/PC2 on the filtered cohort (matches 04_differential)
cat("Computing control-only PC1/PC2...\n")
META_COLS <- c("eid","ms_status","age_at_sampling","age_at_diagnosis",
               "years_to_diagnosis","sex","olink_instance","qc_outlier",
               "UMAP1","UMAP2","mean_npx")
protein_cols <- setdiff(names(ms_qc), META_COLS)
pc_dt <- compute_control_pcs(ms_qc, protein_cols)
ms_qc <- merge(ms_qc, pc_dt, by = "eid", all.x = TRUE)

# Restrict to pre-onset cases + controls only
surv_dt <- ms_qc[ms_status %in% c("pre_onset","control")]
surv_dt[, event    := as.integer(ms_status == "pre_onset")]
surv_dt[, surv_time := fcase(
    ms_status == "pre_onset", as.numeric(years_to_diagnosis),
    default = as.numeric(MAX_FOLLOW_UP)
)]
surv_dt[, sex_num := as.integer(sex)]   # UKB: 0=Female, 1=Male

cat(sprintf("  Cases: %d pre-onset, %d controls\n",
            sum(surv_dt$event), sum(!surv_dt$event)))
cat(sprintf("  Case time range: %.1f – %.1f yr\n",
            min(surv_dt[event==1, surv_time], na.rm=TRUE),
            max(surv_dt[event==1, surv_time], na.rm=TRUE)))

# ── 2. Protein list — all tested proteins ─────────────────────────────────────
deps_all  <- fread(DIFF_FILE, showProgress = FALSE)
all_prots <- tolower(deps_all$protein)
avail_prots <- intersect(all_prots, names(surv_dt))
cat(sprintf("  Proteins to test: %d\n", length(avail_prots)))

# ── 3. Overall Cox — all proteins ─────────────────────────────────────────────
cat("Running Cox regression for each protein...\n")

run_cox <- function(dt, prot) {
    sub <- dt[!is.na(get(prot))]
    if (sum(sub$event) < 5) return(NULL)
    form <- as.formula(sprintf(
        "Surv(surv_time, event) ~ scale(%s) + age_at_sampling + sex_num + PC1 + PC2",
        prot))
    fit <- tryCatch(coxph(form, data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    cf  <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(cf)) return(NULL)
    ci  <- tryCatch(confint(fit), error = function(e) NULL)
    prot_row <- grep(paste0("^scale\\(", prot, "\\)"), rownames(cf))
    if (length(prot_row) == 0) return(NULL)
    data.table(
        protein  = prot,
        HR       = exp(cf[prot_row, "coef"]),
        HR_lo    = if (!is.null(ci)) exp(ci[prot_row, 1]) else NA_real_,
        HR_hi    = if (!is.null(ci)) exp(ci[prot_row, 2]) else NA_real_,
        z        = cf[prot_row, "z"],
        pval     = cf[prot_row, "Pr(>|z|)"],
        n_events = sum(sub$event)
    )
}

cox_list <- lapply(avail_prots, run_cox, dt = surv_dt)
cox_dt   <- rbindlist(Filter(Negate(is.null), cox_list))
cox_dt[, fdr := p.adjust(pval, method = "BH")]
cox_dt   <- cox_dt[order(pval)]
fwrite(cox_dt, file.path(OUT_DIR, "ms_protein_cox_results.csv"))
cat(sprintf("  FDR<0.05: %d proteins\n", sum(cox_dt$fdr < 0.05)))
cat(sprintf("  FDR<0.20: %d proteins\n", sum(cox_dt$fdr < 0.20)))
cat("\nTop 15 proteins by HR (MS risk per 1 SD higher NPX):\n")
print(cox_dt[1:min(15,.N), .(protein=toupper(protein), HR=round(HR,3),
                              HR_lo=round(HR_lo,3), HR_hi=round(HR_hi,3),
                              pval=signif(pval,3), fdr=signif(fdr,3))])

# ── 4. Forest plot — top 30 proteins ─────────────────────────────────────────
sig_df <- cox_dt[1:min(30, .N)]
sig_df[, prot_f := factor(toupper(protein), levels = rev(toupper(protein)))]
sig_df[, direction := fifelse(HR > 1, "Higher → more MS risk", "Lower → less MS risk")]

DIR_COLS <- c("Higher → more MS risk" = "#CC0066", "Lower → less MS risk" = "#56B4E9")

p_forest <- ggplot(sig_df, aes(x = HR, y = prot_f, colour = direction)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
    geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = 0.3, linewidth = 0.4) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = DIR_COLS, name = NULL) +
    scale_x_log10() +
    labs(x = "Hazard ratio (per 1 SD, log scale)",
         y = NULL,
         title = "Per-protein Cox HR for MS incidence",
         subtitle = sprintf("Pre-onset cases n=%d vs controls n=%d | top 30 by p-value",
                            sum(surv_dt$event), sum(!surv_dt$event)),
         caption = "Adjusted for age, sex, PC1, PC2") +
    theme_ukb() +
    theme(axis.text.y = element_text(size = 7), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "ms_protein_cox_forest.pdf"), p_forest,
       width = 6.5, height = 8, device = cairo_pdf)
cat("  Saved: ms_protein_cox_forest.pdf\n")

# ── 5. Rolling HR — time bins ─────────────────────────────────────────────────
cat("\nRunning rolling HR by time-to-diagnosis bin...\n")

# Use top 20 proteins (by overall significance) for rolling analysis
top_prots <- cox_dt[1:min(20, .N), protein]

# Define bins (years BEFORE diagnosis)
time_bins <- list(
    "0–4 yr"   = c(0,  4),
    "4–8 yr"   = c(4,  8),
    "8–12 yr"  = c(8, 12),
    "12+ yr"   = c(12, MAX_FOLLOW_UP + 1)
)

roll_list <- list()
for (bin_name in names(time_bins)) {
    bin_range <- time_bins[[bin_name]]
    # Subset: cases in this time bin + all controls
    dt_bin <- surv_dt[event == 0 |
                      (event == 1 & surv_time >= bin_range[1] & surv_time < bin_range[2])]
    n_cases <- sum(dt_bin$event)
    cat(sprintf("  %s: %d cases\n", bin_name, n_cases))
    if (n_cases < 5) next

    for (prot in top_prots) {
        res <- run_cox(dt_bin, prot)
        if (!is.null(res)) {
            res[, time_bin := bin_name]
            res[, n_bin_cases := n_cases]
            roll_list[[paste(bin_name, prot)]] <- res
        }
    }
}

if (length(roll_list) > 0) {
    roll_dt <- rbindlist(roll_list, fill = TRUE)

    # Order bins and proteins
    roll_dt[, bin_f := factor(time_bin, levels = names(time_bins))]
    roll_dt[, prot_upper := toupper(protein)]
    prot_order <- cox_dt[protein %in% top_prots, toupper(protein)]
    roll_dt[, prot_f := factor(prot_upper, levels = rev(prot_order))]

    # Heatmap of HR (log scale)
    p_roll <- ggplot(roll_dt, aes(x = bin_f, y = prot_f, fill = log2(HR))) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(aes(label = ifelse(pval < 0.05, sprintf("%.2f", HR), "")),
                  size = 2.2, colour = "white") +
        scale_fill_gradient2(low = "#56B4E9", mid = "grey95", high = "#CC0066",
                             midpoint = 0, name = "log₂ HR",
                             limits = c(-0.8, 0.8), oob = scales::squish) +
        labs(x = "Years before MS diagnosis (pre-onset cases)",
             y = NULL,
             title = "Rolling hazard ratio by time window",
             subtitle = "HR per 1 SD protein (log₂ scale); labelled if p<0.05",
             caption = "Controls (n≈49k) included in each window") +
        theme_ukb() +
        theme(axis.text.y = element_text(size = 7),
              axis.text.x = element_text(angle = 30, hjust = 1),
              legend.position = "right")

    ggsave(file.path(OUT_DIR, "ms_protein_cox_rolling.pdf"), p_roll,
           width = 6.5, height = 7, device = cairo_pdf)
    cat("  Saved: ms_protein_cox_rolling.pdf\n")

    # Also: line plot for NEFL specifically
    nefl_roll <- roll_dt[protein == "nefl"]
    if (nrow(nefl_roll) > 0) {
        p_nefl <- ggplot(nefl_roll, aes(x = bin_f, y = HR, group = 1)) +
            geom_ribbon(aes(ymin = HR_lo, ymax = HR_hi), fill = "#CC0066", alpha = 0.2) +
            geom_line(colour = "#CC0066", linewidth = 0.8) +
            geom_point(colour = "#CC0066", size = 2.5) +
            geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
            labs(x = "Years before MS diagnosis",
                 y = "Hazard ratio (95% CI)",
                 title = "NEFL: rolling HR for MS incidence") +
            theme_ukb()
        ggsave(file.path(OUT_DIR, "ms_nefl_rolling_hr.pdf"), p_nefl,
               width = 5, height = 3.5, device = cairo_pdf)
        cat("  Saved: ms_nefl_rolling_hr.pdf\n")
    }

    fwrite(roll_dt[, .(protein=toupper(protein), time_bin, HR, HR_lo, HR_hi, pval, n_bin_cases)],
           file.path(OUT_DIR, "ms_protein_cox_rolling.csv"))
}

cat("\n01_ms_protein_cox.R complete.\n")
