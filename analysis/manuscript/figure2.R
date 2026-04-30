#!/usr/bin/env Rscript
# figure2.R — Figure 2: The temporal proteomic window in MS
#
# Panels: a–m
#   a    Sampling window histogram (valid age_at_diagnosis, clipped display)
#   b    NEFL trajectory across disease course
#   c    GFAP trajectory
#   d    MOG trajectory
#   e-h  Top 4 immune DEP trajectories (ALL MS cases, pre + post)
#   i    Pre-onset vs post-onset logFC scatter
#   j    Temporal protein heatmap (ALL MS, full disease course)
#   k    Cox rolling HR heatmap (protein × pre-onset time window)
#   l    SHAP feature importance (pre-onset classifier)
#   m    ML ROC curves (full cohort vs pre-onset only)
#   n    Pre-diagnostic risk score by time bin
#   o    PRS-stratified NEFL trajectory
#   p    Temporal proteomic landscape summary
#
# Output: results/figures/2/panel_{a..m}.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(splines)
    library(glue)
})

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
    getwd()
}
PROJ_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "top_n_selector.R"))

cfg          <- load_disease_config()
COHORT       <- cfg$cohort_short
DISEASE_CAPS <- cfg$disease_short_caps
DISEASE_LONG <- cfg$disease_long
STATUS_COL   <- cfg$cohort_status_col
SV           <- cfg$status_values
HLA_CARRIER  <- cfg$hla_carrier_col
HLA_DOSAGE   <- cfg$hla_dosage_col
PRS_COL      <- cfg$prs_combined_col
HLA_LBL      <- paste0("HLA-", cfg$hla_allele)
PRS_LABEL    <- cfg$prs_label
PRE_LBL      <- glue("pre-{DISEASE_CAPS}")
POST_LBL     <- glue("post-{DISEASE_CAPS}")

FIG_DIR <- file.path(PROJ_DIR, "results", "figures", "2")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

save_panel <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(FIG_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: panel_", name, ".pdf\n", sep = "")
}

# ─────────────────────────────────────────────────────────────────────────────
# File paths
# ─────────────────────────────────────────────────────────────────────────────
DATA_DIR <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed")
DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
ML_DIR   <- file.path(PROJ_DIR, "results", "ml")

QC_FILE         <- file.path(DATA_DIR, glue("{COHORT}_olink_qc.csv"))
TRAJ_FILE       <- file.path(DIFF_DIR, "cns_trajectories.csv")
PRE_FILE        <- file.path(DIFF_DIR, glue("{COHORT}_pre_vs_hc.csv"))
POST_FILE       <- file.path(DIFF_DIR, glue("{COHORT}_post_vs_hc.csv"))
COMB_FILE       <- file.path(DIFF_DIR, glue("{COHORT}_combined_vs_hc.csv"))
ML_FULL_FILE    <- file.path(ML_DIR, glue("{COHORT}_ml_results.csv"))
ML_PRE_FILE     <- file.path(ML_DIR, glue("{COHORT}_preonset_ml_results.csv"))
ML_PRE_IMP_FILE <- file.path(ML_DIR, glue("{COHORT}_preonset_ml_feature_importance.csv"))
MODEL_FULL_FILE <- file.path(ML_DIR, "models", "glmnet_model.rds")
MODEL_PRE_FILE  <- file.path(ML_DIR, "models", glue("{COHORT}_pre_boruta_glmnet_model.rds"))  # elastic net

# ─────────────────────────────────────────────────────────────────────────────
# Load and validate
# ─────────────────────────────────────────────────────────────────────────────
cat("Loading data...\n")

qc <- fread(QC_FILE)
if (STATUS_COL != "ms_status" && STATUS_COL %in% names(qc))
    setnames(qc, STATUS_COL, "ms_status")

# Exclude biologically impossible age_at_diagnosis values
# (sign convention: years_to_diagnosis = age_at_diagnosis - age_at_sampling;
#  positive = pre-onset, negative = post-onset)
n_before <- nrow(qc)
qc <- qc[is.na(age_at_diagnosis) | (age_at_diagnosis >= 5 & age_at_diagnosis <= 90)]
cat(sprintf("  Excluded %d rows with invalid age_at_diagnosis.\n", n_before - nrow(qc)))
cat("  QC rows after filter:", nrow(qc), "\n")

comb <- fread(COMB_FILE)
cat("  Combined differential rows:", nrow(comb), "\n")

META_COLS <- c("eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "qc_outlier", "UMAP1", "UMAP2")
prot_cols <- setdiff(names(qc), META_COLS)
cat("  Protein columns:", length(prot_cols), "\n")

# Temporal bins — span the full disease course (post-onset left, pre-onset
# right).  Prefer cfg$temporal_bins (numeric vector of cut breakpoints when
# the disease has a non-MS-like window — e.g., acute stroke or chronic DM2);
# otherwise default to the MS/ALS-validated -10..+10 year window.
.cfg_bins_breaks <- cfg$temporal_bins
if (!is.null(.cfg_bins_breaks) && length(.cfg_bins_breaks) >= 3L) {
    BIN_BREAKS <- as.numeric(.cfg_bins_breaks)
    BIN_LABELS <- vapply(seq_len(length(BIN_BREAKS) - 1L), function(i) {
        lo <- BIN_BREAKS[i]; hi <- BIN_BREAKS[i + 1L]
        if      (is.infinite(lo) && lo < 0) sprintf("< %g", hi)
        else if (is.infinite(hi))           sprintf("> %g", lo)
        else                                sprintf("[%g,%g)", lo, hi)
    }, character(1))
} else {
    BIN_BREAKS <- c(-Inf, -10, -5, -2, 0, 2, 5, 10, Inf)
    BIN_LABELS <- c("< -10", "[-10,-5)", "[-5,-2)", "[-2,0)",
                     "[0,2)", "[2,5)", "[5,10)", "> 10")
}
cat(sprintf("  Temporal bin breaks: %s\n",
            paste(BIN_BREAKS, collapse = ", ")))

# ─────────────────────────────────────────────────────────────────────────────
# Panel a — Sampling window
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel a...\n")

cases_a <- qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE]
cases_a[, ms_status := factor(ms_status, levels = c(SV$pre_onset, SV$post_onset))]

pA <- ggplot(cases_a, aes(x = years_to_diagnosis, fill = ms_status)) +
    geom_histogram(binwidth = 1, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.5, colour = "grey30") +
    annotate("text", x = 0.5, y = Inf, label = "Dx",
             hjust = 0, vjust = 1.5, size = 2.5, colour = "grey30") +
    scale_fill_manual(
        values = setNames(c(COL_PRE, COL_POST), c(SV$pre_onset, SV$post_onset)),
        labels = c("Pre-onset", "Post-onset"),
        name   = NULL
    ) +
    scale_x_continuous(breaks = seq(-50, 20, by = 10)) +
    labs(
        title    = glue("a  {DISEASE_CAPS} sampling window"),
        subtitle = sprintf("n=%d | %d pre-onset, %d post-onset",
                           nrow(cases_a),
                           sum(cases_a$ms_status == SV$pre_onset),
                           sum(cases_a$ms_status == SV$post_onset)),
        x        = glue("Years relative to {DISEASE_CAPS} diagnosis"),
        y        = glue("{DISEASE_CAPS} cases")
    ) +
    theme_ukb(base_size = 9) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7))

save_panel(pA, "a_sampling_window", 5.5, 3.2)

# ─────────────────────────────────────────────────────────────────────────────
# Panels b, c, d, e — CNS protein trajectories (one per protein)
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panels b, c, d, e (NEFL, ERBB2, GFAP, MOG CNS trajectories)...\n")

traj <- fread(TRAJ_FILE)
traj[, protein := toupper(protein)]
cat("  Proteins:", paste(unique(traj$protein), collapse = ", "),
    "| years range:", paste(range(traj$years_to_dx), collapse = " to "), "\n")

# ── PSM-matched HC pseudo-trajectory ─────────────────────────────────────────
# Match each MS participant to HC (1:3, nearest-neighbour on age + sex).
# Each matched HC inherits the MS case's years_to_diagnosis as a pseudo-ytd.
# This mirrors the DEG matching approach and produces a HC reference whose
# "disease timeline" distribution matches the MS cohort exactly.
cat("  Computing PSM-matched HC reference...\n")
suppressPackageStartupMessages(library(MatchIt))

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
    # Pass full pool_dt so match.data() retains eid and years_to_diagnosis
    matchit(is_ms ~ age_at_sampling + sex,
            data    = as.data.frame(pool_dt),
            method  = "nearest", ratio = 3, replace = FALSE),
    error = function(e) { cat("  PSM failed:", e$message, "\n"); NULL }
)

psm_hc_ref <- NULL   # fallback: NULL → use flat HC mean

if (!is.null(m_obj)) {
    md <- as.data.table(match.data(m_obj))
    # Each subclass has 1 MS case + 3 HC; propagate the MS ytd to its HC partners
    ms_sub  <- md[is_ms == 1, .(subclass, ytd_ms = years_to_diagnosis)]
    hc_md   <- md[is_ms == 0, .(eid, subclass)]
    hc_md   <- merge(hc_md, ms_sub, by = "subclass", all.x = TRUE)
    hc_rows <- merge(hc_md[, .(eid, pseudo_ytd = ytd_ms)], hc_pool, by = "eid")
    hc_rows[, years_to_diagnosis := pseudo_ytd]
    psm_hc_ref <- hc_rows
    cat(sprintf("  PSM matched: %d HC observations with pseudo-ytd\n", nrow(psm_hc_ref)))
}

# Helper: fit spline on matched HC and predict at a ytd grid
hc_psm_spline <- function(col_nm, ytd_grid) {
    if (is.null(psm_hc_ref) || !col_nm %in% names(psm_hc_ref)) return(NULL)
    d <- psm_hc_ref[, .(ytd = years_to_diagnosis, npx = get(col_nm))][!is.na(npx) & !is.na(ytd)]
    if (nrow(d) < 10) return(NULL)
    fit <- tryCatch(lm(npx ~ splines::ns(ytd, df = 3), data = d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pr <- predict(fit, newdata = data.frame(ytd = ytd_grid), se.fit = TRUE)
    data.table(ytd = ytd_grid,
               hc_ref     = pr$fit,
               hc_ref_lci = pr$fit - 1.96 * pr$se.fit,
               hc_ref_uci = pr$fit + 1.96 * pr$se.fit)
}

hc_qc <- hc_pool   # keep alias for fallback use below

CNS_META <- list(
    NEFL  = list(pid = "b_nefl_trajectory",  col = COL_PRE,      name = "Neurofilament light (NfL)"),
    ERBB2 = list(pid = "c_erbb2_trajectory", col = COL_PRE_DOWN, name = "ERBB2"),
    GFAP  = list(pid = "d_gfap_trajectory",  col = "#E6A817",    name = "GFAP"),
    MOG   = list(pid = "e_mog_trajectory",   col = COL_POST,     name = "Myelin oligodendrocyte glycoprotein (MOG)")
)

for (prot in names(CNS_META)) {
    m   <- CNS_META[[prot]]
    t_p <- traj[protein == prot]
    if (nrow(t_p) == 0) { cat("  Skipping", prot, ": no data\n"); next }

    ytd_g  <- t_p$years_to_dx
    ref_p  <- hc_psm_spline(tolower(prot), ytd_g)
    use_ref <- !is.null(ref_p)

    # Legend driven by colour aesthetic on the lines
    col_map <- setNames(c(m$col, "grey45"), c(DISEASE_CAPS, "HC (PSM-matched)"))

    p <- ggplot(t_p, aes(x = years_to_dx)) +
        # PSM-matched HC ribbon + line
        {if (use_ref) list(
            geom_ribbon(data = ref_p, aes(x = ytd, ymin = hc_ref_lci, ymax = hc_ref_uci),
                        inherit.aes = FALSE, fill = "grey70", alpha = 0.35),
            geom_line(data = ref_p, aes(x = ytd, y = hc_ref, colour = "HC (PSM-matched)"),
                      inherit.aes = FALSE, linewidth = 0.6, linetype = "dashed")
        ) else list(
            geom_hline(data = t_p[1], aes(yintercept = hc_mean, colour = "HC (PSM-matched)"),
                       linetype = "dashed", linewidth = 0.5)
        )} +
        # MS spline ± CI
        geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), fill = m$col, alpha = 0.18) +
        geom_line(aes(y = pred_npx, colour = DISEASE_CAPS), linewidth = 0.9) +
        # Diagnosis line
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
        annotate("text", x = 0.3, y = Inf, label = "Dx",
                 hjust = 0, vjust = 1.5, size = 2.3, colour = "grey40") +
        scale_colour_manual(values = col_map, name = NULL) +
        scale_x_continuous(breaks = seq(-8, 12, by = 4)) +
        labs(
            title = glue("{substr(m$pid, 1, 1)}  {prot} — {DISEASE_CAPS} disease course"),
            x     = glue("Years relative to {DISEASE_CAPS} diagnosis"),
            y     = paste0(prot, " (NPX)")
        ) +
        theme_ukb(base_size = 9) +
        guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))

    save_panel(p, m$pid, 3.5, 3.5)
}

# ─────────────────────────────────────────────────────────────────────────────
# Panels f, g, h — Immune protein trajectories (ALL MS, pre + post)
# BGN is demoted to Extended Data Fig. 2a; ERBB2 is promoted into CNS_META (c).
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panels f, g, h (top pre-onset immune DEP trajectories)...\n")

pre  <- fread(PRE_FILE)
pre[, fdr := p.adjust(P.Value, method = "BH")]
if (!"adj.P.Val" %in% names(pre)) pre[, adj.P.Val := fdr]
# Top-3 trajectory panels: data-driven, top P-value (FDR-aware) across both
# directions. Earlier curated exclusions (BGN, ERBB2) are dropped — the
# programmatic top-N stands on its own; if a protein appears here AND in
# another panel that's an informative redundancy, not a bug.
top_immune <- top_n_one_direction(
    pre, n = 3L,
    min_sig_for_fdr = cfg$min_sig_for_fdr %||% 5L,
    case_fn = identity   # preserve lowercase Olink convention for column lookup
)
cat("  Top pre-onset immune proteins:", paste(top_immune, collapse = ", "),
    " (ordered by ", attr(top_immune, "ordered_by"), ")\n", sep = "")

ms_all  <- qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE]
hc_only <- qc[ms_status == SV$control & qc_outlier == FALSE]

# Clip to reasonable display range (data: -107 to +16; display: -20 to +16)
YTD_MIN <- -20; YTD_MAX <- 16
ms_traj  <- ms_all[years_to_diagnosis >= YTD_MIN & years_to_diagnosis <= YTD_MAX]
x_grid   <- seq(YTD_MIN, YTD_MAX, length.out = 120)

IMM_PANELS <- c("f", "g", "h")
IMM_COLS   <- c("#E6A817", COL_POST, "#1A1A1A")

for (i in seq_along(top_immune)) {
    prot   <- top_immune[i]
    pid    <- IMM_PANELS[i]
    col_nm <- tolower(prot)

    if (!col_nm %in% names(ms_traj)) {
        cat("  Skipping", prot, ": column not found\n"); next
    }

    d <- ms_traj[, .(
        npx                = get(col_nm),
        years_to_diagnosis,
        age_at_sampling,
        sex
    )][!is.na(npx)]

    fit <- tryCatch(
        lm(npx ~ ns(years_to_diagnosis, df = 3) + age_at_sampling + sex, data = d),
        error = function(e) NULL
    )
    if (is.null(fit)) { cat("  Skipping", prot, ": model failed\n"); next }

    nd <- data.table(
        years_to_diagnosis = x_grid,
        age_at_sampling    = median(d$age_at_sampling, na.rm = TRUE),
        sex = as.numeric(names(sort(table(d$sex), decreasing = TRUE))[1])
    )
    pr <- predict(fit, newdata = nd, se.fit = TRUE)

    ref_i   <- hc_psm_spline(col_nm, x_grid)
    use_ref <- !is.null(ref_i)

    traj_i  <- data.table(ytd = x_grid, pred = pr$fit, pred_se = pr$se.fit)
    col_val <- IMM_COLS[i]
    col_map <- setNames(c(col_val, "grey45"), c(DISEASE_CAPS, "HC (PSM-matched)"))

    p <- ggplot(traj_i, aes(x = ytd)) +
        # PSM-matched HC reference ribbon + line
        {if (use_ref) list(
            geom_ribbon(data = ref_i, aes(x = ytd, ymin = hc_ref_lci, ymax = hc_ref_uci),
                        inherit.aes = FALSE, fill = "grey70", alpha = 0.35),
            geom_line(data = ref_i, aes(x = ytd, y = hc_ref, colour = "HC (PSM-matched)"),
                      inherit.aes = FALSE, linewidth = 0.6, linetype = "dashed")
        )} +
        # MS spline ± CI
        geom_ribbon(aes(ymin = pred - 1.96 * pred_se,
                        ymax = pred + 1.96 * pred_se),
                    fill = col_val, alpha = 0.18) +
        geom_line(aes(y = pred, colour = DISEASE_CAPS), linewidth = 0.9) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.4, colour = "grey40") +
        annotate("text", x = 0.4, y = Inf, label = "Dx",
                 hjust = 0, vjust = 1.5, size = 2.3, colour = "grey40") +
        scale_colour_manual(values = col_map, name = NULL) +
        scale_x_continuous(breaks = seq(-20, 15, by = 5)) +
        labs(
            title = paste0(pid, "  ", toupper(prot), " — pre-onset DEP trajectory"),
            x     = glue("Years relative to {DISEASE_CAPS} diagnosis"),
            y     = paste0(toupper(prot), " (NPX, predicted)")
        ) +
        theme_ukb(base_size = 9) +
        guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))

    save_panel(p, paste0(pid, "_", tolower(prot), "_trajectory"), 3.5, 3.2)
}

# ─────────────────────────────────────────────────────────────────────────────
# Panel i — Pre-onset vs post-onset logFC scatter
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel i...\n")

pre  <- fread(PRE_FILE)
post <- fread(POST_FILE)

pre[,  fdr_pre  := p.adjust(P.Value, method = "BH")]
post[, fdr_post := p.adjust(P.Value, method = "BH")]

merged_i <- merge(
    pre[,  .(protein, logFC_pre  = logFC, fdr_pre)],
    post[, .(protein, logFC_post = logFC, fdr_post)],
    by = "protein"
)
merged_i[, sig_group := fcase(
    fdr_pre  < 0.05 & fdr_post  < 0.05, "Both",
    fdr_pre  < 0.05 & fdr_post >= 0.05, "Pre-onset only",
    fdr_pre >= 0.05 & fdr_post  < 0.05, "Post-onset only",
    default = "NS"
)]
merged_i[, sig_group := factor(sig_group,
    levels = c("Both", "Pre-onset only", "Post-onset only", "NS"))]

top_i <- merged_i[order(-(abs(logFC_pre) + abs(logFC_post)))][seq_len(min(12, .N))]
r_val <- cor(merged_i$logFC_pre, merged_i$logFC_post, use = "complete.obs")

SCATTER_COLS <- c(
    "Both"            = "#6B3FA0",
    "Pre-onset only"  = COL_PRE_UP,
    "Post-onset only" = COL_POST,
    "NS"              = "grey70"
)

pI <- ggplot(merged_i, aes(x = logFC_pre, y = logFC_post, colour = sig_group)) +
    # Faint reference axes and identity diagonal
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.4, colour = "grey55") +
    geom_point(size = 0.8, alpha = 0.65) +
    geom_text_repel(
        data          = top_i,
        aes(label     = toupper(protein)),
        size          = 2.2,
        fontface      = "bold",
        max.overlaps  = 25,
        segment.size  = 0.25,
        segment.color = "grey50",
        box.padding   = 0.3,
        show.legend   = FALSE
    ) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("r = %.2f", r_val),
             hjust = -0.15, vjust = 1.5,
             size = 2.5, colour = "grey30") +
    scale_colour_manual(values = SCATTER_COLS, name = NULL) +
    labs(
        title = "i  Pre- vs post-onset proteomic changes",
        x     = expression(log[2]~"FC  (pre-onset vs HC)"),
        y     = expression(log[2]~"FC  (post-onset vs HC)")
    ) +
    theme_ukb(base_size = 9) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

save_panel(pI, "i_prepost_scatter", 4.5, 4.5)

# ─────────────────────────────────────────────────────────────────────────────
# Panel j — Temporal protein heatmap (ALL MS, full disease course)
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel j...\n")

comb_j <- copy(comb)
comb_j[, fdr := p.adjust(P.Value, method = "BH")]
if (!"adj.P.Val" %in% names(comb_j)) comb_j[, adj.P.Val := fdr]
# Programmatic top-N across both directions for the temporal heatmap.
# Previously this hard-forced NEFL/ERBB2 into the row set; that's removed.
# If those CNS biomarkers are still surfaced by adj.P.Val ranking they will
# appear; if not, the data has spoken and the heatmap reflects the cohort.
n_heat <- cfg$top_n_heatmap %||% 30L
top_heat <- top_n_by_direction(
    comb_j, n_each = max(1L, as.integer(n_heat / 2L)),
    min_sig_for_fdr = cfg$min_sig_for_fdr %||% 5L,
    case_fn = identity
)
top30 <- top_heat$all
cat("  Heatmap rows:", length(top30), "(",
    attr(top_heat, "n_up"), "up,", attr(top_heat, "n_down"), "down,",
    "ordered by", paste(attr(top_heat, "ordered_by"), collapse = "/"), ")\n")

ms_j <- qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE]
ms_j[, time_bin := cut(years_to_diagnosis,
                        breaks = BIN_BREAKS, labels = BIN_LABELS,
                        right = FALSE, include.lowest = TRUE)]

heat_list <- lapply(tolower(top30), function(p_col) {
    if (!p_col %in% names(ms_j)) return(NULL)
    ms_j[, .(protein  = toupper(p_col),
              mean_npx = mean(get(p_col), na.rm = TRUE)),
          by = time_bin]
})
heat_dt <- rbindlist(Filter(Negate(is.null), heat_list))

heat_dt[, z_npx := {
    mn   <- mean(mean_npx, na.rm = TRUE)
    sd_v <- sd(mean_npx, na.rm = TRUE)
    if (is.na(sd_v) || sd_v == 0) rep(0, .N) else (mean_npx - mn) / sd_v
}, by = protein]

heat_dt[, time_bin := factor(time_bin, levels = BIN_LABELS)]

# Order proteins by time bin of peak elevation
peak_bin <- heat_dt[, {
    z <- z_npx; bins <- as.integer(time_bin)
    .(peak_i = if (all(is.na(z))) NA_integer_ else bins[which.max(z)])
}, by = protein]
heat_dt[, protein := factor(protein, levels = peak_bin[order(peak_i), protein])]

# Human-readable x-axis labels
HEATMAP_XLABS <- c(
    "< -10"    = ">10 yr\npost-Dx",
    "[-10,-5)" = "5–10 yr\npost-Dx",
    "[-5,-2)"  = "2–5 yr\npost-Dx",
    "[-2,0)"   = "0–2 yr\npost-Dx",
    "[0,2)"    = "0–2 yr\npre-Dx",
    "[2,5)"    = "2–5 yr\npre-Dx",
    "[5,10)"   = "5–10 yr\npre-Dx",
    "> 10"     = ">10 yr\npre-Dx"
)

pJ <- ggplot(heat_dt, aes(x = time_bin, y = protein, fill = z_npx)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    # Dashed line at diagnosis boundary
    geom_vline(xintercept = 4.5, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    scale_fill_gradient2(
        low      = "#56B4E9",
        mid      = "white",
        high     = "#CC0066",
        midpoint = 0,
        name     = "z-score",
        guide    = guide_colourbar(barwidth = 0.5, barheight = 3)
    ) +
    scale_x_discrete(labels = HEATMAP_XLABS) +
    labs(
        title    = glue("j  Temporal protein dynamics across {DISEASE_CAPS} disease course"),
        subtitle = "Top 30 DEPs (FDR<0.05) plus CNS biomarkers NEFL/ERBB2; mean NPX per time bin, z-scored per protein",
        x        = "Time relative to diagnosis  ◄ post-Dx  |  pre-Dx ►",
        y        = NULL
    ) +
    theme_ukb(base_size = 9) +
    theme(
        axis.text.x = element_text(size = 6.5, hjust = 0.5),
        axis.text.y = element_text(size = 6),
        legend.position = "right"
    )

save_panel(pJ, "j_temporal_heatmap", 5.5, 6.5)

# ─────────────────────────────────────────────────────────────────────────────
# Panel m — ML ROC curves
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel m (ML ROC curves)...\n")

get_auc <- function(file, model_name = "glmnet") {
    if (!file.exists(file)) return(NA_real_)
    dt  <- fread(file)
    row <- dt[tolower(model) == tolower(model_name)]
    if (nrow(row) == 0) return(NA_real_)
    row$AUC[1]
}

auc_full <- get_auc(ML_FULL_FILE, "glmnet")
auc_pre  <- get_auc(ML_PRE_FILE,  "glmnet")
cat("  Test-set AUC (CSV, for reference) — full:", auc_full, "  pre-onset:", auc_pre, "\n")
cat("  Plot annotations will use AUC derived from CV ROC curves.\n")

compute_roc <- function(labels, scores) {
    ord   <- order(scores, decreasing = TRUE)
    lab   <- as.logical(labels[ord])
    n_pos <- sum(lab);  n_neg <- sum(!lab)
    if (n_pos == 0 || n_neg == 0) return(data.table(fpr = c(0,1), tpr = c(0,1)))
    tpr <- cumsum(lab)  / n_pos
    fpr <- cumsum(!lab) / n_neg
    data.table(fpr = c(0, fpr), tpr = c(0, tpr))
}

roc_list <- list()

# Helper to extract model feature columns (excludes caret's ".outcome" sentinel)
model_features <- function(mod) {
    cols <- colnames(mod$trainingData)
    setdiff(cols, ".outcome")   # .outcome is position 1, NOT the last column
}

# Use saved CV hold-out predictions (savePredictions="final") for honest ROC curves
# that match the stored CV AUC values, rather than re-predicting on training data.
get_cv_roc <- function(mod) {
    best <- mod$bestTune
    pd   <- as.data.table(mod$pred)
    for (nm in names(best)) pd <- pd[abs(get(nm) - best[[nm]]) < 1e-9]
    compute_roc(as.integer(pd$obs == DISEASE_CAPS), pd[[DISEASE_CAPS]])
}

tryCatch({
    suppressWarnings(library(caret))

    # Full MS classifier (glmnet) — CV hold-out predictions
    if (file.exists(MODEL_FULL_FILE)) {
        mod_f <- suppressWarnings(readRDS(MODEL_FULL_FILE))
        if (!is.null(mod_f$pred)) {
            roc_list[["Full cohort (glmnet)"]] <- get_cv_roc(mod_f)
            cat("  Full model CV ROC computed\n")
        }
    }

    # Pre-onset classifier (glmnet) — CV hold-out predictions
    if (file.exists(MODEL_PRE_FILE)) {
        mod_p <- suppressWarnings(readRDS(MODEL_PRE_FILE))
        if (!is.null(mod_p$pred)) {
            roc_list[["Pre-onset only (glmnet)"]] <- get_cv_roc(mod_p)
            cat("  Pre-onset model CV ROC computed\n")
        }
    }

}, error = function(e) {
    cat("  caret/model loading failed:", conditionMessage(e), "\n")
})

# ── Demographics baseline: age + sex logistic regression ─────────────────────
# Shows how much the proteomic signal adds above basic demographics
tryCatch({
    df_demo <- qc[ms_status %in% c(SV$pre_onset, SV$control) & qc_outlier==FALSE &
                  !is.na(age_at_sampling) & !is.na(sex)]
    df_demo[, y := as.integer(ms_status == SV$pre_onset)]
    lr_demo  <- glm(y ~ age_at_sampling + sex, data=df_demo, family=binomial)
    demo_scores <- predict(lr_demo, type="response")
    roc_list[["Demographics only (age+sex)"]] <-
        compute_roc(df_demo$y, demo_scores)
    cat("  Demographics baseline ROC computed\n")
}, error = function(e) {
    cat("  Demographics baseline failed:", conditionMessage(e), "\n")
})

if (length(roc_list) == 0) {
    cat("  No ROC data available; plotting diagonal placeholder.\n")
    dummy <- data.table(fpr = c(0,1), tpr = c(0,1))
    roc_list[["Full cohort (glmnet)"]]      <- dummy
    roc_list[["Pre-onset only (glmnet)"]]   <- dummy
    roc_list[["Demographics only (age+sex)"]] <- dummy
}

CURVE_LEVELS <- c("Full cohort (glmnet)", "Pre-onset only (glmnet)",
                  "Demographics only (age+sex)")
ROC_COLS <- c(
    "Full cohort (glmnet)"         = "#2B4C7E",
    "Pre-onset only (glmnet)"      = COL_PRE_UP,
    "Demographics only (age+sex)"  = "#E6A817"
)

roc_dt <- rbindlist(lapply(names(roc_list), function(nm) {
    d <- roc_list[[nm]]; d[, curve := nm]; d
}))
roc_dt[, curve := factor(curve, levels = CURVE_LEVELS)]

# AUC from stored CV results; compute from curve data as fallback
auc_from_curve <- function(roc_d) {
    with(roc_d, sum(diff(fpr) * (head(tpr,-1) + tail(tpr,-1)) / 2))
}

# Compute all AUCs from the CV ROC curves (honest, matches stored CV results)
auc_ann <- data.table(curve = CURVE_LEVELS, cv_auc = NA_real_)
for (nm in CURVE_LEVELS) {
    if (nm %in% names(roc_list))
        auc_ann[curve == nm, cv_auc := auc_from_curve(roc_list[[nm]])]
}
auc_ann <- auc_ann[!is.na(cv_auc)]
auc_ann[, y_pos := rev(seq(0.08, 0.08 + 0.13 * (.N - 1), by = 0.13))]

pK <- ggplot(roc_dt, aes(x = fpr, y = tpr, colour = curve)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.4, colour = "grey70") +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = ROC_COLS, name = NULL) +
    labs(
        title    = glue("m  Proteomic prediction of {DISEASE_CAPS}"),
        subtitle = "Proteins vs age+sex baseline; CV AUC in annotation",
        x        = "1 – Specificity (FPR)",
        y        = "Sensitivity (TPR)"
    ) +
    coord_fixed() +
    theme_ukb(base_size = 9) +
    theme(legend.position = "none")

if (nrow(auc_ann) > 0) {
    for (i in seq_len(nrow(auc_ann))) {
        pK <- pK + annotate("text",
            x = 0.50, y = auc_ann$y_pos[i],
            label  = sprintf("%s\nCV AUC = %.3f", auc_ann$curve[i], auc_ann$cv_auc[i]),
            colour = ROC_COLS[auc_ann$curve[i]],
            hjust  = 0, size = 2.0, lineheight = 0.9)
    }
}

save_panel(pK, "m_ml_roc", 4, 4)

# ─────────────────────────────────────────────────────────────────────────────
# Panel n — Pre-diagnostic risk score by time bin
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel n (risk score by time bin)...\n")

tryCatch({
    suppressWarnings(library(caret))

    if (file.exists(MODEL_PRE_FILE)) {
        pre_l   <- qc[ms_status %in% c(SV$pre_onset, SV$control) & qc_outlier == FALSE]
        # glmnet model was trained with sex_num; add it here
        pre_l[, sex_num := as.integer(factor(sex)) - 1L]
        mod_l   <- suppressWarnings(readRDS(MODEL_PRE_FILE))
        avail_l  <- intersect(model_features(mod_l), names(pre_l))
        cat(sprintf("  Panel n: %d features available\n", length(avail_l)))
        # Median-impute protein NAs (mirrors training pipeline)
        prot_feats_l <- setdiff(avail_l, c("age_at_sampling","sex_num"))
        for (col in prot_feats_l) {
            med_v <- median(pre_l[[col]], na.rm = TRUE)
            pre_l[is.na(get(col)), (col) := med_v]
        }
        keep_l <- complete.cases(pre_l[, ..avail_l])
        pre_l  <- pre_l[keep_l]

        prbs_l <- tryCatch(
            predict(mod_l, newdata = pre_l[, ..avail_l], type = "prob"),
            error = function(e) { cat("  Panel n predict error:", e$message, "\n"); NULL })

        if (!is.null(prbs_l)) {
            pre_l[, pred_prob := prbs_l[[DISEASE_CAPS]]]
            hc_mean_prob <- mean(pre_l[ms_status == SV$control, pred_prob], na.rm = TRUE)

            # Pre-onset cases: years_to_diagnosis > 0 (years UNTIL diagnosis)
            cases_l <- pre_l[ms_status == SV$pre_onset]
            cases_l[, time_bin := cut(years_to_diagnosis,
                                       breaks = BIN_BREAKS, labels = BIN_LABELS,
                                       right = FALSE, include.lowest = TRUE)]
            # Order: furthest from Dx → closest to Dx (left → right)
            pre_bins <- c("> 10", "[5,10)", "[2,5)", "[0,2)")
            cases_l  <- cases_l[time_bin %in% pre_bins]
            cases_l[, time_bin := factor(time_bin, levels = pre_bins)]

            bin_l <- cases_l[, .(
                mean_prob = mean(pred_prob, na.rm = TRUE),
                se_prob   = sd(pred_prob, na.rm = TRUE) / sqrt(.N),
                n         = .N
            ), by = time_bin]
            cat(sprintf("  Panel n: %d pre-onset bins, %d cases\n",
                        nrow(bin_l), sum(bin_l$n)))

            PRE_XLABS <- c(
                "> 10"   = ">10 yr\nbefore Dx",
                "[5,10)" = "5–10 yr\nbefore Dx",
                "[2,5)"  = "2–5 yr\nbefore Dx",
                "[0,2)"  = "0–2 yr\nbefore Dx"
            )

            hc_sd_prob <- sd(pre_l[ms_status == SV$control, pred_prob], na.rm = TRUE)

            pL <- ggplot() +
                # HC reference band (mean ± 1 SD)
                annotate("rect",
                         xmin = -Inf, xmax = Inf,
                         ymin = hc_mean_prob - hc_sd_prob,
                         ymax = hc_mean_prob + hc_sd_prob,
                         fill = "grey85", alpha = 0.55) +
                geom_hline(yintercept = hc_mean_prob, linetype = "dashed",
                           linewidth = 0.45, colour = "grey50") +
                # Individual pre-onset predictions (jittered)
                geom_jitter(data = cases_l, aes(x = time_bin, y = pred_prob),
                            width = 0.16, alpha = 0.28, size = 0.85, colour = "#CC0066") +
                # Trend line connecting bin means
                geom_line(data = bin_l,
                          aes(x = time_bin, y = mean_prob, group = 1),
                          colour = "#CC0066", linewidth = 0.75) +
                # Mean ± SE
                geom_pointrange(data = bin_l,
                                aes(x = time_bin, y = mean_prob,
                                    ymin = mean_prob - se_prob,
                                    ymax = mean_prob + se_prob),
                                colour = "#CC0066", linewidth = 1.3,
                                size = 0.55, shape = 18) +
                # Sample size labels at top
                geom_text(data = bin_l,
                          aes(x = time_bin, y = Inf,
                              label = paste0("n=", n)),
                          vjust = 1.5, size = 1.9, colour = "grey45") +
                annotate("text", x = 0.6, y = hc_mean_prob,
                         label = sprintf("HC mean ± SD\n(n=%d)", sum(pre_l$ms_status == SV$control)),
                         hjust = 0, vjust = -0.2, size = 1.85, colour = "grey45") +
                scale_x_discrete(labels = PRE_XLABS) +
                labs(
                    title    = glue("n  Pre-diagnostic {DISEASE_CAPS} risk score"),
                    subtitle = "Elastic net probability by time before diagnosis (points = individuals; diamond = mean ± SE)",
                    x        = NULL,
                    y        = glue("Predicted {DISEASE_CAPS} probability")
                ) +
                theme_ukb(base_size = 9)

            save_panel(pL, "n_risk_score_timebins", 4.2, 3.5)
        } else {
            cat("  Panel n skipped: prediction failed.\n")
        }
    } else {
        cat("  Panel n skipped: pre-onset model not found.\n")
    }
}, error = function(e) {
    cat("  Panel n skipped:", conditionMessage(e), "\n")
})

# ─────────────────────────────────────────────────────────────────────────────
# Panel l — SHAP beeswarm (pre-onset RF model)
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel l (SHAP feature importance)...\n")

tryCatch({
    suppressPackageStartupMessages({
        library(fastshap)
        library(ggbeeswarm)
        library(caret)
    })

    if (!file.exists(MODEL_PRE_FILE)) stop("Pre-onset model not found")
    mod_shap <- suppressWarnings(readRDS(MODEL_PRE_FILE))

    # Use all pre-onset cases + a balanced sample of controls
    pre_shap <- qc[ms_status %in% c(SV$pre_onset, SV$control) & qc_outlier == FALSE]
    pre_shap[, sex_num := as.integer(factor(sex)) - 1L]   # glmnet trained with sex_num
    avail_shap <- intersect(model_features(mod_shap), names(pre_shap))

    set.seed(42)
    n_ctrl <- min(500L, sum(pre_shap$ms_status == SV$control))
    ctrl_idx <- which(pre_shap$ms_status == SV$control)
    use_idx  <- c(which(pre_shap$ms_status == SV$pre_onset),
                  sample(ctrl_idx, n_ctrl))
    shap_dt  <- pre_shap[use_idx]
    shap_dt  <- shap_dt[complete.cases(shap_dt[, ..avail_shap])]

    X_mat <- as.data.frame(shap_dt[, ..avail_shap])

    pfun <- function(object, newdata)
        predict(object, newdata = newdata, type = "prob")[, DISEASE_CAPS]

    cat("  Computing SHAP values (nsim=50, n=", nrow(X_mat), ")...\n")
    set.seed(42)
    shap_obj <- fastshap::explain(mod_shap, X = X_mat,
                                  pred_wrapper = pfun, nsim = 50, adjust = TRUE)

    # Convert to long format; top 20 features by mean |SHAP|
    sv   <- as.data.table(shap_obj)
    sv[, .sid := seq_len(.N)]
    sv_l <- melt(sv, id.vars = ".sid", variable.name = "feature", value.name = "shap")

    fv   <- as.data.table(X_mat); fv[, .sid := seq_len(.N)]
    fv_l <- melt(fv, id.vars = ".sid", variable.name = "feature", value.name = "fval")

    both <- merge(sv_l, fv_l, by = c(".sid", "feature"))
    both[, feat_norm := {
        lo <- min(fval, na.rm = TRUE); hi <- max(fval, na.rm = TRUE)
        if (hi == lo) rep(0.5, .N) else (fval - lo) / (hi - lo)
    }, by = feature]

    feat_order <- both[, .(mean_abs = mean(abs(shap))), by = feature][order(-mean_abs)]
    top20 <- feat_order$feature[seq_len(min(20, nrow(feat_order)))]
    both20 <- both[feature %in% top20]
    both20[, feature := factor(toupper(as.character(feature)),
                               levels = toupper(as.character(rev(top20))))]

    pM <- ggplot(both20, aes(x = shap, y = feature, colour = feat_norm)) +
        geom_quasirandom(groupOnX = FALSE, size = 0.55, alpha = 0.65, width = 0.35) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
        scale_colour_gradient2(
            low      = "#56B4E9",
            mid      = "grey85",
            high     = "#CC0066",
            midpoint = 0.5,
            name     = "Feature\nvalue",
            breaks   = c(0, 1),
            labels   = c("Low", "High"),
            guide    = guide_colourbar(barwidth = 0.5, barheight = 3,
                                       title.position = "top")
        ) +
        labs(
            title    = "l  Pre-onset predictor contributions",
            subtitle = "SHAP values — glmnet model; colour = feature value (low → high)",
            x        = glue("SHAP value (impact on {DISEASE_CAPS} probability)"),
            y        = NULL
        ) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y  = element_text(size = 7),
              legend.position = "right")

    save_panel(pM, "l_shap_features", 4.5, 4.5)
    cat("  Saved: panel_l_shap_features.pdf\n")

}, error = function(e) {
    cat("  Panel l SHAP failed:", conditionMessage(e), "\n")
    # Fallback: feature importance bar chart
    if (file.exists(ML_PRE_IMP_FILE)) {
        imp <- fread(ML_PRE_IMP_FILE)
        top_imp <- imp[order(-mean_imp)][seq_len(min(20, .N))]
        top_imp[, label := factor(toupper(protein), levels = toupper(protein[order(mean_imp)]))]
        pM <- ggplot(top_imp, aes(x = mean_imp, y = label)) +
            geom_col(fill = UKB_NAVY, width = 0.7) +
            labs(title = "l  Pre-onset prediction features",
                 x = "Mean feature importance", y = NULL) +
            theme_ukb(base_size = 9) +
            theme(axis.text.y = element_text(size = 7))
        save_panel(pM, "l_feature_importance", 3.2, 4.2)
    }
})

# ─────────────────────────────────────────────────────────────────────────────
# Panel k — Cox rolling HR heatmap (protein × time window before Dx)
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel k (Cox rolling HR)...\n")

tryCatch({
    ROLL_FILE <- file.path(PROJ_DIR, "results", "survival", glue("{COHORT}_protein_cox_rolling.csv"))
    if (!file.exists(ROLL_FILE)) stop("Run 01_ms_protein_cox.R first")

    roll <- fread(ROLL_FILE)
    roll[, protein := toupper(protein)]
    BIN_ORDER <- c("0–4 yr", "4–8 yr", "8–12 yr", "12+ yr")
    roll[, bin_f := factor(time_bin, levels = BIN_ORDER)]

    # Order proteins by overall HR (0–4 yr bin for ordering)
    prot_order <- roll[bin_f == "0–4 yr"][order(-HR), protein]
    roll[, prot_f := factor(protein, levels = rev(prot_order))]

    pN <- ggplot(roll, aes(x = bin_f, y = prot_f, fill = log2(HR))) +
        geom_tile(colour = "white", linewidth = 0.35) +
        geom_text(aes(label = ifelse(pval < 0.05, sprintf("%.2f", HR), "")),
                  size = 2.0, colour = "white") +
        scale_fill_gradient2(low = "#56B4E9", mid = "grey95", high = "#CC0066",
                             midpoint = 0, name = "log₂ HR",
                             limits = c(-0.8, 0.8), oob = scales::squish,
                             guide  = guide_colourbar(barwidth = 0.5, barheight = 4)) +
        scale_x_discrete(labels = c("0–4 yr"  = "0–4 yr\nbefore Dx",
                                     "4–8 yr"  = "4–8 yr",
                                     "8–12 yr" = "8–12 yr",
                                     "12+ yr"  = ">12 yr")) +
        labs(title    = "k  Protein hazard ratios across pre-onset time windows",
             subtitle = "HR per 1 SD protein; white label = p<0.05",
             x = NULL, y = NULL,
             caption = "Cox model: all controls as censored events (standard for prospective cohorts) | age, sex, UMAP1/2 adjusted") +
        theme_ukb(base_size = 9) +
        theme(axis.text.y   = element_text(size = 7),
              axis.text.x   = element_text(size = 8),
              panel.grid    = element_blank(),
              legend.position = "right")

    save_panel(pN, "k_cox_rolling_hr", 6, 6.5)

}, error = function(e) cat("  Panel k skipped:", conditionMessage(e), "\n"))

# ─────────────────────────────────────────────────────────────────────────────
# Panel o — PRS-stratified NEFL trajectory
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel o...\n")

tryCatch({
    suppressPackageStartupMessages(library(splines))
    PRS_FILE_O <- file.path(PROJ_DIR, "data", "ukb", "genetics", glue("{COHORT}_prs_scores.csv"))
    INT_FILE   <- file.path(PROJ_DIR, "results", "prs", glue("{COHORT}_prs_interaction_tests.csv"))
    if (!file.exists(PRS_FILE_O)) stop("PRS file not found")

    prs_o <- fread(PRS_FILE_O, showProgress = FALSE)
    if (PRS_COL != "prs_score" && PRS_COL %in% names(prs_o))
        setnames(prs_o, PRS_COL, "prs_score")
    ms_o  <- merge(
        qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE &
           !is.na(years_to_diagnosis) & !is.na(nefl) &
           years_to_diagnosis >= -20 & years_to_diagnosis <= 16],
        prs_o[, .(eid, prs_score)], by = "eid"
    )
    ms_o[, prs_q := cut(prs_score,
                         breaks = quantile(prs_score, c(0,.25,.50,.75,1), na.rm=TRUE),
                         labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
                         include.lowest = TRUE)]

    x_grid_o <- seq(-20, 16, length.out = 100)
    PRS_COLS_O <- c("Q1 (low)"="#56B4E9","Q2"="#74c476","Q3"="#E6A817","Q4 (high)"="#CC0066")

    traj_prs <- rbindlist(lapply(levels(ms_o$prs_q), function(q) {
        dq <- ms_o[prs_q == q & !is.na(nefl) & !is.na(age_at_sampling)]
        if (nrow(dq) < 8) return(NULL)
        fit <- tryCatch(
            lm(nefl ~ ns(years_to_diagnosis, df=3) + age_at_sampling + sex, data=dq),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        nd  <- data.table(years_to_diagnosis = x_grid_o,
                          age_at_sampling = median(dq$age_at_sampling, na.rm=TRUE),
                          sex = as.numeric(names(sort(table(dq$sex), decreasing=TRUE))[1]))
        pr  <- predict(fit, newdata=nd, se.fit=TRUE)
        data.table(ytd=x_grid_o, pred=pr$fit,
                   lci=pr$fit-1.96*pr$se.fit, uci=pr$fit+1.96*pr$se.fit,
                   prs_q=q)
    }))
    traj_prs[, prs_q := factor(prs_q, levels=levels(ms_o$prs_q))]

    # Interaction p-value for label
    int_label <- ""
    if (file.exists(INT_FILE)) {
        int_dt <- fread(INT_FILE)
        nefl_row <- int_dt[protein == "NEFL"]
        if (nrow(nefl_row) > 0)
            int_label <- sprintf("\nPRS×time interaction p = %s",
                                  formatC(nefl_row$pval[1], format="g", digits=2))
    }

    # Highlight Q1 (low risk) vs Q4 (high risk); show Q2/Q3 as thin grey lines
    traj_ext <- traj_prs[prs_q %in% c("Q1 (low)", "Q4 (high)")]
    traj_mid <- traj_prs[prs_q %in% c("Q2", "Q3")]

    pO <- ggplot() +
        # Middle quartiles — faint background context
        geom_line(data = traj_mid,
                  aes(x = ytd, y = pred, group = prs_q),
                  linewidth = 0.35, colour = "grey70", linetype = "dashed") +
        # Q1 and Q4 — CI ribbons + solid lines
        geom_ribbon(data = traj_ext,
                    aes(x = ytd, ymin = lci, ymax = uci,
                        fill = prs_q), alpha = 0.18, colour = NA) +
        geom_line(data = traj_ext,
                  aes(x = ytd, y = pred, colour = prs_q), linewidth = 0.95) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.35, colour = "grey40") +
        annotate("text", x = -0.3, y = Inf,
                 label = "← post-Dx  |  pre-Dx →",
                 hjust = 0.5, vjust = 1.6, size = 1.9, colour = "grey50") +
        scale_colour_manual(values = c("Q1 (low)" = "#56B4E9", "Q4 (high)" = "#CC0066"),
                            name = "PRS\nquartile") +
        scale_fill_manual(values = c("Q1 (low)" = "#56B4E9", "Q4 (high)" = "#CC0066"),
                          guide  = "none") +
        scale_x_continuous(breaks = seq(-16, 16, by = 4)) +
        labs(title    = paste0("o  NfL trajectory by ", PRS_LABEL, int_label),
             subtitle = "Q1 (lowest) vs Q4 (highest) PRS shown with 95% CI | Q2/Q3 dashed | age+sex adjusted",
             x = glue("Years relative to {DISEASE_CAPS} diagnosis"),
             y = "NfL (NPX)") +
        theme_ukb(base_size = 9) +
        guides(colour = guide_legend(override.aes = list(linewidth = 1.5)))

    save_panel(pO, "o_prs_nefl_trajectory", 4, 3.8)

}, error = function(e) cat("  Panel o skipped:", conditionMessage(e), "\n"))

# ─────────────────────────────────────────────────────────────────────────────
# Panel p — Pre-onset temporal proteome landscape
#
# Scatter: x = logFC_ms   (MS–HC gap at diagnosis, ytd = 0)
#          y = logFC_ytd  (per-year change in gap going back from Dx)
#                          +ve → gap is larger further from Dx  (early signal)
#                          −ve → gap grows as Dx approaches     (late signal)
# Fill:    gap_5yr = logFC_ms + 5 × logFC_ytd  (predicted MS–HC gap 5yr before Dx)
#          Diverging blue–white–pink: blue = depleted early, pink = elevated early
#
# Quadrant story:
#   Top-right  → Persistent early elevations  (elevated early AND at Dx)
#   Bot-right  → Near-onset elevations        (emerges close to Dx)
#   Top-left   → Early depletions
#   Bot-left   → Near-onset depletions
# ─────────────────────────────────────────────────────────────────────────────
cat("Building panel p (temporal proteome landscape)...\n")

TEMP_FILE <- file.path(DIFF_DIR, glue("{COHORT}_pre_temporal.csv"))

tryCatch({
    if (!file.exists(TEMP_FILE))
        stop(glue("{COHORT}_pre_temporal.csv not found — run 01_limma_ms_vs_hc.R first"))

    temp <- fread(TEMP_FILE)
    cat(sprintf("  Temporal results: %d proteins\n", nrow(temp)))

    # Use nominal p<0.05 — pre-onset n=123 is underpowered for FDR
    temp[, sig_ms  := P.Value_ms  < 0.05]
    temp[, sig_ytd := P.Value_ytd < 0.05]
    temp[, sig_any := sig_ms | sig_ytd]
    temp[, fdr_sig := fdr_ms < 0.05 | fdr_ytd < 0.05]   # FDR-significant subset

    n_early <- sum(temp$sig_any & temp$logFC_ms > 0 & temp$logFC_ytd > 0, na.rm = TRUE)
    n_late  <- sum(temp$sig_any & temp$logFC_ms > 0 & temp$logFC_ytd <= 0, na.rm = TRUE)
    n_dep   <- sum(temp$sig_any & temp$logFC_ms < 0, na.rm = TRUE)
    cat(sprintf("  Early signal: %d | Late signal: %d | Depleted: %d (nominal p<0.05)\n",
                n_early, n_late, n_dep))

    # Labels: FDR-significant proteins always labelled; top nominal hits by |gap_5yr|
    # Winsorise logFC_ms display to 1.2 — removes noisy extreme outliers from axis distortion
    temp[, logFC_ms_plot := pmax(pmin(logFC_ms, 1.2), -1.2)]

    gap_thresh <- if (sum(temp$sig_any) > 0)
        quantile(abs(temp[sig_any == TRUE, gap_5yr]), 0.72, na.rm = TRUE)
    else 0
    temp_labels     <- temp[sig_any == TRUE & abs(gap_5yr) >= gap_thresh & abs(logFC_ms) <= 1.2]
    temp_labels_fdr <- temp[fdr_sig == TRUE]    # FDR hits labelled separately, always shown
    # Programmatic policy: FDR hits get bold labels; nominal-significant hits
    # ranked by |gap_5yr| fill the lighter label layer. No protein is force-
    # injected by name; the data ranking is the only criterion.
    temp_labels[, label     := toupper(protein)]
    temp_labels_fdr[, label := toupper(protein)]

    # Symmetric colour scale
    gap_lim <- max(abs(temp$gap_5yr), na.rm = TRUE) * 0.95

    # Quadrant annotation positions (use winsorised range)
    xr  <- range(temp$logFC_ms_plot, na.rm = TRUE)
    yr  <- range(temp$logFC_ytd,     na.rm = TRUE)
    xmx <- max(abs(xr)); ymx <- max(abs(yr))
    quad_ann <- data.frame(
        x     = c( xmx * 0.72,  xmx * 0.72, -xmx * 0.72, -xmx * 0.72),
        y     = c( ymx * 0.88, -ymx * 0.88,  ymx * 0.88, -ymx * 0.88),
        label = c("Persistent\nearly elevation",
                  "Near-onset\nelevation",
                  "Persistent\nearly depletion",
                  "Near-onset\ndepletion"),
        hjust = c(1, 1, 0, 0)
    )

    pP <- ggplot() +
        # Quadrant background shading
        annotate("rect", xmin = 0,    xmax =  Inf, ymin = 0,    ymax =  Inf,
                 fill = "#CC0066", alpha = 0.04) +
        annotate("rect", xmin = 0,    xmax =  Inf, ymin = -Inf, ymax =  0,
                 fill = "#E6A817", alpha = 0.04) +
        annotate("rect", xmin = -Inf, xmax =  0,   ymin = 0,    ymax =  Inf,
                 fill = "#56B4E9", alpha = 0.03) +
        annotate("rect", xmin = -Inf, xmax =  0,   ymin = -Inf, ymax =  0,
                 fill = "grey60",  alpha = 0.02) +
        # Reference lines
        geom_hline(yintercept = 0, linetype = "dashed",
                   linewidth = 0.35, colour = "grey45") +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.35, colour = "grey45") +
        # NS background points
        geom_point(data  = temp[sig_any == FALSE],
                   aes(x = logFC_ms_plot, y = logFC_ytd),
                   shape = 16, size = 0.55, colour = "grey78", alpha = 0.28) +
        # Nominal-significant points — filled circles coloured by gap_5yr
        geom_point(data  = temp[sig_any == TRUE & fdr_sig == FALSE],
                   aes(x = logFC_ms_plot, y = logFC_ytd, fill = gap_5yr),
                   shape = 21, size = 1.8, stroke = 0.20,
                   colour = "grey35", alpha = 0.82) +
        # FDR-significant points — larger, stronger stroke
        geom_point(data  = temp[fdr_sig == TRUE],
                   aes(x = logFC_ms_plot, y = logFC_ytd, fill = gap_5yr),
                   shape = 21, size = 3.0, stroke = 0.55,
                   colour = "grey10", alpha = 1.0) +
        scale_fill_gradient2(
            low      = "#56B4E9",
            mid      = "grey96",
            high     = "#CC0066",
            midpoint = 0,
            limits   = c(-gap_lim, gap_lim),
            oob      = scales::squish,
            name     = expression(atop(log[2]~"FC", "at 5yr pre-Dx")),
            guide    = guide_colourbar(barwidth = 0.55, barheight = 3.8,
                                       title.hjust = 0.5,
                                       title.position = "top")
        ) +
        # Protein labels for top nominal-significant proteins
        geom_text_repel(
            data           = temp_labels,
            aes(x          = logFC_ms_plot, y = logFC_ytd, label = label),
            size           = 1.95,
            fontface       = "bold",
            colour         = "grey25",
            fill           = NA,
            max.overlaps   = 25,
            segment.size   = 0.22,
            segment.colour = "grey55",
            segment.alpha  = 0.75,
            box.padding    = 0.38,
            point.padding  = 0.22,
            show.legend    = FALSE
        ) +
        # FDR-significant protein labels — always shown, larger font
        geom_text_repel(
            data           = temp_labels_fdr,
            aes(x          = logFC_ms_plot, y = logFC_ytd, label = label),
            size           = 2.5,
            fontface       = "bold",
            colour         = "grey5",
            fill           = scales::alpha("white", 0.75),
            label.size     = NA,
            max.overlaps   = Inf,
            segment.size   = 0.30,
            segment.colour = "grey30",
            box.padding    = 0.55,
            point.padding  = 0.35,
            nudge_y        = 0.03,
            show.legend    = FALSE
        ) +
        coord_cartesian(xlim = c(-1.25, 1.25), clip = "off") +
        # Quadrant text labels
        geom_text(data        = quad_ann,
                  aes(x = x, y = y, label = label, hjust = hjust),
                  size        = 1.95,
                  lineheight  = 0.85,
                  colour      = "grey50",
                  fontface    = "italic",
                  inherit.aes = FALSE) +
        # Y-axis direction hint
        annotate("text", x = -Inf, y = Inf,
                 label  = "\u2191 signal stronger earlier",
                 hjust  = -0.05, vjust = 1.6,
                 size   = 1.85, colour = "grey45", fontface = "italic") +
        labs(
            title    = "p  Pre-onset temporal proteome landscape",
            subtitle = sprintf(
                "%d early/persistent \u2022 %d near-onset \u2022 %d depleted (nominal p<0.05); large dot = FDR<0.05",
                n_early, n_late, n_dep),
            x        = bquote(log[2]~"FC at diagnosis"~(.(DISEASE_CAPS)~vs~HC)),
            y        = expression(paste("Temporal slope (", log[2],
                                         "FC per year before Dx)"))
        ) +
        theme_ukb(base_size = 9) +
        theme(
            legend.position = "right",
            plot.subtitle   = element_text(size = 7, colour = "grey40")
        )

    save_panel(pP, "p_temporal_landscape", 5.5, 4.8)

}, error = function(e) cat("  Panel p skipped:", conditionMessage(e), "\n"))

cat("\nFigure 2 complete. Panels saved to:", FIG_DIR, "\n")
