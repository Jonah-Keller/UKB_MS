#!/usr/bin/env Rscript
# 04_dti_protein_interaction.R
# Protein × disease-diagnosis interaction on DTI white matter metrics
#
# Design (powered via full cohort):
#   Population: participants with both Olink + DTI data
#   Model per protein × tract:
#     FA ~ scale(protein) * is_case + age_at_mri + sex
#   Test: interaction term (scale(protein):is_case)
#   — asks: does the protein-FA slope DIFFER between cases and HC?
#   — FDR correction across all protein × tract × metric combinations
#
# Metrics tested: FA, MD, ICVF (27 tracts each)
# Proteins: DEPs (FDR<0.05 in limma combined contrast)
#
# Output: results/mri/<cohort>_dti_interaction_*.{csv,pdf}

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
DTI_FILE    <- file.path(CADASIL_DIR, "brain_mri", "dti", "dti_tract_protein_data.tsv")
QC_FILE     <- here::here("data", "ukb", "olink", "processed",
                          glue::glue("{cfg$cohort_short}_olink_qc.csv"))
DIFF_FILE   <- here::here("results", "differential",
                          glue::glue("{cfg$cohort_short}_combined_vs_hc.csv"))
OUT_DIR     <- here::here("results", "mri")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR <- 0.05
MIN_N   <- 50   # minimum overlap per protein × tract

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")

dti <- fread(DTI_FILE, sep = "\t", showProgress = FALSE)
cat(sprintf("  DTI file: %d participants\n", nrow(dti)))

# Use age_at_mri and sex already in the DTI file
if (!"age_at_mri" %in% names(dti)) stop("age_at_mri not found in DTI file")
if (!"sex" %in% names(dti))        stop("sex not found in DTI file")

fa_cols   <- grep("^fa_",   names(dti), value = TRUE)
md_cols   <- grep("^md_",   names(dti), value = TRUE)
icvf_cols <- grep("^icvf_", names(dti), value = TRUE)
cat(sprintf("  DTI metrics: %d FA, %d MD, %d ICVF tracts\n",
            length(fa_cols), length(md_cols), length(icvf_cols)))

# Load disease QC (protein levels + status)
case_qc <- fread(QC_FILE, showProgress = FALSE)
case_qc <- case_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]

# DEP list
deps        <- fread(DIFF_FILE, showProgress = FALSE)
dep_prots   <- tolower(deps[adj.P.Val < FDR_THR, protein])
avail_prots <- intersect(dep_prots, names(case_qc))
cat(sprintf("  DEPs available: %d\n", length(avail_prots)))

# Merge: participants with both DTI and Olink
# Use only eid + proteins + status from QC file
prot_cols <- c("eid", STATUS_COL, avail_prots)
case_sub  <- case_qc[, ..prot_cols]

# DTI keeps its own age_at_mri and sex
dti_sub <- dti[, c("eid", "age_at_mri", "sex", fa_cols, md_cols, icvf_cols), with = FALSE]

merged <- merge(dti_sub, case_sub, by = "eid")
merged[, is_case := as.integer(get(STATUS_COL) %in% CASE_VALS)]

cat(sprintf("  Merged: %d participants (%d %s, %d HC)\n",
            nrow(merged), sum(merged$is_case),
            cfg$disease_short_caps, sum(!merged$is_case)))
cat(sprintf("  %s breakdown:\n", cfg$disease_short_caps))
print(table(merged[[STATUS_COL]]))

# ── 2. Interaction model: protein × case_status on each metric ───────────────
cat("\nRunning interaction models...\n")
cat(sprintf("  %d proteins × %d FA + %d MD + %d ICVF tracts\n",
            length(avail_prots),
            length(fa_cols), length(md_cols), length(icvf_cols)))

all_metric_cols <- c(fa_cols, md_cols, icvf_cols)

run_interaction <- function(dt, prot, metric_col) {
    sub <- dt[!is.na(get(prot)) & !is.na(get(metric_col)) &
              !is.na(age_at_mri) & !is.na(sex)]
    if (nrow(sub) < MIN_N) return(NULL)
    if (sum(sub$is_case) < 5) return(NULL)

    form <- as.formula(sprintf(
        "%s ~ scale(%s) * is_case + age_at_mri + sex",
        metric_col, prot))

    fit <- tryCatch(lm(form, data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(cf)) return(NULL)

    # Interaction term: scale(prot):is_case
    int_pat <- paste0("scale\\(", prot, "\\):is_case")
    int_row <- grep(int_pat, rownames(cf))
    if (length(int_row) == 0) return(NULL)

    # Main effect of protein (in HC, i.e. is_case=0)
    main_pat <- paste0("^scale\\(", prot, "\\)$")
    main_row <- grep(main_pat, rownames(cf))

    data.table(
        protein      = prot,
        metric       = metric_col,
        metric_type  = sub("_.*", "", metric_col),   # "fa", "md", "icvf"
        int_beta     = cf[int_row, "Estimate"],
        int_se       = cf[int_row, "Std. Error"],
        int_pval     = cf[int_row, "Pr(>|t|)"],
        main_beta_hc = if (length(main_row)) cf[main_row, "Estimate"] else NA_real_,
        n_total      = nrow(sub),
        n_case       = sum(sub$is_case)
    )
}

res_list <- vector("list", length(avail_prots) * length(all_metric_cols))
idx <- 1L
for (prot in avail_prots) {
    for (mc in all_metric_cols) {
        res_list[[idx]] <- run_interaction(merged, prot, mc)
        idx <- idx + 1L
    }
}

res_dt <- rbindlist(Filter(Negate(is.null), res_list))
res_dt[, fdr := p.adjust(int_pval, method = "BH")]
res_dt <- res_dt[order(int_pval)]

cat(sprintf("  Tests run: %d\n", nrow(res_dt)))
cat(sprintf("  FDR<0.05: %d  |  FDR<0.10: %d  |  FDR<0.20: %d\n",
            sum(res_dt$fdr < 0.05),
            sum(res_dt$fdr < 0.10),
            sum(res_dt$fdr < 0.20)))
cat(sprintf("  p<0.001: %d  |  p<0.01: %d\n",
            sum(res_dt$int_pval < 0.001),
            sum(res_dt$int_pval < 0.01)))

cat("\nTop 20 interactions (by p-value):\n")
print(res_dt[1:min(20, .N),
             .(protein = toupper(protein),
               metric  = sub("fa_|md_|icvf_", "", metric),
               type    = metric_type,
               int_beta = round(int_beta, 4),
               pval    = signif(int_pval, 3),
               fdr     = round(fdr, 3))])

fwrite(res_dt[, .(protein = toupper(protein), metric, metric_type,
                  int_beta, int_se, int_pval, fdr, main_beta_hc,
                  n_total, n_case)],
       file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_dti_interaction_results.csv")))
cat(sprintf("  Saved: %s_dti_interaction_results.csv\n", cfg$cohort_short))

# ── 3. Volcano: interaction beta vs -log10(p) ─────────────────────────────────
cat("\nPlotting interaction volcano...\n")

# Focus on FA for main figure (MD/ICVF in supplement)
fa_res <- res_dt[metric_type == "fa"]
fa_res[, neg_log10p := -log10(int_pval)]
fa_res[, sig_tier := fcase(
    fdr  < 0.05, "FDR<0.05",
    fdr  < 0.20, "FDR<0.20",
    int_pval < 0.01, "p<0.01",
    default = "NS"
)]

# Label top hits
fa_res[, label := NA_character_]
fa_res[fdr < 0.05, label := toupper(protein)]
if (sum(!is.na(fa_res$label)) == 0) {
    # If no FDR hits, label top 10 by p
    fa_res[1:min(10, .N), label := toupper(protein)]
}

TIER_COLS  <- c("FDR<0.05" = "#CC0066", "FDR<0.20" = "#E6A817",
                "p<0.01"   = "#2B4C7E", "NS"       = "grey75")
TIER_SIZE  <- c("FDR<0.05" = 2.0, "FDR<0.20" = 1.5, "p<0.01" = 1.2, "NS" = 0.7)
TIER_ALPHA <- c("FDR<0.05" = 1.0, "FDR<0.20" = 0.8, "p<0.01" = 0.7, "NS" = 0.3)

p_vol <- ggplot(fa_res, aes(x = int_beta, y = neg_log10p,
                             colour = sig_tier, size = sig_tier, alpha = sig_tier)) +
    geom_hline(yintercept = -log10(0.05 / nrow(fa_res)),
               linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point() +
    geom_text_repel(aes(label = label), size = 2.2, na.rm = TRUE,
                    max.overlaps = 25, segment.size = 0.3,
                    segment.colour = "grey50", min.segment.length = 0.2) +
    scale_colour_manual(values = TIER_COLS, name = NULL) +
    scale_size_manual(values = TIER_SIZE, guide = "none") +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none") +
    labs(x = glue::glue("Interaction β (protein × {cfg$disease_short_caps}, FA scale)"),
         y = expression(-log[10](p[interaction])),
         title = glue::glue("Protein × {cfg$disease_short_caps} interaction on white matter FA"),
         subtitle = sprintf("n=%d total (%d %s, %d HC) | %d DEPs × %d FA tracts",
                            nrow(merged), sum(merged$is_case),
                            cfg$disease_short_caps, sum(!merged$is_case),
                            length(avail_prots), length(fa_cols)),
         caption = glue::glue("Dotted line = Bonferroni; positive β = protein-FA slope steeper in {cfg$disease_short_caps}")) +
    theme_ukb() +
    theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_dti_interaction_volcano.pdf")), p_vol,
       width = 7, height = 6, device = cairo_pdf)
cat(sprintf("  Saved: %s_dti_interaction_volcano.pdf\n", cfg$cohort_short))

# ── 4. Heatmap — top proteins across all FA tracts ────────────────────────────
cat("Creating interaction heatmap...\n")

# Select proteins with at least one FA association p<0.01, capped at the top
# 30 by min p per protein.  Without the cap, cohorts with many DEPs (stroke)
# blow past ggsave's 50-inch heatmap-width limit.
HEATMAP_MAX <- 30L
prot_min_p  <- fa_res[, .(min_p = min(int_pval)), by = protein][order(min_p)]
top_prots_fa <- prot_min_p[min_p < 0.01, protein]
if (length(top_prots_fa) < 5) {
    top_prots_fa <- prot_min_p[1:min(15L, .N), protein]
} else if (length(top_prots_fa) > HEATMAP_MAX) {
    cat(sprintf("  %d proteins met p<0.01; capping heatmap at top %d.\n",
                length(top_prots_fa), HEATMAP_MAX))
    top_prots_fa <- prot_min_p[1:HEATMAP_MAX, protein]
}
cat(sprintf("  Proteins for heatmap: %d\n", length(top_prots_fa)))

heat_dt <- fa_res[protein %in% top_prots_fa]

# Tract short labels
tract_short <- sub("^fa_", "", fa_cols)
tract_short <- gsub("_", " ", tract_short)
heat_dt[, tract_label := gsub("^fa_", "", metric)]
heat_dt[, tract_label := gsub("_", " ", tract_label)]

# Order proteins by mean |interaction beta|
prot_order <- heat_dt[, mean(abs(int_beta)), by=protein][order(-V1), protein]
heat_dt[, prot_f   := factor(toupper(protein), levels = rev(toupper(prot_order)))]
heat_dt[, tract_f  := factor(tract_label, levels = unique(tract_label))]

p_heat <- ggplot(heat_dt, aes(x = prot_f, y = tract_f, fill = int_beta)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(data = heat_dt[int_pval < 0.01],
              aes(label = ifelse(fdr < 0.05, "**",
                                 ifelse(int_pval < 0.001, "*", "·"))),
              colour = "white", size = 2.5, vjust = 0.75) +
    scale_fill_gradient2(low = "#56B4E9", mid = "grey95", high = "#CC0066",
                         midpoint = 0, name = "Interaction β",
                         limits = c(-0.025, 0.025), oob = scales::squish) +
    labs(x = NULL, y = NULL,
         title = glue::glue("Protein × {cfg$disease_short_caps} interaction on FA (top proteins)"),
         caption = glue::glue("** FDR<0.05; * p<0.001; · p<0.01 | positive = steeper slope in {cfg$disease_short_caps}")) +
    theme_ukb() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          legend.position = "right")

ggsave(file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_dti_interaction_heatmap.pdf")), p_heat,
       width = max(5, length(top_prots_fa) * 0.5 + 3), height = 8,
       device = cairo_pdf)
cat(sprintf("  Saved: %s_dti_interaction_heatmap.pdf\n", cfg$cohort_short))

# ── 5. Scatter panels for top 6 interaction hits ─────────────────────────────
cat("Creating scatter panels for top hits...\n")

top6 <- res_dt[metric_type == "fa"][1:min(6, .N)]

case_label <- cfg$disease_short_caps
grp_levels <- c("HC", case_label)
grp_colors <- setNames(c("grey60", "#CC0066"), grp_levels)

scatter_list <- lapply(seq_len(nrow(top6)), function(i) {
    prot <- top6$protein[i]
    mc   <- top6$metric[i]
    if (!prot %in% names(merged) || !mc %in% names(merged)) return(NULL)
    sub <- merged[!is.na(get(prot)) & !is.na(get(mc))]
    if (nrow(sub) < 20) return(NULL)

    sub[, grp := factor(fifelse(is_case == 1, case_label, "HC"), levels = grp_levels)]

    # Tract short label
    tract_lab <- gsub("_", " ", sub("^fa_", "", mc))

    # Interaction annotation
    ib   <- round(top6$int_beta[i], 4)
    pval <- signif(top6$int_pval[i], 2)

    ggplot(sub, aes(x = get(prot), y = get(mc), colour = grp)) +
        geom_point(alpha = 0.35, size = 0.8) +
        geom_smooth(aes(group = grp), method = "lm", se = TRUE,
                    linewidth = 0.8, alpha = 0.15) +
        scale_colour_manual(values = grp_colors, name = NULL) +
        annotate("text", x = -Inf, y = Inf,
                 label = sprintf("β_int=%.4f\np=%s", ib, pval),
                 hjust = -0.1, vjust = 1.3, size = 2.2, colour = "grey30") +
        labs(x = sprintf("%s (NPX)", toupper(prot)),
             y = sprintf("FA: %s", tract_lab)) +
        theme_ukb(base_size = 8) +
        theme(legend.position = "top",
              legend.key.size = unit(0.3, "cm"))
})

scatter_list <- Filter(Negate(is.null), scatter_list)
if (length(scatter_list) > 0) {
    p_scatter <- wrap_plots(scatter_list, ncol = 3)
    ggsave(file.path(OUT_DIR, glue::glue("{cfg$cohort_short}_dti_interaction_scatter.pdf")), p_scatter,
           width = 12, height = 4 * ceiling(length(scatter_list) / 3),
           device = cairo_pdf)
    cat(sprintf("  Saved: %s_dti_interaction_scatter.pdf\n", cfg$cohort_short))
}

# ── 6. Summary across metric types ───────────────────────────────────────────
cat("\n=== Summary by metric type ===\n")
for (mt in c("fa", "md", "icvf")) {
    sub <- res_dt[metric_type == mt]
    cat(sprintf("  %s: %d tests | FDR<0.05: %d | FDR<0.20: %d | p<0.01: %d\n",
                toupper(mt), nrow(sub),
                sum(sub$fdr < 0.05),
                sum(sub$fdr < 0.20),
                sum(sub$int_pval < 0.01)))
}

# ── 7. Compare: top interaction proteins vs overall Cox proteins ──────────────
cox_file <- here::here("results", "survival",
                       glue::glue("{cfg$cohort_short}_protein_cox_results.csv"))
if (file.exists(cox_file)) {
    cox <- fread(cox_file)
    cox_sig <- tolower(cox[fdr < 0.05, protein])
    int_sig <- tolower(res_dt[fdr < 0.05, unique(protein)])
    int_nom <- tolower(res_dt[int_pval < 0.01, unique(protein)])
    cat(sprintf("\n=== Overlap with Cox FDR<0.05 proteins (%d) ===\n", length(cox_sig)))
    cat(sprintf("  Cox ∩ interaction FDR<0.05 (%d): %s\n",
                length(intersect(cox_sig, int_sig)),
                paste(toupper(intersect(cox_sig, int_sig)), collapse=", ")))
    cat(sprintf("  Cox ∩ interaction p<0.01 (%d): %s\n",
                length(intersect(cox_sig, int_nom)),
                paste(toupper(intersect(cox_sig, int_nom)), collapse=", ")))
}

cat("\n04_dti_protein_interaction.R complete.\n")
