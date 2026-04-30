#!/usr/bin/env Rscript
# 02_brain_mri_dti_ms.R
# DTI white matter FA analysis in MS — UKB-unique extension
#
# Approach:
#   1. Load DTI FA data (27 tracts, n=5619 UKB participants)
#   2. MS vs HC comparison for each FA tract (Mann-Whitney U + FDR correction)
#   3. Spearman correlation: MS DEPs (NEFL, GFAP, MOG + all 173) vs FA tracts
#   4. Mediation-style analysis: does NEFL mediate MS → lower FA?
#
# Input:
#   data/ukb/brain_mri/dti/dti_tract_protein_data.tsv   (DTI FA, labelled)
#   data/ukb/brain_mri/dti/dti_tract_covariates.tsv     (age, sex at scan)
#   data/ukb/olink/processed/ms_olink_qc.csv
#   results/differential/ms_combined_vs_hc.csv
#
# Output: results/mri/
#   ms_dti_fa_comparison.pdf
#   ms_dti_protein_heatmap.pdf
#   ms_dti_nefl_scatter.pdf
#   ms_dti_results.csv

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
DTI_FILE    <- file.path(CADASIL_DIR, "brain_mri", "dti", "dti_tract_protein_data.tsv")
DTI_COV_FILE <- file.path(CADASIL_DIR, "brain_mri", "dti", "dti_tract_covariates.tsv")
QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DIFF_FILE   <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
OUT_DIR     <- file.path(PROJ_DIR, "results", "mri")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR <- 0.05

# ── 1. Load DTI data ──────────────────────────────────────────────────────────
cat("Loading DTI data...\n")
dti <- fread(DTI_FILE, sep = "\t", showProgress = FALSE)
dti_cov <- fread(DTI_COV_FILE, sep = "\t", showProgress = FALSE)
# dti_cov has age at scan (p21003_i2) and sex (p31)
setnames(dti_cov, "p31", "sex_scan", skip_absent = TRUE)
setnames(dti_cov, "p21003_i2", "age_at_scan", skip_absent = TRUE)

fa_cols <- grep("^fa_", names(dti), value = TRUE)
cat(sprintf("  DTI: %d participants, %d FA tracts\n", nrow(dti), length(fa_cols)))

# Friendly tract names
tract_labels <- c(
    fa_acoustic_radiation_left              = "Acoustic radiation (L)",
    fa_acoustic_radiation_right             = "Acoustic radiation (R)",
    fa_anterior_thalamic_radiation_left     = "Ant. thalamic radiation (L)",
    fa_anterior_thalamic_radiation_right    = "Ant. thalamic radiation (R)",
    fa_cingulate_gyrus_cingulum_left        = "Cingulate cingulum (L)",
    fa_cingulate_gyrus_cingulum_right       = "Cingulate cingulum (R)",
    fa_parahippocampal_cingulum_left        = "Parahippo. cingulum (L)",
    fa_parahippocampal_cingulum_right       = "Parahippo. cingulum (R)",
    fa_corticospinal_tract_left             = "Corticospinal (L)",
    fa_corticospinal_tract_right            = "Corticospinal (R)",
    fa_forceps_major                        = "Forceps major",
    fa_forceps_minor                        = "Forceps minor",
    fa_inferior_fronto_occipital_fasciculus_left  = "IFO (L)",
    fa_inferior_fronto_occipital_fasciculus_right = "IFO (R)",
    fa_inferior_longitudinal_fasciculus_left      = "ILF (L)",
    fa_inferior_longitudinal_fasciculus_right     = "ILF (R)",
    fa_middle_cerebellar_peduncle           = "Mid. cerebellar peduncle",
    fa_medial_lemniscus_left                = "Medial lemniscus (L)",
    fa_medial_lemniscus_right               = "Medial lemniscus (R)",
    fa_posterior_thalamic_radiation_left    = "Post. thalamic radiation (L)",
    fa_posterior_thalamic_radiation_right   = "Post. thalamic radiation (R)",
    fa_superior_longitudinal_fasciculus_left  = "SLF (L)",
    fa_superior_longitudinal_fasciculus_right = "SLF (R)",
    fa_superior_thalamic_radiation_left     = "Sup. thalamic radiation (L)",
    fa_superior_thalamic_radiation_right    = "Sup. thalamic radiation (R)",
    fa_uncinate_fasciculus_left             = "Uncinate fasciculus (L)",
    fa_uncinate_fasciculus_right            = "Uncinate fasciculus (R)"
)

# ── 2. Load MS cohort ─────────────────────────────────────────────────────────
cat("Loading MS cohort...\n")
ms_qc <- fread(QC_FILE, showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & !is.na(age_at_sampling) &
                !is.na(sex) & !is.na(UMAP1)]

# Merge DTI + covariates + MS status
ms_dti <- merge(ms_qc, dti, by = "eid", all.x = FALSE)  # inner join: only those with DTI
# Add scan age/sex
dti_meta <- dti_cov[, .(eid, age_at_scan, sex_scan)]
ms_dti   <- merge(ms_dti, dti_meta, by = "eid", all.x = TRUE)

cat(sprintf("  MS cases with DTI: %d\n", sum(ms_dti$ms_status %in% c("pre_onset","post_onset"))))
cat(sprintf("  Controls with DTI: %d\n", sum(ms_dti$ms_status == "control")))
cat(sprintf("  Status breakdown:\n"))
print(table(ms_dti$ms_status))

ms_dti[, is_ms := as.integer(ms_status %in% c("pre_onset","post_onset"))]

# ── 3. MS vs HC: DTI FA comparison ────────────────────────────────────────────
cat("\nMS vs HC DTI FA comparison (Mann-Whitney U)...\n")
fa_test_list <- lapply(fa_cols, function(fc) {
    sub_dt <- ms_dti[!is.na(get(fc)) & ms_status %in% c("pre_onset","post_onset","control")]
    if (nrow(sub_dt[ms_status != "control"]) < 5) return(NULL)
    wt <- wilcox.test(get(fc) ~ is_ms, data = sub_dt, exact = FALSE)
    ms_mean   <- mean(sub_dt[is_ms == 1, get(fc)], na.rm=TRUE)
    ctrl_mean <- mean(sub_dt[is_ms == 0, get(fc)], na.rm=TRUE)
    data.table(
        tract    = fc,
        label    = tract_labels[fc],
        ms_n     = sum(!is.na(sub_dt[is_ms==1, get(fc)])),
        ctrl_n   = sum(!is.na(sub_dt[is_ms==0, get(fc)])),
        ms_mean  = ms_mean,
        ctrl_mean = ctrl_mean,
        delta_fa = ms_mean - ctrl_mean,
        pval     = wt$p.value
    )
})
fa_test <- rbindlist(Filter(Negate(is.null), fa_test_list))
fa_test[, fdr := p.adjust(pval, method = "BH")]
fa_test[, sig_label := fcase(fdr < 0.01, "**", fdr < 0.05, "*", pval < 0.05, "+", default = "")]
fa_test <- fa_test[order(pval)]
fwrite(fa_test, file.path(OUT_DIR, "ms_dti_results.csv"))
cat(sprintf("  FDR<0.05: %d tracts; p<0.05: %d tracts\n",
            sum(fa_test$fdr < 0.05, na.rm=TRUE),
            sum(fa_test$pval < 0.05, na.rm=TRUE)))

# Plot: lollipop of delta FA
fa_test[, label_f := factor(label, levels = rev(fa_test$label))]
fa_test[, fill_col := fcase(
    fdr < 0.05 & delta_fa < 0, "sig_down",
    fdr < 0.05 & delta_fa > 0, "sig_up",
    default = "ns"
)]
fill_vals <- c(sig_down = "#2B4C7E", sig_up = "#CC0066", ns = "grey70")

p_fa_loll <- ggplot(fa_test, aes(x = delta_fa, y = label_f, colour = fill_col)) +
    geom_segment(aes(xend = 0, yend = label_f), linewidth = 0.5) +
    geom_point(aes(size = -log10(pval)), alpha = 0.9) +
    geom_text(data = fa_test[sig_label != ""],
              aes(x = delta_fa + sign(delta_fa) * 0.0005, label = sig_label),
              size = 3.5, vjust = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    scale_colour_manual(values = fill_vals,
                        labels = c(sig_down="MS<HC (FDR<0.05)",
                                   sig_up="MS>HC (FDR<0.05)", ns="n.s."),
                        name = NULL) +
    scale_size_continuous(name = expression(-log[10](p)), range = c(1.5, 5)) +
    labs(x = "Δ FA (MS − HC)", y = NULL,
         title = "White matter FA: MS vs healthy controls",
         subtitle = sprintf("MS n=%d, HC n=%d | Mann-Whitney U + BH FDR",
                            max(fa_test$ms_n), max(fa_test$ctrl_n))) +
    theme_ukb() +
    theme(axis.text.y = element_text(size = 7), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "ms_dti_fa_comparison.pdf"),
       p_fa_loll, width = 6.5, height = 6, device = cairo_pdf)
cat("  Saved: ms_dti_fa_comparison.pdf\n")

# ── 4. Protein–DTI correlation (MS cases only) ────────────────────────────────
cat("\nProtein–DTI correlations (MS cases)...\n")
deps <- fread(DIFF_FILE)
dep_prots <- tolower(deps[adj.P.Val < FDR_THR, protein])
avail_prots <- intersect(dep_prots, names(ms_dti))

ms_only <- ms_dti[ms_status %in% c("pre_onset","post_onset")]
cat(sprintf("  MS cases for correlation: %d\n", nrow(ms_only)))

corr_list <- list()
for (fc in fa_cols) {
    sub <- ms_only[!is.na(get(fc))]
    if (nrow(sub) < 20) next
    for (prot in avail_prots) {
        if (!prot %in% names(sub)) next
        vals <- sub[!is.na(get(prot)), .(fa = get(fc), npx = get(prot))]
        if (nrow(vals) < 15) next
        ct <- cor.test(vals$fa, vals$npx, method = "spearman", exact = FALSE)
        corr_list[[paste(prot, fc, sep="|")]] <- data.table(
            protein = prot,
            tract   = fc,
            rho     = ct$estimate,
            pval    = ct$p.value,
            n       = nrow(vals)
        )
    }
}

if (length(corr_list) > 0) {
    corr_dt <- rbindlist(corr_list)
    corr_dt[, fdr := p.adjust(pval, method = "BH")]
    corr_dt <- corr_dt[order(pval)]
    cat(sprintf("  %d protein–tract pairs; FDR<0.05: %d\n",
                nrow(corr_dt), sum(corr_dt$fdr < 0.05)))
    fwrite(corr_dt, file.path(OUT_DIR, "ms_dti_protein_correlations.csv"))

    # Top 10 most significant protein–tract pairs
    top_corr <- corr_dt[1:min(10, .N)]
    cat("  Top protein-tract correlations:\n")
    print(top_corr[, .(protein, tract = sub("fa_","",tract), rho, pval, fdr)])

    # Heatmap: key proteins × top tracts (by any significant hit)
    key_prots <- intersect(c("nefl","gfap","mog","uch_l1","nfl","sncb"),
                           unique(corr_dt[fdr < 0.05 | pval < 0.05, protein]))
    if (length(key_prots) == 0) key_prots <- corr_dt[1:min(6,.N), protein]

    heat_dt <- corr_dt[protein %in% key_prots]
    heat_dt[, label_f := factor(tract_labels[tract], levels = unique(tract_labels))]
    heat_dt[, prot_f := factor(toupper(protein), levels = toupper(rev(unique(key_prots))))]

    p_heat <- ggplot(heat_dt[!is.na(label_f)],
                     aes(x = prot_f, y = label_f, fill = rho)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(data = heat_dt[!is.na(label_f) & fdr < 0.05],
                  aes(label = "*"), colour = "white", size = 4, vjust = 0.75) +
        scale_fill_gradient2(low = "#2B4C7E", mid = "white", high = "#CC0066",
                             midpoint = 0, name = "Spearman ρ",
                             limits = c(-0.4, 0.4), oob = scales::squish) +
        labs(x = "Protein", y = NULL,
             title = "Protein–DTI FA correlations (MS cases)",
             caption = "* FDR<0.05") +
        theme_ukb() +
        theme(axis.text.y = element_text(size = 7),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right")

    ggsave(file.path(OUT_DIR, "ms_dti_protein_heatmap.pdf"),
           p_heat, width = max(4, length(key_prots) * 0.8 + 3), height = 7,
           device = cairo_pdf)
    cat("  Saved: ms_dti_protein_heatmap.pdf\n")
}

# ── 5. NEFL vs key tracts — scatter panels ────────────────────────────────────
cat("\nNEFL vs DTI FA scatter panels...\n")
nefl_col <- intersect(c("nefl","nfl"), names(ms_dti))[1]

if (!is.na(nefl_col)) {
    # Top 4 tracts most correlated with NEFL in MS cases
    top_tracts <- if (exists("corr_dt") && nrow(corr_dt) > 0) {
        corr_dt[protein == nefl_col][order(pval)][1:min(4,.N), tract]
    } else {
        fa_cols[1:4]
    }

    scatter_plots <- lapply(top_tracts, function(fc) {
        sub <- ms_dti[ms_status %in% c("pre_onset","post_onset") &
                      !is.na(get(fc)) & !is.na(get(nefl_col))]
        if (nrow(sub) < 10) return(NULL)
        rho  <- cor(sub[[nefl_col]], sub[[fc]], method = "spearman")
        pval <- cor.test(sub[[nefl_col]], sub[[fc]], method = "spearman", exact = FALSE)$p.value
        lab  <- ifelse(fc %in% names(tract_labels), tract_labels[fc], fc)

        col_by_status <- c(pre_onset = "#E6A817", post_onset = "#CC0066")

        ggplot(sub, aes(x = get(nefl_col), y = get(fc),
                        colour = ms_status)) +
            geom_point(alpha = 0.6, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, colour = "#2B4C7E",
                        linewidth = 0.7, fill = "#2B4C7E", alpha = 0.1) +
            scale_colour_manual(values = col_by_status, guide = "none") +
            annotate("text", x = Inf, y = Inf,
                     label = sprintf("ρ=%.2f, p=%s", rho,
                                     ifelse(pval < 0.001, "<0.001",
                                            sprintf("%.3f", pval))),
                     hjust = 1.05, vjust = 1.5, size = 2.8, colour = "grey30") +
            labs(x = "NfL (NPX)", y = sprintf("FA: %s", lab)) +
            theme_ukb(base_size = 8)
    })

    scatter_plots <- Filter(Negate(is.null), scatter_plots)
    if (length(scatter_plots) > 0) {
        p_nefl_scatter <- wrap_plots(scatter_plots, ncol = 2)
        ggsave(file.path(OUT_DIR, "ms_dti_nefl_scatter.pdf"),
               p_nefl_scatter, width = 7, height = max(3, ceiling(length(scatter_plots)/2)*3),
               device = cairo_pdf)
        cat("  Saved: ms_dti_nefl_scatter.pdf\n")
    }
}

cat("\n02_brain_mri_dti_ms.R complete.\n")
