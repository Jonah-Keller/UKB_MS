#!/usr/bin/env Rscript
# 01_ms_prs.R
# MS polygenic risk score analysis
#
# PRS scores: PGS000809 (Barnes 2021, 127 SNPs) + PGS004699 (Loginovic 2024, 307 SNPs Non-HLA)
# Analyses:
#   1. PRS distribution: MS cases vs controls (pre/post split)
#   2. PRS ~ DEP correlation in healthy controls (do high-PRS HC show prodromal elevations?)
#   3. PRS × time-to-diagnosis in pre-onset MS
#   4. PRS quartile × key proteins (NEFL, MOG, GFAP, LILRB4) dose-response
#   5. PRS × proteomic cluster cross-tabulation
#   6. PRS × HLA-DRB1*15:01 interaction
#
# Output: results/prs/

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

PRS_FILE     <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
QC_FILE      <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DIFF_FILE    <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
CLUST_FILE   <- file.path(PROJ_DIR, "results", "clustering", "ms_cluster_assignments.csv")
HLA_FILE     <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
OUT_DIR      <- file.path(PROJ_DIR, "results", "prs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GROUP_COLS <- c("Control"="grey65", "Pre-onset MS"="#E6A817", "Post-onset MS"="#CC0066")

# ── 1. Load and merge ─────────────────────────────────────────────────────────
cat("Loading data...\n")
prs   <- fread(PRS_FILE, showProgress = FALSE)
ms_qc <- fread(QC_FILE,  showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
               !is.na(UMAP1) & !is.na(UMAP2)]

# prs_score is already z-scored mean of the two PGS; use directly
dt <- merge(ms_qc, prs[, .(eid, prs_PGS000809, prs_PGS004699, prs_score)],
            by = "eid", all.x = FALSE)
dt[, ms_group := fcase(
    ms_status == "pre_onset",  "Pre-onset MS",
    ms_status == "post_onset", "Post-onset MS",
    ms_status == "control",    "Control",
    default = NA_character_
)]
dt <- dt[!is.na(ms_group)]

cat(sprintf("  Participants with PRS: %d\n", nrow(dt)))
cat("  Status:\n"); print(table(dt$ms_status))

# DEPs
deps        <- fread(DIFF_FILE, showProgress = FALSE)
dep_prots   <- tolower(deps[adj.P.Val < 0.05, protein])
avail_prots <- intersect(dep_prots, names(dt))
cat(sprintf("  DEPs available: %d\n", length(avail_prots)))

# ── 2. PRS distribution: MS vs HC ─────────────────────────────────────────────
cat("\nPRS distribution...\n")
ms_vec <- dt[ms_status %in% c("pre_onset","post_onset"), prs_score]
hc_vec <- dt[ms_status == "control", prs_score]
wt     <- wilcox.test(ms_vec, hc_vec, exact = FALSE)
d      <- (mean(ms_vec) - mean(hc_vec)) / sd(hc_vec)

cat(sprintf("  MS median=%.3f  HC median=%.3f  p=%.3g  Cohen d=%.3f\n",
            median(ms_vec), median(hc_vec), wt$p.value, d))

# Density plot (subsample HC for display)
dt_plot <- rbind(dt[ms_group != "Control"],
                 dt[ms_group == "Control"][sample(.N, min(5000L, .N))])

p_dist <- ggplot(dt_plot, aes(x = prs_score, fill = ms_group, colour = ms_group)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    scale_fill_manual(values = GROUP_COLS, name = NULL) +
    scale_colour_manual(values = GROUP_COLS, name = NULL) +
    annotate("text", x = 2.5, y = Inf,
             label = sprintf("p=%s\nd=%.2f", format(wt$p.value, digits=2), d),
             vjust = 1.3, hjust = 0, size = 2.8, colour = "grey30") +
    labs(x = "MS PRS (z-score, combined PGS000809+PGS004699)",
         y = "Density",
         title = "MS polygenic risk score distribution",
         subtitle = sprintf("MS n=%d, HC n=%d", length(ms_vec), length(hc_vec))) +
    theme_ukb() +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.key.size = unit(0.35, "cm"))

ggsave(file.path(OUT_DIR, "ms_prs_distribution.pdf"),
       p_dist, width = 5.5, height = 4, device = cairo_pdf)
cat("  Saved: ms_prs_distribution.pdf\n")

# ── 3. PRS quartile × key proteins ───────────────────────────────────────────
cat("\nPRS quartile × protein levels...\n")
dt[, prs_quartile := cut(prs_score, quantile(prs_score, probs=0:4/4),
                          labels=c("Q1 (low)","Q2","Q3","Q4 (high)"),
                          include.lowest=TRUE)]

key_prots <- intersect(c("nefl","mog","gfap","lilrb4","il3","aif1","tnr"),
                        names(dt))
hc_prs <- dt[ms_status == "control" & !is.na(prs_quartile)]

q_plots <- lapply(key_prots, function(prot) {
    sub <- hc_prs[!is.na(get(prot))]
    if (nrow(sub) < 50) return(NULL)
    # Kruskal test
    kw <- kruskal.test(get(prot) ~ prs_quartile, data = sub)
    ggplot(sub, aes(x = prs_quartile, y = get(prot), fill = prs_quartile)) +
        geom_violin(alpha = 0.6, linewidth = 0.3) +
        geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.3, linewidth = 0.35) +
        scale_fill_manual(values = c("#56B4E9","grey65","#E6A817","#CC0066"),
                          guide = "none") +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf("KW p=%s", format(kw$p.value, digits=2)),
                 hjust = 1.1, vjust = 1.5, size = 2.5, colour = "grey30") +
        labs(x = NULL, y = "NPX", title = toupper(prot)) +
        theme_ukb(base_size = 8) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
})
q_plots <- Filter(Negate(is.null), q_plots)

if (length(q_plots) > 0) {
    p_quartile <- wrap_plots(q_plots, ncol = min(4, length(q_plots)))
    ggsave(file.path(OUT_DIR, "ms_prs_quartile_proteins.pdf"), p_quartile,
           width = 14, height = ceiling(length(q_plots)/4) * 3.5,
           device = cairo_pdf)
    cat("  Saved: ms_prs_quartile_proteins.pdf\n")
}

# ── 4. PRS ~ DEP correlation in healthy controls ──────────────────────────────
cat("\nPRS ~ DEP correlation in healthy controls...\n")
hc_only <- dt[ms_status == "control"]
cat(sprintf("  HC with PRS: %d\n", nrow(hc_only)))

corr_list <- lapply(avail_prots, function(prot) {
    vals <- hc_only[!is.na(get(prot)) & !is.na(prs_score),
                    .(prs = prs_score, npx = get(prot))]
    if (nrow(vals) < 20) return(NULL)
    ct <- cor.test(vals$prs, vals$npx, method = "spearman", exact = FALSE)
    data.table(protein=prot, rho=ct$estimate, pval=ct$p.value, n=nrow(vals))
})
corr_dt <- rbindlist(Filter(Negate(is.null), corr_list))
corr_dt[, fdr := p.adjust(pval, method = "BH")]
corr_dt <- corr_dt[order(pval)]

cat(sprintf("  FDR<0.05: %d | FDR<0.20: %d | p<0.01: %d\n",
            sum(corr_dt$fdr < 0.05),
            sum(corr_dt$fdr < 0.20),
            sum(corr_dt$pval < 0.01)))
cat("\nTop 10 PRS-correlated proteins in HC:\n")
print(corr_dt[1:min(10,.N),
              .(protein=toupper(protein), rho=round(rho,3),
                pval=signif(pval,3), fdr=round(fdr,3))])

corr_dt[, label_prot := fifelse(fdr < 0.05 | rank(pval) <= 10,
                                toupper(protein), NA_character_)]
corr_dt[, direction := fifelse(rho > 0, "Higher in high PRS", "Lower in high PRS")]

p_corr <- ggplot(corr_dt, aes(x = rho, y = -log10(pval), colour = direction)) +
    geom_point(aes(alpha = fdr < 0.20, size = fdr < 0.05)) +
    scale_alpha_manual(values = c("TRUE"=0.85, "FALSE"=0.25), guide = "none") +
    scale_size_manual(values  = c("TRUE"=2.5,  "FALSE"=1.2),  guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype="dashed", colour="grey50", linewidth=0.3) +
    geom_text_repel(aes(label = label_prot), size = 2.2, na.rm=TRUE,
                    max.overlaps=20, segment.size=0.3) +
    scale_colour_manual(values = c("Higher in high PRS"="#CC0066",
                                   "Lower in high PRS"="#56B4E9"),
                        name = NULL) +
    labs(x = "Spearman ρ (PRS ~ protein in HC)",
         y = expression(-log[10](p)),
         title = "MS PRS correlates with DEP levels in healthy controls",
         subtitle = sprintf("n≈%d HC | do high-genetic-risk controls show prodromal elevations?",
                            max(corr_dt$n))) +
    theme_ukb() +
    theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "ms_prs_protein_correlation.pdf"),
       p_corr, width = 6, height = 5, device = cairo_pdf)
cat("  Saved: ms_prs_protein_correlation.pdf\n")

# ── 5. PRS × time-to-diagnosis in pre-onset MS ────────────────────────────────
cat("\nPRS × time-to-diagnosis (pre-onset MS)...\n")
pre_dt <- dt[ms_status == "pre_onset" & !is.na(prs_score) & !is.na(years_to_diagnosis)]
cat(sprintf("  Pre-onset with PRS: %d\n", nrow(pre_dt)))

if (nrow(pre_dt) >= 20) {
    lm_fit <- lm(years_to_diagnosis ~ prs_score + age_at_sampling + sex,
                 data = pre_dt)
    cf <- coef(summary(lm_fit))
    cat(sprintf("  PRS ~ years_to_dx: β=%.3f, SE=%.3f, p=%.3g\n",
                cf["prs_score","Estimate"], cf["prs_score","Std. Error"],
                cf["prs_score","Pr(>|t|)"]))

    p_ttd <- ggplot(pre_dt, aes(x = prs_score, y = years_to_diagnosis)) +
        geom_point(colour = "#E6A817", alpha = 0.7, size = 2) +
        geom_smooth(method = "lm", se = TRUE, colour = "#2B4C7E",
                    linewidth = 0.8, fill = "#2B4C7E", alpha = 0.15) +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf("β=%.3f\np=%.3g", cf["prs_score","Estimate"],
                                 cf["prs_score","Pr(>|t|)"]),
                 hjust=1.1, vjust=1.5, size=3, colour="grey30") +
        labs(x = "MS PRS (z-score)", y = "Years to diagnosis",
             title = "MS PRS vs lead time to diagnosis",
             subtitle = sprintf("Pre-onset MS n=%d", nrow(pre_dt)),
             caption = "Adjusted for age at sampling and sex") +
        theme_ukb()

    ggsave(file.path(OUT_DIR, "ms_prs_time_to_diagnosis.pdf"), p_ttd,
           width = 4.5, height = 4, device = cairo_pdf)
    cat("  Saved: ms_prs_time_to_diagnosis.pdf\n")
}

# ── 6. PRS × proteomic cluster ────────────────────────────────────────────────
cat("\nPRS × proteomic cluster...\n")
if (file.exists(CLUST_FILE)) {
    clust <- fread(CLUST_FILE)[, .(eid, cluster)]
    dt_clust <- merge(dt[ms_status %in% c("pre_onset","post_onset")],
                      clust, by = "eid", all.x = TRUE)
    dt_clust <- dt_clust[!is.na(cluster)]
    cat(sprintf("  MS cases with cluster + PRS: %d\n", nrow(dt_clust)))

    wt_clust <- wilcox.test(prs_score ~ cluster, data = dt_clust, exact=FALSE)
    cat(sprintf("  Cluster 1 vs 2 PRS: median %.3f vs %.3f  p=%.3g\n",
                median(dt_clust[cluster==1, prs_score]),
                median(dt_clust[cluster==2, prs_score]),
                wt_clust$p.value))

    p_clust <- ggplot(dt_clust, aes(x = factor(cluster), y = prs_score,
                                     fill = factor(cluster))) +
        geom_violin(alpha = 0.6, linewidth = 0.3) +
        geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5, linewidth=0.35) +
        scale_fill_manual(values = c("1"="#CC0066","2"="#56B4E9"), guide="none") +
        annotate("text", x=1.5, y=Inf, label=sprintf("p=%.3g", wt_clust$p.value),
                 vjust=1.5, size=3, colour="grey30") +
        labs(x = "Proteomic cluster", y = "MS PRS (z-score)",
             title = "PRS by proteomic cluster (MS cases)") +
        theme_ukb()

    ggsave(file.path(OUT_DIR, "ms_prs_cluster.pdf"), p_clust,
           width = 3.5, height = 4, device = cairo_pdf)
    cat("  Saved: ms_prs_cluster.pdf\n")
}

# ── 7. PRS × HLA-DRB1*15:01 ──────────────────────────────────────────────────
cat("\nPRS × HLA-DRB1*15:01 interaction...\n")
if (file.exists(HLA_FILE)) {
    hla <- fread(HLA_FILE, select = c("eid","drb1_1501_carrier"),
                 showProgress = FALSE)
    dt_hla <- merge(dt, hla, by = "eid", all.x = TRUE)
    dt_hla <- dt_hla[!is.na(drb1_1501_carrier)]

    # 2×2 table in MS cases
    ms_hla <- dt_hla[ms_status %in% c("pre_onset","post_onset")]
    cat(sprintf("  MS cases with HLA + PRS: %d\n", nrow(ms_hla)))

    # PRS by HLA carrier status in MS
    wt_hla <- wilcox.test(prs_score ~ drb1_1501_carrier, data=ms_hla, exact=FALSE)
    cat(sprintf("  HLA carrier PRS: %.3f  non-carrier: %.3f  p=%.3g\n",
                median(ms_hla[drb1_1501_carrier==1, prs_score]),
                median(ms_hla[drb1_1501_carrier==0, prs_score]),
                wt_hla$p.value))

    # Also in HC
    hc_hla <- dt_hla[ms_status == "control"]
    wt_hc  <- wilcox.test(prs_score ~ drb1_1501_carrier, data=hc_hla, exact=FALSE)
    cat(sprintf("  HC: HLA carrier PRS: %.3f  non-carrier: %.3f  p=%.3g\n",
                median(hc_hla[drb1_1501_carrier==1, prs_score]),
                median(hc_hla[drb1_1501_carrier==0, prs_score]),
                wt_hc$p.value))

    dt_hla_plot <- dt_hla[ms_status %in% c("pre_onset","post_onset","control")]
    dt_hla_plot[, hla_label := fifelse(drb1_1501_carrier==1, "DRB1*15:01 carrier",
                                                              "Non-carrier")]
    dt_hla_plot[, group_hla := paste(ms_group, hla_label, sep="\n")]

    # Compact: PRS by HLA within MS cases only
    ms_hla[, hla_label := fifelse(drb1_1501_carrier==1, "Carrier", "Non-carrier")]
    p_hla <- ggplot(ms_hla, aes(x=hla_label, y=prs_score, fill=hla_label)) +
        geom_violin(alpha=0.6, linewidth=0.3) +
        geom_boxplot(width=0.12, fill="white", outlier.size=0.5, linewidth=0.35) +
        scale_fill_manual(values=c("Carrier"="#CC0066","Non-carrier"="#56B4E9"),
                          guide="none") +
        annotate("text", x=1.5, y=Inf,
                 label=sprintf("p=%.3g", wt_hla$p.value),
                 vjust=1.5, size=3, colour="grey30") +
        labs(x="HLA-DRB1*15:01", y="MS PRS (z-score)",
             title="PRS by HLA-DRB1*15:01 status (MS cases)",
             subtitle=sprintf("n=%d carriers, %d non-carriers",
                              sum(ms_hla$drb1_1501_carrier), sum(!ms_hla$drb1_1501_carrier))) +
        theme_ukb()

    ggsave(file.path(OUT_DIR, "ms_prs_hla.pdf"), p_hla,
           width=3.5, height=4, device=cairo_pdf)
    cat("  Saved: ms_prs_hla.pdf\n")
}

# ── 8. Save results ───────────────────────────────────────────────────────────
fwrite(corr_dt[, .(protein=toupper(protein), rho, pval, fdr, n)],
       file.path(OUT_DIR, "ms_prs_results.csv"))

cat("\n01_ms_prs.R complete.\n")
