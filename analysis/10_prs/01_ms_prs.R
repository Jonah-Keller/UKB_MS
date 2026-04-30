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
prs_label    <- cfg$prs_label
disease_caps <- cfg$disease_short_caps

PRS_FILE     <- here::here("data", "ukb", "genetics",
                           glue::glue("{cohort}_prs_scores.csv"))
QC_FILE      <- here::here("data", "ukb", "olink", "processed",
                           glue::glue("{cohort}_olink_qc.csv"))
DIFF_FILE    <- here::here("results", "differential",
                           glue::glue("{cohort}_combined_vs_hc.csv"))
CLUST_FILE   <- here::here("results", "clustering",
                           glue::glue("{cohort}_cluster_assignments.csv"))
HLA_FILE     <- here::here("data", "ukb", "genetics", "hla_drb1_imputed.csv")
OUT_DIR      <- here::here("results", "prs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pre_label    <- glue::glue("Pre-onset {disease_caps}")
post_label   <- glue::glue("Post-onset {disease_caps}")
GROUP_COLS   <- setNames(c("grey65", "#E6A817", "#CC0066"),
                         c("Control", pre_label, post_label))

# ── 1. Load and merge ─────────────────────────────────────────────────────────
cat("Loading data...\n")
prs   <- fread(PRS_FILE, showProgress = FALSE)
ms_qc <- fread(QC_FILE,  showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex) &
               !is.na(UMAP1) & !is.na(UMAP2)]

# Combined PRS column is already z-scored mean across PGS IDs; use directly
per_pgs_cols <- intersect(paste0("prs_", cfg$prs_pgs_ids), names(prs))
prs_keep     <- unique(c("eid", per_pgs_cols, prs_col))
dt <- merge(ms_qc, prs[, ..prs_keep], by = "eid", all.x = FALSE)
dt[, ms_group := fcase(
    get(status_col) == sv$pre_onset,  pre_label,
    get(status_col) == sv$post_onset, post_label,
    get(status_col) == sv$control,    "Control",
    default = NA_character_
)]
dt <- dt[!is.na(ms_group)]

cat(sprintf("  Participants with PRS: %d\n", nrow(dt)))
cat("  Status:\n"); print(table(dt[[status_col]]))

# DEPs
deps        <- fread(DIFF_FILE, showProgress = FALSE)
dep_prots   <- tolower(deps[adj.P.Val < 0.05, protein])
avail_prots <- intersect(dep_prots, names(dt))
cat(sprintf("  DEPs available: %d\n", length(avail_prots)))

# ── 2. PRS distribution: cases vs HC ─────────────────────────────────────────
cat("\nPRS distribution...\n")
case_states <- c(sv$pre_onset, sv$post_onset)
ms_vec <- dt[get(status_col) %in% case_states, get(prs_col)]
hc_vec <- dt[get(status_col) == sv$control, get(prs_col)]
wt     <- wilcox.test(ms_vec, hc_vec, exact = FALSE)
d      <- (mean(ms_vec) - mean(hc_vec)) / sd(hc_vec)

cat(sprintf("  %s median=%.3f  HC median=%.3f  p=%.3g  Cohen d=%.3f\n",
            disease_caps, median(ms_vec), median(hc_vec), wt$p.value, d))

# Density plot (subsample HC for display)
dt_plot <- rbind(dt[ms_group != "Control"],
                 dt[ms_group == "Control"][sample(.N, min(5000L, .N))])

p_dist <- ggplot(dt_plot, aes(x = get(prs_col), fill = ms_group, colour = ms_group)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    scale_fill_manual(values = GROUP_COLS, name = NULL) +
    scale_colour_manual(values = GROUP_COLS, name = NULL) +
    annotate("text", x = 2.5, y = Inf,
             label = sprintf("p=%s\nd=%.2f", format(wt$p.value, digits=2), d),
             vjust = 1.3, hjust = 0, size = 2.8, colour = "grey30") +
    labs(x = glue::glue("{prs_label} (z-score, combined {paste(cfg$prs_pgs_ids, collapse='+')})"),
         y = "Density",
         title = glue::glue("{disease_caps} polygenic risk score distribution"),
         subtitle = sprintf("%s n=%d, HC n=%d", disease_caps, length(ms_vec), length(hc_vec))) +
    theme_ukb() +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.key.size = unit(0.35, "cm"))

ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_distribution.pdf")),
       p_dist, width = 5.5, height = 4, device = cairo_pdf)
cat(sprintf("  Saved: %s_prs_distribution.pdf\n", cohort))

# ── 3. PRS quartile × key proteins ───────────────────────────────────────────
cat("\nPRS quartile × protein levels...\n")
dt[, prs_quartile := cut(get(prs_col), quantile(get(prs_col), probs=0:4/4),
                          labels=c("Q1 (low)","Q2","Q3","Q4 (high)"),
                          include.lowest=TRUE)]

key_prots <- intersect(c("nefl","mog","gfap","lilrb4","il3","aif1","tnr"),
                        names(dt))
hc_prs <- dt[get(status_col) == sv$control & !is.na(prs_quartile)]

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
    ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_quartile_proteins.pdf")),
           p_quartile,
           width = 14, height = ceiling(length(q_plots)/4) * 3.5,
           device = cairo_pdf)
    cat(sprintf("  Saved: %s_prs_quartile_proteins.pdf\n", cohort))
}

# ── 4. PRS ~ DEP correlation in healthy controls ──────────────────────────────
cat("\nPRS ~ DEP correlation in healthy controls...\n")
hc_only <- dt[get(status_col) == sv$control]
cat(sprintf("  HC with PRS: %d\n", nrow(hc_only)))

corr_list <- lapply(avail_prots, function(prot) {
    vals <- hc_only[!is.na(get(prot)) & !is.na(get(prs_col)),
                    .(prs = get(prs_col), npx = get(prot))]
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
         title = glue::glue("{prs_label} correlates with DEP levels in healthy controls"),
         subtitle = sprintf("n≈%d HC | do high-genetic-risk controls show prodromal elevations?",
                            max(corr_dt$n))) +
    theme_ukb() +
    theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_protein_correlation.pdf")),
       p_corr, width = 6, height = 5, device = cairo_pdf)
cat(sprintf("  Saved: %s_prs_protein_correlation.pdf\n", cohort))

# ── 5. PRS × time-to-diagnosis in pre-onset cases ─────────────────────────────
cat(sprintf("\nPRS × time-to-diagnosis (pre-onset %s)...\n", disease_caps))
pre_dt <- dt[get(status_col) == sv$pre_onset &
             !is.na(get(prs_col)) & !is.na(years_to_diagnosis)]
cat(sprintf("  Pre-onset with PRS: %d\n", nrow(pre_dt)))

if (nrow(pre_dt) >= 20) {
    lm_formula <- as.formula(paste("years_to_diagnosis ~",
                                   prs_col, "+ age_at_sampling + sex"))
    lm_fit <- lm(lm_formula, data = pre_dt)
    cf <- coef(summary(lm_fit))
    cat(sprintf("  PRS ~ years_to_dx: β=%.3f, SE=%.3f, p=%.3g\n",
                cf[prs_col,"Estimate"], cf[prs_col,"Std. Error"],
                cf[prs_col,"Pr(>|t|)"]))

    p_ttd <- ggplot(pre_dt, aes(x = get(prs_col), y = years_to_diagnosis)) +
        geom_point(colour = "#E6A817", alpha = 0.7, size = 2) +
        geom_smooth(method = "lm", se = TRUE, colour = "#2B4C7E",
                    linewidth = 0.8, fill = "#2B4C7E", alpha = 0.15) +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf("β=%.3f\np=%.3g", cf[prs_col,"Estimate"],
                                 cf[prs_col,"Pr(>|t|)"]),
                 hjust=1.1, vjust=1.5, size=3, colour="grey30") +
        labs(x = glue::glue("{prs_label} (z-score)"), y = "Years to diagnosis",
             title = glue::glue("{prs_label} vs lead time to diagnosis"),
             subtitle = sprintf("Pre-onset %s n=%d", disease_caps, nrow(pre_dt)),
             caption = "Adjusted for age at sampling and sex") +
        theme_ukb()

    ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_time_to_diagnosis.pdf")),
           p_ttd, width = 4.5, height = 4, device = cairo_pdf)
    cat(sprintf("  Saved: %s_prs_time_to_diagnosis.pdf\n", cohort))
}

# ── 6. PRS × proteomic cluster ────────────────────────────────────────────────
cat("\nPRS × proteomic cluster...\n")
if (file.exists(CLUST_FILE)) {
    clust <- fread(CLUST_FILE)[, .(eid, cluster)]
    dt_clust <- merge(dt[get(status_col) %in% case_states],
                      clust, by = "eid", all.x = TRUE)
    dt_clust <- dt_clust[!is.na(cluster)]
    cat(sprintf("  %s cases with cluster + PRS: %d\n", disease_caps, nrow(dt_clust)))

    clust_formula <- as.formula(paste(prs_col, "~ cluster"))
    wt_clust <- wilcox.test(clust_formula, data = dt_clust, exact=FALSE)
    cat(sprintf("  Cluster 1 vs 2 PRS: median %.3f vs %.3f  p=%.3g\n",
                median(dt_clust[cluster==1, get(prs_col)]),
                median(dt_clust[cluster==2, get(prs_col)]),
                wt_clust$p.value))

    p_clust <- ggplot(dt_clust, aes(x = factor(cluster), y = get(prs_col),
                                     fill = factor(cluster))) +
        geom_violin(alpha = 0.6, linewidth = 0.3) +
        geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5, linewidth=0.35) +
        scale_fill_manual(values = c("1"="#CC0066","2"="#56B4E9"), guide="none") +
        annotate("text", x=1.5, y=Inf, label=sprintf("p=%.3g", wt_clust$p.value),
                 vjust=1.5, size=3, colour="grey30") +
        labs(x = "Proteomic cluster", y = glue::glue("{prs_label} (z-score)"),
             title = glue::glue("PRS by proteomic cluster ({disease_caps} cases)")) +
        theme_ukb()

    ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_cluster.pdf")), p_clust,
           width = 3.5, height = 4, device = cairo_pdf)
    cat(sprintf("  Saved: %s_prs_cluster.pdf\n", cohort))
}

# ── 7. PRS × HLA risk allele ─────────────────────────────────────────────────
hla_allele   <- cfg$hla_allele
hla_carrier  <- cfg$hla_carrier_col
cat(sprintf("\nPRS × HLA-%s interaction...\n", hla_allele))
if (file.exists(HLA_FILE)) {
    hla <- fread(HLA_FILE, select = c("eid", hla_carrier),
                 showProgress = FALSE)
    dt_hla <- merge(dt, hla, by = "eid", all.x = TRUE)
    dt_hla <- dt_hla[!is.na(get(hla_carrier))]

    ms_hla <- dt_hla[get(status_col) %in% case_states]
    cat(sprintf("  %s cases with HLA + PRS: %d\n", disease_caps, nrow(ms_hla)))

    hla_formula <- as.formula(paste(prs_col, "~", hla_carrier))
    wt_hla <- wilcox.test(hla_formula, data=ms_hla, exact=FALSE)
    cat(sprintf("  HLA carrier PRS: %.3f  non-carrier: %.3f  p=%.3g\n",
                median(ms_hla[get(hla_carrier)==1, get(prs_col)]),
                median(ms_hla[get(hla_carrier)==0, get(prs_col)]),
                wt_hla$p.value))

    hc_hla <- dt_hla[get(status_col) == sv$control]
    wt_hc  <- wilcox.test(hla_formula, data=hc_hla, exact=FALSE)
    cat(sprintf("  HC: HLA carrier PRS: %.3f  non-carrier: %.3f  p=%.3g\n",
                median(hc_hla[get(hla_carrier)==1, get(prs_col)]),
                median(hc_hla[get(hla_carrier)==0, get(prs_col)]),
                wt_hc$p.value))

    # Compact: PRS by HLA within cases only
    ms_hla[, hla_label := fifelse(get(hla_carrier)==1, "Carrier", "Non-carrier")]
    p_hla <- ggplot(ms_hla, aes(x=hla_label, y=get(prs_col), fill=hla_label)) +
        geom_violin(alpha=0.6, linewidth=0.3) +
        geom_boxplot(width=0.12, fill="white", outlier.size=0.5, linewidth=0.35) +
        scale_fill_manual(values=c("Carrier"="#CC0066","Non-carrier"="#56B4E9"),
                          guide="none") +
        annotate("text", x=1.5, y=Inf,
                 label=sprintf("p=%.3g", wt_hla$p.value),
                 vjust=1.5, size=3, colour="grey30") +
        labs(x=glue::glue("HLA-{hla_allele}"),
             y=glue::glue("{prs_label} (z-score)"),
             title=glue::glue("PRS by HLA-{hla_allele} status ({disease_caps} cases)"),
             subtitle=sprintf("n=%d carriers, %d non-carriers",
                              sum(ms_hla[[hla_carrier]]),
                              sum(!ms_hla[[hla_carrier]]))) +
        theme_ukb()

    ggsave(file.path(OUT_DIR, glue::glue("{cohort}_prs_hla.pdf")), p_hla,
           width=3.5, height=4, device=cairo_pdf)
    cat(sprintf("  Saved: %s_prs_hla.pdf\n", cohort))
}

# ── 8. Save results ───────────────────────────────────────────────────────────
fwrite(corr_dt[, .(protein=toupper(protein), rho, pval, fdr, n)],
       file.path(OUT_DIR, glue::glue("{cohort}_prs_results.csv")))

cat("\n01_ms_prs.R complete.\n")
