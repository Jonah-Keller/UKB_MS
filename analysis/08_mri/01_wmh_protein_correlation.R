#!/usr/bin/env Rscript
# 01_wmh_protein_correlation.R
# Brain MRI: white matter hyperintensity (WMH) correlation with disease DEPs
#
# UKB field p25781_i2 = WMH volume (mm^3) from T2-FLAIR at instance 2
#
# Analyses:
#   1. WMH vs NEFL (Abdelhak key finding: NEFL elevated with WMH burden)
#   2. WMH vs all combined DEPs — Spearman correlation
#   3. WMH vs GFAP, NEFL, MOG (CNS proteins)
#   4. Top correlated proteins
#
# Output: results/mri/
#   wmh_nefl_scatter.pdf
#   wmh_cns_panel.pdf
#   wmh_dep_correlations.csv
#   wmh_top_correlations.pdf

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

QC_FILE   <- here::here("data", "ukb", "olink", "processed",
                        glue::glue("{cfg$cohort_short}_olink_qc.csv"))
WMH_FILE  <- file.path(dirname(here::here()),
                       "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                       "data", "ukb", "misc", "wmh_htn_dm.csv")
DIFF_FILE <- here::here("results", "differential",
                        glue::glue("{cfg$cohort_short}_combined_vs_hc.csv"))
OUT_DIR   <- here::here("results", "mri")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR <- 0.05

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc  <- fread(QC_FILE, showProgress = FALSE)
wmh <- fread(WMH_FILE, showProgress = FALSE)[, .(eid, wmh_vol = p25781_i2)]
dep_list <- fread(DIFF_FILE)[adj.P.Val < FDR_THR, tolower(protein)]
cat(sprintf("  QC file: %d participants, %d DEPs to test\n", nrow(qc), length(dep_list)))

# Merge WMH into QC data
dt <- merge(qc, wmh, by = "eid", all.x = TRUE)

# Restrict to cases with WMH data
case_wmh <- dt[get(STATUS_COL) %in% CASE_VALS & !is.na(wmh_vol)]
case_wmh[, log_wmh := log1p(wmh_vol)]  # log-transform (heavily right-skewed)
case_wmh[, status_f := factor(get(STATUS_COL), levels = CASE_VALS)]

cat(sprintf("  %s cases with WMH: %d (pre=%d, post=%d)\n",
            cfg$disease_short_caps,
            nrow(case_wmh),
            sum(case_wmh[[STATUS_COL]] == PRE),
            sum(case_wmh[[STATUS_COL]] == POST)))

# ── 2. NEFL scatter ────────────────────────────────────────────────────────────
cat("Plotting NEFL vs WMH...\n")
nefl_cor <- cor.test(case_wmh$nefl, case_wmh$log_wmh, method = "spearman", exact = FALSE)

status_colours <- setNames(c("#56B4E9", "#CC0066"), CASE_VALS)
status_labels  <- setNames(c("Pre-onset", "Post-onset"), CASE_VALS)

p_nefl <- ggplot(case_wmh, aes(x = log_wmh, y = nefl, colour = status_f)) +
    geom_point(size = 1.8, alpha = 0.75) +
    geom_smooth(method = "lm", aes(group = 1), colour = UKB_NAVY,
                linewidth = 0.8, se = TRUE, fill = "#56B4E9", alpha = 0.15) +
    scale_colour_manual(values = status_colours,
                        labels = status_labels,
                        name = NULL) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
             label = sprintf("rho = %.2f, p = %.3g",
                             nefl_cor$estimate, nefl_cor$p.value),
             size = 3.2, colour = "grey30") +
    labs(x = "Log WMH volume (log mm\u00b3)",
         y = "NfL NPX (log\u2082)",
         title = glue::glue("NEFL vs white matter hyperintensity burden ({cfg$disease_short_caps} cases)")) +
    theme_ukb()

ggsave(file.path(OUT_DIR, "wmh_nefl_scatter.pdf"),
       p_nefl, width = 5, height = 4.5, device = cairo_pdf)
cat(sprintf("  NEFL rho=%.3f, p=%.3g — Saved: wmh_nefl_scatter.pdf\n",
            nefl_cor$estimate, nefl_cor$p.value))

# ── 3. CNS panel: NEFL / GFAP / MOG ──────────────────────────────────────────
cat("Plotting CNS protein panel...\n")
cns_proteins <- c("nefl","gfap","mog")
cns_labels   <- c("nefl"="NfL","gfap"="GFAP","mog"="MOG")

make_cns_scatter <- function(prot) {
    col <- tolower(prot)
    if (!col %in% names(case_wmh)) return(NULL)
    ct <- cor.test(case_wmh[[col]], case_wmh$log_wmh, method="spearman", exact=FALSE)
    ggplot(case_wmh, aes(x=log_wmh, y=.data[[col]], colour=status_f)) +
        geom_point(size=1.5, alpha=0.7) +
        geom_smooth(method="lm", aes(group=1), colour=UKB_NAVY,
                    linewidth=0.7, se=TRUE, fill="#56B4E9", alpha=0.15) +
        scale_colour_manual(values=status_colours,
                            labels=status_labels,
                            name=NULL) +
        annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.3,
                 label=sprintf("rho=%.2f, p=%.3g", ct$estimate, ct$p.value),
                 size=2.8, colour="grey30") +
        labs(x="Log WMH volume", y=sprintf("%s NPX", cns_labels[prot]),
             title=cns_labels[prot]) +
        theme_ukb() +
        theme(legend.position="none")
}

cns_plots <- Filter(Negate(is.null), lapply(cns_proteins, make_cns_scatter))
if (length(cns_plots) >= 2) {
    p_cns <- wrap_plots(cns_plots, nrow=1) +
        plot_annotation(title=glue::glue("CNS proteins vs WMH burden ({cfg$disease_short_caps} cases, n={nrow(case_wmh)})"),
                        theme=theme(plot.title=element_text(size=10, face="plain")))
    ggsave(file.path(OUT_DIR, "wmh_cns_panel.pdf"),
           p_cns, width=4*length(cns_plots), height=4, device=cairo_pdf)
    cat("  Saved: wmh_cns_panel.pdf\n")
}

# ── 4. All DEPs: Spearman correlation with WMH ───────────────────────────────
cat("Computing WMH correlations for all DEPs...\n")
avail_prots <- intersect(dep_list, names(case_wmh))
cat(sprintf("  Testing %d / %d DEPs\n", length(avail_prots), length(dep_list)))

cor_results <- lapply(avail_prots, function(prot) {
    vals <- case_wmh[[prot]]
    if (sum(!is.na(vals)) < 20) return(NULL)
    ct <- cor.test(vals, case_wmh$log_wmh, method="spearman", exact=FALSE)
    data.table(protein=prot, rho=ct$estimate, p_value=ct$p.value, n=sum(!is.na(vals)))
})
cor_dt <- rbindlist(Filter(Negate(is.null), cor_results))
cor_dt[, fdr := p.adjust(p_value, method="BH")]
cor_dt[, sig := fdr < FDR_THR]
cor_dt <- cor_dt[order(-abs(rho))]

n_sig <- sum(cor_dt$sig)
cat(sprintf("  WMH-correlated DEPs at FDR<0.05: %d\n", n_sig))
print(cor_dt[1:min(15,.N), .(protein, rho=round(rho,3), fdr=round(fdr,4), sig)])

fwrite(cor_dt, file.path(OUT_DIR, "wmh_dep_correlations.csv"))
cat("  Saved: wmh_dep_correlations.csv\n")

# ── 5. Top correlations lollipop ─────────────────────────────────────────────
top_n <- min(25, nrow(cor_dt))
plot_dt <- cor_dt[1:top_n]
plot_dt[, protein_f := factor(protein, levels=rev(protein))]
plot_dt[, dir_col := fifelse(rho > 0, "Positive", "Negative")]

p_top <- ggplot(plot_dt, aes(x=rho, y=protein_f, colour=dir_col, alpha=sig)) +
    geom_segment(aes(x=0, xend=rho, y=protein_f, yend=protein_f),
                 linewidth=0.6) +
    geom_point(size=2.5) +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey50") +
    scale_colour_manual(values=c("Positive"="#CC0066","Negative"="#56B4E9"), name=NULL) +
    scale_alpha_manual(values=c("TRUE"=1.0,"FALSE"=0.4), name="FDR<0.05",
                       labels=c("TRUE"="Yes","FALSE"="No")) +
    labs(x="Spearman rho (protein ~ log WMH)",
         y=NULL,
         title=sprintf("Top %d DEP correlations with WMH burden (n=%d %s cases)",
                       top_n, nrow(case_wmh), cfg$disease_short_caps)) +
    theme_ukb() +
    theme(axis.text.y=element_text(size=7))

ggsave(file.path(OUT_DIR, "wmh_top_correlations.pdf"),
       p_top, width=5.5, height=6, device=cairo_pdf)
cat("  Saved: wmh_top_correlations.pdf\n")

cat("\n01_wmh_protein_correlation.R complete.\n")
