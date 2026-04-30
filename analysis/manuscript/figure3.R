#!/usr/bin/env Rscript
# figure3.R — Figure 3: HLA-stratified molecular heterogeneity in MS
# Panels: a-h
# Output: results/figures/3/panel_{a,b,c,d,e,f,g,h}.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(splines)
    library(glue)
})
HAS_GGSIGNIF <- requireNamespace("ggsignif", quietly = TRUE)
if (HAS_GGSIGNIF) suppressPackageStartupMessages(library(ggsignif))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
    getwd()
}
PROJ_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "celltype_overrep_plot.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "top_n_selector.R"))

cfg          <- load_disease_config()
COHORT       <- cfg$cohort_short
DISEASE_CAPS <- cfg$disease_short_caps
STATUS_COL   <- cfg$cohort_status_col
SV           <- cfg$status_values
HLA_CARRIER  <- cfg$hla_carrier_col
HLA_DOSAGE   <- cfg$hla_dosage_col
PRS_COL      <- cfg$prs_combined_col
HLA_ALLELE   <- cfg$hla_allele
HLA_LBL      <- paste0("HLA-", HLA_ALLELE)
PRS_LABEL    <- cfg$prs_label
PRE_LBL      <- as.character(glue("pre-{DISEASE_CAPS}"))
POST_LBL     <- as.character(glue("post-{DISEASE_CAPS}"))

FIG_DIR  <- file.path(PROJ_DIR, "results", "figures", "3")
SUPP_DIR <- file.path(PROJ_DIR, "results", "figures", "3_supp")
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
save_panel <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(FIG_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: panel_", name, ".pdf\n", sep = "")
}
save_supp <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(SUPP_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: 3_supp/panel_", name, ".pdf\n", sep = "")
}

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------
DATA_DIR  <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed")
GEN_DIR   <- file.path(PROJ_DIR, "data", "ukb", "genetics")
ENDO_DIR  <- file.path(PROJ_DIR, "results", "endophenotype")
PRS_DIR   <- file.path(PROJ_DIR, "results", "prs")
DIFF_DIR  <- file.path(PROJ_DIR, "results", "differential")

QC_FILE       <- file.path(DATA_DIR,  glue("{COHORT}_olink_qc.csv"))
HLA_FILE      <- file.path(GEN_DIR,   "hla_imputed.csv")
PRS_FILE      <- file.path(GEN_DIR,   glue("{COHORT}_prs_scores.csv"))
HLA_DEP_FILE  <- file.path(ENDO_DIR,  glue("{COHORT}_hla_dep_results.csv"))
HLA_INT_FILE      <- file.path(ENDO_DIR, glue("{COHORT}_hla_interaction_all.csv"))
HLA_INT_PRE_FILE  <- file.path(ENDO_DIR, glue("{COHORT}_hla_interaction_pre.csv"))
HLA_INT_POST_FILE <- file.path(ENDO_DIR, glue("{COHORT}_hla_interaction_post.csv"))
HLA_STR_FILE  <- file.path(ENDO_DIR,  glue("{COHORT}_hla_stratified_logfc.csv"))
COMOR_FILE    <- file.path(ENDO_DIR,  glue("{COHORT}_hla_prems_comorbidity.csv"))
PRS_RES_FILE  <- file.path(PRS_DIR,   glue("{COHORT}_prs_results.csv"))
TRAJ_FILE     <- file.path(DIFF_DIR,  "cns_trajectories.csv")
SEX_STRAT_FILE <- file.path(ENDO_DIR, glue("{COHORT}_sex_stratified_genetic.csv"))
PRS_INT_FILE   <- file.path(ENDO_DIR, glue("{COHORT}_prs_interaction_genome.csv"))

# ---------------------------------------------------------------------------
# Load shared data
# ---------------------------------------------------------------------------
cat("Loading data...\n")

qc      <- fread(QC_FILE)
if (STATUS_COL != "ms_status" && STATUS_COL %in% names(qc))
    setnames(qc, STATUS_COL, "ms_status")
hla     <- fread(HLA_FILE)
if (HLA_CARRIER != "drb1_1501_carrier" && HLA_CARRIER %in% names(hla))
    setnames(hla, HLA_CARRIER, "drb1_1501_carrier")
if (HLA_DOSAGE != "drb1_1501_dosage" && HLA_DOSAGE %in% names(hla))
    setnames(hla, HLA_DOSAGE, "drb1_1501_dosage")
hla_dep <- fread(HLA_DEP_FILE)
hla_dep[, fdr := p.adjust(P.Value, method = "BH")]

cat("  QC rows:", nrow(qc), "\n")
cat("  HLA rows:", nrow(hla), "\n")
cat("  HLA DEP rows:", nrow(hla_dep), "\n")

# ---------------------------------------------------------------------------
# Programmatic marker selection
# ---------------------------------------------------------------------------
# HLA_MARKER and PRS_MARKER are the data-driven top-1 protein from each
# stratification analysis.  All downstream panels reference them by the
# sentinel column names `lilrb4` (HLA marker) and `il2ra` (PRS marker), so
# the body of this script reads cleanly across cohorts; the actual protein
# is whatever ranked first in the corresponding stratification CSV for the
# active cohort.  PRS_MARKER is resolved later (after prs_res is loaded).
HLA_MARKER <- {
    sel <- top_n_one_direction(hla_dep, n = 1L,
                                min_sig_for_fdr = cfg$min_sig_for_fdr %||% 5L,
                                case_fn = identity)
    if (length(sel) == 0L) "lilrb4" else sel[1]
}
HLA_MARKER_DISP <- toupper(HLA_MARKER)
cat("  HLA marker (data-driven):", HLA_MARKER_DISP,
    "(ordered by", attr(HLA_MARKER, "ordered_by") %||% "n/a", ")\n")
# Rename in qc so downstream can use the sentinel `lilrb4` column name.
if (HLA_MARKER != "lilrb4" && HLA_MARKER %in% names(qc)) {
    if ("lilrb4" %in% names(qc)) qc[, lilrb4 := NULL]
    setnames(qc, HLA_MARKER, "lilrb4")
}

sex_strat <- if (file.exists(SEX_STRAT_FILE)) fread(SEX_STRAT_FILE) else NULL
prs_int   <- if (file.exists(PRS_INT_FILE))  fread(PRS_INT_FILE)  else NULL
cat("  Sex-stratified rows:", if (!is.null(sex_strat)) nrow(sex_strat) else 0, "\n")
cat("  PRS×MS interaction rows:", if (!is.null(prs_int)) nrow(prs_int) else 0, "\n")

# Metadata columns present in qc
META_COLS <- c("eid", "ms_status", "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "qc_outlier")

# Merge QC + HLA once (shared across panels c, f, k)
qc_hla <- merge(qc, hla[, .(eid, drb1_1501_carrier, drb1_1501_dosage)],
                by = "eid", all.x = TRUE)

# ---------------------------------------------------------------------------
# Panel a — Full HLA allele MS enrichment landscape (80 common alleles)
#   Risk alleles (enriched in MS): warm pink; protective (depleted): blue
#   Gene region inferred from confirmed index anchors:
#     idx 3 = HLA-A*02:01, idx 214 = DQB1*06:02, idx 281 = DRB1*15:01
# ---------------------------------------------------------------------------
cat("Building panel a...\n")

ENRICH_FILE <- file.path(GEN_DIR, "hla_allele_enrichment.csv")
PRUNED_FILE <- file.path(GEN_DIR, "hla_allele_enrichment_pruned.csv")
if (!file.exists(ENRICH_FILE)) {
    cat("  Panel a skipped: hla_allele_enrichment.csv not found\n")
} else {
    enr_src <- if (file.exists(PRUNED_FILE)) PRUNED_FILE else ENRICH_FILE
    enr <- fread(enr_src)
    cat(sprintf("  Panel a: using %s\n", basename(enr_src)))

    # Gene region from confirmed anchor indices
    enr[, gene := fcase(
        idx <= 17,  "HLA-A",
        idx <= 51,  "HLA-B",
        idx <= 65,  "HLA-C",
        idx <= 100, "HLA-DP/DQ\u03b1",
        idx <= 230, "HLA-DQB1",
        idx  > 230, "HLA-DRB1"
    )]

    RISK_LBL <- as.character(glue("Risk (enriched in {DISEASE_CAPS})"))
    PROT_LBL <- as.character(glue("Protective (depleted in {DISEASE_CAPS})"))
    enr[, sig := fdr < 0.05]
    enr[, direction := fcase(
        sig & enrichment > 0, RISK_LBL,
        sig & enrichment < 0, PROT_LBL,
        !sig,                  "NS"
    )]
    enr[, direction := factor(direction,
        levels = c(RISK_LBL, PROT_LBL, "NS"))]

    # Drop LD-redundant alleles (keep_ld == FALSE); NAs are non-sig alleles → keep
    if ("keep_ld" %in% names(enr)) {
        enr <- enr[is.na(keep_ld) | keep_ld == TRUE]
    }

    # Frequency filter + always retain named canonical alleles
    INTEREST_IDX <- c(281L, 3L, 219L)
    enr_plot <- enr[allele_freq >= 0.10 | idx %in% INTEREST_IDX]
    n_plot   <- nrow(enr_plot)
    n_sig    <- sum(enr_plot$fdr < 0.05, na.rm = TRUE)
    cat(sprintf("  %d alleles (freq \u2265 10%%), %d sig after LD pruning\n", n_plot, n_sig))

    enr_plot[, neg_log10_fdr := -log10(pmax(fdr, 1e-40))]

    # Only label the two canonical alleles
    enr_plot[, plot_label := fcase(
        idx == 281L & fdr < 0.05, HLA_ALLELE,
        idx == 3L   & fdr < 0.05, "HLA-A*02:01",
        default = ""
    )]

    DOT_COLS <- setNames(c("#CC0066", COL_POST_DOWN, "grey75"),
                         c(RISK_LBL, PROT_LBL, "NS"))

    fdr_line <- -log10(0.05)

    pA <- ggplot(enr_plot, aes(x = enrichment, y = neg_log10_fdr,
                               colour = direction)) +
        geom_hline(yintercept = fdr_line, linetype = "dashed",
                   linewidth = 0.35, colour = "grey50") +
        geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey60") +
        geom_point(size = 2.2, alpha = 0.85) +
        geom_text_repel(
            data               = enr_plot[plot_label != ""],
            aes(label          = plot_label),
            size               = 2.5,
            fontface           = "bold",
            max.overlaps       = Inf,
            segment.size       = 0.25,
            segment.color      = "grey50",
            box.padding        = 0.5,
            point.padding      = 0.35,
            show.legend        = FALSE,
            min.segment.length = 0.2
        ) +
        scale_colour_manual(values = DOT_COLS, name = NULL) +
        scale_x_continuous(
            labels = function(x) paste0(ifelse(x > 0, "+", ""), round(x), "%"),
            expand = expansion(mult = c(0.08, 0.08))
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        labs(
            title   = glue("a  HLA allele specificity: {DISEASE_CAPS} risk and protection"),
            x       = glue("{DISEASE_CAPS} vs HC carrier rate difference"),
            y       = expression(-log[10]~FDR),
        ) +
        theme_ukb(base_size = 9) +
        theme(
            legend.position  = "bottom",
            legend.key.size  = unit(0.5, "lines")
        )

    save_panel(pA, "a_hla_landscape", w = 5.0, h = 4.5)
}

# ---------------------------------------------------------------------------
# Panel b — Within-MS HLA+ vs HLA− volcano (LILRB4 anchor)
#   Inset: sex-comparison scatter (female vs male within-MS HLA logFC)
#   Nominal p<0.05 shown as coloured; FDR<0.05 as large dot + always labelled
# ---------------------------------------------------------------------------
cat("Building panel b (within-MS HLA+ vs HLA−)...\n")

vol_b <- copy(hla_dep)

vol_b[, tier := fcase(
    fdr   < 0.05 & logFC > 0, "Up FDR",
    P.Value < 0.05 & logFC > 0, "Up nominal",
    fdr   < 0.05 & logFC < 0, "Down FDR",
    P.Value < 0.05 & logFC < 0, "Down nominal",
    default = "NS"
)]
vol_b[, tier := factor(tier,
    levels = c("Up FDR", "Up nominal", "Down FDR", "Down nominal", "NS"))]
vol_b[, neg_log10_p := -log10(P.Value)]

VOL_B_COLS <- c(
    "Up FDR"      = "#E6A817",
    "Up nominal"  = "#F5C97A",
    "Down FDR"    = "#56B4E9",
    "Down nominal"= "#A8D8F0",
    "NS"          = "grey74"
)
VOL_B_SZ <- c("Up FDR"=2.8, "Up nominal"=1.1,
               "Down FDR"=2.8, "Down nominal"=1.1, "NS"=0.45)

# FDR hits always labelled; top nominal by |logFC|
fdr_labs_b <- vol_b[fdr < 0.05][, label := toupper(protein)]
nom_labs_b <- vol_b[P.Value < 0.05 & fdr >= 0.05][order(-abs(logFC))][
    seq_len(min(10, .N))][, label := toupper(protein)]

# FDR threshold line
fdr_thresh_b <- if (any(vol_b$fdr < 0.05)) {
    -log10(max(vol_b[fdr < 0.05, P.Value], na.rm = TRUE))
} else NA_real_

pB_main <- ggplot(vol_b, aes(x = logFC, y = neg_log10_p)) +
    geom_point(data = vol_b[tier == "NS"],
               colour = VOL_B_COLS["NS"], size = VOL_B_SZ["NS"], alpha = 0.35) +
    geom_point(data = vol_b[tier %in% c("Up nominal","Down nominal")],
               aes(colour = tier), size = VOL_B_SZ["Up nominal"], alpha = 0.75) +
    geom_point(data = vol_b[tier %in% c("Up FDR","Down FDR")],
               aes(colour = tier), size = VOL_B_SZ["Up FDR"], alpha = 1.0) +
    {if (!is.na(fdr_thresh_b))
        geom_hline(yintercept = fdr_thresh_b, linetype = "dashed",
                   linewidth = 0.3, colour = "grey45")} +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_text_repel(
        data           = nom_labs_b,
        aes(label      = label),
        size           = 1.95, colour = "grey25",
        max.overlaps   = 20, segment.size = 0.2,
        segment.colour = "grey55", box.padding = 0.28,
        show.legend    = FALSE
    ) +
    geom_text_repel(
        data           = fdr_labs_b,
        aes(label      = label),
        size           = 2.5, fontface = "bold", colour = "grey10",
        max.overlaps   = Inf, segment.size = 0.3,
        segment.colour = "grey30", box.padding = 0.45,
        show.legend    = FALSE
    ) +
    scale_colour_manual(values = VOL_B_COLS, name = NULL) +
    labs(
        title = glue("b  HLA+ vs HLA\u2212 within {DISEASE_CAPS} ({HLA_ALLELE})"),
        x     = expression(log[2]~"FC (HLA+ vs HLA\u2212)"),
        y     = expression(-log[10]~italic(P))
    ) +
    theme_ukb(base_size = 9) +
    theme(legend.position = "none")

# --- Inset: sex-comparison scatter (female vs male within-MS HLA logFC) -----
if (!is.null(sex_strat) && "hla_within_ms" %in% sex_strat$analysis) {
    hw <- sex_strat[analysis == "hla_within_ms",
                    .(protein, logFC, P.Value, fdr, sex_group)]
    hw_wide <- dcast(hw, protein ~ sex_group,
                     value.var = c("logFC", "P.Value", "fdr"))

    hw_wide[, sig_f := get("P.Value_Female") < 0.05]
    hw_wide[, sig_m := get("P.Value_Male")   < 0.05]
    hw_wide[, sig_cat := fcase(
        sig_f & sig_m,  "Both",
        sig_f & !sig_m, "Female only",
        !sig_f & sig_m, "Male only",
        default         = "NS"
    )]
    hw_wide[, sig_cat := factor(sig_cat,
        levels = c("Both", "Female only", "Male only", "NS"))]
    hw_wide[, lbl := fifelse(protein == "lilrb4", HLA_MARKER_DISP, NA_character_)]

    INS_COLS <- c("Both"="#E6A817","Female only"="#CC0066","Male only"="#2B4C7E","NS"="grey82")

    pB_ins <- ggplot(hw_wide, aes(x = logFC_Female, y = logFC_Male,
                                   colour = sig_cat)) +
        geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey70") +
        geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    linewidth = 0.28, colour = "grey55") +
        geom_point(data = hw_wide[sig_cat == "NS"],
                   size = 0.4, alpha = 0.30) +
        geom_point(data = hw_wide[sig_cat != "NS"],
                   size = 1.5, alpha = 0.90) +
        geom_text_repel(
            aes(label = lbl), na.rm = TRUE,
            size = 2.1, fontface = "bold",
            max.overlaps = Inf, segment.size = 0.25,
            nudge_x = 0.18, nudge_y = 0.08,
            min.segment.length = 0, show.legend = FALSE
        ) +
        scale_colour_manual(values = INS_COLS, name = NULL) +
        labs(title = "Sex comparison", x = "Female logFC", y = "Male logFC") +
        theme_ukb(base_size = 7) +
        theme(
            legend.position  = "none",
            plot.title       = element_text(size = 6.5, face = "plain",
                                            margin = margin(b = 1)),
            axis.title       = element_text(size = 6),
            axis.text        = element_text(size = 5.5),
            plot.background  = element_rect(fill = "white", colour = "grey75",
                                            linewidth = 0.35),
            plot.margin      = margin(3, 3, 3, 3)
        )

    # Save sex comparison as its own panel (larger, standalone)
    pB_ins_full <- pB_ins +
        labs(
            title    = glue("c  Within-{DISEASE_CAPS} HLA effect: female vs male concordance"),
            subtitle = "Nominal p<0.05 in female only / male only / both | diagonal = perfect concordance",
            x        = "Female logFC (HLA+ vs HLA\u2212)",
            y        = "Male logFC (HLA+ vs HLA\u2212)"
        ) +
        theme_ukb(base_size = 9) +
        theme(
            plot.title    = element_text(size = 9, face = "plain"),
            plot.subtitle = element_text(size = 6.5, colour = "grey45"),
            axis.title    = element_text(size = 8),
            axis.text     = element_text(size = 7.5),
            legend.position = "bottom",
            plot.background = element_rect(fill = "white", colour = NA),
            plot.margin   = margin(8, 8, 8, 8)
        ) +
        guides(colour = guide_legend(
            title    = NULL,
            override.aes = list(size = 2.5)
        ))
    save_panel(pB_ins_full, "c_sex_comparison", 3.5, 3.5)
}

save_panel(pB_main, "b_hla_within_ms", 3.8, 4.2)

# ---------------------------------------------------------------------------
# Panel c — HLA-stratified violin plots (top 6 HLA DEPs, one panel each)
# ---------------------------------------------------------------------------
cat("Building panel c...\n")

top6_c <- hla_dep[P.Value < 0.05][order(-abs(logFC))][seq_len(min(6, .N)), tolower(protein)]

vio_base <- qc_hla[
    qc_outlier == FALSE &
        ms_status %in% c(SV$pre_onset, SV$post_onset, SV$control) &
        !is.na(drb1_1501_carrier)
]
HC_NEG_LBL <- "HC HLA\u2212"
HC_POS_LBL <- "HC HLA+"
MS_NEG_LBL <- as.character(glue("{DISEASE_CAPS} HLA\u2212"))
MS_POS_LBL <- as.character(glue("{DISEASE_CAPS} HLA+"))
vio_base[, group4 := fcase(
    ms_status == SV$control & drb1_1501_carrier == 0, HC_NEG_LBL,
    ms_status == SV$control & drb1_1501_carrier == 1, HC_POS_LBL,
    ms_status != SV$control & drb1_1501_carrier == 0, MS_NEG_LBL,
    ms_status != SV$control & drb1_1501_carrier == 1, MS_POS_LBL
)]
vio_base[, group4 := factor(group4,
    levels = c(HC_NEG_LBL, HC_POS_LBL, MS_NEG_LBL, MS_POS_LBL))]

VIO_FILLS <- setNames(
    c("grey72", "grey45", "#2B4C7E", "#E6A817"),
    c(HC_NEG_LBL, HC_POS_LBL, MS_NEG_LBL, MS_POS_LBL)
)

# Three biologically meaningful comparisons
COMPS_C <- list(
    c(HC_NEG_LBL, MS_NEG_LBL),
    c(HC_POS_LBL, MS_POS_LBL),
    c(MS_NEG_LBL, MS_POS_LBL)
)

avail_prot_c <- intersect(top6_c, names(vio_base))
cat("  Top HLA proteins:", paste(toupper(avail_prot_c), collapse = ", "), "\n")

supp_vio_counter <- 0L  # tracks supplement violin labelling (a, b, c, …)

for (i in seq_along(avail_prot_c)) {
    prot_lc <- avail_prot_c[i]
    prot_uc <- toupper(prot_lc)
    # Only FDR-significant proteins appear in main figure (panel d); rest → 3_supp
    is_main <- isTRUE(hla_dep[tolower(protein) == prot_lc, fdr[1]] < 0.05)
    if (is_main) {
        sub_ltr <- "d"
    } else {
        supp_vio_counter <- supp_vio_counter + 1L
        sub_ltr <- letters[supp_vio_counter]
    }

    d_v <- vio_base[!is.na(get(prot_lc)),
                    .(group4, npx = get(prot_lc))]

    p_v <- ggplot(d_v, aes(x = group4, y = npx, fill = group4)) +
        geom_violin(trim = TRUE, linewidth = 0.3, alpha = 0.85) +
        geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3,
                     fill = "white", colour = "grey30") +
        scale_fill_manual(values = VIO_FILLS, guide = "none") +
        labs(
            title = sprintf("%s  %s", sub_ltr, prot_uc),
            x     = NULL,
            y     = "NPX"
        ) +
        theme_ukb(base_size = 9) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

    if (HAS_GGSIGNIF) {
        sig_comps_v <- Filter(function(comp) {
            d1 <- d_v$npx[d_v$group4 == comp[1]]
            d2 <- d_v$npx[d_v$group4 == comp[2]]
            if (length(d1) < 3 || length(d2) < 3) return(FALSE)
            tryCatch(wilcox.test(d1, d2, exact = FALSE)$p.value < 0.05,
                     error = function(e) FALSE)
        }, COMPS_C)
        if (length(sig_comps_v) > 0) {
            p_v <- p_v + geom_signif(
                comparisons      = sig_comps_v,
                test             = "wilcox.test",
                map_signif_level = TRUE,
                step_increase    = 0.09,
                tip_length       = 0.01,
                textsize         = 2.8,
                colour           = "grey25"
            )
        }
    }

    if (is_main) {
        save_panel(p_v, sprintf("%s_%s", sub_ltr, prot_lc), 2.8, 3.8)
    } else {
        save_supp(p_v, sprintf("%s_%s", sub_ltr, prot_lc), 2.8, 3.8)
    }
}

# ---------------------------------------------------------------------------
# Panel i — MS-effect logFC scatter: HLA− vs HLA+ background
# ---------------------------------------------------------------------------
cat("Building panel e...\n")

strat_d <- fread(HLA_STR_FILE)
cat("  Stratified logFC rows:", nrow(strat_d), "\n")

strat_d[, direction_d := fcase(
    fdr_hlaneg < 0.05 & fdr_hlapos < 0.05, "DEP in both",
    fdr_hlaneg < 0.05 | fdr_hlapos < 0.05, "DEP in one",
    default = "NS"
)]
strat_d[, direction_d := factor(direction_d,
    levels = c("DEP in both", "DEP in one", "NS"))]

int_sig_d  <- strat_d[fdr_hlaneg < 0.05 | fdr_hlapos < 0.05][
    order(-(abs(logFC_hlaneg) + abs(logFC_hlapos)))][seq_len(min(15, .N))]
r_val_d    <- cor(strat_d$logFC_hlaneg, strat_d$logFC_hlapos, use = "complete.obs")
r_lab_d    <- sprintf("r = %.2f", r_val_d)

# Count how many are sig in each single stratum for annotation
n_hlaneg_sig <- sum(strat_d$fdr_hlaneg < 0.05, na.rm = TRUE)
n_hlapos_sig <- sum(strat_d$fdr_hlapos < 0.05, na.rm = TRUE)

x_range_d <- range(strat_d$logFC_hlaneg, na.rm = TRUE)
y_range_d <- range(strat_d$logFC_hlapos, na.rm = TRUE)

pD <- ggplot(strat_d, aes(x = logFC_hlaneg, y = logFC_hlapos,
                           colour = direction_d)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.4, colour = "grey50") +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_text_repel(
        data          = int_sig_d,
        aes(label     = toupper(protein)),
        size          = 2.2,
        fontface      = "bold",
        max.overlaps  = 25,
        segment.size  = 0.25,
        segment.color = "grey50",
        box.padding   = 0.3,
        show.legend   = FALSE
    ) +
    annotate("text",
             x = x_range_d[2] * 0.95, y = y_range_d[1] * 0.95,
             label = paste0(n_hlaneg_sig, " DEPs"),
             size = 2.2, hjust = 1, vjust = 0, colour = "grey30") +
    annotate("text",
             x = x_range_d[1] * 0.95, y = y_range_d[2] * 0.95,
             label = paste0(n_hlapos_sig, " DEPs"),
             size = 2.2, hjust = 0, vjust = 1, colour = "grey30") +
    annotate("text",
             x = x_range_d[1] * 0.9, y = y_range_d[1] * 0.9,
             label = r_lab_d,
             size = 2.5, hjust = 0, colour = "grey30") +
    scale_colour_manual(
        values = c("DEP in both" = "#E6A817", "DEP in one" = "#56B4E9", "NS" = "grey65"),
        name   = NULL
    ) +
    labs(
        title = glue("e  {DISEASE_CAPS} effect: HLA\u2212 vs HLA+ background"),
        x     = expression(log[2]~"FC (HLA\u2212 background)"),
        y     = expression(log[2]~"FC (HLA+ background)")
    ) +
    theme_ukb(base_size = 9)

save_panel(pD, "e_hla_logfc_scatter", 3.5, 3.5)

# ---------------------------------------------------------------------------
# Panel j — Pathway enrichment on HLA-specific DEPs
# ---------------------------------------------------------------------------
cat("Building panel f...\n")

tryCatch({
    suppressPackageStartupMessages({
        library(clusterProfiler)
        library(org.Hs.eg.db)
    })

    calc_fe <- function(gr, br) {
        p  <- strsplit(gr, "/"); num <- as.numeric(sapply(p, `[[`, 1)); den <- as.numeric(sapply(p, `[[`, 2))
        bp <- strsplit(br, "/"); bn  <- as.numeric(sapply(bp,`[[`, 1)); bd  <- as.numeric(sapply(bp,`[[`, 2))
        (num / den) / (bn / bd)
    }

    hla_up   <- hla_dep[P.Value < 0.05 & logFC > 0, toupper(protein)]
    hla_down <- hla_dep[P.Value < 0.05 & logFC < 0, toupper(protein)]
    cat(sprintf("  HLA GO: %d up in HLA+,  %d down in HLA+\n", length(hla_up), length(hla_down)))

    run_ego <- function(genes, label) {
        if (length(genes) < 5) return(NULL)
        ego <- tryCatch(
            enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     readable = TRUE),
            error = function(e) NULL
        )
        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
        dt <- as.data.table(as.data.frame(ego))
        dt[, FoldEnrichment := calc_fe(GeneRatio, BgRatio)]
        dt[, direction := label]
        dt[order(p.adjust)][seq_len(min(8, .N))]
    }

    go_rows <- list(
        run_ego(hla_up,   "Up in HLA+"),
        run_ego(hla_down, "Down in HLA+")
    )
    go_f <- rbindlist(Filter(Negate(is.null), go_rows), fill = TRUE)

    if (nrow(go_f) == 0) stop("No significant GO terms for either direction")

    go_f[, direction := factor(direction, levels = c("Up in HLA+", "Down in HLA+"))]
    go_f[, Description := factor(Description,
        levels = unique(go_f[order(direction, FoldEnrichment), Description]))]

    GO_F_COLS <- c("Up in HLA+" = "#E6A817", "Down in HLA+" = "#2B4C7E")

    pE <- ggplot(go_f, aes(x = FoldEnrichment, y = Description,
                            size = Count, colour = direction)) +
        geom_point() +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
        scale_colour_manual(values = GO_F_COLS,
                             guide = guide_legend(override.aes = list(size = 3.5))) +
        scale_size_continuous(name = "Gene count", range = c(2, 6),
                               guide = guide_legend(override.aes = list(colour = "grey40"))) +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        labs(title    = "f  GO:BP enrichment \u2014 HLA-stratified DEPs",
             subtitle = sprintf("Nominal p<0.05 HLA\u00b1 within %s: %d up in HLA+, %d down in HLA+",
                                DISEASE_CAPS, length(hla_up), length(hla_down)),
             x = "Fold enrichment", y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y    = element_text(size = 7),
              strip.text     = element_text(size = 8, face = "bold"),
              legend.position = "right",
              plot.subtitle  = element_text(size = 7, colour = "grey40"))

    save_panel(pE, "f_hla_go", 5.5, 6)

}, error = function(e) {
    cat("  Panel f GO: clusterProfiler/org.Hs.eg.db not available or no enrichment —",
        conditionMessage(e), "\n")
    pE_ph <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "Requires clusterProfiler\nand org.Hs.eg.db",
                 size = 4, hjust = 0.5, vjust = 0.5, colour = "grey40") +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = "f  GO:BP enrichment \u2014 HLA-stratified DEPs") +
        theme_void() + theme(plot.title = element_text(size = 9, face = "plain"))
    save_panel(pE_ph, "f_hla_go", 5.5, 6)
})

# ---------------------------------------------------------------------------
# Panels l, m, n — HLA-stratified full disease-course trajectories
#   Top within-MS HLA DEPs across the full disease course (pre + post).
#   ytd_display = -years_to_diagnosis for all MS cases:
#     pre_onset ytd > 0  → display < 0 (years before Dx)
#     post_onset ytd < 0 → display > 0 (years after Dx)
# ---------------------------------------------------------------------------
cat("Building panels l, m, n (HLA-stratified disease-course trajectories)...\n")

# Programmatic: top-3 HLA-stratified DEPs by adj.P.Val (FDR-aware fallback to
# nominal P). The first protein gets the main-figure trajectory panel; the
# remaining two go to Extended Data Figure 3 supp.
hla_top <- top_n_one_direction(
    hla_dep, n = 3L,
    min_sig_for_fdr = cfg$min_sig_for_fdr %||% 5L,
    case_fn = toupper
)
TRAJ_HLA_PROTS  <- as.character(hla_top)
# First protein → main figure panel "g"; remainder → ED supp panels f, g (or
# whatever letters the layout assigns).  Keep panel-letter mapping defensive
# in case fewer than 3 hits surface.
.traj_letters <- c("g", "f", "g")
TRAJ_HLA_LABELS <- setNames(.traj_letters[seq_along(TRAJ_HLA_PROTS)],
                             TRAJ_HLA_PROTS)
TRAJ_HLA_MAIN   <- setNames(c(TRUE, rep(FALSE, length(TRAJ_HLA_PROTS) - 1L)),
                             TRAJ_HLA_PROTS)
cat("  HLA trajectory panels (data-driven, ordered by ",
    attr(hla_top, "ordered_by"), "): ",
    paste(TRAJ_HLA_PROTS, collapse = ", "), "\n", sep = "")
HLA_TRAJ_COLS   <- c("HLA+" = "#E6A817", "HLA\u2212" = "#2B4C7E")
YTD_LMN_RANGE   <- c(-15, 15)

ms_traj <- qc_hla[
    ms_status %in% c(SV$pre_onset, SV$post_onset) &
    qc_outlier == FALSE &
    !is.na(drb1_1501_carrier) &
    !is.na(years_to_diagnosis)
]
ms_traj[, ytd_display := -years_to_diagnosis]
ms_traj <- ms_traj[ytd_display >= YTD_LMN_RANGE[1] & ytd_display <= YTD_LMN_RANGE[2]]
cat(sprintf("  MS for trajectory: n=%d (pre=%d, post=%d)\n",
    nrow(ms_traj),
    sum(ms_traj$ms_status == SV$pre_onset),
    sum(ms_traj$ms_status == SV$post_onset)))

x_grid_lmn <- seq(YTD_LMN_RANGE[1], YTD_LMN_RANGE[2], length.out = 200)

for (prot in TRAJ_HLA_PROTS) {
    col_nm  <- tolower(prot)
    pan_ltr <- TRAJ_HLA_LABELS[[prot]]
    if (!col_nm %in% names(ms_traj)) {
        cat("  Skipping", prot, "(not in data)\n"); next
    }

    traj_rows <- list()
    for (hla_val in c(0, 1)) {
        hla_label <- if (hla_val == 1) "HLA+" else "HLA\u2212"
        d <- ms_traj[drb1_1501_carrier == hla_val,
                     .(npx = get(col_nm), ytd_display, age_at_sampling, sex)]
        d <- d[!is.na(npx)]
        if (nrow(d) < 15) next
        fit <- tryCatch(
            lm(npx ~ ns(ytd_display, df = 4) + age_at_sampling + sex, data = d),
            error = function(e) NULL
        )
        if (is.null(fit)) next
        nd <- data.table(
            ytd_display     = x_grid_lmn,
            age_at_sampling = median(d$age_at_sampling, na.rm = TRUE),
            sex             = as.numeric(names(sort(table(d$sex), decreasing = TRUE))[1])
        )
        pr <- tryCatch(predict(fit, newdata = nd, se.fit = TRUE), error = function(e) NULL)
        if (is.null(pr)) next
        traj_rows[[hla_label]] <- data.table(
            hla_group   = hla_label,
            ytd_display = x_grid_lmn,
            pred_npx    = pr$fit,
            pred_lci    = pr$fit - 1.96 * pr$se.fit,
            pred_uci    = pr$fit + 1.96 * pr$se.fit
        )
    }
    if (length(traj_rows) == 0) { cat("  Skipping", prot, "(insufficient data)\n"); next }
    traj_dt <- rbindlist(traj_rows)
    traj_dt[, hla_group := factor(hla_group, levels = c("HLA\u2212", "HLA+"))]

    # Global HLA offset t-test at each timebin (pre vs post separately)
    n_hlapos <- ms_traj[drb1_1501_carrier == 1 & !is.na(get(col_nm)), .N]
    n_hlaneg <- ms_traj[drb1_1501_carrier == 0 & !is.na(get(col_nm)), .N]

    raw_p <- ms_traj[!is.na(get(col_nm)),
                     .(ytd_display, npx = get(col_nm),
                       hla_group = fifelse(drb1_1501_carrier == 1, "HLA+", "HLA\u2212"))]
    raw_p[, hla_group := factor(hla_group, levels = c("HLA\u2212", "HLA+"))]

    p_traj <- ggplot() +
        geom_point(data = raw_p, aes(x = ytd_display, y = npx, colour = hla_group),
                   size = 0.4, alpha = 0.12) +
        geom_ribbon(data = traj_dt,
                    aes(x = ytd_display, ymin = pred_lci, ymax = pred_uci, fill = hla_group),
                    alpha = 0.18) +
        geom_line(data = traj_dt,
                  aes(x = ytd_display, y = pred_npx, colour = hla_group),
                  linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
        annotate("text", x = 0.5, y = Inf, label = "Dx",
                 hjust = 0, vjust = 1.5, size = 2.2, colour = "grey40") +
        annotate("text", x = YTD_LMN_RANGE[1] + 0.5, y = -Inf,
                 label = sprintf("HLA+ n=%d  |  HLA\u2212 n=%d", n_hlapos, n_hlaneg),
                 hjust = 0, vjust = -0.4, size = 1.8, colour = "grey45") +
        scale_colour_manual(values = HLA_TRAJ_COLS, name = NULL) +
        scale_fill_manual(values = HLA_TRAJ_COLS, name = NULL, guide = "none") +
        scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
        labs(
            title = sprintf("%s  %s: disease-course trajectory by HLA status (%s)", pan_ltr, prot, HLA_ALLELE),
            x     = glue("Years relative to {DISEASE_CAPS} diagnosis"),
            y     = paste0(prot, " (NPX)")
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom")

    if (isTRUE(TRAJ_HLA_MAIN[[prot]])) {
        save_panel(p_traj, sprintf("%s_%s_hla_trajectory", pan_ltr, tolower(prot)), 4.0, 3.5)
    } else {
        save_supp(p_traj, sprintf("%s_%s_hla_trajectory", pan_ltr, tolower(prot)), 4.0, 3.5)
    }
}

# ---------------------------------------------------------------------------
# Panel h — LILRB4 trajectory sex-stratified by HLA (Female | Male facets)
# ---------------------------------------------------------------------------
cat("Building panel h (LILRB4 sex-stratified trajectory)...\n")

if ("lilrb4" %in% names(ms_traj)) {
    SEX_LABELS_H <- c("0" = "Female", "1" = "Male")

    traj_sh_rows <- list()
    for (sex_val in c(0, 1)) {
        sex_label <- SEX_LABELS_H[as.character(sex_val)]
        for (hla_val in c(0, 1)) {
            hla_label <- if (hla_val == 1) "HLA+" else "HLA\u2212"
            d <- ms_traj[sex == sex_val & drb1_1501_carrier == hla_val,
                         .(npx = lilrb4, ytd_display, age_at_sampling)]
            d <- d[!is.na(npx)]
            if (nrow(d) < 10) next
            fit <- tryCatch(
                lm(npx ~ ns(ytd_display, df = 3) + age_at_sampling, data = d),
                error = function(e) NULL
            )
            if (is.null(fit)) next
            nd  <- data.table(ytd_display     = x_grid_lmn,
                              age_at_sampling = median(d$age_at_sampling, na.rm = TRUE))
            pr  <- tryCatch(predict(fit, newdata = nd, se.fit = TRUE), error = function(e) NULL)
            if (is.null(pr)) next
            key <- paste(sex_label, hla_label)
            traj_sh_rows[[key]] <- data.table(
                sex_group   = sex_label,
                hla_group   = hla_label,
                ytd_display = x_grid_lmn,
                pred_npx    = pr$fit,
                pred_lci    = pr$fit - 1.96 * pr$se.fit,
                pred_uci    = pr$fit + 1.96 * pr$se.fit,
                n_grp       = nrow(d)
            )
        }
    }

    if (length(traj_sh_rows) > 0) {
        traj_sh <- rbindlist(traj_sh_rows)
        traj_sh[, sex_group := factor(sex_group, levels = c("Female", "Male"))]
        traj_sh[, hla_group := factor(hla_group, levels = c("HLA\u2212", "HLA+"))]

        raw_sh <- ms_traj[!is.na(lilrb4),
                          .(ytd_display, npx = lilrb4,
                            sex_group   = SEX_LABELS_H[as.character(sex)],
                            hla_group   = fifelse(drb1_1501_carrier == 1, "HLA+", "HLA\u2212"))]
        raw_sh[, sex_group := factor(sex_group, levels = c("Female", "Male"))]
        raw_sh[, hla_group := factor(hla_group, levels = c("HLA\u2212", "HLA+"))]

        # n labels per facet × HLA group
        n_labs_h <- traj_sh[ytd_display == x_grid_lmn[1],
                             .(label = sprintf("n=%d", n_grp)), by = .(sex_group, hla_group)]

        pHh <- ggplot() +
            geom_point(data = raw_sh,
                       aes(x = ytd_display, y = npx, colour = hla_group),
                       size = 0.4, alpha = 0.12) +
            geom_ribbon(data = traj_sh,
                        aes(x = ytd_display, ymin = pred_lci, ymax = pred_uci, fill = hla_group),
                        alpha = 0.18) +
            geom_line(data = traj_sh,
                      aes(x = ytd_display, y = pred_npx, colour = hla_group),
                      linewidth = 0.85) +
            geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
            annotate("text", x = 0.5, y = Inf, label = "Dx",
                     hjust = 0, vjust = 1.5, size = 2.2, colour = "grey40") +
            facet_wrap(~ sex_group, nrow = 1) +
            scale_colour_manual(values = HLA_TRAJ_COLS, name = NULL) +
            scale_fill_manual(values = HLA_TRAJ_COLS, name = NULL, guide = "none") +
            scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
            labs(
                title = glue("h  {HLA_MARKER_DISP} trajectory by HLA status: female vs male"),
                x     = glue("Years relative to {DISEASE_CAPS} diagnosis"),
                y     = glue("{HLA_MARKER_DISP} (NPX)")
            ) +
            theme_ukb(base_size = 9) +
            theme(legend.position = "bottom")

        save_panel(pHh, glue("h_{tolower(HLA_MARKER_DISP)}_sex_hla_trajectory"), 6.0, 3.5)
    } else {
        cat("  Insufficient data for sex-stratified LILRB4 trajectory\n")
    }
} else {
    cat("  lilrb4 not in ms_traj data\n")
}

# ---------------------------------------------------------------------------
# Panel i — PRS-proteome association
# ---------------------------------------------------------------------------
cat("Building panel i...\n")

prs_res <- fread(PRS_RES_FILE)
cat("  PRS result rows:", nrow(prs_res), "\n")

prs_res[, neg_log10_p := -log10(pval)]
prs_res[, sig_prs := fdr < 0.05]
n_sig_prs <- sum(prs_res$sig_prs, na.rm = TRUE)
annot_prs <- sprintf("%d protein%s FDR<0.05 (HLA-adjusted)",
                     n_sig_prs, if (n_sig_prs == 1) "" else "s")

prs_res[, col_group := fifelse(sig_prs, "FDR<0.05", "NS")]
prs_res[, col_group := factor(col_group, levels = c("FDR<0.05", "NS"))]

PRS_COLS <- c("FDR<0.05" = "#CC0066", "NS" = "grey65")

# Label sig hits
sig_prs_labs <- prs_res[sig_prs == TRUE]
sig_prs_labs[, label := toupper(protein)]

pG <- ggplot(prs_res, aes(x = rho, y = neg_log10_p, colour = col_group)) +
    geom_point(size = 0.6, alpha = 0.7) +
    {
        if (nrow(sig_prs_labs) > 0) {
            geom_text_repel(
                data          = sig_prs_labs,
                aes(label     = label),
                size          = 2.2,
                fontface      = "bold",
                max.overlaps  = 20,
                segment.size  = 0.25,
                segment.color = "grey50",
                box.padding   = 0.3,
                show.legend   = FALSE
            )
        } else {
            list()
        }
    } +
    annotate("text",
             x = max(prs_res$rho, na.rm = TRUE) * 0.95,
             y = max(prs_res$neg_log10_p, na.rm = TRUE) * 0.95,
             label = annot_prs,
             size = 2.2, hjust = 1, colour = "grey30") +
    geom_hline(
        yintercept = -log10(max(prs_res[fdr < 0.05, pval], na.rm = TRUE)),
        linetype = "dashed", linewidth = 0.3, colour = "grey50"
    ) +
    scale_colour_manual(values = PRS_COLS, name = NULL) +
    labs(
        title = "i  PRS-proteome association (HLA-adjusted)",
        x     = "Spearman rho (PRS vs protein)",
        y     = expression(-log[10](p))
    ) +
    theme_ukb(base_size = 9)

save_panel(pG, "i_prs_proteome", 3.5, 3.5)

# ---------------------------------------------------------------------------
# Panel j — IL2RA disease-course trajectory by PRS quartile (Q1 vs Q4)
#   Top hit from panel i (PRS-proteome: rho=0.063, FDR<10^-39).
#   Illustrates that high-PRS MS cases show elevated IL2RA across the course.
# ---------------------------------------------------------------------------
cat("Building panel j (IL2RA trajectory by PRS quartile)...\n")

# Programmatic: top-1 PRS-correlated protein. Renamed to sentinel `il2ra`
# in qc so downstream panels reference the correct cohort-specific protein
# without further edits.
PRS_MARKER <- {
    # prs_res uses lowercase `pval` column; harmonise with the helper expectation
    prs_res_h <- copy(prs_res)
    if ("pval" %in% names(prs_res_h) && !"P.Value" %in% names(prs_res_h))
        prs_res_h[, P.Value := pval]
    if ("fdr"  %in% names(prs_res_h) && !"adj.P.Val" %in% names(prs_res_h))
        prs_res_h[, adj.P.Val := fdr]
    sel <- top_n_one_direction(prs_res_h, n = 1L,
                                min_sig_for_fdr = cfg$min_sig_for_fdr %||% 5L,
                                case_fn = identity)
    if (length(sel) == 0L) "il2ra" else sel[1]
}
PRS_MARKER_DISP <- toupper(PRS_MARKER)
cat("  PRS marker (data-driven):", PRS_MARKER_DISP,
    "(ordered by", attr(PRS_MARKER, "ordered_by") %||% "n/a", ")\n")
if (PRS_MARKER != "il2ra" && PRS_MARKER %in% names(qc)) {
    if ("il2ra" %in% names(qc)) qc[, il2ra := NULL]
    setnames(qc, PRS_MARKER, "il2ra")
}
# qc_hla was built before PRS_MARKER was resolved; mirror the rename there
if (exists("qc_hla") && PRS_MARKER != "il2ra" && PRS_MARKER %in% names(qc_hla)) {
    if ("il2ra" %in% names(qc_hla)) qc_hla[, il2ra := NULL]
    setnames(qc_hla, PRS_MARKER, "il2ra")
}
PRS_TRAJ_PROT <- "il2ra"   # sentinel column name (see PRS_MARKER above)

# Display labels — programmatic counterpart to the previously hardcoded
# c("LILRB4 (HLA marker)", "IL2RA (PRS marker)") factor levels.
HLA_MKR_LABEL <- paste0(HLA_MARKER_DISP, " (HLA marker)")
PRS_MKR_LABEL <- paste0(PRS_MARKER_DISP, " (PRS marker)")
MKR_LEVELS    <- c(HLA_MKR_LABEL, PRS_MKR_LABEL)

if (file.exists(PRS_FILE) && PRS_TRAJ_PROT %in% names(qc)) {
    prs_scores_j <- fread(PRS_FILE, showProgress = FALSE)
    if (PRS_COL != "prs_score" && PRS_COL %in% names(prs_scores_j))
        setnames(prs_scores_j, PRS_COL, "prs_score")
    ms_prs_traj  <- merge(
        qc[qc_outlier == FALSE &
           ms_status %in% c(SV$pre_onset, SV$post_onset) &
           !is.na(years_to_diagnosis)],
        prs_scores_j[, .(eid, prs_score)],
        by = "eid"
    )
    ms_prs_traj[, ytd_display := -years_to_diagnosis]
    ms_prs_traj <- ms_prs_traj[ytd_display >= YTD_LMN_RANGE[1] &
                                ytd_display <= YTD_LMN_RANGE[2] &
                                !is.na(get(PRS_TRAJ_PROT)) & !is.na(prs_score)]

    if (nrow(ms_prs_traj) >= 30) {
        ms_prs_traj[, prs_q := cut(prs_score,
            breaks = quantile(prs_score, c(0, .25, .5, .75, 1), na.rm = TRUE),
            labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
            include.lowest = TRUE)]

        PRS_Q_TRAJ <- c("Q1 (low)" = "#56B4E9", "Q4 (high)" = "#CC0066")

        traj_prs_rows <- list()
        for (q_val in c("Q1 (low)", "Q4 (high)")) {
            d <- ms_prs_traj[prs_q == q_val,
                             .(npx = get(PRS_TRAJ_PROT), ytd_display, age_at_sampling, sex)]
            d <- d[!is.na(npx)]
            if (nrow(d) < 15) next
            fit <- tryCatch(
                lm(npx ~ ns(ytd_display, df = 4) + age_at_sampling + sex, data = d),
                error = function(e) NULL
            )
            if (is.null(fit)) next
            nd <- data.table(
                ytd_display     = x_grid_lmn,
                age_at_sampling = median(d$age_at_sampling, na.rm = TRUE),
                sex             = as.numeric(names(sort(table(d$sex), decreasing = TRUE))[1])
            )
            pr <- tryCatch(predict(fit, newdata = nd, se.fit = TRUE), error = function(e) NULL)
            if (is.null(pr)) next
            traj_prs_rows[[q_val]] <- data.table(
                prs_q = q_val, n_grp = nrow(d),
                ytd_display = x_grid_lmn,
                pred_npx    = pr$fit,
                pred_lci    = pr$fit - 1.96 * pr$se.fit,
                pred_uci    = pr$fit + 1.96 * pr$se.fit
            )
        }

        if (length(traj_prs_rows) > 0) {
            traj_prs <- rbindlist(traj_prs_rows)
            traj_prs[, prs_q := factor(prs_q, levels = c("Q1 (low)", "Q4 (high)"))]

            pJ_new <- ggplot() +
                geom_ribbon(data = traj_prs,
                            aes(x = ytd_display, ymin = pred_lci, ymax = pred_uci, fill = prs_q),
                            alpha = 0.18) +
                geom_line(data = traj_prs,
                          aes(x = ytd_display, y = pred_npx, colour = prs_q),
                          linewidth = 0.85) +
                geom_vline(xintercept = 0, linetype = "dashed",
                           linewidth = 0.4, colour = "grey40") +
                annotate("text", x = 0.5, y = Inf, label = "Dx",
                         hjust = 0, vjust = 1.5, size = 2.2, colour = "grey40") +
                scale_colour_manual(values = PRS_Q_TRAJ, name = "PRS quartile") +
                scale_fill_manual(values = PRS_Q_TRAJ, name = "PRS quartile", guide = "none") +
                scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
                labs(
                    title    = paste0("j  ", toupper(PRS_TRAJ_PROT),
                                      ": disease-course trajectory by PRS quartile"),
                    subtitle = "Top PRS-proteome hit (panel i, FDR<10\u207b\u00b3\u2079) | Q1=lowest, Q4=highest genetic risk | age + sex adjusted",
                    x        = glue("Years relative to {DISEASE_CAPS} diagnosis"),
                    y        = paste0(toupper(PRS_TRAJ_PROT), " (NPX)")
                ) +
                theme_ukb(base_size = 9) +
                theme(legend.position  = "bottom",
                      plot.subtitle    = element_text(size = 6.5, colour = "grey45"))

            save_panel(pJ_new, sprintf("j_%s_prs_trajectory", PRS_TRAJ_PROT), 4.0, 3.5)
        }
    } else {
        cat("  Insufficient data for", toupper(PRS_TRAJ_PROT), PRS_LABEL, "trajectory\n")
    }
} else {
    cat("  PRS file or", toupper(PRS_TRAJ_PROT), "not available, skipping panel j\n")
}

# ---------------------------------------------------------------------------
# Panel k — PRS×MS interaction volcano
# ---------------------------------------------------------------------------
cat("Building panel k (PRS\u00d7MS interaction)...\n")

comor <- fread(COMOR_FILE)
cat("  Comorbidity rows:", nrow(comor), "\n")

comor[, log2OR   := log2(OR)]
comor[, log2OR_lo := log2(OR_lo)]
comor[, log2OR_hi := log2(OR_hi)]
comor[, direction_h := fifelse(OR >= 1, "Enriched in HLA+", "Depleted in HLA+")]
comor[, direction_h := factor(direction_h,
    levels = c("Enriched in HLA+", "Depleted in HLA+"))]
comor[, sig_star := fcase(
    fdr < 0.001, "***",
    fdr < 0.01,  "**",
    fdr < 0.05,  "*",
    fdr < 0.2,   "~",
    default      = ""
)]
comor[, label := factor(label, levels = label[order(log2OR)])]

COMOR_COLS <- c(
    "Enriched in HLA+" = "#E6A817",
    "Depleted in HLA+" = "#56B4E9"
)

# x placement for stars: just past the CI edge
comor[, star_x := fifelse(log2OR >= 0, log2OR_hi + 0.05, log2OR_lo - 0.05)]
comor[, star_hjust := fifelse(log2OR >= 0, 0, 1)]

pH <- ggplot(comor, aes(y = label, colour = direction_h)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.4, colour = "grey40") +
    geom_errorbar(aes(xmin = log2OR_lo, xmax = log2OR_hi),
                  orientation = "y", width = 0.35, linewidth = 0.5) +
    geom_point(aes(x = log2OR), size = 2) +
    geom_text(
        aes(x = star_x, label = sig_star, hjust = star_hjust),
        size   = 3,
        colour = "grey20"
    ) +
    scale_colour_manual(values = COMOR_COLS, name = NULL) +
    labs(
        title = glue("h  {PRE_LBL} comorbidity enrichment by HLA status"),
        x     = expression(log[2]~"OR (HLA+ vs HLA\u2212)"),
        y     = NULL
    ) +
    theme_ukb(base_size = 9) +
    theme(
        axis.text.y     = element_text(size = 7),
        legend.position = "bottom"
    )

save_supp(pH, "h_comorbidity", 3.5, 4.5)

# ---------------------------------------------------------------------------
# Panel n — PRS×MS interaction: does higher genetic risk amplify the MS
#   protein signature? Volcano of interaction coefficient (ms_case:prs_z)
#   across all 2,911 proteins.
# ---------------------------------------------------------------------------
cat("  (PRS\u00d7MS volcano code follows)\n")

if (!is.null(prs_int) && nrow(prs_int) > 0) {
    pi <- copy(prs_int)
    if (!"fdr" %in% names(pi)) pi[, fdr := p.adjust(P.Value, method = "BH")]

    pi[, tier := fcase(
        fdr   < 0.05 & logFC > 0, "Up FDR",
        P.Value < 0.05 & logFC > 0, "Up nominal",
        fdr   < 0.05 & logFC < 0, "Down FDR",
        P.Value < 0.05 & logFC < 0, "Down nominal",
        default = "NS"
    )]
    pi[, tier := factor(tier,
        levels = c("Up FDR", "Up nominal", "Down FDR", "Down nominal", "NS"))]
    pi[, neg_log10_p := -log10(P.Value)]

    PRS_VOL_COLS <- c(
        "Up FDR"       = "#CC0066",
        "Up nominal"   = "#F5A0C5",
        "Down FDR"     = "#56B4E9",
        "Down nominal" = "#A8D8F0",
        "NS"           = "grey74"
    )
    PRS_VOL_SZ <- c("Up FDR" = 2.8, "Up nominal" = 1.1,
                    "Down FDR" = 2.8, "Down nominal" = 1.1, "NS" = 0.45)

    # Label FDR hits (if any) + top nominal by |logFC|
    fdr_labs_n <- pi[fdr < 0.05][, label := toupper(protein)]
    nom_labs_n <- pi[P.Value < 0.05 & fdr >= 0.05][order(-abs(logFC))][
        seq_len(min(8, .N))][, label := toupper(protein)]

    n_fdr_n <- nrow(fdr_labs_n)
    n_nom_n <- sum(pi$P.Value < 0.05)

    pN <- ggplot(pi, aes(x = logFC, y = neg_log10_p)) +
        geom_point(data = pi[tier == "NS"],
                   colour = PRS_VOL_COLS["NS"], size = PRS_VOL_SZ["NS"], alpha = 0.35) +
        geom_point(data = pi[tier %in% c("Up nominal", "Down nominal")],
                   aes(colour = tier), size = PRS_VOL_SZ["Up nominal"], alpha = 0.75) +
        geom_point(data = pi[tier %in% c("Up FDR", "Down FDR")],
                   aes(colour = tier), size = PRS_VOL_SZ["Up FDR"], alpha = 1.0) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
        {if (n_fdr_n > 0)
            geom_hline(yintercept = -log10(max(pi[fdr < 0.05, P.Value], na.rm = TRUE)),
                       linetype = "dashed", linewidth = 0.3, colour = "grey45")} +
        {if (nrow(nom_labs_n) > 0)
            geom_text_repel(
                data = nom_labs_n,
                aes(label = label),
                size = 1.95, colour = "grey25",
                max.overlaps = 20, segment.size = 0.2,
                segment.colour = "grey55", box.padding = 0.28,
                show.legend = FALSE
            )} +
        {if (n_fdr_n > 0)
            geom_text_repel(
                data = fdr_labs_n,
                aes(label = label),
                size = 2.5, fontface = "bold", colour = "grey10",
                max.overlaps = Inf, segment.size = 0.3,
                segment.colour = "grey30", box.padding = 0.45,
                show.legend = FALSE
            )} +
        annotate("text",
                 x = max(pi$logFC, na.rm = TRUE) * 0.95,
                 y = max(pi$neg_log10_p, na.rm = TRUE) * 0.04,
                 label = sprintf("%d FDR<0.05  |  %d nominal p<0.05", n_fdr_n, n_nom_n),
                 size = 1.95, hjust = 1, colour = "grey40") +
        scale_colour_manual(values = PRS_VOL_COLS, name = NULL) +
        labs(
            title    = glue("k  PRS\u00d7{DISEASE_CAPS} interaction: amplification of proteomic signature"),
            subtitle = glue("Interaction coefficient ms_case:prs_z | positive = HLA-adjusted PRS amplifies {DISEASE_CAPS} signal"),
            x        = bquote(beta[PRS%*%.(DISEASE_CAPS)]~"(interaction, per SD PRS)"),
            y        = expression(-log[10]~italic(P))
        ) +
        theme_ukb(base_size = 9) +
        theme(
            legend.position = "none",
            plot.subtitle   = element_text(size = 6.5, colour = "grey45")
        )

    save_panel(pN, "k_prs_ms_interaction", 4.0, 3.8)
} else {
    cat("  Panel k skipped: ms_prs_interaction_genome.csv not found\n")
}

# ---------------------------------------------------------------------------
# Panel l — GO enrichment on nominal PRS×MS interaction proteins
# ---------------------------------------------------------------------------
cat("Building panel l (PRS\u00d7MS GO enrichment)...\n")

tryCatch({
    suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })

    if (!exists("calc_fe")) {
        calc_fe <- function(gr, br) {
            p  <- strsplit(gr, "/"); num <- as.numeric(sapply(p, `[[`, 1)); den <- as.numeric(sapply(p, `[[`, 2))
            bp <- strsplit(br, "/"); bn  <- as.numeric(sapply(bp,`[[`, 1)); bd  <- as.numeric(sapply(bp,`[[`, 2))
            (num / den) / (bn / bd)
        }
    }

    prs_up   <- if (!is.null(prs_int)) prs_int[P.Value < 0.05 & logFC > 0, toupper(protein)] else character(0)
    prs_down <- if (!is.null(prs_int)) prs_int[P.Value < 0.05 & logFC < 0, toupper(protein)] else character(0)
    cat(sprintf("  PRS\u00d7MS GO: %d amplified, %d attenuated (nominal p<0.05)\n",
                length(prs_up), length(prs_down)))

    run_ego_l <- function(genes, label) {
        if (length(genes) < 5) return(NULL)
        ego <- tryCatch(
            enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     readable = TRUE),
            error = function(e) NULL
        )
        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
        dt <- as.data.table(as.data.frame(ego))
        dt[, FoldEnrichment := calc_fe(GeneRatio, BgRatio)]
        dt[, direction := label]
        dt[order(p.adjust)][seq_len(min(8, .N))]
    }

    go_l_rows <- list(
        run_ego_l(prs_up,   "PRS amplifies"),
        run_ego_l(prs_down, "PRS attenuates")
    )
    go_l <- rbindlist(Filter(Negate(is.null), go_l_rows), fill = TRUE)
    if (nrow(go_l) == 0) stop("No significant GO terms for PRS\u00d7MS interaction")

    go_l[, direction := factor(direction, levels = c("PRS amplifies", "PRS attenuates"))]
    go_l[, Description := factor(Description,
        levels = unique(go_l[order(direction, FoldEnrichment), Description]))]

    GO_L_COLS <- c("PRS amplifies" = "#CC0066", "PRS attenuates" = "#56B4E9")

    pL <- ggplot(go_l, aes(x = FoldEnrichment, y = Description,
                            size = Count, colour = direction)) +
        geom_point() +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
        scale_colour_manual(values = GO_L_COLS,
                             guide = guide_legend(override.aes = list(size = 3.5))) +
        scale_size_continuous(name = "Gene count", range = c(2, 6),
                               guide = guide_legend(override.aes = list(colour = "grey40"))) +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        labs(title    = glue("l  GO:BP enrichment \u2014 PRS\u00d7{DISEASE_CAPS} interaction proteins"),
             subtitle = sprintf("Nominal p<0.05 PRS\u00d7%s proteins: %d amplified (logFC>0), %d attenuated (logFC<0)",
                                DISEASE_CAPS, length(prs_up), length(prs_down)),
             x = "Fold enrichment", y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y    = element_text(size = 7),
              strip.text     = element_text(size = 8, face = "bold"),
              legend.position = "right",
              plot.subtitle  = element_text(size = 7, colour = "grey40"))

    save_panel(pL, "l_prs_go", 5.5, 6)

}, error = function(e) {
    cat("  Panel l GO:", conditionMessage(e), "\n")
    pL_ph <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = sprintf("PRS\u00d7%s GO enrichment\n(%d nominal p<0.05 proteins)",
                                 DISEASE_CAPS,
                                 if (!is.null(prs_int)) sum(prs_int$P.Value < 0.05) else 0),
                 size = 3.5, hjust = 0.5, colour = "grey40") +
        xlim(0,1) + ylim(0,1) +
        labs(title = glue("l  GO:BP enrichment \u2014 PRS\u00d7{DISEASE_CAPS} interaction proteins")) +
        theme_void() + theme(plot.title = element_text(size = 9, face = "plain"))
    save_panel(pL_ph, "l_prs_go", 5.5, 6)
})

# ---------------------------------------------------------------------------
# Panels m & n — HLA × PRS × sex integration
#
#   m: 4-group violin (HLA−/+ × PRS-low/high, faceted by sex) for LILRB4 and
#      IL2RA — shows that LILRB4 is driven by HLA while IL2RA is lifted by PRS
#      specifically in the HLA+ background.
#
#   n: LILRB4 (y) × IL2RA (x) scatter per MS case.  HLA-coloured, sex-faceted.
#      LILRB4 and IL2RA are both immune proteins that co-vary (r≈0.35) but HLA
#      status provides the floor level for LILRB4, independently of PRS.
#
# Supp j: Panel e re-plotted with PRS-sig proteins (IL2RA, TNFRSF6B, CD5L,
#      KLRK1) explicitly labelled — confirms they are NOT HLA-stratified DEPs.
# ---------------------------------------------------------------------------
cat("Building panels m, n and supp j (HLA \u00d7 PRS \u00d7 sex integration)...\n")

if (file.exists(PRS_FILE)) {
    prs_mn <- fread(PRS_FILE, showProgress = FALSE)
    if (PRS_COL != "prs_score" && PRS_COL %in% names(prs_mn))
        setnames(prs_mn, PRS_COL, "prs_score")

    ms_mn <- merge(
        qc_hla[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier == FALSE &
               !is.na(drb1_1501_carrier) & !is.na(lilrb4) & !is.na(il2ra)],
        prs_mn[, .(eid, prs_score)],
        by = "eid"
    )
    ms_mn[, prs_q2     := fifelse(prs_score >= median(prs_score, na.rm=TRUE), "PRS-high", "PRS-low")]
    ms_mn[, hla_label  := fifelse(drb1_1501_carrier == 1, "HLA+", "HLA\u2212")]
    ms_mn[, sex_label  := factor(fifelse(sex == 0, "Female", "Male"), levels = c("Female","Male"))]

    GRP_LEVELS <- c("HLA\u2212 PRS-low", "HLA\u2212 PRS-high", "HLA+ PRS-low", "HLA+ PRS-high")
    ms_mn[, geno_grp := factor(paste(hla_label, prs_q2), levels = GRP_LEVELS)]

    GENO_FILLS <- c(
        "HLA\u2212 PRS-low"  = "#A8D8F0",
        "HLA\u2212 PRS-high" = "#2B4C7E",
        "HLA+ PRS-low"       = "#F5C97A",
        "HLA+ PRS-high"      = "#E6A817"
    )

    cat(sprintf("  Integration cohort: n=%d  (F=%d, M=%d)\n",
                nrow(ms_mn), sum(ms_mn$sex==0), sum(ms_mn$sex==1)))

    # ── Panel m: 4-group violin for LILRB4 (HLA marker) + IL2RA (PRS marker) ──
    m_long <- melt(ms_mn[, .(eid, geno_grp, sex_label, lilrb4, il2ra)],
                   id.vars    = c("eid","geno_grp","sex_label"),
                   measure.vars = c("lilrb4","il2ra"),
                   variable.name = "protein", value.name = "npx")
    m_long[, prot_label := fifelse(protein == "lilrb4",
                                   HLA_MKR_LABEL,
                                   PRS_MKR_LABEL)]
    m_long[, prot_label := factor(prot_label,
                                   levels = MKR_LEVELS)]
    m_long[, fill_col := geno_grp]

    pM <- ggplot(m_long[!is.na(npx)],
                 aes(x = geno_grp, y = npx, fill = fill_col)) +
        geom_violin(trim = TRUE, linewidth = 0.25, alpha = 0.85) +
        geom_boxplot(width = 0.10, outlier.shape = NA, linewidth = 0.25,
                     fill = "white", colour = "grey30") +
        facet_grid(prot_label ~ sex_label, scales = "free_y") +
        scale_fill_manual(values = GENO_FILLS, guide = "none") +
        labs(
            title    = "m  Genetic stratification: HLA \u00d7 PRS \u00d7 sex",
            subtitle = glue("{DISEASE_CAPS} cases: HLA\u00b1 \u00d7 PRS median split. ",
                              "{HLA_MARKER_DISP} driven by HLA; {PRS_MARKER_DISP} highest in HLA+ \u00d7 PRS-high"),
            x = NULL, y = "NPX"
        ) +
        theme_ukb(base_size = 9) +
        theme(axis.text.x    = element_text(angle = 35, hjust = 1, size = 6.5),
              strip.text     = element_text(size = 7.5, face = "bold"),
              plot.subtitle  = element_text(size = 7, colour = "grey40"))

    save_panel(pM, "m_hla_prs_sex_stratification", 6.5, 5.5)

    # ── Panel n: LILRB4 × IL2RA scatter — within-individual co-variation ──────
    n_dat <- ms_mn[!is.na(lilrb4) & !is.na(il2ra)]

    r_pos <- round(cor(n_dat[drb1_1501_carrier==1, lilrb4],
                       n_dat[drb1_1501_carrier==1, il2ra], use="complete.obs"), 2)
    r_neg <- round(cor(n_dat[drb1_1501_carrier==0, lilrb4],
                       n_dat[drb1_1501_carrier==0, il2ra], use="complete.obs"), 2)

    HLA_N_COLS <- c("HLA+" = "#E6A817", "HLA\u2212" = "#2B4C7E")

    pN_main <- ggplot(n_dat, aes(x = il2ra, y = lilrb4, colour = hla_label)) +
        geom_point(size = 0.65, alpha = 0.40) +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                    linewidth = 0.75, alpha = 0.14) +
        scale_colour_manual(values = HLA_N_COLS, name = NULL) +
        annotate("text",
                 x = Inf, y = Inf,
                 label = sprintf("HLA+ r = %.2f\nHLA\u2212 r = %.2f", r_pos, r_neg),
                 hjust = 1.1, vjust = 1.5, size = 2.5, colour = "grey30") +
        facet_wrap(~ sex_label, nrow = 1) +
        scale_x_continuous(expand = expansion(mult = c(0.04, 0.12))) +
        labs(
            title    = glue("n  {HLA_MARKER_DISP} vs {PRS_MARKER_DISP}: co-variation within {DISEASE_CAPS}"),
            subtitle = "HLA and PRS protein markers co-vary (shared immune activation) but are independently driven",
            x        = glue("{PRS_MARKER_DISP} (NPX, PRS-correlated)"),
            y        = glue("{HLA_MARKER_DISP} (NPX, HLA marker)")
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom",
              plot.subtitle   = element_text(size = 7, colour = "grey40"))

    save_panel(pN_main,
               glue("n_{tolower(HLA_MARKER_DISP)}_{tolower(PRS_MARKER_DISP)}_scatter"),
               5.5, 3.8)

} else {
    cat("  PRS file not found; skipping panels m and n\n")
}

# ── Supp panel j — Panel e rescored: HLA logFC scatter with PRS-sig overlay ──
cat("Building supp panel j (HLA\u00d7PRS protein overlap)...\n")

if (exists("strat_d")) {
    PRS_SIG_PROTS <- c("il2ra", "tnfrsf6b", "cd5l", "klrk1")

    strat_j <- copy(strat_d)
    strat_j[, overlap_cat := fcase(
        tolower(protein) %in% PRS_SIG_PROTS & (fdr_hlaneg < 0.05 | fdr_hlapos < 0.05),
            "PRS-sig AND HLA-dep",
        tolower(protein) %in% PRS_SIG_PROTS,
            "PRS-sig (HLA-independent)",
        fdr_hlaneg < 0.05 | fdr_hlapos < 0.05,
            "HLA-dep only",
        default = "NS"
    )]
    strat_j[, overlap_cat := factor(overlap_cat,
        levels = c("PRS-sig AND HLA-dep", "PRS-sig (HLA-independent)",
                   "HLA-dep only", "NS"))]

    OVL_COLS <- c(
        "PRS-sig AND HLA-dep"       = "#CC0066",
        "PRS-sig (HLA-independent)" = "#56B4E9",
        "HLA-dep only"               = "#E6A817",
        "NS"                         = "grey75"
    )
    OVL_SZ   <- c("PRS-sig AND HLA-dep"=3.0, "PRS-sig (HLA-independent)"=2.8,
                  "HLA-dep only"=2.0, "NS"=0.55)
    prs_lbl_j <- strat_j[tolower(protein) %in% PRS_SIG_PROTS,
                          .(protein = toupper(protein), logFC_hlaneg, logFC_hlapos, overlap_cat)]

    pSJ <- ggplot(strat_j, aes(x = logFC_hlaneg, y = logFC_hlapos,
                                colour = overlap_cat, size = overlap_cat)) +
        geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey70") +
        geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    linewidth = 0.35, colour = "grey50") +
        geom_point(data = strat_j[overlap_cat == "NS"], alpha = 0.30) +
        geom_point(data = strat_j[overlap_cat != "NS"], alpha = 0.92) +
        geom_text_repel(
            data          = prs_lbl_j,
            aes(label     = protein),
            size          = 2.5, fontface = "bold",
            max.overlaps  = Inf, segment.size = 0.25,
            segment.color = "grey50", box.padding = 0.4,
            show.legend   = FALSE
        ) +
        scale_colour_manual(values = OVL_COLS, name = NULL) +
        scale_size_manual(values  = OVL_SZ,   guide = "none") +
        annotate("text",
                 x = range(strat_j$logFC_hlaneg, na.rm=TRUE)[2] * 0.9,
                 y = range(strat_j$logFC_hlapos, na.rm=TRUE)[1] * 0.9,
                 label = "PRS-sig proteins\nnot HLA-stratified",
                 size = 2.2, hjust = 1, colour = "grey35", fontface = "italic") +
        labs(
            title    = "j  HLA-logFC scatter: PRS-sig proteins are HLA-independent",
            subtitle = "PRS FDR<0.05 proteins (IL2RA, TNFRSF6B, CD5L, KLRK1) plotted on panel e axes",
            x = expression(log[2]~"FC (HLA\u2212 background)"),
            y = expression(log[2]~"FC (HLA+ background)")
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom",
              plot.subtitle   = element_text(size = 7, colour = "grey40"))

    save_supp(pSJ, "j_hla_prs_overlap", 4.0, 4.2)
} else {
    cat("  strat_d not available; skipping supp panel j\n")
}

# ---------------------------------------------------------------------------
# Supp panel i — Pre-MS comorbidity enrichment by LILRB4 level
#   Parallel to the HLA comorbidity (supp h), but stratified by LILRB4
#   expression level (high vs low, median split within pre-MS cases).
#   Uses the pre-computed composite feature matrix from
#   analysis/15_endophenotype/02_ms_prems_comorbidity_endophenotype.R.
# ---------------------------------------------------------------------------
cat("Building supp panel i (LILRB4 pre-MS comorbidity)...\n")

FEAT_FILE <- file.path(PROJ_DIR, "results", "endophenotype",
                        glue("{COHORT}_prems_composite_features.csv"))

if (file.exists(FEAT_FILE) && "lilrb4" %in% names(qc)) {
    feat_all <- fread(FEAT_FILE, showProgress = FALSE)
    if (STATUS_COL != "ms_status" && STATUS_COL %in% names(feat_all))
        setnames(feat_all, STATUS_COL, "ms_status")
    # Restrict to pre-onset cases; merge with QC to get LILRB4 NPX
    feat_pre <- feat_all[ms_status == SV$pre_onset]
    feat_pre <- merge(feat_pre,
                      qc[qc_outlier == FALSE, .(eid, lilrb4)],
                      by = "eid", all.x = FALSE)
    feat_pre <- feat_pre[!is.na(lilrb4)]
    feat_pre[, lilrb4_high := as.integer(lilrb4 >= median(lilrb4, na.rm = TRUE))]
    n_high_lr <- sum(feat_pre$lilrb4_high == 1)
    n_low_lr  <- sum(feat_pre$lilrb4_high == 0)
    cat(sprintf("  %s with LILRB4: n=%d (high=%d, low=%d)\n",
                PRE_LBL, nrow(feat_pre), n_high_lr, n_low_lr))

    COMPOSITE_GROUPS <- c(
        "Migraine", "Optic_neuritis", "Demyelinating_NOS", "Neurological",
        "Eye", "Hypertension", "Cardiovascular", "Atopic",
        "Resp_infections", "Digestive", "Skin", "Musculoskeletal",
        "Gynaecological", "Urological"
    )
    COMP_LABELS <- c(
        Migraine="Migraine", Optic_neuritis="Optic neuritis",
        Demyelinating_NOS="Demyelinating NOS", Neurological="Neurological",
        Eye="Eye disorders", Hypertension="Hypertension",
        Cardiovascular="Cardiovascular", Atopic="Atopic/allergic",
        Resp_infections="Resp. infections", Digestive="Digestive",
        Skin="Skin", Musculoskeletal="Musculoskeletal",
        Gynaecological="Gynaecological", Urological="Urological"
    )

    avail_groups <- intersect(COMPOSITE_GROUPS, names(feat_pre))

    lilrb4_comor <- rbindlist(lapply(avail_groups, function(grp) {
        d <- feat_pre[!is.na(get(grp)),
                      .(y = as.integer(get(grp)), lilrb4_high,
                        age = age_at_diagnosis, sex)]
        if (sum(d$y == 1, na.rm = TRUE) < 3) return(NULL)
        fit <- tryCatch(
            glm(y ~ lilrb4_high + age + sex, data = d, family = binomial),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        cf <- summary(fit)$coefficients
        rn <- rownames(cf)
        if (!"lilrb4_high" %in% rn) return(NULL)
        row_lr <- cf["lilrb4_high", ]
        est    <- row_lr["Estimate"]
        se     <- row_lr["Std. Error"]
        OR     <- exp(est);  OR_lo <- exp(est - 1.96 * se);  OR_hi <- exp(est + 1.96 * se)
        pval   <- row_lr["Pr(>|z|)"]
        n_h <- sum(d$lilrb4_high == 1 & d$y == 1, na.rm = TRUE)
        n_l <- sum(d$lilrb4_high == 0 & d$y == 1, na.rm = TRUE)
        data.table(group = grp, label = COMP_LABELS[grp],
                   OR, OR_lo, OR_hi, pval, n_high = n_h, n_low = n_l)
    }))

    if (!is.null(lilrb4_comor) && nrow(lilrb4_comor) > 0) {
        lilrb4_comor[, fdr       := p.adjust(pval, method = "BH")]
        lilrb4_comor[, log2OR    := log2(OR)]
        lilrb4_comor[, log2OR_lo := log2(OR_lo)]
        lilrb4_comor[, log2OR_hi := log2(OR_hi)]
        lilrb4_comor[, direction := fifelse(OR >= 1,
                                            glue("Enriched in {HLA_MARKER_DISP}-high"),
                                            glue("Depleted in {HLA_MARKER_DISP}-high"))]
        lilrb4_comor[, direction := factor(direction,
            levels = c(glue("Enriched in {HLA_MARKER_DISP}-high"),
                       glue("Depleted in {HLA_MARKER_DISP}-high")))]
        lilrb4_comor[, sig_star := fcase(
            fdr < 0.001, "***", fdr < 0.01, "**", fdr < 0.05, "*", fdr < 0.2, "~", default = ""
        )]
        lilrb4_comor[, label := factor(label, levels = label[order(log2OR)])]
        lilrb4_comor[, star_x    := fifelse(log2OR >= 0, log2OR_hi + 0.05, log2OR_lo - 0.05)]
        lilrb4_comor[, star_hjust := fifelse(log2OR >= 0, 0, 1)]

        fwrite(lilrb4_comor,
               file.path(PROJ_DIR, "results", "endophenotype",
                         glue("{COHORT}_{tolower(HLA_MARKER_DISP)}_prems_comorbidity.csv")))

        LILRB4_COMOR_COLS <- setNames(
            c("#E6A817", "#2B4C7E"),
            c(glue("Enriched in {HLA_MARKER_DISP}-high"),
              glue("Depleted in {HLA_MARKER_DISP}-high")))

        pI_supp <- ggplot(lilrb4_comor, aes(y = label, colour = direction)) +
            geom_vline(xintercept = 0, linetype = "dashed",
                       linewidth = 0.4, colour = "grey40") +
            geom_errorbar(aes(xmin = log2OR_lo, xmax = log2OR_hi),
                          orientation = "y", width = 0.35, linewidth = 0.5) +
            geom_point(aes(x = log2OR), size = 2) +
            geom_text(
                aes(x = star_x, label = sig_star, hjust = star_hjust),
                size = 3, colour = "grey20"
            ) +
            scale_colour_manual(values = LILRB4_COMOR_COLS, name = NULL) +
            labs(
                title    = glue("i  {PRE_LBL} comorbidity enrichment by LILRB4 level"),
                subtitle = sprintf(
                    "High LILRB4 (top 50%%, n=%d) vs Low (n=%d) within %s cases | age + sex adjusted",
                    n_high_lr, n_low_lr, PRE_LBL),
                x        = expression(log[2]~"OR (LILRB4-high vs LILRB4-low)"),
                y        = NULL
            ) +
            theme_ukb(base_size = 9) +
            theme(
                axis.text.y     = element_text(size = 7),
                legend.position = "bottom",
                plot.subtitle   = element_text(size = 6.5, colour = "grey45")
            )

        save_supp(pI_supp, "i_lilrb4_comorbidity", 4.0, 4.5)

        n_sig_lc <- sum(lilrb4_comor$fdr < 0.05, na.rm = TRUE)
        n_nom_lc <- sum(lilrb4_comor$pval < 0.05, na.rm = TRUE)
        cat(sprintf("  LILRB4 comorbidity: FDR<0.05=%d  nominal p<0.05=%d\n",
                    n_sig_lc, n_nom_lc))
        top_lc <- lilrb4_comor[order(pval)][seq_len(min(5, .N))]
        cat("  Top associations:\n")
        print(top_lc[, .(label, OR=round(OR,2), pval=signif(pval,3), fdr=signif(fdr,3))],
              row.names=FALSE)
    } else {
        cat("  No LILRB4 comorbidity results — skipping panel i\n")
    }
} else {
    cat("  Skipping supp i: ms_prems_composite_features.csv or lilrb4 not found\n")
}

# ---------------------------------------------------------------------------
# Panel o — HLA × PRS genetic groups: MS cases vs age/sex-matched HC
#   Panel m re-drawn with PSM-matched HC (faded) alongside MS (coloured).
#   Left violin per group = HC; right = MS. Shows which genetic effects
#   are MS-specific vs present at baseline in healthy controls.
# ---------------------------------------------------------------------------
cat("Building panel o (HLA \u00d7 PRS: MS vs matched HC)...\n")

HAS_MATCHIT <- requireNamespace("MatchIt", quietly = TRUE)

if (file.exists(PRS_FILE) && exists("ms_mn") && exists("prs_mn")) {

    # --- HC pool: controls with HLA + PRS + at least one protein ----------------
    hc_pool_o <- merge(
        qc_hla[ms_status == SV$control & qc_outlier == FALSE &
               !is.na(drb1_1501_carrier) & !is.na(age_at_sampling) & !is.na(sex)],
        prs_mn[, .(eid, prs_score)],
        by = "eid"
    )
    hc_pool_o <- hc_pool_o[!is.na(lilrb4) | !is.na(il2ra)]
    cat(sprintf("  HC with HLA+PRS: %d  (HLA+=%d  HLA\u2212=%d)\n",
                nrow(hc_pool_o),
                sum(hc_pool_o$drb1_1501_carrier == 1),
                sum(hc_pool_o$drb1_1501_carrier == 0)))

    # --- PSM: match HC 1:3 to MS cases (nearest-neighbour on age + sex) ----------
    prs_med_ms_o <- median(ms_mn$prs_score, na.rm = TRUE)

    if (HAS_MATCHIT) {
        suppressPackageStartupMessages(library(MatchIt))
        pool_o <- rbind(
            ms_mn[, .(eid, is_ms = 1L, age_at_sampling, sex)],
            hc_pool_o[!eid %in% ms_mn$eid, .(eid, is_ms = 0L, age_at_sampling, sex)]
        )
        set.seed(42)
        m_obj_o <- tryCatch(
            matchit(is_ms ~ age_at_sampling + sex,
                    data   = as.data.frame(pool_o),
                    method = "nearest", ratio = 3, replace = FALSE),
            error = function(e) {
                cat("  MatchIt failed:", conditionMessage(e), "\n"); NULL
            }
        )
        if (!is.null(m_obj_o)) {
            md_o         <- match.data(m_obj_o)
            matched_eids <- as.integer(md_o[md_o$is_ms == 0, "eid"])
            hc_matched_o <- hc_pool_o[eid %in% matched_eids]
            cat(sprintf("  PSM-matched HC: n=%d (3:1 ratio to %d MS cases)\n",
                        nrow(hc_matched_o), nrow(ms_mn)))
        } else {
            hc_matched_o <- hc_pool_o
        }
    } else {
        cat("  MatchIt not available \u2014 using all HC (no PSM)\n")
        hc_matched_o <- hc_pool_o
    }

    # --- Assign HC to 4 genetic groups (MS-derived PRS median) ------------------
    hc_matched_o[, prs_q2_hc  := fifelse(prs_score >= prs_med_ms_o, "PRS-high", "PRS-low")]
    hc_matched_o[, hla_lbl_hc := fifelse(drb1_1501_carrier == 1, "HLA+", "HLA\u2212")]
    hc_matched_o[, geno_grp   := factor(paste(hla_lbl_hc, prs_q2_hc), levels = GRP_LEVELS)]
    hc_matched_o[, disease_grp := "HC"]

    ms_o <- copy(ms_mn)
    ms_o[, disease_grp := DISEASE_CAPS]

    # --- Stack and melt ---------------------------------------------------------
    shared_o <- c("eid", "geno_grp", "disease_grp", "lilrb4", "il2ra")
    o_stack  <- rbind(ms_o[, ..shared_o], hc_matched_o[, ..shared_o])
    o_stack[, disease_grp := factor(disease_grp, levels = c("HC", DISEASE_CAPS))]

    o_long <- melt(o_stack,
                   id.vars      = c("eid", "geno_grp", "disease_grp"),
                   measure.vars = c("lilrb4", "il2ra"),
                   variable.name = "protein", value.name = "npx")
    o_long[, prot_label := fifelse(protein == "lilrb4",
                                   HLA_MKR_LABEL, PRS_MKR_LABEL)]
    o_long[, prot_label := factor(prot_label,
                                   levels = MKR_LEVELS)]
    # group = HC or MS within each geno_grp → drives position_dodge split
    o_long[, fill_key := paste(as.character(geno_grp), as.character(disease_grp))]

    # --- Palette: MS = GENO_FILLS; HC = same colour faded (alpha ≈ 0.35) --------
    HC_FILLS_O <- setNames(
        sapply(GENO_FILLS, function(col) adjustcolor(col, alpha.f = 0.35)),
        paste(names(GENO_FILLS), "HC")
    )
    MS_FILLS_O  <- setNames(GENO_FILLS, paste(names(GENO_FILLS), DISEASE_CAPS))
    ALL_FILLS_O <- c(HC_FILLS_O, MS_FILLS_O)

    n_o <- o_long[!is.na(npx) & protein == "lilrb4",
                  .(n = .N), by = .(geno_grp, disease_grp)]
    cat("  LILRB4 n per group:\n"); print(n_o, row.names = FALSE)

    pO <- ggplot(o_long[!is.na(npx)],
                 aes(x     = geno_grp,
                     y     = npx,
                     fill  = fill_key,
                     group = fill_key)) +
        geom_violin(trim = TRUE, linewidth = 0.22,
                    position = position_dodge(width = 0.9)) +
        geom_boxplot(width = 0.08, outlier.shape = NA, linewidth = 0.22,
                     fill = "white", colour = "grey30",
                     position = position_dodge(width = 0.9)) +
        annotate("text", x = Inf, y = Inf,
                 label = glue("Left per group: HC (faded)\nRight per group: {DISEASE_CAPS} (solid)"),
                 hjust = 1.08, vjust = 1.5, size = 2.0, colour = "grey40",
                 fontface = "italic") +
        facet_grid(prot_label ~ ., scales = "free_y") +
        scale_fill_manual(values = ALL_FILLS_O, guide = "none") +
        labs(
            title    = glue("o  Genetic groups: {DISEASE_CAPS} cases vs age/sex-matched controls"),
            subtitle = paste0(
                "HC PSM-matched 3:1 (age+sex). Faded = HC baseline; solid = ", DISEASE_CAPS, ". ",
                "Elevation above HC reveals disease-specific genetic amplification."),
            x = NULL, y = "NPX"
        ) +
        theme_ukb(base_size = 9) +
        theme(
            axis.text.x   = element_text(angle = 35, hjust = 1, size = 6.5),
            strip.text    = element_text(size = 7.5, face = "bold"),
            plot.subtitle = element_text(size = 7, colour = "grey40")
        )

    save_panel(pO, "o_hla_prs_ms_vs_hc", 6.5, 5.5)

} else {
    cat("  Panel o skipped: PRS file or ms_mn / prs_mn not available\n")
}

# ---------------------------------------------------------------------------
# Supp panel k — Formal regression: IL2RA & LILRB4 ~ HLA × PRS_z + age + sex
#   Coefficient forest plot (95% CI) for key terms from within-MS lm.
#   Tests whether IL2RA elevation in HLA+ × PRS-high is statistically supported.
# ---------------------------------------------------------------------------
cat("Building supp panel k (IL2RA/LILRB4 regression forest)...\n")

if (file.exists(PRS_FILE) && exists("ms_mn")) {

    ms_reg <- copy(ms_mn)
    ms_reg[, prs_z := as.numeric(scale(prs_score))]
    ms_reg[, hla_f := factor(drb1_1501_carrier, levels = c(0, 1),
                              labels = c("HLA\u2212", "HLA+"))]
    ms_reg[, sex_f := factor(sex, levels = c(0, 1), labels = c("Female", "Male"))]

    TERM_LABELS_K <- c(
        "hla_fHLA+"         = "HLA+ (main)",
        "prs_z"             = "PRS-z (main)",
        "hla_fHLA+:prs_z"  = "HLA+ \u00d7 PRS-z (interaction)",
        "age_at_sampling"   = "Age at sampling",
        "sex_fMale"         = "Sex: Male"
    )

    reg_rows <- rbindlist(lapply(c("il2ra", "lilrb4"), function(prot) {
        d <- ms_reg[!is.na(get(prot)) & !is.na(prs_z)]
        if (nrow(d) < 50) return(NULL)
        fit <- tryCatch(
            lm(as.formula(sprintf(
                "%s ~ hla_f * prs_z + age_at_sampling + sex_f", prot)),
               data = d),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        cf   <- summary(fit)$coefficients
        ci   <- confint(fit)
        rn   <- rownames(cf)
        keep <- rn %in% names(TERM_LABELS_K)
        if (!any(keep)) return(NULL)
        data.table(
            protein  = toupper(prot),
            term     = rn[keep],
            estimate = cf[keep, "Estimate"],
            lci      = ci[keep, 1],
            uci      = ci[keep, 2],
            pval     = cf[keep, "Pr(>|t|)"]
        )
    }))

    if (!is.null(reg_rows) && nrow(reg_rows) > 0) {
        reg_rows[, term_label := factor(TERM_LABELS_K[term],
                                         levels = rev(unname(TERM_LABELS_K)))]
        reg_rows[, sig_star := fcase(
            pval < 0.001, "***", pval < 0.01, "**",
            pval < 0.05,  "*",   pval < 0.1,  "~",
            default = ""
        )]
        reg_rows[, is_int := term == "hla_fHLA+:prs_z"]

        cat("  Interaction + main terms:\n")
        print(
            reg_rows[term %in% c("hla_fHLA+", "prs_z", "hla_fHLA+:prs_z"),
                     .(protein, term_label, estimate = round(estimate, 4),
                       lci = round(lci, 4), uci = round(uci, 4),
                       pval = signif(pval, 3))],
            row.names = FALSE
        )

        PROT_COLS_K <- setNames(c("#CC0066", "#E6A817"),
                                 c(PRS_MARKER_DISP, HLA_MARKER_DISP))

        pSK <- ggplot(reg_rows,
                      aes(y = term_label, x = estimate,
                          colour = protein, shape = is_int)) +
            geom_vline(xintercept = 0, linetype = "dashed",
                       linewidth = 0.4, colour = "grey40") +
            geom_errorbar(aes(xmin = lci, xmax = uci),
                          orientation = "y", width = 0.25, linewidth = 0.6,
                          position = position_dodge(width = 0.7)) +
            geom_point(size = 2.4,
                       position = position_dodge(width = 0.7)) +
            geom_text(aes(x = uci, label = sig_star),
                      hjust = -0.3, size = 3.0, colour = "grey20",
                      position = position_dodge(width = 0.7)) +
            scale_colour_manual(values = PROT_COLS_K, name = NULL) +
            scale_shape_manual(values = c("FALSE" = 16L, "TRUE" = 18L), guide = "none") +
            labs(
                title    = "k  Formal regression: IL2RA & LILRB4 ~ HLA \u00d7 PRS",
                subtitle = paste0(
                    "Within-", DISEASE_CAPS, " cases | protein ~ HLA + PRS-z + HLA:PRS-z + age + sex",
                    " | diamond = interaction term | 95% CI"),
                x = "\u0394NPX per unit (regression coefficient)",
                y = NULL
            ) +
            theme_ukb(base_size = 9) +
            theme(
                legend.position = "bottom",
                plot.subtitle   = element_text(size = 7, colour = "grey40"),
                axis.text.y     = element_text(size = 8)
            )

        save_supp(pSK, "k_il2ra_lilrb4_regression", 5.5, 4.0)

    } else {
        cat("  Regression failed; skipping supp k\n")
    }
} else {
    cat("  Supp k skipped: PRS file or ms_mn not available\n")
}

# ---------------------------------------------------------------------------
# Supp panel l — MS vs matched-HC delta per HLA×PRS group (formalises panel o)
#   Forest plot: mean difference (MS − HC) ± 95% CI for LILRB4 and IL2RA.
#   Wilcoxon FDR corrected within each protein across the 4 genetic groups.
#   Key finding: LILRB4 = baseline HLA-carrier effect (no MS-specific signal);
#                IL2RA = MS-specific signal, significant in HLA+ × PRS-high.
# ---------------------------------------------------------------------------
cat("Building supp panel l (MS vs HC delta forest)...\n")

if (file.exists(PRS_FILE) && exists("o_stack") && nrow(o_stack) > 0) {

    o_long_l <- melt(o_stack,
                     id.vars      = c("eid", "geno_grp", "disease_grp"),
                     measure.vars = c("lilrb4", "il2ra"),
                     variable.name = "protein", value.name = "npx")
    o_long_l[, prot_label := fifelse(protein == "lilrb4",
                                     HLA_MKR_LABEL, PRS_MKR_LABEL)]
    o_long_l[, prot_label := factor(prot_label,
                                     levels = MKR_LEVELS)]

    # Mean, SE, Wilcoxon p per group × protein
    delta_l <- rbindlist(lapply(levels(o_stack$geno_grp), function(grp) {
        rbindlist(lapply(c("lilrb4", "il2ra"), function(prot) {
            ms_v <- o_long_l[geno_grp == grp & disease_grp == DISEASE_CAPS & protein == prot &
                             !is.na(npx), npx]
            hc_v <- o_long_l[geno_grp == grp & disease_grp == "HC" & protein == prot &
                             !is.na(npx), npx]
            if (length(ms_v) < 3 || length(hc_v) < 3) return(NULL)
            wt    <- wilcox.test(ms_v, hc_v, exact = FALSE)
            delta <- mean(ms_v) - mean(hc_v)
            se    <- sqrt(var(ms_v) / length(ms_v) + var(hc_v) / length(hc_v))
            data.table(geno_grp = grp, protein = prot,
                       prot_label = fifelse(prot == "lilrb4",
                                            "LILRB4 (HLA marker)", "IL2RA (PRS marker)"),
                       delta, lci = delta - 1.96 * se, uci = delta + 1.96 * se,
                       p = wt$p.value, n_ms = length(ms_v), n_hc = length(hc_v))
        }))
    }))
    delta_l[, fdr := p.adjust(p, method = "BH"), by = protein]
    delta_l[, sig_star := fcase(
        fdr < 0.001, "***", fdr < 0.01, "**",
        fdr < 0.05,  "*",   fdr < 0.1,  "~",
        default = ""
    )]
    delta_l[, prot_label  := factor(prot_label,
                                     levels = MKR_LEVELS)]
    delta_l[, geno_grp    := factor(geno_grp, levels = rev(GRP_LEVELS))]
    delta_l[, sig_cat := fcase(
        fdr < 0.05,  "FDR<0.05",
        fdr < 0.1,   "Trending (FDR<0.1)",
        default      = "NS"
    )]
    delta_l[, sig_cat := factor(sig_cat,
        levels = c("FDR<0.05", "Trending (FDR<0.1)", "NS"))]

    DELTA_COLS <- c("FDR<0.05" = "#CC0066", "Trending (FDR<0.1)" = "#F5A0C5", "NS" = "grey60")

    pSL <- ggplot(delta_l, aes(y = geno_grp, x = delta, colour = sig_cat)) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.4, colour = "grey40") +
        geom_errorbar(aes(xmin = lci, xmax = uci),
                      orientation = "y", width = 0.3, linewidth = 0.6) +
        geom_point(size = 2.4) +
        geom_text(aes(x = uci, label = sig_star),
                  hjust = -0.3, size = 3.0, colour = "grey20") +
        facet_wrap(~ prot_label, nrow = 1, scales = "free_x") +
        scale_colour_manual(values = DELTA_COLS, name = NULL) +
        labs(
            title    = glue("l  {DISEASE_CAPS} vs matched HC: protein elevation per genetic group"),
            subtitle = paste0(
                "Mean difference (", DISEASE_CAPS, " \u2212 HC) \u00b1 95% CI | PSM-matched HC (3:1 age+sex) | ",
                "Wilcoxon FDR corrected per protein"),
            x = glue("\u0394NPX ({DISEASE_CAPS} \u2212 HC mean)"),
            y = NULL
        ) +
        theme_ukb(base_size = 9) +
        theme(
            legend.position = "bottom",
            strip.text      = element_text(size = 8, face = "bold"),
            plot.subtitle   = element_text(size = 7, colour = "grey40")
        )

    save_supp(pSL, "l_ms_hc_delta_per_group", 6.0, 3.5)

} else {
    cat("  Supp l skipped: o_stack not available (panel o must run first)\n")
}

# ---------------------------------------------------------------------------
# Supp m — Walchli+HPA cell-type heatmap: HLA-interaction DEPs
#   Mirrors figure 1 panel g for the HLA-stratified DEP list.
# ---------------------------------------------------------------------------
cat("Building supp panel m (cell-type heatmap for HLA DEPs)...\n")
tryCatch({
    suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer) })

    WH_FILE_F3 <- file.path(PROJ_DIR, "results", "networks",
                             glue("{COHORT}_combined_walchli_hpa_matrix.csv"))
    if (!file.exists(WH_FILE_F3)) stop("Walchli+HPA matrix not found")

    wh_f3 <- fread(WH_FILE_F3)
    ct_cols_f3 <- setdiff(names(wh_f3), c("gene","logFC","direction","abs_logFC",
                                            "vm","cm","im","compartment"))

    # Use Walchli proteins ranked by HLA P.Value — matrix is a curated subset
    hla_in_wh  <- hla_dep[toupper(protein) %in% toupper(wh_f3$gene)]
    hla_sig    <- hla_in_wh[order(P.Value)][seq_len(min(40, .N))]
    hla_genes_uc <- toupper(hla_sig$protein)

    wh_hla <- wh_f3[toupper(gene) %in% hla_genes_uc]
    if (nrow(wh_hla) < 3) stop("Too few HLA DEPs in matrix")

    mat_m <- as.matrix(wh_hla[, ..ct_cols_f3])
    rownames(mat_m) <- toupper(wh_hla$gene)
    mat_m[is.na(mat_m)] <- 0

    ct_order_m <- c("Endothelial","Pericytes","Smooth muscle","Fibroblasts",
                    "Oligodendrocytes","Microglia","Astrocytes","Neurons",
                    "T cells","B cells","NK cells","Monocytes",
                    "Dendritic cells","Macrophages","Plasma cells","Granulocytes")
    ct_pres_m <- intersect(ct_order_m, colnames(mat_m))
    mat_disp_m <- t(mat_m[, ct_pres_m, drop=FALSE])

    ct_grp_m <- ifelse(ct_pres_m %in% c("Endothelial","Pericytes","Smooth muscle","Fibroblasts"),
                       "Vascular",
                ifelse(ct_pres_m %in% c("Oligodendrocytes","Microglia","Astrocytes","Neurons"),
                       "CNS", "Immune"))
    ann_row_m  <- data.frame(Compartment = ct_grp_m, row.names = ct_pres_m)
    ann_col_m  <- data.frame(Direction = hla_sig[match(colnames(mat_disp_m), toupper(protein)),
                                                  ifelse(logFC>0,"Up HLA+","Down HLA+")],
                              row.names = colnames(mat_disp_m))
    gaps_m <- c(sum(ct_grp_m=="Vascular"), sum(ct_grp_m=="Vascular")+sum(ct_grp_m=="CNS"))
    gaps_m <- gaps_m[gaps_m > 0 & gaps_m < nrow(mat_disp_m)]

    ann_cols_m <- list(
        Direction   = c("Up HLA+"="#E6A817", "Down HLA+"="#56B4E9"),
        Compartment = c("Vascular"="#E67E22","CNS"="#27AE60","Immune"="#2980B9")
    )

    out_m <- file.path(SUPP_DIR, "panel_m_hla_celltype_heatmap.pdf")
    pdf(out_m, width=12, height=5, onefile=FALSE)
    pheatmap(mat_disp_m,
             color             = colorRampPalette(c("white","#E6A817"))(100),
             breaks            = seq(0, 1, length.out=101),
             cluster_rows      = FALSE, cluster_cols = TRUE,
             fontsize_row      = 8, fontsize_col = 6,
             annotation_col    = ann_col_m, annotation_row = ann_row_m,
             annotation_colors = ann_cols_m,
             annotation_names_col = FALSE, annotation_names_row = FALSE,
             border_color      = NA,
             gaps_row          = if (length(gaps_m)>0) gaps_m else NULL,
             main = glue("m  HLA-interaction {DISEASE_CAPS} DEPs: Walchli+HPA cell-type expression"))
    dev.off()
    cat("  Saved: 3_supp/panel_m_hla_celltype_heatmap.pdf\n")
}, error = function(e) cat("  Supp m skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Supp n — Cell-type enrichment: HLA+ vs HLA- proteins  [REFACTORED 2025]
#   Preranked fgsea on limma moderated-t from the HLA main-effect test
#   (`ms_hla_dep_results.csv`), using specificity-filtered markers from the
#   shared helper. Replaces the previous Walchli+HPA matrix Fisher test which
#   used an un-normalised `expr > 0.5` threshold and was confounded by marker
#   set size (see analysis/helpers/celltype_overrep_plot.R header).
# ---------------------------------------------------------------------------
cat("Building supp panel n (cell-type fgsea for HLA DEPs)...\n")
tryCatch({
    WALCHLI_CSV_N <- file.path(PROJ_DIR, "data", "reference",
                                "walchli_avg_expr_by_celltype.csv")
    HPA_FILE_N    <- file.path(PROJ_DIR, "data", "reference",
                                "hpa_rna_single_cell_type.tsv.zip")
    marker_sets_n <- load_celltype_markers(
        walchli_csv = WALCHLI_CSV_N, hpa_file = HPA_FILE_N,
        method = "specificity", verbose = FALSE
    )
    if (!"t" %in% names(hla_dep))
        stop("hla_dep missing `t` column required for fgsea")
    rank_stat_n <- build_rank_stat(hla_dep, gene_col = "protein", stat_col = "t")
    gsea_n <- run_celltype_gsea(rank_stat_n, marker_sets_n,
                                  min_size = 5, max_size = 500, seed = 42)
    if (is.null(gsea_n) || nrow(gsea_n) == 0L)
        stop("No GSEA results for HLA main-effect contrast")
    gsea_n[, analysis := "hla_main"]

    pN_supp <- make_celltype_gsea_plot(
        gsea_n,
        contrast_col    = "analysis",
        contrast_order  = "hla_main",
        contrast_labels = c(hla_main = "HLA+ vs HLA-"),
        title_str       = "n  Cell-type enrichment \u2014 HLA-main-effect DEPs",
        subtitle_str    = paste0(
            "fgsea on limma moderated-t (HLA+ vs HLA-) | ",
            "NES>0: markers higher in HLA+ | NES<0: lower | * FDR<0.05"
        ),
        high_col        = "#E6A817",
        low_col         = "#56B4E9"
    )
    if (is.null(pN_supp)) stop("make_celltype_gsea_plot returned NULL")
    save_supp(pN_supp, "n_hla_celltype_gsea", 4.5, 6.5)
}, error = function(e) cat("  Supp n skipped:", conditionMessage(e), "\n"))

cat("\nFigure 3 complete. Panels saved to:", FIG_DIR, "\n")
cat("Supplement panels saved to:", SUPP_DIR, "\n")
