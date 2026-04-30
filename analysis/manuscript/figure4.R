#!/usr/bin/env Rscript
# figure4.R — Figure 4: Sex differences in the MS proteomic signature
#
# Panels:
#   a  Female vs male logFC scatter — full proteome, labelled by category + known markers
#   b  DEP landscape: female-specific / concordant / discordant / male-specific
#   c  Female-MS DEPs: Walchli+HPA cell-type expression heatmap
#   d  GO:BP — female DEPs (up + down, direction-faceted)
#   e  Male-MS DEPs: Walchli+HPA cell-type expression heatmap
#   f  GO:BP — male DEPs (up + down, direction-faceted)
#   g  Top female-specific marker trajectory by sex (excl. GFAP/NEFL/MOG already in Fig 2)
#   h  Top male-specific marker trajectory by sex
#   i  Pre-onset sex differences: female vs male logFC scatter (coloured by full-cohort category)
#   j  Sex-stratified ML ROC curves: female/male × combined/pre-onset/post-onset (6 curves)
# Supp:
#   a  Sex × MS interaction volcano
#   b  All sex DEPs: cell-type heatmap (top 50 by |logFC|)
#   c  Cell-type enrichment (fgsea NES): female MS vs HC
#   d  Cell-type enrichment (fgsea NES): male MS vs HC
#   e  SHAP paired importance: female vs male combined classifier (top 35)
# Extended Data Fig. 4:
#   a  Pre-onset: female vs male logFC concordance
#   b  Post-onset: female vs male logFC concordance
#
# Output: results/figures/4/

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(splines)
    library(glue)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "celltype_overrep_plot.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))

cfg          <- load_disease_config()
COHORT       <- cfg$cohort_short
DISEASE_CAPS <- cfg$disease_short_caps
STATUS_COL   <- cfg$cohort_status_col
SV           <- cfg$status_values
PRE_LBL      <- glue("pre-{DISEASE_CAPS}")
POST_LBL     <- glue("post-{DISEASE_CAPS}")

FIG_DIR  <- file.path(PROJ_DIR, "results", "figures", "4")
SUPP_DIR <- file.path(PROJ_DIR, "results", "figures", "4_supp")
SEX_DIR  <- file.path(PROJ_DIR, "results", "sex_stratified")
DATA_DIR <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed")
NET_DIR  <- file.path(PROJ_DIR, "results", "networks")
WH_FILE  <- file.path(NET_DIR, glue("{COHORT}_combined_walchli_hpa_matrix.csv"))
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

# Compartment constants (consistent with figure 1)
COMP_FILL  <- c("CNS" = "#27AE60", "Vascular" = "#E67E22", "Immune" = "#2980B9",
                "Unknown" = "grey78")
COMP_ORDER <- c("CNS", "Vascular", "Immune", "Unknown")
COMP_ALPHA <- 0.11

build_comp_map <- function(wh_file) {
    if (!file.exists(wh_file))
        return(data.table(gene = character(0), compartment = character(0)))
    wh <- fread(wh_file)
    vcts <- intersect(c("Endothelial","Pericytes","Smooth muscle","Fibroblasts"), names(wh))
    ccts <- intersect(c("Oligodendrocytes","Microglia","Astrocytes","Neurons"),   names(wh))
    icts <- intersect(c("T cells","B cells","NK cells","Monocytes","Dendritic cells",
                         "Macrophages","Plasma cells","Granulocytes"),              names(wh))
    wh[, vm := if (length(vcts) > 0) rowMeans(.SD, na.rm=TRUE) else 0, .SDcols = vcts]
    wh[, cm := if (length(ccts) > 0) rowMeans(.SD, na.rm=TRUE) else 0, .SDcols = ccts]
    wh[, im := if (length(icts) > 0) rowMeans(.SD, na.rm=TRUE) else 0, .SDcols = icts]
    wh[, compartment := fcase(
        cm >= vm & cm >= im, "CNS",
        vm >= cm & vm >= im, "Vascular",
        default              = "Immune"
    )]
    wh[, .(gene = tolower(gene), compartment)]
}
comp_map <- build_comp_map(WH_FILE)

save_panel <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(FIG_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: panel_", name, ".pdf\n", sep = "")
}
save_supp <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(SUPP_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: 4_supp/panel_", name, ".pdf\n", sep = "")
}

# ---------------------------------------------------------------------------
# Load sex-stratified results
# ---------------------------------------------------------------------------
cat("Loading sex-stratified DEP results...\n")

fa_file <- file.path(SEX_DIR, glue("{COHORT}_female_vs_hc_all.csv"))
ma_file <- file.path(SEX_DIR, glue("{COHORT}_male_vs_hc_all.csv"))
si_file <- file.path(SEX_DIR, glue("{COHORT}_sex_interaction.csv"))

if (!file.exists(fa_file) || !file.exists(ma_file)) {
    stop("Sex-stratified results not found in ", SEX_DIR)
}

fa <- fread(fa_file)
ma <- fread(ma_file)
si <- if (file.exists(si_file)) {
    s <- fread(si_file); s[, fdr := p.adjust(P.Value, method = "BH")]; s
} else NULL

cat(sprintf("  Female all: %d  Male all: %d  Interaction: %d\n",
            nrow(fa), nrow(ma), if (!is.null(si)) nrow(si) else 0))

# Merge into wide table for concordance analysis
wide <- merge(
    fa[, .(protein, logFC_f = logFC, fdr_f = adj.P.Val, t_f = t)],
    ma[, .(protein, logFC_m = logFC, fdr_m = adj.P.Val, t_m = t)],
    by = "protein"
)
wide[, dep_cat := fcase(
    fdr_f < 0.05 & fdr_m < 0.05 & sign(logFC_f) == sign(logFC_m), "Both concordant",
    fdr_f < 0.05 & fdr_m < 0.05 & sign(logFC_f) != sign(logFC_m), "Both discordant",
    fdr_f < 0.05 & fdr_m >= 0.05,                                  "Female-specific",
    fdr_f >= 0.05 & fdr_m < 0.05,                                  "Male-specific",
    default                                                         = "NS"
)]
wide[, dep_cat := factor(dep_cat,
    levels = c("Female-specific", "Both concordant", "Both discordant",
               "Male-specific", "NS"))]

n_cat <- wide[dep_cat != "NS", .N, by = dep_cat][order(-N)]
cat("DEP categories:\n"); print(n_cat, row.names = FALSE)

n_f_dep <- sum(wide$fdr_f < 0.05)
n_m_dep <- sum(wide$fdr_m < 0.05)

CAT_COLS <- c(
    "Female-specific" = "#CC0066",
    "Both concordant" = "#7B4F8E",
    "Both discordant" = "#E6A817",
    "Male-specific"   = "#2B4C7E",
    "NS"              = "grey80"
)
CAT_SZ <- c("Female-specific"=1.6,"Both concordant"=2.2,"Both discordant"=2.2,
            "Male-specific"=1.6,"NS"=0.35)

# ---------------------------------------------------------------------------
# Panel a — Full proteome scatter: female logFC vs male logFC
#   Leading panel. Labels: known MS markers + all discordant + top 5 each sex.
# ---------------------------------------------------------------------------
cat("Building panel a (F vs M logFC scatter)...\n")

# Programmatic label strategy (no by-name marker injection):
#   1. All sex-discordant proteins (both FDR<0.05, opposite direction)
#   2. Top 5 female-specific by |logFC_f|
#   3. Top 5 male-specific by |logFC_m|
disc_prot    <- wide[dep_cat == "Both discordant", protein]
f_top5       <- wide[dep_cat == "Female-specific"][order(-abs(logFC_f))][
    seq_len(min(5,.N)), protein]
m_top5       <- wide[dep_cat == "Male-specific"][order(-abs(logFC_m))][
    seq_len(min(5,.N)), protein]

lab_prots <- unique(c(disc_prot, f_top5, m_top5))
lab_a <- wide[protein %in% lab_prots, .(protein = toupper(protein), logFC_f, logFC_m, dep_cat)]
cat(sprintf("  Labelling %d proteins in scatter\n", nrow(lab_a)))

r_all  <- round(cor(wide$logFC_f, wide$logFC_m, use = "complete.obs"), 3)
r_deps <- round(cor(wide[fdr_f < 0.05 | fdr_m < 0.05, logFC_f],
                    wide[fdr_f < 0.05 | fdr_m < 0.05, logFC_m],
                    use = "complete.obs"), 3)

pA <- ggplot(wide, aes(x = logFC_m, y = logFC_f,
                        colour = dep_cat, size = dep_cat)) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey55") +
    geom_point(data = wide[dep_cat == "NS"], alpha = 0.22) +
    geom_point(data = wide[dep_cat != "NS"], alpha = 0.90) +
    geom_text_repel(
        data         = lab_a,
        aes(label    = protein),
        size         = 2.2, fontface = "bold",
        max.overlaps = Inf, segment.size = 0.22,
        segment.color = "grey50", box.padding = 0.35,
        show.legend  = FALSE,
        min.segment.length = 0.2
    ) +
    annotate("text",
             x = max(wide$logFC_m, na.rm=TRUE) * 0.92,
             y = min(wide$logFC_f, na.rm=TRUE) * 0.85,
             label = sprintf("r (all) = %.3f\nr (DEPs) = %.3f", r_all, r_deps),
             size = 2.3, hjust = 1, colour = "grey30") +
    scale_colour_manual(values = CAT_COLS, name = NULL,
                        guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    scale_size_manual(values = CAT_SZ, guide = "none") +
    labs(
        title    = glue("a  Female vs male {DISEASE_CAPS} proteome: concordance and sex-specificity"),
        subtitle = sprintf(
            "Female FDR<0.05: %d | Male FDR<0.05: %d | Discordant (orange): FDR<0.05 in both, opposite direction",
            n_f_dep, n_m_dep),
        x = bquote(log[2]~"FC (Male"~.(DISEASE_CAPS)~"vs HC)"),
        y = bquote(log[2]~"FC (Female"~.(DISEASE_CAPS)~"vs HC)")
    ) +
    theme_ukb(base_size = 9) +
    theme(legend.position = "bottom",
          plot.subtitle   = element_text(size = 7, colour = "grey40"))

save_panel(pA, "a_sex_logfc_scatter", 4.5, 4.5)

# ---------------------------------------------------------------------------
# Panel b — DEP landscape count bar
# ---------------------------------------------------------------------------
cat("Building panel b (DEP landscape)...\n")

bar_dat <- n_cat[order(N)]
bar_dat[, dep_cat := factor(dep_cat, levels = bar_dat$dep_cat)]

pB <- ggplot(bar_dat, aes(y = dep_cat, x = N, colour = dep_cat)) +
    geom_segment(aes(xend = 0, yend = dep_cat), linewidth = 1.2) +
    geom_point(size = 5.5) +
    scale_colour_manual(values = CAT_COLS, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)), limits = c(0, NA)) +
    labs(
        title    = glue("b  Sex-stratified {DISEASE_CAPS} DEP landscape"),
        subtitle = sprintf("Female FDR<0.05: %d | Male FDR<0.05: %d | Discordant = opposite direction in each sex",
                           n_f_dep, n_m_dep),
        x = "Number of proteins", y = NULL
    ) +
    theme_ukb(base_size = 9) +
    theme(
        axis.text.y   = element_text(size = 8),
        plot.subtitle = element_text(size = 7, colour = "grey40")
    )

save_panel(pB, "b_dep_landscape", 4.0, 2.8)

# ---------------------------------------------------------------------------
# Panel i — Pre-onset sex differences: female vs male logFC scatter
#   Asks: is sex dimorphism already present before diagnosis?
#   Points coloured by full-cohort dep_cat from panels a/b.
# ---------------------------------------------------------------------------
cat("Building panel i (pre-onset sex differences scatter)...\n")

PRE_F_FILE <- file.path(SEX_DIR, glue("{COHORT}_pre_female_vs_hc_all.csv"))
PRE_M_FILE <- file.path(SEX_DIR, glue("{COHORT}_pre_male_vs_hc_all.csv"))

tryCatch({
    if (!file.exists(PRE_F_FILE) || !file.exists(PRE_M_FILE))
        stop("Pre-onset sex-stratified results not found — run 01_sex_stratified_ms_deps.R first")

    pre_f <- fread(PRE_F_FILE)
    pre_m <- fread(PRE_M_FILE)

    pre_wide <- merge(
        pre_f[, .(protein, pre_logFC_f = logFC, pre_fdr_f = adj.P.Val)],
        pre_m[, .(protein, pre_logFC_m = logFC, pre_fdr_m = adj.P.Val)],
        by = "protein"
    )
    pre_wide <- merge(pre_wide, wide[, .(protein, dep_cat)], by = "protein", all.x = TRUE)
    pre_wide[is.na(dep_cat), dep_cat := "NS"]
    pre_wide[, dep_cat := factor(dep_cat,
        levels = c("Female-specific","Both concordant","Both discordant","Male-specific","NS"))]

    # Correlation annotation
    r_pre <- round(cor(pre_wide$pre_logFC_f, pre_wide$pre_logFC_m, use = "complete.obs"), 3)
    r_pre_deps <- {
        sub <- pre_wide[pre_fdr_f < 0.05 | pre_fdr_m < 0.05]
        if (nrow(sub) > 5) round(cor(sub$pre_logFC_f, sub$pre_logFC_m, use = "complete.obs"), 3)
        else NA
    }

    # Label strategy: programmatic top-5 sex-specific by |pre_logFC|
    f_top_c <- pre_wide[pre_fdr_f < 0.05][order(-abs(pre_logFC_f))][
                   seq_len(min(5, .N)), protein]
    m_top_c <- pre_wide[pre_fdr_m < 0.05][order(-abs(pre_logFC_m))][
                   seq_len(min(5, .N)), protein]
    lab_c   <- unique(c(f_top_c, m_top_c))
    lab_c_dt <- pre_wide[protein %in% lab_c,
                          .(protein = toupper(protein), pre_logFC_f, pre_logFC_m, dep_cat)]

    r_label <- if (!is.na(r_pre_deps))
        sprintf("r (all) = %.3f\nr (pre-onset DEPs) = %.3f", r_pre, r_pre_deps)
    else
        sprintf("r (all) = %.3f", r_pre)

    n_pre_f <- sum(pre_wide$pre_fdr_f < 0.05, na.rm = TRUE)
    n_pre_m <- sum(pre_wide$pre_fdr_m < 0.05, na.rm = TRUE)

    pC <- ggplot(pre_wide, aes(x = pre_logFC_m, y = pre_logFC_f,
                                colour = dep_cat, size = dep_cat)) +
        geom_hline(yintercept = 0, linewidth = 0.22, colour = "grey70") +
        geom_vline(xintercept = 0, linewidth = 0.22, colour = "grey70") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    linewidth = 0.3, colour = "grey55") +
        geom_point(data = pre_wide[dep_cat == "NS"], alpha = 0.22) +
        geom_point(data = pre_wide[dep_cat != "NS"], alpha = 0.90) +
        geom_text_repel(
            data          = lab_c_dt,
            aes(label     = protein),
            size          = 2.2, fontface = "bold",
            max.overlaps  = Inf, segment.size = 0.22,
            segment.color = "grey50", box.padding = 0.35,
            show.legend   = FALSE, min.segment.length = 0.2
        ) +
        annotate("text",
                 x = max(pre_wide$pre_logFC_m, na.rm=TRUE) * 0.92,
                 y = min(pre_wide$pre_logFC_f, na.rm=TRUE) * 0.85,
                 label = r_label,
                 size = 2.3, hjust = 1, colour = "grey30") +
        scale_colour_manual(values = CAT_COLS, name = NULL,
                            guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
        scale_size_manual(values = CAT_SZ, guide = "none") +
        labs(
            title    = glue("i  Pre-onset sex differences in the {DISEASE_CAPS} proteome"),
            subtitle = sprintf(
                "Pre-onset: Female FDR<0.05: %d | Male FDR<0.05: %d | Coloured by full-cohort category",
                n_pre_f, n_pre_m),
            x = bquote(log[2]~"FC (Pre-onset Male"~.(DISEASE_CAPS)~"vs HC)"),
            y = bquote(log[2]~"FC (Pre-onset Female"~.(DISEASE_CAPS)~"vs HC)")
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom",
              plot.subtitle   = element_text(size = 7, colour = "grey40"))

    save_panel(pC, "i_preonset_sex_scatter", 4.5, 4.5)

}, error = function(e) cat("  Panel i skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Panels c & d — Cell-type expression heatmaps (female and male DEPs)
# ---------------------------------------------------------------------------
cat("Building panels c and d (cell-type heatmaps)...\n")

tryCatch({
    suppressPackageStartupMessages({ library(pheatmap) })
    if (!file.exists(WH_FILE)) stop("Walchli+HPA matrix not found")
    wh_c  <- fread(WH_FILE)
    ct_cols_hm <- setdiff(names(wh_c), c("gene","logFC","direction","abs_logFC",
                                           "vm","cm","im","compartment"))

    # Helper: build pheatmap for a set of DEPs
    make_dep_heatmap <- function(dep_dt, logfc_col, fdr_col, title_str, pan_name,
                                  up_lbl, dn_lbl, up_col, dn_col,
                                  heat_col = "#CC0066", w = 11, h = 5.5) {
        # Prefer FDR<0.05 genes in Walchli matrix; fall back to top by p-value
        sig <- dep_dt[get(fdr_col) < 0.05][toupper(protein) %in% toupper(wh_c$gene)]
        if (nrow(sig) < 5)
            sig <- dep_dt[toupper(protein) %in% toupper(wh_c$gene)][order(get(fdr_col))][
                seq_len(min(40, .N))]
        sig <- sig[seq_len(min(40, .N))]

        wh_sub <- wh_c[toupper(gene) %in% toupper(sig$protein)]
        if (nrow(wh_sub) < 3) stop(sprintf("Too few DEPs in Walchli matrix for %s", pan_name))

        mat <- as.matrix(wh_sub[, ..ct_cols_hm])
        rownames(mat) <- toupper(wh_sub$gene)
        mat[is.na(mat)] <- 0

        ct_ord <- c("Endothelial","Pericytes","Smooth muscle","Fibroblasts",
                    "Oligodendrocytes","Microglia","Astrocytes","Neurons",
                    "T cells","B cells","NK cells","Monocytes",
                    "Dendritic cells","Macrophages","Plasma cells","Granulocytes")
        ct_pres <- intersect(ct_ord, colnames(mat))
        mat_d   <- t(mat[, ct_pres, drop = FALSE])

        ct_grp <- ifelse(ct_pres %in% c("Endothelial","Pericytes","Smooth muscle","Fibroblasts"),
                         "Vascular",
                  ifelse(ct_pres %in% c("Oligodendrocytes","Microglia","Astrocytes","Neurons"),
                         "CNS", "Immune"))
        ann_row <- data.frame(Compartment = ct_grp, row.names = ct_pres)

        dir_vec <- sig[match(colnames(mat_d), toupper(protein)),
                       ifelse(get(logfc_col) > 0, up_lbl, dn_lbl)]
        ann_col  <- data.frame(Direction = dir_vec, row.names = colnames(mat_d))

        gaps <- c(sum(ct_grp == "Vascular"),
                  sum(ct_grp == "Vascular") + sum(ct_grp == "CNS"))
        gaps <- gaps[gaps > 0 & gaps < nrow(mat_d)]

        ann_colors <- list(
            Direction   = setNames(c(up_col, dn_col), c(up_lbl, dn_lbl)),
            Compartment = c("Vascular" = "#E67E22","CNS" = "#27AE60","Immune" = "#2980B9")
        )

        out <- file.path(FIG_DIR, paste0("panel_", pan_name, ".pdf"))
        pdf(out, width = w, height = h, onefile = FALSE)
        pheatmap(mat_d,
                 color             = colorRampPalette(c("white", heat_col))(100),
                 breaks            = seq(0, 1, length.out = 101),
                 cluster_rows      = FALSE, cluster_cols = TRUE,
                 fontsize_row      = 8, fontsize_col = 6,
                 annotation_col    = ann_col, annotation_row = ann_row,
                 annotation_colors = ann_colors,
                 annotation_names_col = FALSE, annotation_names_row = FALSE,
                 border_color      = NA,
                 gaps_row          = if (length(gaps) > 0) gaps else NULL,
                 main              = title_str)
        dev.off()
        cat("  Saved: panel_", pan_name, ".pdf\n", sep = "")
    }

    make_dep_heatmap(wide, "logFC_f", "fdr_f",
                     title_str = glue("c  Female-{DISEASE_CAPS} DEPs: Walchli+HPA cell-type expression"),
                     pan_name  = "c_female_celltype_heatmap",
                     up_lbl = glue("Up in F-{DISEASE_CAPS}"), dn_lbl = glue("Down in F-{DISEASE_CAPS}"),
                     up_col = "#CC0066", dn_col = "#2B4C7E", heat_col = "#CC0066")

    make_dep_heatmap(wide[fdr_m < 0.05], "logFC_m", "fdr_m",
                     title_str = glue("e  Male-{DISEASE_CAPS} DEPs: Walchli+HPA cell-type expression"),
                     pan_name  = "e_male_celltype_heatmap",
                     up_lbl = glue("Up in M-{DISEASE_CAPS}"), dn_lbl = glue("Down in M-{DISEASE_CAPS}"),
                     up_col = "#2B4C7E", dn_col = "#56B4E9", heat_col = "#2B4C7E")

}, error = function(e) cat("  Panels c/d skipped:", conditionMessage(e), "\n"))


# ---------------------------------------------------------------------------
# Panels e & f — GO enrichment (clusterProfiler)
# ---------------------------------------------------------------------------
cat("Building panels e and f (GO enrichment)...\n")

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

    run_ego <- function(genes, label, n_top = 8) {
        if (length(genes) < 5) return(NULL)
        ego <- tryCatch(
            enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 1.0, qvalueCutoff = 1.0,
                     readable = TRUE),
            error = function(e) NULL
        )
        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
        dt <- as.data.table(as.data.frame(ego))
        dt[, FoldEnrichment := calc_fe(GeneRatio, BgRatio)]
        dt[, direction := label]
        # Prefer FDR-significant; fall back to top N by nominal p-value
        sig <- dt[p.adjust < 0.05][order(p.adjust)][seq_len(min(n_top, .N))]
        if (nrow(sig) > 0) return(sig)
        dt[order(pvalue)][seq_len(min(5, .N))]
    }

    make_go_panel <- function(up_genes, down_genes, title_str, up_lbl, dn_lbl,
                              up_col, dn_col, fname) {
        cat(sprintf("  %s GO: %d up, %d down\n", title_str, length(up_genes), length(down_genes)))
        go_rows <- rbindlist(Filter(Negate(is.null), list(
            run_ego(up_genes, up_lbl),
            run_ego(down_genes, dn_lbl)
        )), fill = TRUE)
        if (nrow(go_rows) == 0) { cat("  No significant GO terms\n"); return(invisible(NULL)) }
        go_rows[, Description := stringr::str_wrap(Description, width = 42)]
        go_rows[, direction := factor(direction, levels = c(up_lbl, dn_lbl))]
        go_rows[, Description := factor(Description,
            levels = unique(go_rows[order(direction, FoldEnrichment), Description]))]
        go_cols <- setNames(c(up_col, dn_col), c(up_lbl, dn_lbl))
        p <- ggplot(go_rows, aes(x = FoldEnrichment, y = Description,
                                  size = Count, colour = direction)) +
            geom_point() +
            facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
            scale_colour_manual(values = go_cols,
                                 guide = guide_legend(override.aes = list(size = 3.5))) +
            scale_size_continuous(name = "Gene count", range = c(2, 6),
                                   guide = guide_legend(override.aes = list(colour = "grey40"))) +
            scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
            labs(
                title    = sprintf("%s  GO:BP \u2014 %s", substr(fname,1,1), title_str),
                subtitle = sprintf("%d up, %d down (FDR<0.05 in this sex)",
                                   length(up_genes), length(down_genes)),
                x = "Fold enrichment", y = NULL
            ) +
            theme_ukb(base_size = 9) +
            theme(axis.text.y    = element_text(size = 7, lineheight = 0.88),
                  strip.text     = element_text(size = 8, face = "bold"),
                  legend.position = "right",
                  plot.subtitle  = element_text(size = 7, colour = "grey40"))
        save_panel(p, fname, 6.5, 6.5)
    }

    # Panel d: female DEPs
    make_go_panel(
        up_genes  = wide[fdr_f < 0.05 & logFC_f > 0, toupper(protein)],
        down_genes = wide[fdr_f < 0.05 & logFC_f < 0, toupper(protein)],
        title_str = glue("female-{DISEASE_CAPS} DEPs"),
        up_lbl = glue("Up in F-{DISEASE_CAPS}"), dn_lbl = glue("Down in F-{DISEASE_CAPS}"),
        up_col = "#CC0066", dn_col = "#2B4C7E",
        fname = "d_female_go"
    )

    # Panel f: male DEPs
    make_go_panel(
        up_genes  = wide[fdr_m < 0.05 & logFC_m > 0, toupper(protein)],
        down_genes = wide[fdr_m < 0.05 & logFC_m < 0, toupper(protein)],
        title_str = glue("male-{DISEASE_CAPS} DEPs"),
        up_lbl = glue("Up in M-{DISEASE_CAPS}"), dn_lbl = glue("Down in M-{DISEASE_CAPS}"),
        up_col = "#2B4C7E", dn_col = "#56B4E9",
        fname = "f_male_go"
    )

}, error = function(e) {
    cat("  GO (clusterProfiler):", conditionMessage(e), "\n")
    for (pan in c("e_female_go", "f_male_go")) {
        ph <- ggplot() +
            annotate("text", x=0.5, y=0.5, label="Requires clusterProfiler",
                     size=4, hjust=0.5, colour="grey40") +
            xlim(0,1) + ylim(0,1) + theme_void()
        save_panel(ph, pan, 5.5, 6)
    }
})

# ---------------------------------------------------------------------------
# Panels g & h — Sex-stratified disease-course trajectories
#   g: top female-specific DEP (auto-selected by |logFC_f|, FDR<0.05 female)
#   h: top male-specific DEP (auto-selected by |logFC_m|, FDR<0.05 male)
#   Side-by-side facets (Female | Male) make the sex-specificity vs. shared
#   biology immediately legible.
# ---------------------------------------------------------------------------
cat("Building panels g and h (sex-stratified trajectories)...\n")

QC_FILE       <- file.path(DATA_DIR, glue("{COHORT}_olink_qc.csv"))
SEX_TRAJ_COLS <- c("Female" = "#CC0066", "Male" = "#2B4C7E")
YTD_RANGE     <- c(-15, 15)

make_sex_traj <- function(qc, prot_lc, pan_letter, dep_note, save_name,
                           hc_ref_dt = NULL,
                           save_w = 5.5, save_h = 3.5) {
    prot_uc <- toupper(prot_lc)
    if (!prot_lc %in% names(qc)) {
        cat(sprintf("  %s not in QC file; skipping panel %s\n", prot_uc, pan_letter))
        return(invisible(NULL))
    }

    traj_d <- qc[
        ms_status %in% c(SV$pre_onset, SV$post_onset) &
        qc_outlier == FALSE &
        !is.na(get(prot_lc)) & !is.na(years_to_diagnosis) & !is.na(sex)
    ]
    traj_d[, ytd_display := -years_to_diagnosis]
    traj_d <- traj_d[ytd_display >= YTD_RANGE[1] & ytd_display <= YTD_RANGE[2]]
    traj_d[, sex_label := factor(fifelse(sex == 0, "Female", "Male"),
                                  levels = c("Female", "Male"))]

    x_grid   <- seq(YTD_RANGE[1], YTD_RANGE[2], length.out = 200)
    n_sex    <- integer(2); names(n_sex) <- c("Female", "Male")
    ms_rows  <- list()
    hc_rows  <- list()

    for (s_lbl in c("Female", "Male")) {
        s_val <- if (s_lbl == "Female") 0L else 1L

        # MS spline
        d_ms <- traj_d[sex == s_val, .(npx = get(prot_lc), ytd_display, age_at_sampling)]
        d_ms <- d_ms[!is.na(npx)]
        n_sex[s_lbl] <- nrow(d_ms)
        if (nrow(d_ms) >= 15) {
            fit_ms <- tryCatch(
                lm(npx ~ ns(ytd_display, df = 4) + age_at_sampling, data = d_ms),
                error = function(e) NULL
            )
            if (!is.null(fit_ms)) {
                nd_ms <- data.table(ytd_display = x_grid,
                                    age_at_sampling = median(d_ms$age_at_sampling, na.rm=TRUE))
                pr_ms <- tryCatch(predict(fit_ms, newdata=nd_ms, se.fit=TRUE), error=function(e) NULL)
                if (!is.null(pr_ms))
                    ms_rows[[s_lbl]] <- data.table(
                        sex_label   = s_lbl, ytd_display = x_grid,
                        pred_npx    = pr_ms$fit,
                        pred_lci    = pr_ms$fit - 1.96 * pr_ms$se.fit,
                        pred_uci    = pr_ms$fit + 1.96 * pr_ms$se.fit)
            }
        }

        # HC reference spline (PSM-matched, per sex)
        if (!is.null(hc_ref_dt) && prot_lc %in% names(hc_ref_dt)) {
            d_hc <- hc_ref_dt[sex == s_val & !is.na(get(prot_lc)) & !is.na(ytd_display),
                               .(npx = get(prot_lc), ytd_display, age_at_sampling)]
            d_hc <- d_hc[ytd_display >= YTD_RANGE[1] & ytd_display <= YTD_RANGE[2]]
            if (nrow(d_hc) >= 10) {
                fit_hc <- tryCatch(
                    lm(npx ~ ns(ytd_display, df = 3) + age_at_sampling, data = d_hc),
                    error = function(e) NULL
                )
                if (!is.null(fit_hc)) {
                    nd_hc <- data.table(ytd_display = x_grid,
                                        age_at_sampling = median(d_hc$age_at_sampling, na.rm=TRUE))
                    pr_hc <- tryCatch(predict(fit_hc, newdata=nd_hc, se.fit=TRUE), error=function(e) NULL)
                    if (!is.null(pr_hc))
                        hc_rows[[s_lbl]] <- data.table(
                            sex_label   = s_lbl, ytd_display = x_grid,
                            hc_pred     = pr_hc$fit,
                            hc_lci      = pr_hc$fit - 1.96 * pr_hc$se.fit,
                            hc_uci      = pr_hc$fit + 1.96 * pr_hc$se.fit)
                }
            }
        }
    }

    if (length(ms_rows) == 0) {
        cat(sprintf("  Insufficient data for %s; skipping panel %s\n", prot_uc, pan_letter))
        return(invisible(NULL))
    }

    traj    <- rbindlist(ms_rows)
    traj[,  sex_label := factor(sex_label, levels = c("Female", "Male"))]
    hc_traj <- if (length(hc_rows) > 0) {
        ht <- rbindlist(hc_rows)
        ht[, sex_label := factor(sex_label, levels = c("Female", "Male"))]
        ht
    } else NULL

    n_lbl <- data.table(sex_label = factor(c("Female","Male"), levels=c("Female","Male")),
                         n = c(n_sex["Female"], n_sex["Male"]))

    prot_row <- wide[tolower(protein) == prot_lc]
    lfc_note <- if (nrow(prot_row) > 0)
        sprintf("F logFC=%.2f (FDR=%s)  |  M logFC=%.2f (FDR=%s)",
                prot_row$logFC_f, signif(prot_row$fdr_f, 2),
                prot_row$logFC_m, signif(prot_row$fdr_m, 2))
    else ""
    hc_note <- if (!is.null(hc_traj)) " | grey dashed = HC (PSM-matched)" else ""

    pP <- ggplot() +
        # HC ribbon + dashed line (rendered behind MS)
        {if (!is.null(hc_traj)) list(
            geom_ribbon(data = hc_traj,
                        aes(x = ytd_display, ymin = hc_lci, ymax = hc_uci,
                            group = sex_label),
                        fill = "grey72", alpha = 0.35, inherit.aes = FALSE),
            geom_line(data = hc_traj,
                      aes(x = ytd_display, y = hc_pred, group = sex_label),
                      colour = "grey45", linewidth = 0.6, linetype = "dashed",
                      inherit.aes = FALSE)
        )} +
        geom_point(data = traj_d,
                   aes(x = ytd_display, y = get(prot_lc), colour = sex_label),
                   size = 0.4, alpha = 0.10) +
        geom_ribbon(data = traj,
                    aes(x = ytd_display, ymin = pred_lci, ymax = pred_uci,
                        fill = sex_label),
                    alpha = 0.20) +
        geom_line(data = traj,
                  aes(x = ytd_display, y = pred_npx, colour = sex_label),
                  linewidth = 0.85) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.4, colour = "grey40") +
        geom_text(data = n_lbl,
                  aes(label = paste0("n=", n)),
                  x = YTD_RANGE[1] + 0.5, y = -Inf, hjust = 0, vjust = -0.5,
                  size = 1.8, colour = "grey45") +
        annotate("text", x = 0.5, y = Inf, label = "Dx",
                 hjust = 0, vjust = 1.5, size = 2.2, colour = "grey40") +
        facet_wrap(~sex_label, nrow = 1) +
        scale_colour_manual(values = SEX_TRAJ_COLS, name = NULL, guide = "none") +
        scale_fill_manual(values   = SEX_TRAJ_COLS, name = NULL, guide = "none") +
        scale_x_continuous(breaks = seq(-12, 12, by = 4)) +
        labs(
            title    = sprintf("%s  %s disease-course trajectory by sex", pan_letter, prot_uc),
            subtitle = paste0(dep_note, "\n", lfc_note, hc_note),
            x = glue("Years relative to {DISEASE_CAPS} diagnosis"),
            y = sprintf("%s (NPX)", prot_uc)
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position  = "none",
              strip.text       = element_text(size = 9, face = "bold"),
              plot.subtitle    = element_text(size = 6.0, colour = "grey45"))

    save_panel(pP, save_name, save_w, save_h)
}

if (file.exists(QC_FILE)) {
    qc <- fread(QC_FILE)
    if (STATUS_COL != "ms_status" && STATUS_COL %in% names(qc))
        setnames(qc, STATUS_COL, "ms_status")

    # PSM-matched HC reference: match MS cases to HC on age + sex (1:3),
    # propagate each MS case's years_to_diagnosis as pseudo-ytd to its matched HCs
    suppressPackageStartupMessages(library(MatchIt))
    hc_ref_traj <- tryCatch({
        ms_pool_t <- qc[ms_status %in% c(SV$pre_onset, SV$post_onset) & qc_outlier==FALSE &
                        !is.na(age_at_sampling) & !is.na(sex) & !is.na(years_to_diagnosis)]
        hc_pool_t <- qc[ms_status == SV$control & qc_outlier==FALSE &
                        !is.na(age_at_sampling) & !is.na(sex)]
        pool_t <- rbind(
            ms_pool_t[, .(eid, is_ms=1L, age_at_sampling, sex, years_to_diagnosis)],
            hc_pool_t[, .(eid, is_ms=0L, age_at_sampling, sex, years_to_diagnosis=NA_real_)]
        )
        set.seed(42)
        m_obj <- matchit(is_ms ~ age_at_sampling + sex, data=as.data.frame(pool_t),
                         method="nearest", ratio=3, replace=FALSE)
        md     <- as.data.table(match.data(m_obj))
        ms_sub <- md[is_ms==1, .(subclass, ytd_ms=years_to_diagnosis)]
        hc_md  <- merge(md[is_ms==0, .(eid, subclass)], ms_sub, by="subclass", all.x=TRUE)
        hc_rows_t <- merge(hc_md[, .(eid, pseudo_ytd=ytd_ms)], qc, by="eid")
        hc_rows_t[, ytd_display := -pseudo_ytd]
        cat(sprintf("  HC PSM reference: %d observations\n", nrow(hc_rows_t)))
        hc_rows_t
    }, error=function(e) { cat("  HC PSM failed:", e$message, "\n"); NULL })

    # Panel h: top female-specific DEP by |logFC_f|
    # Exclude proteins already featuring as trajectory panels in figures 2/3.
    EXCL_TRAJ <- c("gfap", "nefl", "mog", "lilrb4")
    h_prot_f <- tolower(
        wide[dep_cat == "Female-specific" & !tolower(protein) %in% EXCL_TRAJ
             ][order(-abs(logFC_f)), protein[1]]
    )
    cat(sprintf("  Panel g protein: %s (top female-specific by |logFC|)\n", toupper(h_prot_f)))
    make_sex_traj(qc, h_prot_f, "g",
                  "Top female-specific DEP by |logFC| (FDR<0.05 female, NS male)",
                  sprintf("g_%s_sex_trajectory", h_prot_f),
                  hc_ref_dt = hc_ref_traj)

    # Panel i: top male-specific DEP by |logFC_m|
    i_prot_m <- tolower(
        wide[dep_cat == "Male-specific"][order(-abs(logFC_m)), protein[1]]
    )
    if (!is.na(i_prot_m) && length(i_prot_m) > 0) {
        cat(sprintf("  Panel h protein: %s (top male-specific by |logFC|)\n", toupper(i_prot_m)))
        make_sex_traj(qc, i_prot_m, "h",
                      "Top male-specific DEP by |logFC| (FDR<0.05 male, NS female)",
                      sprintf("h_%s_sex_trajectory", i_prot_m),
                      hc_ref_dt = hc_ref_traj)
    } else {
        cat("  No male-specific proteins; skipping panel h\n")
    }
} else {
    cat("  QC file not found; skipping panels g and h\n")
}

# ---------------------------------------------------------------------------
# Supp b — Walchli+HPA cell-type expression heatmap: female + male DEPs
# ---------------------------------------------------------------------------
cat("Building supp b (cell-type heatmap for sex DEPs)...\n")
tryCatch({
    suppressPackageStartupMessages({ library(pheatmap); library(RColorBrewer) })

    if (!file.exists(WH_FILE)) stop("Walchli+HPA matrix not found")
    wh_all <- fread(WH_FILE)
    ct_cols <- setdiff(names(wh_all), c("gene","logFC","direction","abs_logFC",
                                         "vm","cm","im","compartment"))

    dep_all <- wide[fdr_f < 0.05 | fdr_m < 0.05]
    dep_all[, cat_col := fcase(
        dep_cat == "Female-specific", "#CC0066",
        dep_cat == "Male-specific",   "#2B4C7E",
        dep_cat == "Both concordant", "#7B4F8E",
        default = "#E6A817"
    )]
    top_hm_genes <- dep_all[order(-pmax(abs(logFC_f), abs(logFC_m)))][
        seq_len(min(50, .N)), toupper(protein)]

    wh_sub <- wh_all[toupper(gene) %in% top_hm_genes]
    if (nrow(wh_sub) < 3) stop("Insufficient proteins in Walchli matrix")

    mat_h <- as.matrix(wh_sub[, ..ct_cols])
    rownames(mat_h) <- toupper(wh_sub$gene)
    mat_h[is.na(mat_h)] <- 0

    ct_order <- c("Endothelial","Pericytes","Smooth muscle","Fibroblasts",
                  "Oligodendrocytes","Microglia","Astrocytes","Neurons",
                  "T cells","B cells","NK cells","Monocytes",
                  "Dendritic cells","Macrophages","Plasma cells","Granulocytes")
    ct_present <- intersect(ct_order, colnames(mat_h))
    mat_disp   <- t(mat_h[, ct_present, drop=FALSE])

    ct_grp <- ifelse(ct_present %in% c("Endothelial","Pericytes","Smooth muscle","Fibroblasts"),
                     "Vascular",
              ifelse(ct_present %in% c("Oligodendrocytes","Microglia","Astrocytes","Neurons"),
                     "CNS", "Immune"))
    ann_row_b <- data.frame(Compartment = ct_grp, row.names = ct_present)
    ann_col_b <- data.frame(
        Category = dep_all[match(colnames(mat_disp), toupper(protein)), dep_cat],
        row.names = colnames(mat_disp)
    )
    gaps_b <- c(sum(ct_grp=="Vascular"), sum(ct_grp=="Vascular")+sum(ct_grp=="CNS"))
    gaps_b <- gaps_b[gaps_b>0 & gaps_b<nrow(mat_disp)]
    ann_cols_b <- list(
        Compartment = c("Vascular"="#E67E22","CNS"="#27AE60","Immune"="#2980B9"),
        Category    = c("Female-specific"="#CC0066","Male-specific"="#2B4C7E",
                        "Both concordant"="#7B4F8E","Both discordant"="#E6A817","NS"="grey80")
    )

    out_b <- file.path(SUPP_DIR, "panel_b_celltype_heatmap.pdf")
    pdf(out_b, width=12, height=5, onefile=FALSE)
    pheatmap(mat_disp,
             color             = colorRampPalette(c("white","#CC0066"))(100),
             breaks            = seq(0, 1, length.out=101),
             cluster_rows      = FALSE,
             cluster_cols      = TRUE,
             fontsize_row      = 8, fontsize_col = 5.5,
             annotation_col    = ann_col_b,
             annotation_row    = ann_row_b,
             annotation_colors = ann_cols_b,
             annotation_names_col = FALSE, annotation_names_row = FALSE,
             border_color      = NA,
             gaps_row          = if (length(gaps_b)>0) gaps_b else NULL,
             main = glue("b  Sex-stratified {DISEASE_CAPS} DEPs: Walchli+HPA cell-type expression"))
    dev.off()
    cat("  Saved: 4_supp/panel_b_celltype_heatmap.pdf\n")
}, error = function(e) cat("  Supp b skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Supp c — Cell-type enrichment: female MS DEPs   [REFACTORED 2025]
# Supp d — Cell-type enrichment: male   MS DEPs   [REFACTORED 2025]
#   Preranked fgsea on per-sex limma moderated-t (fa / ma), against
#   specificity-filtered markers from the shared helper. Replaces the previous
#   Walchli+HPA matrix Fisher test (expr > 0.5 threshold) which measured
#   marker-set size bias rather than cell-type specificity. See helper header
#   for full rationale.
# ---------------------------------------------------------------------------
cat("Building supp c/d (cell-type fgsea for sex DEPs)...\n")
tryCatch({
    WALCHLI_CSV_S <- file.path(PROJ_DIR, "data", "reference",
                                "walchli_avg_expr_by_celltype.csv")
    HPA_FILE_S    <- file.path(PROJ_DIR, "data", "reference",
                                "hpa_rna_single_cell_type.tsv.zip")
    marker_sets_s <- load_celltype_markers(
        walchli_csv = WALCHLI_CSV_S, hpa_file = HPA_FILE_S,
        method = "specificity", verbose = FALSE
    )

    make_sex_gsea_panel <- function(dep_dt, pan_letter, dep_label, fname,
                                      high_col = "#CC0066") {
        if (!"t" %in% names(dep_dt)) {
            cat(sprintf("  %s: no `t` column in sex DEP table\n", dep_label))
            return(invisible(NULL))
        }
        rank_stat <- build_rank_stat(dep_dt, gene_col = "protein", stat_col = "t")
        gs <- run_celltype_gsea(rank_stat, marker_sets_s,
                                  min_size = 5, max_size = 500, seed = 42)
        if (is.null(gs) || nrow(gs) == 0L) {
            cat(sprintf("  %s: no GSEA results\n", dep_label))
            return(invisible(NULL))
        }
        gs[, analysis := dep_label]
        p <- make_celltype_gsea_plot(
            gs,
            contrast_col    = "analysis",
            contrast_order  = dep_label,
            contrast_labels = setNames(dep_label, dep_label),
            title_str       = sprintf("%s  Cell-type enrichment \u2014 %s vs HC",
                                        pan_letter, dep_label),
            subtitle_str    = paste0(
                "fgsea on limma moderated-t (", dep_label, " vs sex-matched HC) | ",
                "NES>0: markers higher in ", dep_label, " | NES<0: lower | * FDR<0.05"
            ),
            high_col        = high_col,
            low_col         = "#56B4E9"
        )
        if (!is.null(p)) save_supp(p, fname, 4.5, 6.5)
    }

    # fa / ma loaded earlier in figure4.R from ms_female_vs_hc_all.csv /
    # ms_male_vs_hc_all.csv; both have limma's `t` column.
    make_sex_gsea_panel(fa, "c", glue("F-{DISEASE_CAPS}"), "c_female_celltype_gsea",
                         high_col = "#CC0066")
    make_sex_gsea_panel(ma, "d", glue("M-{DISEASE_CAPS}"), "d_male_celltype_gsea",
                         high_col = "#2B4C7E")
}, error = function(e) cat("  Supp c/d skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Supp a — Sex interaction volcano (0 FDR-sig, 175 nominal)
# ---------------------------------------------------------------------------
cat("Building supp a (sex interaction volcano)...\n")

if (!is.null(si)) {
    si[, neg_log10_p := -log10(P.Value)]
    si[, tier := fcase(
        fdr   < 0.05 & logFC > 0, "Up FDR",
        P.Value < 0.05 & logFC > 0, "Up nominal",
        fdr   < 0.05 & logFC < 0, "Down FDR",
        P.Value < 0.05 & logFC < 0, "Down nominal",
        default = "NS"
    )]
    si[, tier := factor(tier, levels = c("Up FDR","Up nominal","Down FDR","Down nominal","NS"))]

    SI_COLS <- c("Up FDR"="#CC0066","Up nominal"="#F5A0C5",
                 "Down FDR"="#56B4E9","Down nominal"="#A8D8F0","NS"="grey74")
    SI_SZ   <- c("Up FDR"=2.8,"Up nominal"=1.1,"Down FDR"=2.8,"Down nominal"=1.1,"NS"=0.45)

    nom_si <- si[P.Value < 0.05][order(-abs(logFC))][seq_len(min(12, .N))]
    nom_si[, label := toupper(protein)]

    n_nom_si <- sum(si$P.Value < 0.05)

    pSA <- ggplot(si, aes(x = logFC, y = neg_log10_p)) +
        geom_point(data = si[tier == "NS"],
                   colour = SI_COLS["NS"], size = SI_SZ["NS"], alpha = 0.35) +
        geom_point(data = si[tier %in% c("Up nominal","Down nominal")],
                   aes(colour = tier), size = SI_SZ["Up nominal"], alpha = 0.75) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
        geom_text_repel(
            data = nom_si, aes(label = label),
            size = 2.0, colour = "grey25",
            max.overlaps = 20, segment.size = 0.2,
            segment.colour = "grey55", box.padding = 0.28,
            show.legend = FALSE
        ) +
        annotate("text",
                 x = max(si$logFC, na.rm=TRUE) * 0.95,
                 y = max(si$neg_log10_p, na.rm=TRUE) * 0.04,
                 label = sprintf("%d FDR<0.05 | %d nominal p<0.05",
                                 sum(si$fdr < 0.05), n_nom_si),
                 size = 1.9, hjust = 1, colour = "grey40") +
        scale_colour_manual(values = SI_COLS, name = NULL) +
        labs(
            title    = glue("a  Sex \u00d7 {DISEASE_CAPS} interaction volcano"),
            subtitle = sprintf("Interaction: sex:%s | %d FDR<0.05 | %d nominal p<0.05",
                               STATUS_COL, sum(si$fdr < 0.05), n_nom_si),
            x        = bquote(beta[sex%*%.(DISEASE_CAPS)]~"(interaction coefficient)"),
            y        = expression(-log[10]~italic(P))
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "none",
              plot.subtitle   = element_text(size = 6.5, colour = "grey45"))

    save_supp(pSA, "a_sex_interaction_volcano", 4.0, 3.8)
}

# ---------------------------------------------------------------------------
# Panel j — Sex-stratified ML: 6-curve ROC (female/male × combined/pre/post)
# ---------------------------------------------------------------------------
cat("Building panel j (6-curve ROC)...\n")

ML_ROC_FILE  <- file.path(PROJ_DIR, "results", "ml", glue("{COHORT}_sex_ml_roc_data.csv"))
ML_SHAP_FILE <- file.path(PROJ_DIR, "results", "ml", glue("{COHORT}_sex_ml_shap.csv"))

tryCatch({
    if (!file.exists(ML_ROC_FILE))
        stop("ROC data not found — run 05_ms_sex_stratified_ml.R first")

    roc_all <- fread(ML_ROC_FILE)

    ROC_COLS_J <- c(
        "Female Combined"   = "#CC0066",
        "Female Post-onset" = "#E8779D",
        "Female Pre-onset"  = "#F2B8CE",
        "Male Combined"     = "#2B4C7E",
        "Male Post-onset"   = "#5B7FA6",
        "Male Pre-onset"    = "#A8BFCF"
    )
    ROC_LTY_J <- c(
        "Female Combined"   = "solid",
        "Female Post-onset" = "dashed",
        "Female Pre-onset"  = "dotted",
        "Male Combined"     = "solid",
        "Male Post-onset"   = "dashed",
        "Male Pre-onset"    = "dotted"
    )
    lvl_order_j <- names(ROC_COLS_J)

    # Faceted: one panel per stage, female vs male within each
    STAGE_LEVELS <- c("Combined", "Post-onset", "Pre-onset")
    SEX_COLS_J   <- c("Female" = "#CC0066", "Male" = "#2B4C7E")

    roc_all[, stage_f := factor(stage, levels = STAGE_LEVELS)]
    roc_all[, sex_f   := factor(sex_lbl, levels = c("Female", "Male"))]
    roc_all <- roc_all[!is.na(stage_f)]

    # Per-facet AUC labels placed inside each panel
    auc_lbl <- unique(roc_all[, .(sex_f, stage_f, auc)])
    auc_lbl[, x_pos := 0.52]
    auc_lbl[, y_pos := fifelse(sex_f == "Female", 0.18, 0.09)]
    auc_lbl[, label := sprintf("%s: %.2f", as.character(sex_f), auc)]
    auc_lbl[, colour := SEX_COLS_J[as.character(sex_f)]]

    pJ <- ggplot(roc_all, aes(x = fpr, y = tpr, colour = sex_f)) +
        geom_abline(slope = 1, intercept = 0, colour = "grey80", linewidth = 0.3) +
        geom_line(linewidth = 0.85) +
        geom_text(data = auc_lbl,
                  aes(x = x_pos, y = y_pos, label = label, colour = sex_f),
                  size = 2.4, hjust = 0, inherit.aes = FALSE, show.legend = FALSE) +
        facet_wrap(~stage_f, nrow = 1) +
        scale_colour_manual(values = SEX_COLS_J, name = NULL) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(
            title    = glue("j  Sex-stratified {DISEASE_CAPS} classifier ROC curves"),
            subtitle = "glmnet elastic net | held-out test set",
            x = "1 \u2013 Specificity (FPR)",
            y = "Sensitivity (TPR)"
        ) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom",
              strip.text      = element_text(size = 9, face = "bold"),
              plot.subtitle   = element_text(size = 7, colour = "grey40"))

    save_panel(pJ, "j_ml_roc_6curves", 8.0, 3.5)

}, error = function(e) cat("  Panel j skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Supp e — SHAP paired bar chart (combined female vs male models)
# ---------------------------------------------------------------------------
cat("Building supp e (SHAP paired importance)...\n")
tryCatch({
    if (!file.exists(ML_SHAP_FILE))
        stop("SHAP file not found — run 05_ms_sex_stratified_ml.R first")

    shap <- fread(ML_SHAP_FILE)
    shap <- merge(shap, wide[, .(protein, dep_cat)], by = "protein", all.x = TRUE)
    shap[is.na(dep_cat), dep_cat := "NS"]
    shap[, max_shap := pmax(female_shap, male_shap)]
    top_shap <- shap[order(-max_shap)][seq_len(min(35, .N))]
    top_shap[, prot_f := factor(toupper(protein), levels = rev(toupper(protein)))]

    shap_long <- melt(top_shap,
                      id.vars     = c("protein","prot_f","dep_cat","max_shap"),
                      measure.vars = c("female_shap","male_shap"),
                      variable.name = "sex", value.name = "shap_val")
    shap_long[, sex_lbl := factor(
        fifelse(sex == "female_shap", "Female", "Male"),
        levels = c("Female","Male")
    )]

    pSE <- ggplot(shap_long,
                  aes(x = shap_val, y = prot_f, fill = sex_lbl)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.9) +
        scale_fill_manual(values = c("Female" = "#CC0066", "Male" = "#2B4C7E"), name = NULL) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
        labs(
            title    = glue("e  Sex-stratified {DISEASE_CAPS} classifier: SHAP feature importance"),
            subtitle = "Linear SHAP = |coef| \u00d7 sd(feature) | combined models | top 35 | empty bar = absent in that sex",
            x = "Mean |SHAP| (log-odds scale)", y = NULL
        ) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y     = element_text(size = 6.5),
              legend.position = "top",
              plot.subtitle   = element_text(size = 6.5, colour = "grey40"))

    save_supp(pSE, "e_shap_paired", 5.5, 8.0)

}, error = function(e) cat("  Supp e skipped:", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------------
# Supplemental Extended Data Fig. 4 — Pre- and post-onset sex logFC concordance
#   2-panel scatter: F-pre vs M-pre (left), F-post vs M-post (right)
#   Demonstrates whether sex dimorphism is stable or shifts with disease stage.
# ---------------------------------------------------------------------------
cat("Building Extended Data Fig. 4 (pre/post sex concordance)...\n")

tryCatch({
    PRE_F  <- file.path(SEX_DIR, glue("{COHORT}_pre_female_vs_hc_all.csv"))
    PRE_M  <- file.path(SEX_DIR, glue("{COHORT}_pre_male_vs_hc_all.csv"))
    POST_F <- file.path(SEX_DIR, glue("{COHORT}_post_female_vs_hc_all.csv"))
    POST_M <- file.path(SEX_DIR, glue("{COHORT}_post_male_vs_hc_all.csv"))

    missing <- !file.exists(c(PRE_F, PRE_M, POST_F, POST_M))
    if (any(missing)) stop("Stage-stratified sex results missing — run 01_sex_stratified_ms_deps.R")

    pre_f  <- fread(PRE_F);  pre_m  <- fread(PRE_M)
    post_f <- fread(POST_F); post_m <- fread(POST_M)

    make_stage_scatter <- function(res_f, res_m, stage_lbl, pan_letter) {
        d <- merge(
            res_f[, .(protein, logFC_f = logFC, fdr_f = adj.P.Val)],
            res_m[, .(protein, logFC_m = logFC, fdr_m = adj.P.Val)],
            by = "protein"
        )
        d <- merge(d, wide[, .(protein, dep_cat)], by = "protein", all.x = TRUE)
        d[is.na(dep_cat), dep_cat := "NS"]
        d[, dep_cat := factor(dep_cat,
            levels = c("Female-specific","Both concordant","Both discordant","Male-specific","NS"))]

        r_all  <- round(cor(d$logFC_f, d$logFC_m, use = "complete.obs"), 3)
        n_f    <- sum(d$fdr_f < 0.05, na.rm = TRUE)
        n_m    <- sum(d$fdr_m < 0.05, na.rm = TRUE)

        lab_prot <- unique(c(
            d[fdr_f < 0.05][order(-abs(logFC_f))][seq_len(min(4,.N)), protein],
            d[fdr_m < 0.05][order(-abs(logFC_m))][seq_len(min(4,.N)), protein]
        ))
        lab_dt <- d[protein %in% lab_prot, .(protein=toupper(protein), logFC_f, logFC_m, dep_cat)]

        ggplot(d, aes(x=logFC_m, y=logFC_f, colour=dep_cat, size=dep_cat)) +
            geom_hline(yintercept=0, linewidth=0.22, colour="grey70") +
            geom_vline(xintercept=0, linewidth=0.22, colour="grey70") +
            geom_abline(slope=1, intercept=0, linetype="dashed",
                        linewidth=0.3, colour="grey55") +
            geom_point(data=d[dep_cat=="NS"], alpha=0.22) +
            geom_point(data=d[dep_cat!="NS"], alpha=0.88) +
            geom_text_repel(data=lab_dt, aes(label=protein),
                            size=2.1, fontface="bold", max.overlaps=Inf,
                            segment.size=0.22, segment.color="grey50",
                            box.padding=0.32, show.legend=FALSE,
                            min.segment.length=0.2) +
            annotate("text",
                     x=max(d$logFC_m, na.rm=TRUE)*0.9,
                     y=min(d$logFC_f, na.rm=TRUE)*0.85,
                     label=sprintf("r = %.3f", r_all),
                     size=2.2, hjust=1, colour="grey35") +
            scale_colour_manual(values=CAT_COLS, name=NULL,
                                guide=guide_legend(override.aes=list(size=2.5,alpha=1))) +
            scale_size_manual(values=CAT_SZ, guide="none") +
            labs(
                title    = sprintf("%s  %s: Female vs Male logFC (full cohort categories)",
                                   pan_letter, stage_lbl),
                subtitle = sprintf("F FDR<0.05: %d | M FDR<0.05: %d | r=%.3f",
                                   n_f, n_m, r_all),
                x = bquote(log[2]~"FC (Male"~.(DISEASE_CAPS)~"vs HC)"),
                y = bquote(log[2]~"FC (Female"~.(DISEASE_CAPS)~"vs HC)")
            ) +
            theme_ukb(base_size=9) +
            theme(legend.position="bottom",
                  plot.subtitle=element_text(size=7, colour="grey40"))
    }

    pE4a <- make_stage_scatter(pre_f,  pre_m,  "Pre-onset",  "a")
    pE4b <- make_stage_scatter(post_f, post_m, "Post-onset", "b")

    out_e4a <- file.path(SUPP_DIR, "efig4_a_preonset_sex_concordance.pdf")
    out_e4b <- file.path(SUPP_DIR, "efig4_b_postonset_sex_concordance.pdf")
    ggsave(out_e4a, pE4a, width=4.5, height=4.5, device=cairo_pdf)
    ggsave(out_e4b, pE4b, width=4.5, height=4.5, device=cairo_pdf)
    cat("  Saved: efig4_a_preonset_sex_concordance.pdf\n")
    cat("  Saved: efig4_b_postonset_sex_concordance.pdf\n")

}, error = function(e) cat("  Extended Data Fig. 4 skipped:", conditionMessage(e), "\n"))

cat("\nFigure 4 complete. Panels saved to:", FIG_DIR, "\n")
cat("Supplement panels saved to:", SUPP_DIR, "\n")
