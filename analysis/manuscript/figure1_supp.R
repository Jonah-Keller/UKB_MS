#!/usr/bin/env Rscript
# figure1_supp.R — Supplementary Figure 1: Full 4-analysis × 5-viztype matrix
#
# For each of the 4 analyses runs the same 5 visualization types:
#   V1. Volcano plot
#   V2. Ranked barplot (top 20 by nominal p)
#   V3. GO:BP enrichment (up + down)
#   V4. HPA cell-type expression heatmap (pheatmap)
#   V5. Cell-type enrichment dot plot (preranked fgsea NES; legacy Fisher fallback)
#
# Panel naming: panel_{analysis}_{viz}.pdf
#   e.g., panel_combined_volcano.pdf, panel_pre_go.pdf, panel_post_hpa.pdf
#
# Output: results/figures/1_supp/

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(scales)
    library(glue)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
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

FIG_DIR  <- file.path(PROJ_DIR, "results", "figures", "1_supp")
DATA_DIR <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed")
DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
NET_DIR  <- file.path(PROJ_DIR, "results", "networks")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# All colour constants come from ukb_theme.R (sourced above).
# Local aliases for readability only — do not redefine values here.

# ---------------------------------------------------------------------------
# Helper: save panel
# ---------------------------------------------------------------------------
save_panel <- function(p, name, w = 3.5, h = 3.5) {
    out <- file.path(FIG_DIR, paste0("panel_", name, ".pdf"))
    ggsave(out, p, width = w, height = h, device = cairo_pdf)
    cat("  Saved: panel_", name, ".pdf\n", sep = "")
}

# ---------------------------------------------------------------------------
# Helper: volcano
# ---------------------------------------------------------------------------
make_volcano <- function(df, title_str, n_label = 15, fdr_line = 0.05,
                         up_col = COL_PRE, down_col = COL_DOWN,
                         up_label = "Up", down_label = "Down") {
    df <- copy(df)
    if ("adj.P.Val" %in% names(df) && !"fdr" %in% names(df)) df[, fdr := adj.P.Val]
    df[, direction := fcase(fdr < fdr_line & logFC > 0, "Up",
                             fdr < fdr_line & logFC < 0, "Down", default = "NS")]
    df[, direction := factor(direction, levels = c("Up", "Down", "NS"))]
    df[, nlp := -log10(P.Value)]
    bonf <- -log10(0.05 / nrow(df))
    top  <- df[fdr < fdr_line][order(-abs(logFC))][seq_len(min(n_label, .N))]
    n_up <- sum(df$direction == "Up", na.rm = TRUE)
    n_dn <- sum(df$direction == "Down", na.rm = TRUE)
    ggplot(df, aes(x = logFC, y = nlp, colour = direction)) +
        geom_point(size = 0.55, alpha = 0.55) +
        geom_hline(yintercept = bonf, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
        geom_text_repel(data = top, aes(label = toupper(protein)),
                        size = 2.0, fontface = "bold", max.overlaps = 20,
                        segment.size = 0.2, box.padding = 0.25, show.legend = FALSE) +
        scale_colour_manual(
            values = c("Up" = up_col, "Down" = down_col, "NS" = "grey65"),
            labels = c("Up" = sprintf("%s (n=%d)", up_label, n_up),
                       "Down" = sprintf("%s (n=%d)", down_label, n_dn), "NS" = "NS"),
            name = NULL
        ) +
        annotate("text", x = max(df$logFC, na.rm=TRUE) * 0.92, y = bonf + 0.4,
                 label = "Bonferroni", size = 2.1, hjust = 1, colour = "grey45") +
        labs(title = title_str, x = expression(log[2]~"fold change"),
             y = expression(-log[10]~italic(P))) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

# ---------------------------------------------------------------------------
# Helper: ranked barplot
# ---------------------------------------------------------------------------
make_ranked_barplot <- function(df, title_str, n_each = 20,
                                up_col = COL_PRE, down_col = COL_DOWN,
                                up_label = glue("Up in {DISEASE_CAPS}"),
                                down_label = glue("Down in {DISEASE_CAPS}")) {
    df <- copy(df)
    if ("adj.P.Val" %in% names(df) && !"fdr" %in% names(df)) df[, fdr := adj.P.Val]
    df[, direction := fifelse(logFC > 0, up_label, down_label)]
    df[, protein_upper := toupper(protein)]
    top_up   <- df[logFC > 0][order(P.Value)][seq_len(min(n_each, .N))]
    top_down <- df[logFC < 0][order(P.Value)][seq_len(min(n_each, .N))]
    bar <- rbind(top_up, top_down)
    bar[, protein_f := factor(protein_upper, levels = protein_upper[order(logFC)])]
    bar[, sig_str := fcase(fdr < 0.001, "***", fdr < 0.01, "**", fdr < 0.05, "*",
                            P.Value < 0.05, ".", default = "")]
    ggplot(bar, aes(x = logFC, y = protein_f, fill = direction)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = sig_str,
                      x = ifelse(logFC > 0, logFC + 0.01, logFC - 0.01),
                      hjust = ifelse(logFC > 0, 0, 1)),
                  size = 2.5, colour = "grey30") +
        geom_vline(xintercept = 0, linewidth = 0.3) +
        scale_fill_manual(values = c(setNames(up_col, up_label), setNames(down_col, down_label)),
                          name = NULL) +
        labs(title = title_str,
             subtitle = "Top 20 per direction by nominal p | * FDR<0.05",
             x = expression(log[2]~"fold change"), y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y = element_text(size = 6), legend.position = "bottom",
              plot.subtitle = element_text(size = 6.5, colour = "grey40"))
}

# ---------------------------------------------------------------------------
# Helper: GO dotplot
# ---------------------------------------------------------------------------
make_go_dotplot <- function(go_dt, title_str, up_col = COL_PRE, down_col = COL_DOWN,
                             up_dir_label = "Up", down_dir_label = "Down", n_top = 10) {
    go_dt <- copy(go_dt)
    up   <- go_dt[gene_set == "up_DEPs"   & p.adjust < 0.05]
    down <- go_dt[gene_set == "down_DEPs" & p.adjust < 0.05]
    if (nrow(up) > 0)   up[, direction   := up_dir_label]
    if (nrow(down) > 0) down[, direction := down_dir_label]
    top <- rbindlist(Filter(Negate(is.null), list(
        if (nrow(up)   > 0) up[order(p.adjust)][seq_len(min(n_top, .N))]   else NULL,
        if (nrow(down) > 0) down[order(p.adjust)][seq_len(min(n_top, .N))] else NULL
    )))
    if (nrow(top) == 0) return(NULL)
    top[, direction := factor(direction, levels = c(up_dir_label, down_dir_label))]
    top[, Description := factor(Description, levels = top[order(direction, FoldEnrichment), Description])]
    dir_vals <- c(setNames(up_col, up_dir_label), setNames(down_col, down_dir_label))
    ggplot(top, aes(x = FoldEnrichment, y = Description, size = Count, colour = direction)) +
        geom_point() +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
        scale_colour_manual(values = dir_vals,
                            guide = guide_legend(override.aes = list(size = 3.5))) +
        scale_size_continuous(name = "Gene count", range = c(2, 6),
                              guide = guide_legend(override.aes = list(colour = "grey40"))) +
        scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        labs(title = title_str, x = "Fold enrichment", y = NULL) +
        theme_ukb(base_size = 9) +
        theme(axis.text.y = element_text(size = 7), strip.text = element_text(size = 8, face = "bold"),
              legend.position = "right")
}

# ---------------------------------------------------------------------------
# Helper: cell-type enrichment dot plot (single analysis)
# Uses the new preranked-fgsea pipeline (signed NES, specificity-filtered
# markers). Falls back to the legacy Fisher CSV + make_celltype_overrep_plot
# if the GSEA CSV is absent, so pre-refactor runs still render.
# ---------------------------------------------------------------------------
make_overrep_plot <- function(gsea_dt, analysis_name, title_str,
                               up_col = COL_PRE, down_col = COL_DOWN,
                               overrep_dt = NULL) {
    dt_gsea <- if (!is.null(gsea_dt) && nrow(gsea_dt) > 0L)
                    gsea_dt[analysis == analysis_name] else NULL
    if (!is.null(dt_gsea) && nrow(dt_gsea) > 0L) {
        return(make_celltype_gsea_plot(
            dt_gsea,
            contrast_col    = "analysis",
            contrast_order  = analysis_name,
            contrast_labels = setNames(analysis_name, analysis_name),
            title_str       = title_str,
            subtitle_str    = paste0(
                "Preranked fgsea on limma moderated-t | ",
                "NES>0: markers higher in cases vs HC | NES<0: lower | * FDR<0.05"
            ),
            high_col        = up_col,
            low_col         = down_col
        ))
    }
    # Legacy fallback (retained for reproducibility of pre-refactor runs).
    if (is.null(overrep_dt) || nrow(overrep_dt) == 0L) return(NULL)
    dt <- overrep_dt[analysis == analysis_name]
    if (nrow(dt) == 0L) return(NULL)
    make_celltype_overrep_plot(
        dt,
        title_str    = title_str,
        subtitle_str = "Fisher's exact (legacy) | nTPM >= 10 | * p_adj<0.05",
        high_col     = up_col,
        low_col      = down_col
    )
}

# ---------------------------------------------------------------------------
# Load shared data
# ---------------------------------------------------------------------------
cat("Loading data...\n")

qc       <- fread(file.path(DATA_DIR, glue("{COHORT}_olink_qc.csv")))
if (STATUS_COL != "ms_status" && STATUS_COL %in% names(qc))
    setnames(qc, STATUS_COL, "ms_status")
# Cell-type enrichment: primary source is <cohort>_celltype_gsea.csv (new
# preranked fgsea + specificity markers). <cohort>_celltype_overrep.csv is the
# legacy Fisher output kept for fallback in pre-refactor runs.
gsea     <- if (file.exists(file.path(NET_DIR, glue("{COHORT}_celltype_gsea.csv"))))
                fread(file.path(NET_DIR, glue("{COHORT}_celltype_gsea.csv"))) else data.table()
overrep  <- if (file.exists(file.path(NET_DIR, glue("{COHORT}_celltype_overrep.csv"))))
                fread(file.path(NET_DIR, glue("{COHORT}_celltype_overrep.csv"))) else data.table()

# Pre vs post (computed inline in figure1.R → read from saved CSV)
PREVSPOST_FILE <- file.path(DIFF_DIR, glue("{COHORT}_pre_vs_post.csv"))
prevspost_diff <- if (file.exists(PREVSPOST_FILE)) fread(PREVSPOST_FILE) else data.table()

# GO results for all 3 analyses that have them
go_comb  <- if (file.exists(file.path(NET_DIR, glue("{COHORT}_combined_go_results.csv"))))
                fread(file.path(NET_DIR, glue("{COHORT}_combined_go_results.csv"))) else data.table()
go_pre   <- if (file.exists(file.path(NET_DIR, glue("{COHORT}_pre_go_results.csv"))))
                fread(file.path(NET_DIR, glue("{COHORT}_pre_go_results.csv")))  else data.table()
go_post  <- if (file.exists(file.path(NET_DIR, glue("{COHORT}_post_go_results.csv"))))
                fread(file.path(NET_DIR, glue("{COHORT}_post_go_results.csv"))) else data.table()

# GO for pre vs post — run inline at nominal p < 0.05 (only 3 FDR-sig DEPs)
go_prevspost <- tryCatch({
    suppressPackageStartupMessages({
        library(clusterProfiler)
        library(org.Hs.eg.db)
    })
    pvp <- fread(file.path(DIFF_DIR, glue("{COHORT}_pre_vs_post.csv")))
    up_genes   <- toupper(pvp[P.Value < 0.05 & logFC > 0, protein])
    down_genes <- toupper(pvp[P.Value < 0.05 & logFC < 0, protein])
    run_go <- function(genes, set_label) {
        if (length(genes) < 5) return(data.table())
        ids <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID",
                      keytype = "SYMBOL", multiVals = "first")
        ids <- ids[!is.na(ids)]
        if (length(ids) < 5) return(data.table())
        res <- enrichGO(gene = ids, OrgDb = org.Hs.eg.db, ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2, readable = TRUE)
        if (is.null(res) || nrow(as.data.frame(res)) == 0) return(data.table())
        dt <- as.data.table(as.data.frame(res))
        dt[, gene_set := set_label]
        dt
    }
    rbindlist(list(run_go(up_genes, "up_DEPs"), run_go(down_genes, "down_DEPs")), fill = TRUE)
}, error = function(e) { cat("  Pre-vs-post GO failed:", conditionMessage(e), "\n"); data.table() })
cat(sprintf("  Pre-vs-post GO: %d terms (up %d / down %d)\n",
            nrow(go_prevspost),
            sum(go_prevspost$gene_set == "up_DEPs",  na.rm = TRUE),
            sum(go_prevspost$gene_set == "down_DEPs", na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Define the 4 analyses
# ---------------------------------------------------------------------------
analyses <- list(
    list(
        name       = "combined",
        label      = glue("Combined {DISEASE_CAPS} vs HC"),
        diff_file  = glue("{COHORT}_combined_vs_hc.csv"),
        fdr_thr    = 0.05,
        up_col     = COL_COMB_UP,    # golden yellow
        down_col   = COL_COMB_DOWN,  # blue-violet
        up_label   = glue("Up in {DISEASE_CAPS}"),
        down_label = glue("Down in {DISEASE_CAPS}"),
        go_data    = "combined"
    ),
    # Pre-onset panels are in Figure 1 (c–h) — omitted here to avoid duplication
    list(
        name       = "post",
        label      = "Post-onset vs HC",
        diff_file  = glue("{COHORT}_post_vs_hc.csv"),
        fdr_thr    = 0.05,
        up_col     = COL_POST_UP,    # vivid spring green
        down_col   = COL_POST_DOWN,  # vermillion
        up_label   = glue("Up in {DISEASE_CAPS}"),
        down_label = glue("Down in {DISEASE_CAPS}"),
        go_data    = "post"
    ),
    list(
        name       = "preVsPost",
        label      = "Pre vs post-onset",
        diff_file  = glue("{COHORT}_pre_vs_post.csv"),
        fdr_thr    = 0.05,
        nominal    = TRUE,           # use P.Value < 0.05, not FDR
        up_col     = COL_PRE,        # hot pink = pre-onset identity
        down_col   = COL_POST,       # navy    = post-onset identity
        up_label   = "Higher pre-onset",
        down_label = "Higher post-onset",
        go_data    = "preVsPost"
    )
)

# ---------------------------------------------------------------------------
# Main loop: 4 analyses × 5 viz types, panels named a, b, c... sequentially
#
# Panel order (manuscript order):
#   Combined  (A-E): volcano, barplot, GO, HPA, over-rep
#   Post-onset (F-J): same 5
#   Pre vs post (K-N): volcano, barplot, GO (nominal p<0.05), over-rep  (HPA skipped)
#   Pre-onset panels are in Figure 1 (c–h) — omitted here to avoid duplication
# ---------------------------------------------------------------------------
panel_idx <- 0L   # increments before each successful save

next_panel <- function() {
    panel_idx <<- panel_idx + 1L
    LETTERS[panel_idx]
}

for (an in analyses) {
    cat(sprintf("\n=== %s ===\n", an$label))

    diff_path <- file.path(DIFF_DIR, an$diff_file)
    if (!file.exists(diff_path)) {
        cat(sprintf("  Skipping: %s not found\n", an$diff_file))
        next
    }
    diff <- fread(diff_path)

    n_cases <- if (an$name == "combined") {
        sum(qc$ms_status %in% c(SV$pre_onset, SV$post_onset), na.rm = TRUE)
    } else if (an$name == "pre") {
        sum(qc$ms_status == SV$pre_onset, na.rm = TRUE)
    } else if (an$name == "post") {
        sum(qc$ms_status == SV$post_onset, na.rm = TRUE)
    } else {
        sum(qc$ms_status %in% c(SV$pre_onset, SV$post_onset), na.rm = TRUE)
    }

    an_tag   <- an$name   # used in descriptive filenames
    use_nom  <- isTRUE(an$nominal)
    sig_lbl  <- if (use_nom) "nominal p<0.05" else sprintf("FDR<%.2f", an$fdr_thr)

    # For nominal analyses, override adj.P.Val with P.Value so volcano/barplot
    # colour and count by nominal p rather than BH-adjusted p
    if (use_nom) diff[, adj.P.Val := P.Value]

    # V1: Volcano
    cat("  V1 volcano...\n")
    tryCatch({
        lbl <- next_panel()
        pV1 <- make_volcano(
            diff,
            sprintf("S1%s  %s — volcano (n=%d, %s)", lbl, an$label, n_cases, sig_lbl),
            fdr_line = an$fdr_thr, up_col = an$up_col, down_col = an$down_col,
            up_label = an$up_label, down_label = an$down_label
        )
        save_panel(pV1, paste0(lbl, "_volcano_", an_tag), 3.5, 4)
    }, error = function(e) { panel_idx <<- panel_idx - 1L; cat("  V1 skipped:", conditionMessage(e), "\n") })

    # V2: Ranked barplot
    cat("  V2 barplot...\n")
    tryCatch({
        lbl <- next_panel()
        pV2 <- make_ranked_barplot(
            diff,
            sprintf("S1%s  %s — ranked barplot (top 20, %s)", lbl, an$label, sig_lbl),
            up_col = an$up_col, down_col = an$down_col,
            up_label = an$up_label, down_label = an$down_label
        )
        save_panel(pV2, paste0(lbl, "_barplot_", an_tag), 5.5, 5.5)
    }, error = function(e) { panel_idx <<- panel_idx - 1L; cat("  V2 skipped:", conditionMessage(e), "\n") })

    # V3: GO enrichment (skip if no data)
    cat("  V3 GO...\n")
    tryCatch({
        go_dt <- if (!is.null(an$go_data)) {
            if (an$go_data == "combined")   go_comb
            else if (an$go_data == "pre")   go_pre
            else if (an$go_data == "post")  go_post
            else if (an$go_data == "preVsPost") go_prevspost
            else NULL
        } else NULL
        if (is.null(go_dt) || nrow(go_dt) == 0) stop("No GO data for this analysis")
        lbl <- next_panel()
        go_suffix <- if (use_nom) sprintf(" (%s)", sig_lbl) else ""
        pV3 <- make_go_dotplot(
            go_dt,
            sprintf("S1%s  %s — GO:BP enrichment%s", lbl, an$label, go_suffix),
            up_col = an$up_col, down_col = an$down_col,
            up_dir_label = an$up_label, down_dir_label = an$down_label
        )
        if (!is.null(pV3)) save_panel(pV3, paste0(lbl, "_go_", an_tag), 5.5, 6)
        else { panel_idx <<- panel_idx - 1L; cat("  V3 no terms to plot\n") }
    }, error = function(e) cat("  V3 skipped:", conditionMessage(e), "\n"))

    # V4: HPA cell-type heatmap (skip if file absent)
    cat("  V4 HPA heatmap...\n")
    tryCatch({
        suppressPackageStartupMessages(library(pheatmap))
        hpa_file <- file.path(NET_DIR, sprintf("%s_%s_hpa_celltype_matrix.csv", COHORT, an$name))
        if (!file.exists(hpa_file)) stop(sprintf("%s not found", basename(hpa_file)))
        lbl <- next_panel()

        hpa_mat  <- fread(hpa_file)
        ct_cols  <- setdiff(names(hpa_mat), c("gene", "logFC", "direction", "abs_logFC"))
        top40    <- hpa_mat[!is.na(logFC)][order(-abs(logFC))][seq_len(min(40, .N)), gene]
        mat_sub  <- hpa_mat[gene %in% top40, ..ct_cols]
        genes_v4 <- hpa_mat[gene %in% top40, gene]
        mat_sub  <- mat_sub[match(top40, genes_v4)]
        expr_v4  <- as.matrix(mat_sub); rownames(expr_v4) <- top40

        dep_dir_v4 <- hpa_mat[gene %in% top40, .(gene, direction)]
        dep_dir_v4 <- dep_dir_v4[match(top40, dep_dir_v4$gene)]
        ann_v4 <- data.frame(Direction = dep_dir_v4$direction, row.names = top40)
        unique_dirs_v4 <- unique(dep_dir_v4$direction[!is.na(dep_dir_v4$direction)])
        dir_col_v4 <- c()
        for (d in unique_dirs_v4)
            dir_col_v4[d] <- if (grepl("Up|Higher pre", d, ignore.case = TRUE)) an$up_col else an$down_col

        mat_v4 <- t(expr_v4); mat_v4[is.na(mat_v4)] <- 0
        ct_rename <- c("Microglial cells" = "Microglia",
                       "Oligodendrocyte precursor cells" = "OPCs",
                       "Excitatory neurons" = "Exc. neurons",
                       "Inhibitory neurons" = "Inh. neurons",
                       "dendritic cells" = "Dendritic cells",
                       "monocytes" = "Monocytes", "granulocytes" = "Granulocytes")
        rownames(mat_v4) <- ifelse(rownames(mat_v4) %in% names(ct_rename),
                                    ct_rename[rownames(mat_v4)], rownames(mat_v4))

        out_v4 <- file.path(FIG_DIR, sprintf("panel_%s_celltype_%s.pdf", lbl, an_tag))
        pdf(out_v4, width = 10, height = 5, onefile = FALSE)
        pheatmap(mat_v4,
            color = colorRampPalette(c("white", an$up_col))(100),
            breaks = seq(0, 1, length.out = 101),
            cluster_rows = TRUE, cluster_cols = TRUE,
            show_rownames = TRUE, show_colnames = TRUE,
            fontsize_row = 9, fontsize_col = 6,
            annotation_col = ann_v4, annotation_colors = list(Direction = dir_col_v4),
            annotation_names_col = FALSE, border_color = NA,
            main = sprintf("S1%s  %s — HPA cell-type expression (top 40 DEPs)", lbl, an$label)
        )
        dev.off()
        cat(sprintf("  Saved: panel_%s_celltype_%s.pdf\n", lbl, an_tag))
    }, error = function(e) cat("  V4 skipped:", conditionMessage(e), "\n"))

    # V5: Cell-type enrichment (new fgsea pipeline; Fisher fallback)
    cat("  V5 cell-type enrichment...\n")
    tryCatch({
        if (nrow(gsea) == 0 && nrow(overrep) == 0)
            stop(glue("neither {COHORT}_celltype_gsea.csv nor {COHORT}_celltype_overrep.csv present"))
        lbl <- next_panel()
        pV5 <- make_overrep_plot(
            gsea_dt      = gsea,
            analysis_name = an$name,
            title_str    = sprintf("S1%s  %s \u2014 cell-type enrichment (fgsea NES)", lbl, an$label),
            up_col       = an$up_col,
            down_col     = an$down_col,
            overrep_dt   = overrep
        )
        if (!is.null(pV5)) save_panel(pV5, paste0(lbl, "_celltype_gsea_", an_tag), 4.5, 7.5)
        else { panel_idx <<- panel_idx - 1L; cat("  V5 no data\n") }
    }, error = function(e) cat("  V5 skipped:", conditionMessage(e), "\n"))
}

cat(sprintf("\nSupplementary Figure 1 complete: %d panels (a-%s) saved to: %s\n",
            panel_idx, LETTERS[panel_idx], FIG_DIR))
