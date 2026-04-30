#!/usr/bin/env Rscript
# go_dotplot.R — reusable GO:BP enrichment dot plot
#
# Extracted from figure1.R panel e. Call make_go_dotplot() to produce a dot
# plot identical in style to figure 1 panel e (panel_e_go_pre.pdf).
#
# Required columns in go_dt (same schema as ms_pre_go_results.csv):
#   Description     — GO term name
#   FoldEnrichment  — fold enrichment vs background
#   Count           — number of DEP genes in GO term
#   p.adjust        — BH-adjusted p-value
#   gene_set        — "up_DEPs" or "down_DEPs"
#
# Caller must have sourced ukb_theme.R and loaded ggplot2 + data.table.
# clusterProfiler + org.Hs.eg.db are required for run_go() / collect_go().

# ── GO:BP ORA helpers (shared by all cluster GO analyses) ────────────────────

run_go <- function(gene_symbols, label, min_genes = 5L) {
    if (length(gene_symbols) < min_genes) {
        message(sprintf("  Skipping %s: only %d genes (need \u2265%d)",
                        label, length(gene_symbols), min_genes))
        return(NULL)
    }
    fg <- suppressMessages(
        clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                              toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
    if (nrow(fg) == 0L) { message("  No Entrez mappings for ", label); return(NULL) }
    ego <- tryCatch(
        clusterProfiler::enrichGO(
            gene          = fg$ENTREZID,
            OrgDb         = org.Hs.eg.db,
            ont           = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff  = 0.05,
            qvalueCutoff  = 0.20,
            readable      = TRUE
        ),
        error = function(e) { message("  enrichGO failed for ", label, ": ", e$message); NULL }
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0L) {
        message(sprintf("  No significant GO BP terms for %s", label))
        return(NULL)
    }
    cat(sprintf("  %s: %d significant GO BP terms\n", label, nrow(as.data.frame(ego))))
    ego
}

collect_go <- function(ego, gene_set_label) {
    if (is.null(ego)) return(NULL)
    dt <- data.table::as.data.table(ego@result)[p.adjust < 0.05]
    dt[, gene_set := gene_set_label]
    dt
}

make_go_dotplot <- function(
    go_dt,
    title_str,
    subtitle_str = NULL,
    up_col,
    down_col   = COL_DOWN,
    up_label   = "Up-regulated DEPs",
    down_label = "Down-regulated DEPs",
    max_n      = 10L
) {
    dt <- copy(go_dt)
    dt[gene_set == "up_DEPs",   direction := up_label]
    dt[gene_set == "down_DEPs", direction := down_label]
    dt <- dt[!is.na(direction) & !is.na(p.adjust)]

    pre_up   <- dt[direction == up_label   & p.adjust < 0.05][order(p.adjust)][seq_len(min(max_n, .N))]
    pre_down <- dt[direction == down_label & p.adjust < 0.05][order(p.adjust)][seq_len(min(max_n, .N))]

    top <- rbindlist(Filter(Negate(is.null), list(
        if (nrow(pre_up)   > 0L) pre_up   else NULL,
        if (nrow(pre_down) > 0L) pre_down else NULL
    )))
    if (nrow(top) == 0L) {
        message(sprintf("  make_go_dotplot: no FDR<0.05 GO terms for '%s'", title_str))
        return(NULL)
    }

    top[, direction := factor(direction, levels = c(up_label, down_label))]
    top[, Description := factor(Description,
                                levels = unique(top[order(direction, FoldEnrichment), Description]))]

    p <- ggplot(top, aes(x = FoldEnrichment, y = Description,
                          size = Count, colour = direction)) +
        geom_point() +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
        scale_colour_manual(
            values = setNames(c(up_col, down_col), c(up_label, down_label)),
            guide  = guide_legend(override.aes = list(size = 3.5))
        ) +
        scale_size_continuous(
            name  = "Gene count",
            range = c(2, 6),
            guide = guide_legend(override.aes = list(colour = "grey40"))
        ) +
        scale_x_continuous(limits = c(0, NA),
                           expand  = expansion(mult = c(0, 0.05))) +
        scale_y_discrete(labels = function(x)
            ifelse(nchar(x) > 38, paste0(substr(x, 1, 35), "\u2026"), x)) +
        labs(title    = title_str,
             subtitle = subtitle_str,
             x = "Fold enrichment", y = NULL) +
        theme_ukb(base_size = 9) +
        theme(
            axis.text.y     = element_text(size = 7, lineheight = 1.1),
            strip.text      = element_text(size = 8, face = "bold"),
            legend.position = "right",
            plot.subtitle   = element_text(size = 7, colour = "grey40")
        )
    p
}

# ── make_cross_cluster_go_dotplot() ──────────────────────────────────────────
# Cross-cluster GO dotplot (figure 5 panel m / 5S panel g style).
# go_all must have columns: cluster, gene_set, ID, Description,
#   FoldEnrichment, Count, p.adjust.
make_cross_cluster_go_dotplot <- function(
    go_all,
    title_str,
    subtitle_str      = NULL,
    max_n_per_cluster = 8L,
    high_col          = COL_PRE
) {
    top_terms <- go_all[, .SD[order(p.adjust)][seq_len(min(max_n_per_cluster, .N))],
                         by = .(cluster, gene_set)]
    keep_ids  <- unique(top_terms$ID)

    plot_dt <- go_all[ID %in% keep_ids,
                       .(cluster, gene_set, Description, FoldEnrichment, Count, p.adjust)]
    plot_dt[, cluster  := factor(cluster, levels = c("C0", "C1", "C2"))]
    plot_dt[, gene_set := factor(gene_set, levels = c("up_DEPs", "down_DEPs"),
                                  labels = c("Up-regulated", "Down-regulated"))]

    desc_order <- plot_dt[, .(mean_fe = mean(FoldEnrichment, na.rm = TRUE)), by = Description
                           ][order(-mean_fe), Description]
    plot_dt[, Description := factor(Description, levels = desc_order)]

    ggplot(plot_dt,
           aes(x = cluster, y = Description,
               size = pmin(-log10(p.adjust), 4),
               colour = FoldEnrichment)) +
        geom_point(alpha = 0.88) +
        facet_wrap(~ gene_set, nrow = 1) +
        scale_colour_gradient(
            low = "grey90", high = high_col,
            limits = c(0, NA),
            name = "Fold\nenrichment",
            guide = guide_colourbar(barwidth = 0.5, barheight = 3.5)
        ) +
        scale_size_continuous(
            name   = expression(-log[10] ~ italic(p[adj])),
            range  = c(1.5, 7), breaks = c(1, 2, 3, 4),
            labels = c("1", "2", "3", "\u22654")
        ) +
        labs(x = NULL, y = NULL, title = title_str, subtitle = subtitle_str) +
        scale_y_discrete(labels = function(x)
            ifelse(nchar(x) > 38, paste0(substr(x, 1, 35), "\u2026"), x)) +
        theme_ukb(base_size = 9) +
        theme(
            axis.text.x      = element_text(size = 8.5, face = "bold"),
            axis.text.y      = element_text(size = 7.5, lineheight = 1.1),
            strip.text       = element_text(size = 8.5, face = "bold"),
            panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
            legend.position  = "right"
        )
}
