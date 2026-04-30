#!/usr/bin/env Rscript
# 01_ms_als_comparison.R
# Cross-disease comparison: MS vs ALS differential proteomics
#
# Panels produced:
#   A) Venn diagram — shared/unique DEPs (combined, pre-onset, post-onset)
#   B) logFC scatter — all 2911 proteins, coloured by DEP membership
#   C) Pre-onset focused scatter (presymptomatic signals)
#   D) Direction-concordance bar chart for shared proteins
#
# Outputs:  results/cross_disease/
#   ms_als_venn_combined.pdf
#   ms_als_logfc_scatter.pdf
#   ms_als_preonset_scatter.pdf
#   ms_als_shared_proteins.csv

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

# ── Path resolution ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
    getwd()
}
PROJ_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

DIFF_DIR   <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR    <- file.path(PROJ_DIR, "results", "cross_disease")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR <- 0.05

# ── Load all six DEP result files ─────────────────────────────────────────────
load_deps <- function(disease, contrast) {
    fname <- file.path(DIFF_DIR, sprintf("%s_%s_vs_hc.csv", disease, contrast))
    dt <- fread(fname)
    dt[, disease   := disease]
    dt[, contrast  := contrast]
    dt[, sig       := adj.P.Val < FDR_THR]
    dt[, direction := fifelse(sig & logFC > 0, "up",
                      fifelse(sig & logFC < 0, "down", "ns"))]
    dt
}

ms_comb  <- load_deps("ms",  "combined")
ms_pre   <- load_deps("ms",  "pre")
ms_post  <- load_deps("ms",  "post")
als_comb <- load_deps("als", "combined")
als_pre  <- load_deps("als", "pre")
als_post <- load_deps("als", "post")

cat(sprintf("MS  combined DEPs: %d  |  ALS combined DEPs: %d\n",
            sum(ms_comb$sig),  sum(als_comb$sig)))
cat(sprintf("MS  pre-onset DEPs: %d |  ALS pre-onset DEPs: %d\n",
            sum(ms_pre$sig),   sum(als_pre$sig)))
cat(sprintf("MS  post-onset DEPs: %d | ALS post-onset DEPs: %d\n",
            sum(ms_post$sig),  sum(als_post$sig)))


# ── Shared / unique sets ───────────────────────────────────────────────────────
shared_combined <- function(ms_dt, als_dt) {
    ms_sig  <- ms_dt[sig == TRUE,  protein]
    als_sig <- als_dt[sig == TRUE, protein]
    list(
        ms_only   = setdiff(ms_sig,  als_sig),
        als_only  = setdiff(als_sig, ms_sig),
        shared    = intersect(ms_sig, als_sig)
    )
}

sets_comb <- shared_combined(ms_comb,  als_comb)
sets_pre  <- shared_combined(ms_pre,   als_pre)
sets_post <- shared_combined(ms_post,  als_post)

cat("\n=== Combined contrast ===\n")
cat(sprintf("  MS-only: %d  |  ALS-only: %d  |  Shared: %d\n",
            length(sets_comb$ms_only), length(sets_comb$als_only), length(sets_comb$shared)))
cat("  Shared proteins:", paste(sort(sets_comb$shared), collapse=", "), "\n")

cat("\n=== Pre-onset contrast ===\n")
cat(sprintf("  MS-only: %d  |  ALS-only: %d  |  Shared: %d\n",
            length(sets_pre$ms_only), length(sets_pre$als_only), length(sets_pre$shared)))
cat("  Shared proteins:", paste(sort(sets_pre$shared), collapse=", "), "\n")

cat("\n=== Post-onset contrast ===\n")
cat(sprintf("  MS-only: %d  |  ALS-only: %d  |  Shared: %d\n",
            length(sets_post$ms_only), length(sets_post$als_only), length(sets_post$shared)))
cat("  Shared proteins:", paste(sort(sets_post$shared), collapse=", "), "\n")


# ── Shared-protein detail table ────────────────────────────────────────────────
build_shared_table <- function(ms_dt, als_dt, contrast_label) {
    ms_sig  <- ms_dt[sig == TRUE,  .(protein, ms_logFC  = logFC, ms_fdr   = adj.P.Val)]
    als_sig <- als_dt[sig == TRUE, .(protein, als_logFC = logFC, als_fdr  = adj.P.Val)]
    shared  <- merge(ms_sig, als_sig, by = "protein")
    shared[, contrast   := contrast_label]
    shared[, same_dir   := (ms_logFC > 0) == (als_logFC > 0)]
    shared[, ms_dir     := fifelse(ms_logFC  > 0, "up", "down")]
    shared[, als_dir    := fifelse(als_logFC > 0, "up", "down")]
    shared[order(-abs(ms_logFC))]
}

shared_comb_dt  <- build_shared_table(ms_comb,  als_comb,  "combined")
shared_pre_dt   <- build_shared_table(ms_pre,   als_pre,   "pre_onset")
shared_post_dt  <- build_shared_table(ms_post,  als_post,  "post_onset")

all_shared <- rbind(shared_comb_dt, shared_pre_dt, shared_post_dt, fill = TRUE)
fwrite(all_shared, file.path(OUT_DIR, "ms_als_shared_proteins.csv"))
cat(sprintf("\nShared protein table written: %d rows\n", nrow(all_shared)))

# Print combined table
cat("\n--- Combined shared proteins detail ---\n")
print(shared_comb_dt[, .(protein, ms_logFC = round(ms_logFC,3),
                          als_logFC = round(als_logFC,3), same_dir, ms_dir, als_dir)])


# ── Panel B: logFC scatter (combined) ─────────────────────────────────────────
# Merge all proteins present in both
scatter_dt <- merge(
    ms_comb[,  .(protein, ms_logFC  = logFC, ms_sig  = sig,  ms_dir  = direction)],
    als_comb[, .(protein, als_logFC = logFC, als_sig = sig, als_dir = direction)],
    by = "protein"
)

# Colour category
scatter_dt[, group := fcase(
    ms_sig == TRUE  & als_sig == TRUE  & (ms_logFC > 0) == (als_logFC > 0), "Shared (concordant)",
    ms_sig == TRUE  & als_sig == TRUE  & (ms_logFC > 0) != (als_logFC > 0), "Shared (discordant)",
    ms_sig == TRUE  & als_sig == FALSE, "MS only",
    ms_sig == FALSE & als_sig == TRUE,  "ALS only",
    default = "NS"
)]

scatter_dt[, label := fifelse(
    group %in% c("Shared (concordant)", "Shared (discordant)"), toupper(protein), NA_character_
)]

SCATTER_COLS <- c(
    "Shared (concordant)"  = "#CC0066",
    "Shared (discordant)"  = "#E6A817",
    "MS only"              = "#56B4E9",
    "ALS only"             = "#2B4C7E",
    "NS"                   = "grey80"
)
SCATTER_ALPHA <- c(
    "Shared (concordant)"  = 0.95,
    "Shared (discordant)"  = 0.95,
    "MS only"              = 0.7,
    "ALS only"             = 0.7,
    "NS"                   = 0.25
)
SCATTER_SIZE <- c(
    "Shared (concordant)"  = 2.0,
    "Shared (discordant)"  = 2.0,
    "MS only"              = 1.4,
    "ALS only"             = 1.4,
    "NS"                   = 0.7
)

# Force NS to bottom layer
scatter_dt[, group_f := factor(group, levels = c("NS", "MS only", "ALS only",
                                                   "Shared (discordant)",
                                                   "Shared (concordant)"))]

n_labels <- c(
    sprintf("Shared concordant n=%d",  sum(scatter_dt$group == "Shared (concordant)")),
    sprintf("Shared discordant n=%d",  sum(scatter_dt$group == "Shared (discordant)")),
    sprintf("MS only n=%d",            sum(scatter_dt$group == "MS only")),
    sprintf("ALS only n=%d",           sum(scatter_dt$group == "ALS only")),
    sprintf("NS n=%d",                 sum(scatter_dt$group == "NS"))
)

p_scatter <- ggplot(scatter_dt, aes(x = ms_logFC, y = als_logFC,
                                     colour = group_f, alpha = group_f, size = group_f)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50", linetype = "dashed") +
    geom_point(shape = 16) +
    geom_text_repel(
        aes(label = label),
        size = 2.5, colour = "black",
        max.overlaps = 25,
        box.padding  = 0.3,
        segment.size = 0.25,
        segment.colour = "grey40",
        na.rm = TRUE
    ) +
    scale_colour_manual(values = SCATTER_COLS,
                        labels = n_labels,
                        name   = NULL,
                        guide  = guide_legend(override.aes = list(size=2, alpha=1))) +
    scale_alpha_manual(values = SCATTER_ALPHA, guide = "none") +
    scale_size_manual(values  = SCATTER_SIZE,  guide = "none") +
    labs(
        x     = "MS logFC (combined vs HC)",
        y     = "ALS logFC (combined vs HC)",
        title = "Cross-disease logFC: MS vs ALS (all 2,911 proteins)"
    ) +
    theme_ukb()

ggsave(file.path(OUT_DIR, "ms_als_logfc_scatter.pdf"),
       p_scatter, width = 6, height = 5.5, device = cairo_pdf)
cat("  Saved: ms_als_logfc_scatter.pdf\n")


# ── Panel C: Pre-onset scatter ─────────────────────────────────────────────────
pre_scatter <- merge(
    ms_pre[,  .(protein, ms_logFC  = logFC, ms_sig  = sig)],
    als_pre[, .(protein, als_logFC = logFC, als_sig = sig)],
    by = "protein"
)

pre_scatter[, group := fcase(
    ms_sig & als_sig  & (ms_logFC > 0) == (als_logFC > 0), "Shared (concordant)",
    ms_sig & als_sig  & (ms_logFC > 0) != (als_logFC > 0), "Shared (discordant)",
    ms_sig & !als_sig, "MS pre-onset only",
    !ms_sig & als_sig, "ALS pre-onset only",
    default = "NS"
)]

pre_scatter[, label := fifelse(
    ms_sig | als_sig, toupper(protein), NA_character_
)]

PRE_COLS <- c(
    "Shared (concordant)"  = "#CC0066",
    "Shared (discordant)"  = "#E6A817",
    "MS pre-onset only"    = "#56B4E9",
    "ALS pre-onset only"   = "#2B4C7E",
    "NS"                   = "grey80"
)

pre_n_labels <- c(
    sprintf("Shared concordant n=%d",    sum(pre_scatter$group == "Shared (concordant)")),
    sprintf("Shared discordant n=%d",    sum(pre_scatter$group == "Shared (discordant)")),
    sprintf("MS pre-onset only n=%d",    sum(pre_scatter$group == "MS pre-onset only")),
    sprintf("ALS pre-onset only n=%d",   sum(pre_scatter$group == "ALS pre-onset only")),
    sprintf("NS n=%d",                   sum(pre_scatter$group == "NS"))
)

pre_scatter[, group_f := factor(group, levels = c("NS", "MS pre-onset only", "ALS pre-onset only",
                                                    "Shared (discordant)",
                                                    "Shared (concordant)"))]

PRE_SIZE  <- c("Shared (concordant)"=2.2, "Shared (discordant)"=2.2,
               "MS pre-onset only"=1.6,   "ALS pre-onset only"=1.6, "NS"=0.7)
PRE_ALPHA <- c("Shared (concordant)"=1.0, "Shared (discordant)"=1.0,
               "MS pre-onset only"=0.8,   "ALS pre-onset only"=0.8, "NS"=0.2)

p_pre <- ggplot(pre_scatter, aes(x = ms_logFC, y = als_logFC,
                                  colour = group_f, alpha = group_f, size = group_f)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50", linetype = "dashed") +
    geom_point(shape = 16) +
    geom_text_repel(
        aes(label = label),
        size = 2.8, colour = "black",
        max.overlaps = 40,
        box.padding  = 0.4,
        segment.size = 0.25,
        segment.colour = "grey40",
        na.rm = TRUE
    ) +
    scale_colour_manual(values = PRE_COLS,
                        labels = pre_n_labels,
                        name   = NULL,
                        guide  = guide_legend(override.aes = list(size=2.2, alpha=1))) +
    scale_alpha_manual(values = PRE_ALPHA, guide = "none") +
    scale_size_manual(values  = PRE_SIZE,  guide = "none") +
    labs(
        x     = "MS logFC (pre-onset vs HC)",
        y     = "ALS logFC (pre-onset vs HC)",
        title = "Pre-onset presymptomatic signals: MS vs ALS"
    ) +
    theme_ukb()

ggsave(file.path(OUT_DIR, "ms_als_preonset_scatter.pdf"),
       p_pre, width = 5.5, height = 5, device = cairo_pdf)
cat("  Saved: ms_als_preonset_scatter.pdf\n")


# ── Panel D: Venn-style upset counts bar chart ────────────────────────────────
# Simple count bars — MS-only / ALS-only / Shared, per contrast
venn_counts <- rbind(
    data.table(contrast="Combined",  group="MS only",  n=length(sets_comb$ms_only)),
    data.table(contrast="Combined",  group="ALS only", n=length(sets_comb$als_only)),
    data.table(contrast="Combined",  group="Shared",   n=length(sets_comb$shared)),
    data.table(contrast="Pre-onset", group="MS only",  n=length(sets_pre$ms_only)),
    data.table(contrast="Pre-onset", group="ALS only", n=length(sets_pre$als_only)),
    data.table(contrast="Pre-onset", group="Shared",   n=length(sets_pre$shared)),
    data.table(contrast="Post-onset",group="MS only",  n=length(sets_post$ms_only)),
    data.table(contrast="Post-onset",group="ALS only", n=length(sets_post$als_only)),
    data.table(contrast="Post-onset",group="Shared",   n=length(sets_post$shared))
)

venn_counts[, contrast_f := factor(contrast, levels = c("Pre-onset","Combined","Post-onset"))]
venn_counts[, group_f    := factor(group,    levels = c("ALS only","Shared","MS only"))]

VENN_COLS <- c("MS only" = "#56B4E9", "ALS only" = "#2B4C7E", "Shared" = "#CC0066")

p_venn <- ggplot(venn_counts, aes(x = contrast_f, y = n, fill = group_f)) +
    geom_col(position = "stack", width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = fifelse(n > 0, as.character(n), "")),
              position = position_stack(vjust = 0.5),
              size = 3, colour = "white", fontface = "bold") +
    scale_fill_manual(values = VENN_COLS, name = NULL) +
    labs(
        x     = NULL,
        y     = "Number of DEPs (FDR < 0.05)",
        title = "DEP overlap: MS vs ALS (by contrast)"
    ) +
    theme_ukb() +
    theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "ms_als_dep_counts.pdf"),
       p_venn, width = 4.5, height = 4, device = cairo_pdf)
cat("  Saved: ms_als_dep_counts.pdf\n")


# ── Combined figure (3-panel): counts + combined scatter + pre-onset scatter ──
p_combined <- (p_venn | p_scatter | p_pre) +
    plot_annotation(
        title   = "Cross-disease proteomics: Multiple Sclerosis vs ALS",
        subtitle = sprintf("MS: %d combined DEPs | ALS: %d combined DEPs | %d shared",
                           sum(ms_comb$sig), sum(als_comb$sig), length(sets_comb$shared)),
        theme   = theme(
            plot.title    = element_text(size = 11, face = "plain"),
            plot.subtitle = element_text(size = 9,  colour = "grey40")
        )
    )

ggsave(file.path(OUT_DIR, "ms_als_combined_figure.pdf"),
       p_combined, width = 15, height = 5.5, device = cairo_pdf)
cat("  Saved: ms_als_combined_figure.pdf\n")


# ── Direction concordance summary ─────────────────────────────────────────────
cat("\n=== Direction concordance (combined DEPs) ===\n")
concordant  <- shared_comb_dt[same_dir == TRUE]
discordant  <- shared_comb_dt[same_dir == FALSE]
cat(sprintf("  Concordant (same direction): %d\n", nrow(concordant)))
if (nrow(concordant) > 0) {
    cat("    Up in both:   ", paste(concordant[ms_dir=="up",  protein], collapse=", "), "\n")
    cat("    Down in both: ", paste(concordant[ms_dir=="down",protein], collapse=", "), "\n")
}
cat(sprintf("  Discordant (opposite direction): %d\n", nrow(discordant)))
if (nrow(discordant) > 0) {
    cat("    Proteins: ", paste(discordant$protein, collapse=", "), "\n")
    print(discordant[, .(protein, ms_logFC=round(ms_logFC,3), als_logFC=round(als_logFC,3), ms_dir, als_dir)])
}

cat("\n01_ms_als_comparison.R complete.\n")
