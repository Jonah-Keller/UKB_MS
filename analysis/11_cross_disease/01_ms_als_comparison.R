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
    library(glue)
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
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

PRIMARY_SHORT <- cfg$cohort_short
SECOND_SHORT  <- cfg$comparison_cohort_short
if (is.null(SECOND_SHORT) || !nzchar(SECOND_SHORT)) {
    stop("comparison_cohort_short is empty in configs/disease.yaml; ",
         "set it (e.g., \"als\") to enable cross-disease comparison.")
}
PRIMARY_CAPS <- cfg$disease_short_caps
SECOND_CAPS  <- toupper(SECOND_SHORT)

DIFF_DIR   <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR    <- file.path(PROJ_DIR, "results", "cross_disease")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR <- 0.05

# ── Load all six DEP result files ─────────────────────────────────────────────
load_deps <- function(disease, contrast) {
    fname <- file.path(DIFF_DIR, glue("{disease}_{contrast}_vs_hc.csv"))
    dt <- fread(fname)
    dt[, disease   := disease]
    dt[, contrast  := contrast]
    dt[, sig       := adj.P.Val < FDR_THR]
    dt[, direction := fifelse(sig & logFC > 0, "up",
                      fifelse(sig & logFC < 0, "down", "ns"))]
    dt
}

ms_comb  <- load_deps(PRIMARY_SHORT, "combined")
ms_pre   <- load_deps(PRIMARY_SHORT, "pre")
ms_post  <- load_deps(PRIMARY_SHORT, "post")
als_comb <- load_deps(SECOND_SHORT,  "combined")
als_pre  <- load_deps(SECOND_SHORT,  "pre")
als_post <- load_deps(SECOND_SHORT,  "post")

cat(sprintf("%s combined DEPs: %d  |  %s combined DEPs: %d\n",
            PRIMARY_CAPS, sum(ms_comb$sig),  SECOND_CAPS, sum(als_comb$sig)))
cat(sprintf("%s pre-onset DEPs: %d |  %s pre-onset DEPs: %d\n",
            PRIMARY_CAPS, sum(ms_pre$sig),   SECOND_CAPS, sum(als_pre$sig)))
cat(sprintf("%s post-onset DEPs: %d | %s post-onset DEPs: %d\n",
            PRIMARY_CAPS, sum(ms_post$sig),  SECOND_CAPS, sum(als_post$sig)))


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

print_overlap <- function(label, sets) {
    cat(sprintf("\n=== %s contrast ===\n", label))
    cat(sprintf("  %s-only: %d  |  %s-only: %d  |  Shared: %d\n",
                PRIMARY_CAPS, length(sets$ms_only),
                SECOND_CAPS,  length(sets$als_only),
                length(sets$shared)))
    cat("  Shared proteins:", paste(sort(sets$shared), collapse=", "), "\n")
}
print_overlap("Combined",   sets_comb)
print_overlap("Pre-onset",  sets_pre)
print_overlap("Post-onset", sets_post)


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
PAIR_PREFIX <- glue("{PRIMARY_SHORT}_{SECOND_SHORT}")
fwrite(all_shared, file.path(OUT_DIR, glue("{PAIR_PREFIX}_shared_proteins.csv")))
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
LBL_PRIM_ONLY   <- glue("{PRIMARY_CAPS} only")
LBL_SECOND_ONLY <- glue("{SECOND_CAPS} only")
scatter_dt[, group := fcase(
    ms_sig == TRUE  & als_sig == TRUE  & (ms_logFC > 0) == (als_logFC > 0), "Shared (concordant)",
    ms_sig == TRUE  & als_sig == TRUE  & (ms_logFC > 0) != (als_logFC > 0), "Shared (discordant)",
    ms_sig == TRUE  & als_sig == FALSE, LBL_PRIM_ONLY,
    ms_sig == FALSE & als_sig == TRUE,  LBL_SECOND_ONLY,
    default = "NS"
)]

scatter_dt[, label := fifelse(
    group %in% c("Shared (concordant)", "Shared (discordant)"), toupper(protein), NA_character_
)]

SCATTER_COLS <- setNames(
    c("#CC0066", "#E6A817", "#56B4E9", "#2B4C7E", "grey80"),
    c("Shared (concordant)", "Shared (discordant)",
      LBL_PRIM_ONLY, LBL_SECOND_ONLY, "NS")
)
SCATTER_ALPHA <- setNames(
    c(0.95, 0.95, 0.7, 0.7, 0.25),
    c("Shared (concordant)", "Shared (discordant)",
      LBL_PRIM_ONLY, LBL_SECOND_ONLY, "NS")
)
SCATTER_SIZE <- setNames(
    c(2.0, 2.0, 1.4, 1.4, 0.7),
    c("Shared (concordant)", "Shared (discordant)",
      LBL_PRIM_ONLY, LBL_SECOND_ONLY, "NS")
)

# Force NS to bottom layer
scatter_dt[, group_f := factor(group, levels = c("NS", LBL_PRIM_ONLY, LBL_SECOND_ONLY,
                                                   "Shared (discordant)",
                                                   "Shared (concordant)"))]

n_labels <- c(
    sprintf("Shared concordant n=%d",  sum(scatter_dt$group == "Shared (concordant)")),
    sprintf("Shared discordant n=%d",  sum(scatter_dt$group == "Shared (discordant)")),
    sprintf("%s n=%d", LBL_PRIM_ONLY,   sum(scatter_dt$group == LBL_PRIM_ONLY)),
    sprintf("%s n=%d", LBL_SECOND_ONLY, sum(scatter_dt$group == LBL_SECOND_ONLY)),
    sprintf("NS n=%d",                  sum(scatter_dt$group == "NS"))
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
        x     = glue("{PRIMARY_CAPS} logFC (combined vs HC)"),
        y     = glue("{SECOND_CAPS} logFC (combined vs HC)"),
        title = glue("Cross-disease logFC: {PRIMARY_CAPS} vs {SECOND_CAPS} (all 2,911 proteins)")
    ) +
    theme_ukb()

SCATTER_PDF <- glue("{PAIR_PREFIX}_logfc_scatter.pdf")
ggsave(file.path(OUT_DIR, SCATTER_PDF),
       p_scatter, width = 6, height = 5.5, device = cairo_pdf)
cat(sprintf("  Saved: %s\n", SCATTER_PDF))


# ── Panel C: Pre-onset scatter ─────────────────────────────────────────────────
pre_scatter <- merge(
    ms_pre[,  .(protein, ms_logFC  = logFC, ms_sig  = sig)],
    als_pre[, .(protein, als_logFC = logFC, als_sig = sig)],
    by = "protein"
)

PRE_PRIM   <- glue("{PRIMARY_CAPS} pre-onset only")
PRE_SECOND <- glue("{SECOND_CAPS} pre-onset only")
pre_scatter[, group := fcase(
    ms_sig & als_sig  & (ms_logFC > 0) == (als_logFC > 0), "Shared (concordant)",
    ms_sig & als_sig  & (ms_logFC > 0) != (als_logFC > 0), "Shared (discordant)",
    ms_sig & !als_sig, PRE_PRIM,
    !ms_sig & als_sig, PRE_SECOND,
    default = "NS"
)]

pre_scatter[, label := fifelse(
    ms_sig | als_sig, toupper(protein), NA_character_
)]

PRE_COLS <- setNames(
    c("#CC0066", "#E6A817", "#56B4E9", "#2B4C7E", "grey80"),
    c("Shared (concordant)", "Shared (discordant)", PRE_PRIM, PRE_SECOND, "NS")
)

pre_n_labels <- c(
    sprintf("Shared concordant n=%d",  sum(pre_scatter$group == "Shared (concordant)")),
    sprintf("Shared discordant n=%d",  sum(pre_scatter$group == "Shared (discordant)")),
    sprintf("%s n=%d", PRE_PRIM,        sum(pre_scatter$group == PRE_PRIM)),
    sprintf("%s n=%d", PRE_SECOND,      sum(pre_scatter$group == PRE_SECOND)),
    sprintf("NS n=%d",                  sum(pre_scatter$group == "NS"))
)

pre_scatter[, group_f := factor(group, levels = c("NS", PRE_PRIM, PRE_SECOND,
                                                    "Shared (discordant)",
                                                    "Shared (concordant)"))]

PRE_SIZE  <- setNames(c(2.2, 2.2, 1.6, 1.6, 0.7),
                      c("Shared (concordant)","Shared (discordant)", PRE_PRIM, PRE_SECOND, "NS"))
PRE_ALPHA <- setNames(c(1.0, 1.0, 0.8, 0.8, 0.2),
                      c("Shared (concordant)","Shared (discordant)", PRE_PRIM, PRE_SECOND, "NS"))

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
        x     = glue("{PRIMARY_CAPS} logFC (pre-onset vs HC)"),
        y     = glue("{SECOND_CAPS} logFC (pre-onset vs HC)"),
        title = glue("Pre-onset presymptomatic signals: {PRIMARY_CAPS} vs {SECOND_CAPS}")
    ) +
    theme_ukb()

PRE_PDF <- glue("{PAIR_PREFIX}_preonset_scatter.pdf")
ggsave(file.path(OUT_DIR, PRE_PDF),
       p_pre, width = 5.5, height = 5, device = cairo_pdf)
cat(sprintf("  Saved: %s\n", PRE_PDF))


# ── Panel D: Venn-style upset counts bar chart ────────────────────────────────
contrast_sets <- list(
    Combined     = sets_comb,
    `Pre-onset`  = sets_pre,
    `Post-onset` = sets_post
)
venn_counts <- rbindlist(lapply(names(contrast_sets), function(label) {
    s <- contrast_sets[[label]]
    data.table(
        contrast = label,
        group    = c(LBL_PRIM_ONLY, LBL_SECOND_ONLY, "Shared"),
        n        = c(length(s$ms_only), length(s$als_only), length(s$shared))
    )
}))

venn_counts[, contrast_f := factor(contrast, levels = c("Pre-onset","Combined","Post-onset"))]
venn_counts[, group_f    := factor(group,    levels = c(LBL_SECOND_ONLY, "Shared", LBL_PRIM_ONLY))]

VENN_COLS <- setNames(c("#56B4E9", "#2B4C7E", "#CC0066"),
                      c(LBL_PRIM_ONLY, LBL_SECOND_ONLY, "Shared"))

p_venn <- ggplot(venn_counts, aes(x = contrast_f, y = n, fill = group_f)) +
    geom_col(position = "stack", width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = fifelse(n > 0, as.character(n), "")),
              position = position_stack(vjust = 0.5),
              size = 3, colour = "white", fontface = "bold") +
    scale_fill_manual(values = VENN_COLS, name = NULL) +
    labs(
        x     = NULL,
        y     = "Number of DEPs (FDR < 0.05)",
        title = glue("DEP overlap: {PRIMARY_CAPS} vs {SECOND_CAPS} (by contrast)")
    ) +
    theme_ukb() +
    theme(legend.position = "right")

VENN_PDF <- glue("{PAIR_PREFIX}_dep_counts.pdf")
ggsave(file.path(OUT_DIR, VENN_PDF),
       p_venn, width = 4.5, height = 4, device = cairo_pdf)
cat(sprintf("  Saved: %s\n", VENN_PDF))


# ── Combined figure (3-panel): counts + combined scatter + pre-onset scatter ──
p_combined <- (p_venn | p_scatter | p_pre) +
    plot_annotation(
        title   = glue("Cross-disease proteomics: {cfg$disease_long} vs {SECOND_CAPS}"),
        subtitle = sprintf("%s: %d combined DEPs | %s: %d combined DEPs | %d shared",
                           PRIMARY_CAPS, sum(ms_comb$sig),
                           SECOND_CAPS,  sum(als_comb$sig),
                           length(sets_comb$shared)),
        theme   = theme(
            plot.title    = element_text(size = 11, face = "plain"),
            plot.subtitle = element_text(size = 9,  colour = "grey40")
        )
    )

COMBINED_PDF <- glue("{PAIR_PREFIX}_combined_figure.pdf")
ggsave(file.path(OUT_DIR, COMBINED_PDF),
       p_combined, width = 15, height = 5.5, device = cairo_pdf)
cat(sprintf("  Saved: %s\n", COMBINED_PDF))


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
