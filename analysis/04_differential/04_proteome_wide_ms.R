#!/usr/bin/env Rscript
# 04_proteome_wide_ms.R
# Proteome-wide visualisation of full MS vs HC results (all 2,911 proteins)
#
# Context:
#   01_limma_ms_vs_hc.R already tested ALL proteins and saved results to
#   ms_combined_vs_hc.csv (2,911 rows, 173 at FDR<0.05).
#   This script produces the full proteome-wide plots:
#     - Enhanced volcano: all proteins, FDR<0.05 and Bonferroni highlighted
#     - Effect-size rank plot (Manhattan-style, sorted by logFC magnitude)
#     - Proteins significant at FDR<0.05 vs 0.05–0.20 ("suggestive") comparison
#     - Full results table with Bonferroni flag
#
# Output: results/differential/
#   ms_proteome_volcano_full.pdf
#   ms_proteome_rankplot.pdf
#   ms_proteome_wide_results.csv   (all 2,911 proteins, annotated)

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

DIFF_FILE <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
OUT_DIR   <- file.path(PROJ_DIR, "results", "differential")

# ── 1. Load full proteome results ─────────────────────────────────────────────
cat("Loading proteome-wide results...\n")
dt <- fread(DIFF_FILE, showProgress = FALSE)
setnames(dt, "adj.P.Val", "fdr", skip_absent = TRUE)
cat(sprintf("  %d proteins tested\n", nrow(dt)))

# Bonferroni threshold
BONF_THR <- 0.05 / nrow(dt)
cat(sprintf("  Bonferroni threshold: p < %.2e\n", BONF_THR))
cat(sprintf("  FDR<0.05:     %d\n", sum(dt$fdr < 0.05,    na.rm = TRUE)))
cat(sprintf("  FDR 0.05–0.20: %d (suggestive)\n",
            sum(dt$fdr >= 0.05 & dt$fdr < 0.20, na.rm = TRUE)))
cat(sprintf("  Bonferroni:   %d\n", sum(dt$P.Value < BONF_THR, na.rm = TRUE)))

# Significance tier
dt[, sig_tier := fcase(
    P.Value < BONF_THR,              "Bonferroni",
    fdr     < 0.05,                  "FDR<0.05",
    fdr     < 0.20,                  "Suggestive (FDR<0.20)",
    default =                         "NS"
)]
dt[, neg_log10p := -log10(P.Value)]

# Labels: top proteins per tier
dt[, label := NA_character_]
dt[sig_tier == "Bonferroni", label := toupper(protein)]
n_fdr <- sum(dt$sig_tier == "FDR<0.05")
if (n_fdr <= 30) {
    dt[sig_tier == "FDR<0.05", label := toupper(protein)]
} else {
    # label top 20 by effect size within FDR tier
    top_fdr <- dt[sig_tier == "FDR<0.05"][order(-abs(logFC))][1:20, protein]
    dt[protein %in% top_fdr, label := toupper(protein)]
}

# Colour scale
TIER_COLS <- c(
    "Bonferroni"              = "#1A1A1A",
    "FDR<0.05"                = "#CC0066",
    "Suggestive (FDR<0.20)"  = "#E6A817",
    "NS"                      = "grey75"
)
TIER_ALPHA <- c("Bonferroni"=1, "FDR<0.05"=0.8, "Suggestive (FDR<0.20)"=0.6, "NS"=0.25)
TIER_SIZE  <- c("Bonferroni"=2, "FDR<0.05"=1.5, "Suggestive (FDR<0.20)"=1.2, "NS"=0.8)

# ── 2. Full volcano ────────────────────────────────────────────────────────────
cat("Plotting full proteome volcano...\n")

p_vol <- ggplot(dt, aes(x = logFC, y = neg_log10p, colour = sig_tier,
                         alpha = sig_tier, size = sig_tier)) +
    geom_hline(yintercept = -log10(BONF_THR), linetype = "dotted",
               colour = "#1A1A1A", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05 / 1),   # arbitrary FDR line position
               linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    annotate("text", x = max(dt$logFC, na.rm=TRUE) * 0.9,
             y = -log10(BONF_THR) + 0.3, label = "Bonferroni",
             size = 2.5, colour = "#1A1A1A", hjust = 1) +
    geom_point() +
    geom_text_repel(aes(label = label), size = 2.0, na.rm = TRUE,
                    max.overlaps = 30, segment.colour = "grey50",
                    segment.size = 0.3, min.segment.length = 0.2) +
    scale_colour_manual(values = TIER_COLS, name = NULL) +
    scale_alpha_manual(values = TIER_ALPHA, guide = "none") +
    scale_size_manual(values = TIER_SIZE, guide = "none") +
    labs(x = "logFC (MS vs HC)",
         y = expression(-log[10](p)),
         title = sprintf("Proteome-wide MS vs HC: %d proteins", nrow(dt)),
         subtitle = sprintf("Bonferroni: %d | FDR<0.05: %d | Suggestive FDR<0.20: %d",
                            sum(dt$P.Value < BONF_THR, na.rm=TRUE),
                            sum(dt$fdr < 0.05, na.rm=TRUE),
                            sum(dt$fdr >= 0.05 & dt$fdr < 0.20, na.rm=TRUE))) +
    theme_ukb() +
    theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "ms_proteome_volcano_full.pdf"), p_vol,
       width = 7, height = 6, device = cairo_pdf)
cat("  Saved: ms_proteome_volcano_full.pdf\n")

# ── 3. Effect-size rank plot ───────────────────────────────────────────────────
dt_rank <- copy(dt)
dt_rank <- dt_rank[order(-logFC)]
dt_rank[, rank := seq_len(.N)]
dt_rank[, rank_label := fcase(
    rank <= 10 | rank >= (.N - 9), toupper(protein),
    sig_tier == "Bonferroni",       toupper(protein),
    default = NA_character_
)]

p_rank <- ggplot(dt_rank, aes(x = rank, y = logFC, colour = sig_tier, size = sig_tier)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_text_repel(aes(label = rank_label), size = 2.0, na.rm = TRUE,
                    max.overlaps = 20, segment.size = 0.3) +
    scale_colour_manual(values = TIER_COLS, name = NULL) +
    scale_size_manual(values = TIER_SIZE, guide = "none") +
    labs(x = "Rank (by logFC, high → low)",
         y = "logFC (MS vs HC)",
         title = "Proteome-wide effect size rank plot",
         subtitle = "All 2,911 tested proteins") +
    theme_ukb() +
    theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "ms_proteome_rankplot.pdf"), p_rank,
       width = 7, height = 5, device = cairo_pdf)
cat("  Saved: ms_proteome_rankplot.pdf\n")

# ── 4. Suggestive hits table ───────────────────────────────────────────────────
suggestive <- dt[sig_tier == "Suggestive (FDR<0.20)"][order(fdr)]
cat(sprintf("\nSuggestive hits (FDR 0.05–0.20), n=%d:\n", nrow(suggestive)))
print(suggestive[1:min(20, .N),
                  .(protein = toupper(protein), logFC = round(logFC,3),
                    P.Value = signif(P.Value,3), fdr = round(fdr,3))])

# ── 5. Save annotated full results ─────────────────────────────────────────────
dt[, bonferroni_sig := P.Value < BONF_THR]
dt[, bonferroni_threshold := BONF_THR]
fwrite(dt[, .(protein, logFC, AveExpr, t, P.Value, fdr, B,
              sig_tier, bonferroni_sig, bonferroni_threshold)],
       file.path(OUT_DIR, "ms_proteome_wide_results.csv"))
cat(sprintf("\nFull results written: %d proteins\n", nrow(dt)))

cat("\n04_proteome_wide_ms.R complete.\n")
