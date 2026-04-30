#!/usr/bin/env Rscript
# 04_prepost_go.R
# GO BP enrichment for pre-onset and post-onset MS DEPs separately.
# Used by figure1.R panels e (pre-onset GO dotplot) and i (pre vs post pathway comparison).
#
# Pre-onset uses FDR < 0.20 (liberal threshold appropriate for low-power pre-diagnosis window).
# Post-onset uses FDR < 0.05 (standard threshold; well-powered).
#
# Outputs:
#   results/networks/ms_pre_go_results.csv
#   results/networks/ms_post_go_results.csv

suppressPackageStartupMessages({
    library(data.table)
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "go_dotplot.R"))

DIFF_DIR <- file.path(PROJ_DIR, "results", "differential")
OUT_DIR  <- file.path(PROJ_DIR, "results", "networks")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Thresholds
PRE_FDR  <- 0.20   # pre-onset: low power, use liberal FDR
POST_FDR <- 0.05   # post-onset: well powered

# ---------------------------------------------------------------------------
# Pre-onset DEPs
# ---------------------------------------------------------------------------
cat("=== Pre-onset DEPs (FDR <", PRE_FDR, ") ===\n")
pre_all <- fread(file.path(DIFF_DIR, "ms_pre_vs_hc.csv"))
pre_sig <- pre_all[adj.P.Val < PRE_FDR]
cat(sprintf("  Pre-onset DEPs: %d total (%d up, %d down)\n",
            nrow(pre_sig),
            sum(pre_sig$logFC > 0), sum(pre_sig$logFC < 0)))

bg_pre    <- toupper(pre_all$protein)
pre_up    <- toupper(pre_sig[logFC > 0, protein])
pre_down  <- toupper(pre_sig[logFC < 0, protein])

ego_pre_up   <- run_go(pre_up,   "Pre-onset up")
ego_pre_down <- run_go(pre_down, "Pre-onset down")

go_pre_dt <- rbindlist(Filter(Negate(is.null), list(
    collect_go(ego_pre_up,   "up_DEPs"),
    collect_go(ego_pre_down, "down_DEPs")
)))
if (nrow(go_pre_dt) > 0) {
    fwrite(go_pre_dt, file.path(OUT_DIR, "ms_pre_go_results.csv"))
    cat(sprintf("  Saved ms_pre_go_results.csv (%d rows)\n", nrow(go_pre_dt)))
} else {
    cat("  No significant GO terms — writing empty file\n")
    fwrite(data.table(), file.path(OUT_DIR, "ms_pre_go_results.csv"))
}

# ---------------------------------------------------------------------------
# Post-onset DEPs
# ---------------------------------------------------------------------------
cat("\n=== Post-onset DEPs (FDR <", POST_FDR, ") ===\n")
post_all <- fread(file.path(DIFF_DIR, "ms_post_vs_hc.csv"))
post_sig <- post_all[adj.P.Val < POST_FDR]
cat(sprintf("  Post-onset DEPs: %d total (%d up, %d down)\n",
            nrow(post_sig),
            sum(post_sig$logFC > 0), sum(post_sig$logFC < 0)))

bg_post    <- toupper(post_all$protein)
post_up    <- toupper(post_sig[logFC > 0, protein])
post_down  <- toupper(post_sig[logFC < 0, protein])

ego_post_up   <- run_go(post_up,   "Post-onset up")
ego_post_down <- run_go(post_down, "Post-onset down")

go_post_dt <- rbindlist(Filter(Negate(is.null), list(
    collect_go(ego_post_up,   "up_DEPs"),
    collect_go(ego_post_down, "down_DEPs")
)))
if (nrow(go_post_dt) > 0) {
    fwrite(go_post_dt, file.path(OUT_DIR, "ms_post_go_results.csv"))
    cat(sprintf("  Saved ms_post_go_results.csv (%d rows)\n", nrow(go_post_dt)))
} else {
    cat("  No significant GO terms — writing empty file\n")
    fwrite(data.table(), file.path(OUT_DIR, "ms_post_go_results.csv"))
}

# ---------------------------------------------------------------------------
# Combined MS vs HC DEPs (for supplementary figure)
# ---------------------------------------------------------------------------
cat("\n=== Combined MS vs HC DEPs (FDR < 0.05) ===\n")
comb_all <- fread(file.path(DIFF_DIR, "ms_combined_vs_hc.csv"))
comb_sig <- comb_all[adj.P.Val < 0.05]
cat(sprintf("  Combined DEPs: %d total (%d up, %d down)\n",
            nrow(comb_sig),
            sum(comb_sig$logFC > 0), sum(comb_sig$logFC < 0)))

bg_comb    <- toupper(comb_all$protein)
comb_up    <- toupper(comb_sig[logFC > 0, protein])
comb_down  <- toupper(comb_sig[logFC < 0, protein])

ego_comb_up   <- run_go(comb_up,   bg_comb, "Combined up")
ego_comb_down <- run_go(comb_down, bg_comb, "Combined down")

go_comb_dt <- rbindlist(Filter(Negate(is.null), list(
    collect_go(ego_comb_up,   "up_DEPs"),
    collect_go(ego_comb_down, "down_DEPs")
)))
if (nrow(go_comb_dt) > 0) {
    fwrite(go_comb_dt, file.path(OUT_DIR, "ms_combined_go_results.csv"))
    cat(sprintf("  Saved ms_combined_go_results.csv (%d rows)\n", nrow(go_comb_dt)))
} else {
    cat("  No significant GO terms — writing empty file\n")
    fwrite(data.table(), file.path(OUT_DIR, "ms_combined_go_results.csv"))
}

cat("\nDone.\n")
