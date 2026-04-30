#!/usr/bin/env Rscript
# 02_temporal_go_enrichment.R
# Temporal GO enrichment: which biological processes activate at each
# time window before MS diagnosis?
#
# Method:
#   1. Read rolling Cox HR table (produced by 01_ms_protein_cox.R)
#   2. For each time bin, select proteins with HR > 1 and p < 0.05
#      (significant RISK-increasing associations)
#   3. Submit each protein list to gprofiler2 (GO:BP + KEGG + REAC)
#   4. Visualise as a heatmap: term × time bin, coloured by -log10(p_value)
#
# Output: results/survival/
#   ms_temporal_go_enrichment.csv
#   ms_temporal_go_heatmap.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(gprofiler2)
    library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "disease_config.R"))
suppressPackageStartupMessages({
    library(glue)
    library(here)
})

cfg     <- load_disease_config()
COHORT  <- cfg$cohort_short
DISEASE <- cfg$disease_short_caps

ROLLING_FILE <- file.path(PROJ_DIR, "results", "survival",
                          glue("{COHORT}_protein_cox_rolling.csv"))
OUT_DIR      <- file.path(PROJ_DIR, "results", "survival")

# Bail gracefully when upstream Cox stage hasn't run for this cohort.
if (!file.exists(ROLLING_FILE)) {
    cat(sprintf("Temporal GO stage skipped: %s not found.\n",
                basename(ROLLING_FILE)))
    cat("  Run analysis/13_survival/01_ms_protein_cox.R first.\n")
    quit(save = "no", status = 0)
}

roll <- fread(ROLLING_FILE)
if (nrow(roll) == 0L) {
    cat(sprintf("Temporal GO stage skipped: %s is empty (no time bins\n",
                basename(ROLLING_FILE)))
    cat("  reached the survival_min_events_per_bin threshold).\n")
    quit(save = "no", status = 0)
}
# Standardise column names
if ("protein" %in% names(roll)) roll[, protein := toupper(protein)]

cat("Rolling HR table:\n")
print(roll[, .N, by = time_bin])

# Bin order discovered from the rolling Cox output (so cohort-specific
# bins from cfg$ytd_bins flow through automatically).
BIN_ORDER <- unique(roll$time_bin)
roll[, bin_f := factor(time_bin, levels = BIN_ORDER)]

# ── For each bin: risk-associated proteins (HR>1, p<0.05) ────────────────────
cat("\nProteins per time bin (HR>1, p<0.05):\n")
go_results <- list()

# Use full Cox protein list as background (all proteins tested in Cox analysis)
COX_ALL_FILE <- file.path(PROJ_DIR, "results", "survival",
                          glue("{COHORT}_protein_cox_results.csv"))
bg_prots <- if (file.exists(COX_ALL_FILE)) {
    toupper(fread(COX_ALL_FILE)$protein)
} else NULL
cat(sprintf("Background: %d proteins\n", length(bg_prots)))

for (bn in BIN_ORDER) {
    # Broader threshold (p<0.10, HR>1) to capture pre-onset signal
    sub <- roll[bin_f == bn & HR > 1 & pval < 0.10, protein]
    cat(sprintf("  %s: %d proteins — %s\n", bn, length(sub),
                paste(head(sub, 8), collapse=", ")))
    if (length(sub) < 2) next

    res <- tryCatch(
        gost(query             = sub,
             custom_bg         = bg_prots,   # background = full Cox protein list
             organism          = "hsapiens",
             sources           = c("GO:BP"),
             significant       = FALSE,
             evcodes           = FALSE,
             user_threshold    = 0.5,         # permissive; filter below
             correction_method = "fdr"),      # standard FDR over g_SCS for small lists
        error = function(e) { message("  gprofiler2 failed: ", e$message); NULL }
    )
    if (is.null(res) || is.null(res$result)) next
    dt <- as.data.table(res$result)
    dt <- dt[p_value < 0.2]   # keep FDR<0.2 (exploratory for small lists)
    if (nrow(dt) == 0) { cat("    No terms at FDR<0.20\n"); next }
    dt[, time_bin := bn]
    go_results[[bn]] <- dt
}

if (length(go_results) == 0) {
    cat("No enrichment results — protein lists too small for GO ORA.\n")
    cat("Writing empty placeholder and producing Cox-only heatmap.\n")
    # Fall back to just the rolling HR heatmap (no GO annotation)
    fwrite(data.table(),
           file.path(OUT_DIR, glue("{COHORT}_temporal_go_enrichment.csv")))
    cat("\n02_temporal_go_enrichment.R complete (no GO terms found).\n")
    quit(save = "no", status = 0)
}

all_go <- rbindlist(go_results, fill = TRUE)
fwrite(all_go,
       file.path(OUT_DIR, glue("{COHORT}_temporal_go_enrichment.csv")))
cat(sprintf("\nTotal enriched terms: %d\n", nrow(all_go)))

# ── Select top terms for display ──────────────────────────────────────────────
# Top 5 GO:BP terms per bin, ordered by significance
top_terms <- all_go[source == "GO:BP"
                   ][order(p_value),
                   ][, .SD[seq_len(min(5, .N))], by = time_bin]

# Shorten long GO term names
shorten <- function(s, n = 45) ifelse(nchar(s) > n, paste0(substr(s, 1, n), "…"), s)
top_terms[, term_short := shorten(term_name)]

# All unique terms across bins
all_terms <- unique(top_terms$term_short)

# For heatmap: fill in -log10(p_adj) for every term × bin combination
heat_grid <- CJ(term_short = all_terms, time_bin = BIN_ORDER)
heat_grid <- merge(heat_grid,
                   top_terms[, .(term_short, time_bin, p_value)],
                   by = c("term_short","time_bin"), all.x = TRUE)
heat_grid[, neg_log10_p := ifelse(is.na(p_value), 0, -log10(p_value))]
heat_grid[, bin_f := factor(time_bin, levels = BIN_ORDER)]

# Order terms by earliest bin of peak significance
term_order <- top_terms[, .(min_bin = min(as.integer(factor(time_bin, levels=BIN_ORDER)))),
                          by = term_short][order(-min_bin), term_short]
heat_grid[, term_f := factor(term_short, levels = term_order)]

p_heat <- ggplot(heat_grid, aes(x = bin_f, y = term_f, fill = neg_log10_p)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(neg_log10_p > 0,
                                  sprintf("%.1f", neg_log10_p), "")),
              size = 2.3, colour = "white") +
    scale_fill_gradient(low = "grey92", high = "#CC0066",
                        name = expression(-log[10]~italic(p)[adj]),
                        guide = guide_colourbar(barwidth = 0.6, barheight = 4)) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(title    = glue("Temporal GO enrichment of pre-{DISEASE} protein changes"),
         subtitle = "GO:BP terms enriched in risk-associated proteins (HR>1, p<0.05) per time window",
         x = NULL, y = NULL) +
    theme_ukb(base_size = 9) +
    theme(axis.text.y    = element_text(size = 7.5),
          axis.text.x    = element_text(size = 8),
          panel.grid     = element_blank(),
          legend.position = "right")

h_fig <- max(5, length(all_terms) * 0.32 + 2)
ggsave(file.path(OUT_DIR, glue("{COHORT}_temporal_go_heatmap.pdf")),
       p_heat, width = 7.5, height = h_fig, device = cairo_pdf)
cat(sprintf("  Saved: %s_temporal_go_heatmap.pdf\n", COHORT))

cat("\n02_temporal_go_enrichment.R complete.\n")
