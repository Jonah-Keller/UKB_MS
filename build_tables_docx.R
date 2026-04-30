#!/usr/bin/env Rscript
# ============================================================================
# build_tables_docx.R — Emit each registered MAIN table as a standalone .docx
# ============================================================================
# Stages source CSVs from results/ into _build_intermediates/ (per
# `source_csv:` in figures.yaml), then re-renders each main table with house
# styling into a standalone .docx (main_tables/Table_N_<slug>.docx).
#
# Supplementary tables are NOT emitted here — they are packed into
# Supplementary_Tables.xlsx by build_supp_xlsx.R (called from
# build_supplementary.R). One sheet per supp table, plus a Contents sheet.
#
# Called by build_submission.R as the penultimate step. Standalone usage:
# `Rscript build_tables_docx.R`.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(flextable)
  library(officer)
})

source(here("code", "config", "00_config.R"))
source(here("code", "registry", "figure_registry.R"))
source(here("code", "tables", "style_flextable.R"))
source(here("code", "tables", "stage_table_csvs.R"))

registry_load(FIGURES_YAML)

cat("Staging table CSVs...\n")
stage_table_csvs()

cat("Building table .docx files...\n")
for (t in .REG$tables) {
  if (!identical(t$section, "main")) next   # supp tables ship in the xlsx
  csv  <- here(t$output$csv)
  docx <- here(t$output$docx)
  if (!file.exists(csv)) {
    cat("  \u26a0  skipping", t$id, "\u2014 CSV not found:", csv, "\n")
    next
  }
  dir.create(dirname(docx), showWarnings = FALSE, recursive = TRUE)
  d <- read.csv(csv, check.names = FALSE)
  ft <- flextable(d) |>
    style_flextable(caption = tbl_caption(t$id)) |>
    align(j = 2:ncol(d), align = "center", part = "body") |>
    align(j = 2:ncol(d), align = "center", part = "header")
  save_flextable_docx(ft, docx, orient = "landscape")
  cat("  \u2192", basename(docx), "\n")
}
cat("Done.\n")
