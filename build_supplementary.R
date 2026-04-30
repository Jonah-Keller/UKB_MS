#!/usr/bin/env Rscript
# ============================================================================
# build_supplementary.R — Build Supplementary_Material.docx from registry
# ============================================================================
# Contents: title page + contents listing of Supplementary Tables (with a
# pointer to Supplementary_Tables.xlsx, where the actual table data lives)
# + Extended Data Figures (PDFs converted to PNGs and embedded).
#
# Supplementary tables are NOT inlined as flextables — flextable rendering of
# 3k+-row tables is unworkable, and a 17k-row Word table is unreadable anyway.
# All supp tables are packed into a single .xlsx (one sheet per table) by
# build_supp_xlsx().
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(officer)
  library(flextable)
})

source(here("code", "config", "00_config.R"))
source(here("code", "registry", "figure_registry.R"))
source(here("code", "tables", "style_flextable.R"))
source(here("code", "tables", "stage_table_csvs.R"))
source(here("code", "tables", "build_supp_xlsx.R"))

registry_load(FIGURES_YAML)

cat("Staging supplementary table CSVs...\n")
stage_table_csvs()

cat("Building Supplementary_Material.docx...\n")

fp_title <- fp_text(bold = TRUE,  font.size = 14, font.family = "Arial")
fp_bold  <- fp_text(bold = TRUE,  font.size = 12, font.family = "Arial")
fp_norm  <- fp_text(bold = FALSE, font.size = 12, font.family = "Arial")

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size    = page_size(orient = "portrait"),
  page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1)))

# Title page
doc <- body_add_fpar(doc, fpar(ftext("Supplementary Material", prop = fp_title)))
doc <- body_add_par(doc, "")
doc <- body_add_fpar(doc, fpar(ftext(PROJECT_TITLE, prop = fp_bold)))
doc <- body_add_par(doc, "")
doc <- body_add_fpar(doc, fpar(ftext(paste(PROJECT_AUTHORS, collapse = ", "),
                                     prop = fp_norm)))
doc <- body_add_par(doc, "")

# Contents listing
doc <- body_add_fpar(doc, fpar(ftext("Contents", prop = fp_bold)))
doc <- body_add_par(doc, "")
for (t in .REG$tables) {
  if (t$section != "supplement") next
  doc <- body_add_fpar(doc, fpar(
    ftext(paste0(tbl(t$id), ". "), prop = fp_bold),
    ftext(trimws(t$caption), prop = fp_norm)))
}
doc <- body_add_par(doc, "")
for (f in .REG$figures) {
  if (f$section != "supplement") next
  doc <- add_legend_fpar(doc, paste0(fig(f$id), ". "),
                         split_legend_panels(f$caption), fp_bold, fp_norm)
}
doc <- body_add_break(doc, pos = "after")

# Supplementary Tables — caption-only listing; data lives in the xlsx
doc <- body_add_fpar(doc, fpar(ftext("Supplementary Tables", prop = fp_bold)))
doc <- body_add_par(doc, "")
doc <- body_add_fpar(doc, fpar(ftext(
  paste0("All supplementary tables are provided as sheets in ",
         "Supplementary_Tables.xlsx (one sheet per table; see the Contents ",
         "sheet for the index). Captions are listed below."),
  prop = fp_norm)))
doc <- body_add_par(doc, "")
for (t in .REG$tables) {
  if (t$section != "supplement") next
  doc <- body_add_fpar(doc, fpar(
    ftext(paste0(tbl(t$id), ". "), prop = fp_bold),
    ftext(trimws(t$caption), prop = fp_norm)))
  doc <- body_add_par(doc, "")
}
doc <- body_add_break(doc, pos = "after")

# Extended Data Figures
for (f in .REG$figures) {
  if (f$section != "supplement") next
  doc <- add_legend_fpar(doc, paste0(fig(f$id), ". "),
                         split_legend_panels(f$caption), fp_bold, fp_norm)
  doc <- body_add_par(doc, "")
  png_p <- here(f$output$png)
  pdf_p <- here(f$output$pdf)
  if (file.exists(png_p) && requireNamespace("magick", quietly = TRUE)) {
    img  <- magick::image_read(png_p)
    info <- magick::image_info(img)
    doc  <- body_add_img(doc, png_p, width = 6.5,
                         height = 6.5 * info$height / info$width)
  } else if (file.exists(pdf_p) && requireNamespace("magick", quietly = TRUE)) {
    tmp_png <- tempfile(fileext = ".png")
    img <- magick::image_read_pdf(pdf_p, density = 300)
    magick::image_write(img, tmp_png, format = "png")
    info <- magick::image_info(img)
    doc <- body_add_img(doc, tmp_png, width = 6.5,
                        height = 6.5 * info$height / info$width)
    unlink(tmp_png)
  } else {
    doc <- body_add_fpar(doc, fpar(ftext(
      paste0("[Figure not found: ", basename(png_p), "]"), prop = fp_norm)))
  }
  doc <- body_add_break(doc, pos = "after")
}

supp_path <- file.path(DIR_SUPP, "Supplementary_Material.docx")
print(doc, target = supp_path)
cat("  →", basename(supp_path), "\n")

cat("Building Supplementary_Tables.xlsx...\n")
build_supp_xlsx()
