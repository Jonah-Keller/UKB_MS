# ============================================================================
# style_flextable.R — Canonical flextable styling for manuscript tables
# ============================================================================
# Ligh Lab house style (derived from NfL_Cardiovascular_UKB):
#   Arial, body 9 pt, header 9.5 pt bold on #f5f5f5,
#   1.5 pt top-header / bottom-body rules, 0.5 pt under header,
#   padding 2 pt, autofit add_w = 0.1
# ============================================================================

suppressPackageStartupMessages({
  library(flextable)
  library(officer)
})

# Set font family + base size + padding globally for all flextables created
# in this R session. This SIDESTEPS a compatibility bug in flextable 0.9.11
# paired with ggplot2 4.x where calling `font(fontname = "Arial", part = ...)`
# internally forwards fontname/part into ggplot2::element_text and crashes
# with "EXPR must be a length 1 vector". Using set_flextable_defaults
# establishes the same typography without triggering the buggy code path.
flextable::set_flextable_defaults(
  font.family  = "Arial",
  font.size    = 9,
  padding      = 2,
  border.color = "black"
)

#' Apply Ligh Lab style to a flextable.
#'
#' Font family ("Arial") and body size (9 pt) come from
#' `set_flextable_defaults()` above; this function only applies per-part
#' overrides (header typography/background) and border/alignment rules.
style_flextable <- function(ft, caption = NULL,
                            header_size = 9.5) {
  if (!is.null(caption)) {
    ft <- set_caption(ft,
      caption = as_paragraph(as_chunk(caption,
        props = fp_text_default(bold = TRUE, font.size = 11,
                                font.family = "Arial"))),
      fp_p = fp_par(padding.bottom = 6))
  }
  ft |>
    fontsize(size = header_size, part = "header") |>
    bold(part = "header") |>
    bg(bg = "#f5f5f5", part = "header") |>
    align(align = "left", part = "body") |>
    align(align = "left", part = "header") |>
    border_remove() |>
    hline_top(border    = fp_border(width = 1.5), part = "header") |>
    hline_bottom(border = fp_border(width = 0.5), part = "header") |>
    hline_bottom(border = fp_border(width = 1.5), part = "body") |>
    autofit(add_w = 0.1)
}

#' Save a flextable as a standalone .docx with landscape orientation
save_flextable_docx <- function(ft, path, footnote = NULL,
                                orient = c("landscape", "portrait"),
                                margin_in = 0.8) {
  orient <- match.arg(orient)
  doc <- read_docx()
  doc <- body_set_default_section(doc, prop_section(
    page_size    = page_size(orient = orient),
    page_margins = page_mar(top = margin_in, bottom = margin_in,
                            left = margin_in, right = margin_in)))
  doc <- body_add_flextable(doc, ft)
  if (!is.null(footnote)) {
    doc <- body_add_fpar(doc,
      fpar(ftext(footnote,
        prop = fp_text(font.size = 8, font.family = "Arial"))),
      pos = "after")
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  print(doc, target = path)
  invisible(path)
}
