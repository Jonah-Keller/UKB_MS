#!/usr/bin/env Rscript
# ============================================================================
# compose_figure.R — Build a single-page composite image for ONE figure
# ============================================================================
# Usage:   Rscript compose_figure.R <figure_id>
# Example: Rscript compose_figure.R fig1_preonset_signature
#
# Reads the panel PDFs listed under the figure's `panels_dir` in figures.yaml,
# arranges them in a grid at 300 DPI, adds bold panel letters, and writes a
# single composite PNG to the figure's registered `output$png` path.
# Override the grid layout per-figure in figures.yaml via:
#   composite:
#     ncol: 3
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(magick)
})

source(here("code", "config", "00_config.R"))
source(here("code", "registry", "figure_registry.R"))
registry_load(FIGURES_YAML)

# ── Resolve target figure ────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
fig_id <- if (length(args) >= 1) args[[1]] else "fig1_preonset_signature"

fig <- NULL
for (f in .REG$figures) if (identical(f$id, fig_id)) { fig <- f; break }
if (is.null(fig)) stop("Figure id not found in registry: ", fig_id)

# ── Locate panel PDFs (same logic as build_submission.R::assemble_figure) ────
panels_dir <- here(fig$panels_dir)
if (!dir.exists(panels_dir)) stop("panels_dir missing: ", panels_dir)

panel_specs <- fig$panels %||% list()
if (length(panel_specs) == 0) {
  found <- sort(list.files(panels_dir,
                           pattern = "\\.(pdf|png|jpe?g|tiff?)$",
                           ignore.case = TRUE, full.names = TRUE))
  panel_specs <- lapply(found, function(p)
    list(letter = NA_character_, colspan = 1L, rowspan = 1L, pdf = p))
} else {
  panel_specs <- lapply(panel_specs, function(p) {
    # Optional `source` field: if it contains "/" treat as a path relative to
    # the project root; otherwise use as a letter pattern within panels_dir.
    if (!is.null(p$source) && grepl("/", p$source, fixed = TRUE)) {
      full_path <- here(p$source)
      fs <- if (file.exists(full_path)) full_path else character(0)
    } else {
      lookup <- p$source %||% p$letter
      pat <- paste0("^panel_[", tolower(lookup), toupper(lookup), "]([._]|$)")
      fs  <- list.files(panels_dir, pattern = pat, full.names = TRUE)
      # Filter to supported raster/vector extensions; prefer PDF when multiple
      # files exist for the same panel (e.g. a .pdf and a .png side by side).
      fs  <- fs[grepl("\\.(pdf|png|jpe?g|tiff?)$", fs, ignore.case = TRUE)]
      if (length(fs) > 1) {
        ext_rank <- match(tolower(tools::file_ext(fs)),
                          c("pdf", "svg", "png", "tif", "tiff", "jpg", "jpeg"),
                          nomatch = 99L)
        fs <- fs[order(ext_rank)]
      }
    }
    list(letter  = p$letter,
         colspan = as.integer(p$colspan %||% p$span %||% 1L),
         rowspan = as.integer(p$rowspan %||% 1L),
         pdf     = if (length(fs) == 0) NA_character_ else fs[[1]])
  })
  missing <- vapply(panel_specs, function(p) is.na(p$pdf), logical(1))
  if (any(missing))
    warning(fig_id, ": missing panels ",
            paste(vapply(panel_specs[missing], function(p) p$letter, character(1)),
                  collapse = ", "))
  panel_specs <- panel_specs[!missing]
}
if (length(panel_specs) == 0) stop("No panel PDFs resolved for ", fig_id)

cat(sprintf("Composing %s (%d panels): %s\n", fig_id, length(panel_specs),
            paste(vapply(panel_specs, function(p) {
              suffix <- ""
              if (p$colspan > 1 || p$rowspan > 1)
                suffix <- sprintf("(%dc\u00d7%dr)", p$colspan, p$rowspan)
              sprintf("%s%s", p$letter, suffix)
            }, character(1)), collapse = ", ")))

# ── Compose ─────────────────────────────────────────────────────────────────
# 2D bin-packed composite with col + row spans. Each panel is scaled to fit
# (colspan * cell_w); if the scaled height exceeds its allocated rowspan
# height, the panel is rescaled by height. Panels are first placed via a
# first-fit 2D bin-packer, then composited onto a single canvas at 300 DPI.
compose_figure <- function(fig, panel_specs, dpi = 300,
                           cell_width_in = 3.0,   # inches per column
                           gutter_px     = 20) {  # white space between cells
  ncol_grid <- fig$composite$ncol %||% 3L
  cell_w_px <- round(dpi * cell_width_in)
  col_w     <- cell_w_px + gutter_px       # per-column allocation incl. gutter

  colspans <- vapply(panel_specs, function(p) p$colspan, integer(1))
  rowspans <- vapply(panel_specs, function(p) p$rowspan, integer(1))

  # Guard: colspan must fit within the grid. Without this check, the
  # row-major packer below wraps forever (cur_c + cs - 1 > ncol_grid
  # is permanently true), silently hanging the build.
  bad <- which(colspans > ncol_grid)
  if (length(bad) > 0) {
    letters_bad <- vapply(panel_specs[bad], function(p) p$letter %||% "?",
                          character(1))
    stop(sprintf(
      "compose_figure: %s has panel(s) %s with colspan=%d > composite.ncol=%d. Reduce colspan or increase ncol in figures.yaml.",
      fig$id,
      paste(sprintf("'%s'(span=%d)", letters_bad, colspans[bad]), collapse = ", "),
      max(colspans[bad]), ncol_grid),
      call. = FALSE)
  }

  # Read each panel (PDF rendered at target DPI; raster formats read verbatim
  # then scaled) and rescale to its (colspan * cell_w) width
  read_panel <- function(path, dpi) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "pdf") magick::image_read_pdf(path, density = dpi, pages = 1)
    else              magick::image_read(path)
  }
  imgs <- lapply(seq_along(panel_specs), function(i) {
    img <- read_panel(panel_specs[[i]]$pdf, dpi)
    target_w <- cell_w_px * colspans[i] + gutter_px * (colspans[i] - 1L)
    magick::image_scale(img, as.character(target_w))
  })

  # ── Order-preserving row-major packing ────────────────────────────────────
  # A cursor (cur_r, cur_c) moves left-to-right, top-to-bottom. Each panel is
  # placed at the cursor if its (colspan x rowspan) block fits; otherwise the
  # cursor advances (wrapping rows) and any skipped cells remain empty. This
  # keeps the visual layout in the same order as the YAML panel list.
  grid <- matrix(FALSE, nrow = 0, ncol = ncol_grid)
  positions <- vector("list", length(panel_specs))
  cur_r <- 1L; cur_c <- 1L
  for (i in seq_along(panel_specs)) {
    cs <- colspans[i]; rs <- rowspans[i]
    placed <- FALSE
    while (!placed) {
      # Wrap to next row if current column can't fit the span
      if (cur_c + cs - 1L > ncol_grid) {
        cur_r <- cur_r + 1L; cur_c <- 1L
        next
      }
      # Ensure grid has enough rows
      while (nrow(grid) < cur_r + rs - 1L)
        grid <- rbind(grid, matrix(FALSE, nrow = 1, ncol = ncol_grid))
      block <- grid[cur_r:(cur_r + rs - 1L),
                    cur_c:(cur_c + cs - 1L), drop = FALSE]
      if (!any(block)) {
        grid[cur_r:(cur_r + rs - 1L), cur_c:(cur_c + cs - 1L)] <- TRUE
        positions[[i]] <- list(row = cur_r, col = cur_c,
                               colspan = cs, rowspan = rs)
        cur_c  <- cur_c + cs
        placed <- TRUE
      } else {
        cur_c <- cur_c + 1L  # skip past the blocker
      }
    }
  }
  nrow_total <- nrow(grid)

  # ── Row heights: max of (panel_h / rowspan) across panels in each row ─────
  row_h <- integer(nrow_total)
  for (i in seq_along(panel_specs)) {
    p     <- positions[[i]]
    img_h <- magick::image_info(imgs[[i]])$height
    per_r <- ceiling(img_h / p$rowspan)
    for (r in p$row:(p$row + p$rowspan - 1L))
      row_h[r] <- max(row_h[r], per_r)
  }
  row_h <- row_h + gutter_px  # add gutter

  # ── Render onto single canvas ─────────────────────────────────────────────
  total_w     <- col_w * ncol_grid
  total_h     <- sum(row_h)
  row_offset  <- c(0L, cumsum(row_h))
  canvas      <- magick::image_blank(total_w, total_h, color = "white")

  for (i in seq_along(panel_specs)) {
    p   <- positions[[i]]
    img <- imgs[[i]]

    cell_x <- (p$col - 1L) * col_w
    cell_y <- row_offset[p$row]
    cell_w_full <- col_w * p$colspan
    cell_h_full <- sum(row_h[p$row:(p$row + p$rowspan - 1L)])

    # Re-scale by height if panel is taller than its allocated cell
    ii <- magick::image_info(img)
    if (ii$height > cell_h_full - gutter_px) {
      img <- magick::image_scale(img, sprintf("x%d", cell_h_full - gutter_px))
      ii  <- magick::image_info(img)
    }

    pad_x <- cell_x + (cell_w_full - ii$width)  %/% 2L
    pad_y <- cell_y + (cell_h_full - ii$height) %/% 2L
    canvas <- magick::image_composite(canvas, img,
                                      offset = sprintf("+%d+%d", pad_x, pad_y))
  }

  png_path <- here(fig$output$png)
  dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
  magick::image_write(canvas, png_path, format = "png",
                      density = as.character(dpi))

  info <- magick::image_info(canvas)
  cat(sprintf("  \u2713 %s \u2192 %s\n    %d\u00d7%d px @ %d DPI  (%d rows \u00d7 %d cols, cell %.1f\")\n",
              fig$id, png_path, info$width, info$height, dpi,
              nrow_total, ncol_grid, cell_width_in))
  invisible(png_path)
}

compose_figure(fig, panel_specs)
