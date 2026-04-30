#!/usr/bin/env Rscript
# ============================================================================
# build_submission.R — UKB MS proteomics manuscript build pipeline
# ============================================================================
#
# OUTPUTS (all under manuscript/submission/):
#
#   main_text/
#     manuscript.docx           — Main text: title + authors + intro + results
#                                  + discussion + methods + references + legends
#                                  Arial 12pt, 1" margins, continuous line numbers, page numbers
#
#   main_figures/
#     Figure_1_Preonset_Signature.{pdf,png}
#     Figure_2_Temporal_Window.{pdf,png}
#     Figure_3_Genetic_Heterogeneity.{pdf,png}
#     Figure_4_Sex_Differences.{pdf,png}
#     Figure_5_Endophenotypes.{pdf,png}
#     (written by analysis scripts during knit, not by this file)
#
#   main_tables/
#     Table_1_Cohort.docx       — Cohort characteristics (standalone landscape)
#
#   main_materials.docx         — Cover sheet + all 5 figures + Table 1 (landscape)
#                                  for reviewer/submission bundle
#
#   supplementary/
#     Supplementary_Material.docx  — Supplementary Tables 1–33 + Extended Data Figures 1–5 in one doc
#     Supplementary Tables/Supplementary Table_*.docx        — 33 standalone Supplementary Tables (landscape)
#     Extended Data Figures/Extended Data Fig._*.{pdf,png} — 5 Extended Data Figures
#                                    (written by analysis scripts during knit)
#
#   _build_intermediates/
#     _knitted.md               — Raw knitr output
#     manuscript.md             — Formatted markdown fed to pandoc
#     *.csv                     — Table data for each registered table
#
# BUILD STEPS:
#   0. Preflight: verify pandoc, python3, python-docx, reference.docx
#   1. Load config + registry (figures.yaml → 10 figures, 34 tables)
#   2. Locate manuscript.Rmd
#   3. knit() → _knitted.md  [HARD FAIL if knitr errors]
#   4. transform_markdown() + legends_block() → manuscript.md
#   5. registry_lint() on manuscript.md  [HARD FAIL on citation-order violation]
#   6. pandoc → manuscript.docx  [HARD FAIL if docx not produced]
#   7. python3 docx_formatter → margins + line numbers + page numbers + [REF] hl
#      [HARD FAIL if python-docx not installed]
#   8. build_tables_docx.R → standalone Table_N.docx  [HARD FAIL on error]
#   9. build_supplementary.R → Supplementary_Material.docx  [HARD FAIL on error]
#  10. Assemble main_materials.docx  [HARD FAIL on error]
#
# USAGE:
#   Rscript build_submission.R
#   MANUSCRIPT_RMD=path/to.Rmd Rscript build_submission.R
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(knitr)
  library(rmarkdown)
  library(officer)
  library(flextable)
})

source(here("code", "config", "00_config.R"))
source(here("code", "registry", "figure_registry.R"))
source(here("code", "tables", "style_flextable.R"))

# ── Debug/timing helper ──────────────────────────────────────────────────────
# Set DEBUG_BUILD=0 to silence these. Every line is flushed immediately so the
# user sees exactly which step is running (or hung) in real time.
.BUILD_T0 <- Sys.time()
.dbg <- function(...) {
  if (!identical(Sys.getenv("DEBUG_BUILD", "1"), "0")) {
    dt <- sprintf("%6.1fs", as.numeric(difftime(Sys.time(), .BUILD_T0, units = "secs")))
    cat("  [", dt, "] ", ..., "\n", sep = "")
    flush.console()
  }
}

cat("\n══ UKB MS submission build ══\n")
cat("Project:", PROJECT_TITLE, "\n")
cat("Journal:", TARGET_JOURNAL, "\n")

# ── [0] Preflight checks ─────────────────────────────────────────────────────
cat("\n[0/7] preflight checks\n")

# pandoc
pandoc_ver <- tryCatch(system("pandoc --version", intern = TRUE)[1], error = function(e) "")
if (!nzchar(pandoc_ver)) stop("pandoc not found. Install from https://pandoc.org/installing.html")
cat("  pandoc:", pandoc_ver, "\n")

# reference.docx
if (!file.exists(REFERENCE_DOCX))
  stop("Reference .docx not found: ", REFERENCE_DOCX,
       "\nPlace a styled reference document at manuscript/templates/reference.docx")
cat("  reference.docx: OK\n")

# python3 + python-docx
py_check <- system(
  'python3 -c "import docx; print(docx.__version__)"',
  intern = TRUE, ignore.stderr = TRUE)
if (!length(py_check) || grepl("Error|error", paste(py_check, collapse = "")))
  stop("python-docx not importable. Run: pip install python-docx")
cat("  python-docx:", py_check[1], "\n")

# magick (needed for PDF → PNG in main_materials and supplementary)
if (!requireNamespace("magick", quietly = TRUE))
  stop("R package 'magick' required. Install with: install.packages('magick')")
cat("  magick: OK\n")

# figures.yaml
if (!file.exists(FIGURES_YAML))
  stop("figures.yaml not found: ", FIGURES_YAML)
cat("  figures.yaml: OK\n")

# ── [1/7] Load registry ───────────────────────────────────────────────────────
cat("\n[1/7] registry\n")
registry_load(FIGURES_YAML)
s <- registry_summary()
cat(sprintf("  %d figures (%d main, %d supplement)\n",
    nrow(s[s$kind == "figure", ]),
    sum(s$kind == "figure" & s$section == "main"),
    sum(s$kind == "figure" & s$section == "supplement")))
cat(sprintf("  %d tables  (%d main, %d supplement)\n",
    nrow(s[s$kind == "table", ]),
    sum(s$kind == "table"  & s$section == "main"),
    sum(s$kind == "table"  & s$section == "supplement")))

# ── [2/7] Locate Rmd ─────────────────────────────────────────────────────────
cat("\n[2/7] locate Rmd\n")
rmd_path <- Sys.getenv("MANUSCRIPT_RMD", unset = "")
if (!nzchar(rmd_path)) {
  hits <- c(
    Sys.glob(here("analysis", "manuscript", "manuscript.Rmd")),
    Sys.glob(here("analysis", "02_*", "manuscript.Rmd"))
  )
  if (length(hits) == 0)
    stop("No manuscript.Rmd found. Set MANUSCRIPT_RMD env var or place at ",
         "analysis/manuscript/manuscript.Rmd")
  rmd_path <- hits[1]
}
if (!file.exists(rmd_path)) stop("Rmd not found: ", rmd_path)
cat("  Rmd:", rmd_path, "\n")

# ── [3/7] Knit Rmd → markdown ────────────────────────────────────────────────
cat("\n[3/7] knit\n")
.dbg("knit start \u2192 ", rmd_path)
md_intermediate <- file.path(DIR_BUILD_INT, "_knitted.md")
knit(rmd_path, output = md_intermediate, quiet = FALSE)
.dbg("knit \u2713")
if (!file.exists(md_intermediate) || file.info(md_intermediate)$size == 0)
  stop("Knit produced no output at: ", md_intermediate)
cat("  →", md_intermediate, sprintf("(%.1f KB)\n", file.info(md_intermediate)$size / 1024))

# ── [4/7] Format markdown + front matter + legends ───────────────────────────
cat("\n[4/7] format + lint\n")

#' Collapse runs of 3+ consecutive same-figure panel citations into a range.
#' e.g. "Figure 2a, Figure 2b, Figure 2c, Figure 2e" -> "Figure 2a\u2013c, Figure 2e"
#' Partial runs (<3 consecutive) are preserved. Mixed runs within a single
#' comma-separated list are handled independently.
collapse_panel_ranges <- function(txt) {
  cite_core <- "(?:Supplementary Table|Extended Data Fig\\.|Figure|Table)\\s+\\d+[a-zA-Z]"
  list_pat  <- sprintf("%s(?:,\\s+%s)+", cite_core, cite_core)
  cite_pat  <- "(Supplementary Table|Extended Data Fig\\.|Figure|Table)\\s+(\\d+)([a-zA-Z])"

  collapse_one_list <- function(s) {
    matches <- regmatches(s, gregexpr(cite_pat, s, perl = TRUE))[[1]]
    if (length(matches) < 2) return(s)
    parsed <- do.call(rbind, lapply(matches, function(x) {
      y <- regmatches(x, regexec(cite_pat, x, perl = TRUE))[[1]]
      data.frame(prefix = y[2], number = y[3], letter = y[4],
                 stringsAsFactors = FALSE)
    }))
    n <- nrow(parsed)
    is_consec <- logical(n)
    if (n > 1) {
      for (i in 2:n) {
        is_consec[i] <- parsed$prefix[i] == parsed$prefix[i - 1] &&
                        parsed$number[i] == parsed$number[i - 1] &&
                        utf8ToInt(parsed$letter[i]) ==
                          utf8ToInt(parsed$letter[i - 1]) + 1L
      }
    }
    run_id <- cumsum(!is_consec)
    out <- character()
    for (rid in unique(run_id)) {
      idx <- which(run_id == rid)
      if (length(idx) >= 3) {
        first <- parsed[idx[1], ]
        last  <- parsed[idx[length(idx)], ]
        out <- c(out, sprintf("%s %s%s\u2013%s",
                              first$prefix, first$number,
                              first$letter, last$letter))
      } else {
        for (i in idx) {
          p <- parsed[i, ]
          out <- c(out, sprintf("%s %s%s", p$prefix, p$number, p$letter))
        }
      }
    }
    paste(out, collapse = ", ")
  }

  m <- regmatches(txt, gregexpr(list_pat, txt, perl = TRUE))[[1]]
  if (length(m) == 0) return(txt)
  m <- unique(m)
  # Process longest matches first to avoid partial substring collisions
  m <- m[order(-nchar(m))]
  for (ms in m) {
    collapsed <- collapse_one_list(ms)
    if (!identical(collapsed, ms))
      txt <- gsub(ms, collapsed, txt, fixed = TRUE)
  }
  txt
}

transform_markdown <- function(md_lines) {
  yaml_end <- which(md_lines == "---")
  if (length(yaml_end) >= 2) md_lines <- md_lines[-(1:yaml_end[2])]

  in_code <- FALSE
  keep <- logical(length(md_lines))
  for (i in seq_along(md_lines)) {
    if (grepl("^```", md_lines[i])) { in_code <- !in_code; keep[i] <- FALSE }
    else keep[i] <- !in_code
  }
  md_lines <- md_lines[keep]
  md_lines <- md_lines[!grepl("<[a-zA-Z/!][^>]*>", md_lines)]
  md_lines <- md_lines[!grepl("^---$|^!\\[", md_lines)]

  txt <- paste(md_lines, collapse = "\n")
  txt <- gsub("\n{4,}", "\n\n\n", txt)

  txt <- gsub("(?m)^(#{2,4} .+)$", "**\\1**", txt, perl = TRUE)
  txt <- gsub("\\*\\*#+ ", "**", txt)

  txt <- gsub("(?<![*\\w])([PN])( [=<>\u2264\u2265])", "*\\1*\\2", txt, perl = TRUE)
  txt <- gsub("P~int~", "*P*~int~", txt, fixed = TRUE)

  for (phrase in c("et al","in vivo","in vitro","de novo","a priori",
                   "in situ","per se","ad hoc")) {
    pat <- sprintf("(?<!\\*)\\b(%s)\\b(?!\\*)", phrase)
    txt <- gsub(pat, "*\\1*", txt, perl = TRUE)
  }

  txt <- gsub(">=", "\u2265", txt)
  txt <- gsub("<=", "\u2264", txt)
  txt <- gsub("log2", "log~2~", txt)

  # Collapse 3+ consecutive same-figure panel citations into a hyphenated range
  # (must run before the bolding step so the range regex below matches cleanly)
  txt <- collapse_panel_ranges(txt)

  txt <- gsub(
    "(?<!\\*)((?:Supplementary Table|Extended Data Fig\\.|e?Figure|e?Table)\\s+\\d+[a-zA-Z]?(?:[\u2013\\-][a-zA-Z])?)(?!\\*)",
    "**\\1**", txt, perl = TRUE)

  txt <- gsub("p=<0.001",   "*P* < .001", txt, fixed = TRUE)
  txt <- gsub("p=<\\.001",  "*P* < .001", txt)
  txt <- gsub("p= ?([0-9])", "*P* = \\1", txt, perl = TRUE)

  txt <- gsub("(?m)^\\*\\*Session Information\\*\\*.*$", "", txt, perl = TRUE)
  txt
}

# Parse "Name1,2" into list(name="Name", sup="1,2")
parse_author <- function(a) {
  m <- regmatches(a, regexpr("[0-9][0-9,]*[*]?$", a))
  if (length(m) > 0 && nchar(m) > 0)
    list(name = sub("[0-9][0-9,]*[*]?$", "", a), sup = m)
  else
    list(name = a, sup = "")
}

# Build an fpar with author names and superscript affiliation numbers
author_fpar <- function(authors, fp_norm, fp_sup) {
  parts <- list()
  for (i in seq_along(authors)) {
    if (i > 1) parts <- c(parts, list(ftext(", ", prop = fp_norm)))
    p <- parse_author(authors[[i]])
    parts <- c(parts, list(ftext(p$name, prop = fp_norm)))
    if (nchar(p$sup) > 0)
      parts <- c(parts, list(ftext(p$sup, prop = fp_sup)))
  }
  do.call(fpar, parts)
}

# Helper: insert line breaks before each panel letter in a figure legend
split_legend_panels <- function(caption) {
  gsub("([.)])\\s+([a-zA-Z](?:[\u2013\\-][a-zA-Z])?,\\s)",
       "\\1\n\\2", trimws(caption))
}

# Helper: add legend paragraph with line breaks between panels (officer docx)
add_legend_fpar <- function(doc, label, text, fp_bold, fp_norm) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  parts <- list(ftext(label, prop = fp_bold))
  for (i in seq_along(lines)) {
    if (i > 1) parts <- c(parts, list(run_linebreak()))
    parts <- c(parts, list(ftext(lines[i], prop = fp_norm)))
  }
  body_add_fpar(doc, do.call(fpar, parts))
}

legends_block <- function() {
  # Only include entries that were actually cited during knit (have assigned numbers).
  # Iterating uncited entries here would create spurious new assignments.
  out <- c("\n\n**Table Legends**\n")
  for (t in .REG$tables) {
    if (t$section != "main") next
    if (is.null(.FIG_NUM$tbls[[t$id]])) next
    out <- c(out, sprintf("**%s.** %s\n", tbl(t$id), trimws(t$caption)))
  }
  out <- c(out, "\n**Figure Legends**\n")
  for (f in .REG$figures) {
    if (f$section != "main") next
    if (is.null(.FIG_NUM$figs[[f$id]])) next
    out <- c(out, sprintf("**%s.** %s\n", fig(f$id), split_legend_panels(f$caption)))
  }
  paste(out, collapse = "\n")
}

raw_lines  <- readLines(md_intermediate, warn = FALSE)
body_md    <- transform_markdown(raw_lines)
final_md   <- paste0(body_md, legends_block())
final_md_path <- file.path(DIR_BUILD_INT, "manuscript.md")
writeLines(final_md, final_md_path)

# Lint — hard fail on citation-order violations
.dbg("lint start \u2192 ", final_md_path)
registry_lint(final_md_path, strict = TRUE)
.dbg("lint \u2713")
cat("  → lint OK,", final_md_path, "\n")

# ── [5/7] Pandoc → docx ──────────────────────────────────────────────────────
cat("\n[5/7] pandoc → docx\n")
docx_path <- file.path(DIR_MAIN_TEXT, "manuscript.docx")
pandoc_cmd <- sprintf('pandoc -s "%s" -o "%s" --reference-doc "%s"',
                      final_md_path, docx_path, REFERENCE_DOCX)
status <- system(pandoc_cmd)
if (status != 0) stop("pandoc exited with status ", status)
if (!file.exists(docx_path)) stop("pandoc ran but manuscript.docx not created at: ", docx_path)
cat("  →", docx_path, sprintf("(%.0f KB)\n", file.info(docx_path)$size / 1024))

# ── [6/7] Python post-processing ─────────────────────────────────────────────
cat("\n[6/7] python post-process (margins, line numbers, page numbers, [REF] highlights)\n")
py_script <- tempfile(fileext = ".py")
writeLines(c(
  paste0('import sys; sys.path.insert(0, r"', here(), '")'),
  "from analysis.helpers.docx_formatter import postprocess_docx",
  "from pathlib import Path",
  paste0('postprocess_docx(Path(r"', docx_path, '"))')
), py_script)
status <- system(paste("python3", shQuote(py_script)))
if (status != 0) stop("docx_formatter.py exited non-zero (status ", status, "). ",
                      "Check: pip install python-docx")
cat("  → post-processing applied\n")

# ── [6.5/7] Assemble multi-panel figures from individual panel PDFs ───────────
cat("\n[6.5/7] assemble figures from panels\n")
if (!requireNamespace("pdftools", quietly = TRUE))
  stop("R package 'pdftools' required. Install with: install.packages('pdftools')")

# Returns named list of panel PDF paths (in panel order), or NULL on failure.
# Also writes the combined multi-page PDF to the registered output path.
assemble_figure <- function(fig) {
  pd <- fig$panels_dir
  if (is.null(pd)) { cat("  skip", fig$id, "(no panels_dir)\n"); return(invisible(NULL)) }
  panels_dir <- here(pd)
  if (!dir.exists(panels_dir)) {
    warning("panels_dir not found: ", panels_dir); return(invisible(NULL))
  }
  panel_ext_pat <- "\\.(pdf|png|jpe?g|tiff?)$"
  pick_panel_file <- function(fs) {
    fs <- fs[grepl(panel_ext_pat, fs, ignore.case = TRUE)]
    if (length(fs) == 0) return(NA_character_)
    # Prefer vector (PDF) over raster when both exist for the same panel
    ext_rank <- match(tolower(tools::file_ext(fs)),
                      c("pdf", "svg", "png", "tif", "tiff", "jpg", "jpeg"),
                      nomatch = 99L)
    fs[order(ext_rank)][[1]]
  }
  letters_wanted <- vapply(fig$panels %||% list(), function(p) p$letter, character(1))
  if (length(letters_wanted) == 0) {
    panel_pdfs <- sort(list.files(panels_dir, pattern = panel_ext_pat,
                                  ignore.case = TRUE, full.names = TRUE))
    letters_wanted <- rep(NA_character_, length(panel_pdfs))
  } else {
    panel_pdfs <- unlist(lapply(letters_wanted, function(l) {
      pat <- paste0("^panel_[", tolower(l), toupper(l), "]([._]|$)")
      fs  <- list.files(panels_dir, pattern = pat, full.names = TRUE)
      pick_panel_file(fs)
    }))
    missing_panels <- letters_wanted[is.na(panel_pdfs)]
    if (length(missing_panels) > 0)
      warning(fig$id, ": missing panels ", paste(missing_panels, collapse = ", "))
    keep <- !is.na(panel_pdfs)
    panel_pdfs <- panel_pdfs[keep]
    letters_wanted <- letters_wanted[keep]
  }
  if (length(panel_pdfs) == 0) {
    warning("No panel files found for ", fig$id); return(invisible(NULL))
  }
  out_pdf <- here(fig$output$pdf)
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

  # pdftools::pdf_combine requires PDF inputs. Convert raster panels (PNG/JPG/
  # TIFF) to a temporary single-page PDF via magick so they can be combined.
  to_pdf_for_combine <- function(path) {
    if (tolower(tools::file_ext(path)) == "pdf") return(path)
    tmp <- tempfile(fileext = ".pdf")
    img <- magick::image_read(path)
    magick::image_write(img, tmp, format = "pdf", density = "300")
    tmp
  }
  combine_inputs <- vapply(panel_pdfs, to_pdf_for_combine, character(1))
  pdftools::pdf_combine(combine_inputs, out_pdf)

  # Copy individual panel files (original format) to submission panels dir
  # for Illustrator editing
  panels_out_dir <- file.path(dirname(out_pdf), "panels", fig$id)
  unlink(panels_out_dir, recursive = TRUE)
  dir.create(panels_out_dir, showWarnings = FALSE, recursive = TRUE)
  for (pdf in panel_pdfs)
    file.copy(pdf, file.path(panels_out_dir, basename(pdf)), overwrite = TRUE)

  cat(sprintf("  ✓ %s (%d panels → %s)\n", fig$id, length(panel_pdfs), basename(out_pdf)))
  invisible(setNames(as.list(panel_pdfs), letters_wanted))
}

fig_panels <- list()  # fig$id -> named list of letter -> panel_pdf path
.dbg("assemble: ", length(.REG$figures), " figures to process")
for (f in .REG$figures) {
  .dbg("assemble \u2192 ", f$id, " (", f$section, ")")
  t0 <- Sys.time()
  res <- assemble_figure(f)
  .dbg("assemble \u2713 ", f$id, " in ",
       sprintf("%.2fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))),
       if (is.null(res)) " [NULL]" else sprintf(" [%d panels]", length(res)))
  if (!is.null(res)) fig_panels[[f$id]] <- res
}

# Compose single-page composite PNGs at 300 DPI for each assembled figure.
# Uses magick to grid-arrange panels, add bold letter labels, and output a
# single high-resolution image suitable for submission.
compose_figure <- function(fig, panel_pdfs, letters_wanted, dpi = 300) {
  ncol_grid <- fig$composite$ncol %||% min(4L, length(panel_pdfs))

  # Read each panel at target DPI
  imgs <- lapply(panel_pdfs, function(pdf) {
    magick::image_read_pdf(pdf, density = dpi, pages = 1)
  })

  # Add panel letter labels (bold, top-left)
  label_size <- round(dpi * 0.16)
  for (i in seq_along(imgs)) {
    lbl <- letters_wanted[i]
    if (!is.na(lbl) && nzchar(lbl)) {
      imgs[[i]] <- magick::image_annotate(
        imgs[[i]], lbl, size = label_size, weight = 700,
        location = "+10+5", color = "black")
    }
  }

  # Build grid row by row
  nrow_grid <- ceiling(length(imgs) / ncol_grid)
  rows <- list()
  for (r in seq_len(nrow_grid)) {
    idx <- seq((r - 1L) * ncol_grid + 1L, min(r * ncol_grid, length(imgs)))
    row_imgs <- imgs[idx]
    # Normalise heights within each row
    max_h <- max(vapply(row_imgs, function(x) magick::image_info(x)$height, numeric(1)))
    row_imgs <- lapply(row_imgs, function(x) {
      info <- magick::image_info(x)
      if (info$height < max_h)
        magick::image_extent(x, geometry = sprintf("%dx%d", info$width, max_h),
                             gravity = "north", color = "white")
      else x
    })
    rows[[r]] <- magick::image_append(Reduce(c, row_imgs), stack = FALSE)
  }

  # Normalise widths across rows and stack
  max_w <- max(vapply(rows, function(x) magick::image_info(x)$width, numeric(1)))
  rows <- lapply(rows, function(x) {
    info <- magick::image_info(x)
    if (info$width < max_w)
      magick::image_extent(x, geometry = sprintf("%dx%d", max_w, info$height),
                           gravity = "west", color = "white")
    else x
  })
  composite <- magick::image_append(Reduce(c, rows), stack = TRUE)

  # Write composite PNG
  png_path <- here(fig$output$png)
  dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
  magick::image_write(composite, png_path, format = "png", density = as.character(dpi))

  info <- magick::image_info(composite)
  cat(sprintf("  \u2713 composite %s \u2192 %s (%d panels, %d\u00d7%d px @ %d DPI)\n",
              fig$id, basename(png_path), length(panel_pdfs),
              info$width, info$height, dpi))
  invisible(png_path)
}

cat("\n[compose] generating composite PNGs via compose_figure.R\n")
.dbg("compose: ", length(.REG$figures), " figures to compose")
fig_composites  <- list()
compose_script  <- here("compose_figure.R")
compose_log_dir <- file.path(DIR_BUILD_INT, "compose_logs")
dir.create(compose_log_dir, showWarnings = FALSE, recursive = TRUE)
for (f in .REG$figures) {
  if (is.null(f$panels_dir)) {
    .dbg("compose skip ", f$id, " (no panels_dir)")
    next
  }
  n_panels <- length(f$panels %||% list())
  .dbg("compose \u2192 ", f$id, " (", n_panels, " panels, dir=", f$panels_dir, ")")
  t0 <- Sys.time()
  log_p <- file.path(compose_log_dir, paste0(f$id, ".log"))
  status <- tryCatch(
    system2("Rscript", args = c(shQuote(compose_script), f$id),
            stdout = log_p, stderr = log_p),
    error = function(e) 1L
  )
  png_p <- here(f$output$png)
  dt <- sprintf("%.2fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (status == 0L && file.exists(png_p)) {
    fig_composites[[f$id]] <- png_p
    .dbg("compose \u2713 ", f$id, " in ", dt)
  } else {
    .dbg("compose \u2717 ", f$id, " FAILED in ", dt, " (status=", status, ", log=", log_p, ")")
    warning("Composite failed for ", f$id, " \u2014 see ", log_p)
  }
}

# ── [7/7] Downstream artifacts ───────────────────────────────────────────────
cat("\n[7/7] tables + supplementary + main_materials\n")

# Standalone main tables
source(here("build_tables_docx.R"), local = TRUE)

# Supplementary_Material.docx
source(here("build_supplementary.R"), local = TRUE)

# main_materials.docx ─────────────────────────────────────────────────────────
cat("Building main_materials.docx...\n")

mm <- read_docx()
mm <- body_set_default_section(mm, prop_section(
  page_size    = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.8, bottom = 0.8, left = 0.8, right = 0.8)))

fp_title <- fp_text(bold = TRUE,  font.size = 14, font.family = "Arial")
fp_bold  <- fp_text(bold = TRUE,  font.size = 12, font.family = "Arial")
fp_norm  <- fp_text(bold = FALSE, font.size = 12, font.family = "Arial")

mm <- body_add_fpar(mm, fpar(ftext("Main Figures and Tables", prop = fp_title)))
mm <- body_add_par(mm, "")
mm <- body_add_fpar(mm, fpar(ftext(PROJECT_TITLE, prop = fp_bold)))
mm <- body_add_par(mm, "")
fp_sup <- fp_text(bold = FALSE, font.size = 10, font.family = "Arial",
                  vertical.align = "superscript")
mm <- body_add_fpar(mm, author_fpar(PROJECT_AUTHORS, fp_norm, fp_sup))
mm <- body_add_par(mm, "")

mm <- body_add_fpar(mm, fpar(ftext("Contents", prop = fp_bold)))
mm <- body_add_par(mm, "")

full_legend_text <- function(f) {
  split_legend_panels(f$caption)
}

for (f in .REG$figures)
  if (f$section == "main") {
    mm <- add_legend_fpar(mm, paste0(fig(f$id), ". "),
                          full_legend_text(f), fp_bold, fp_norm)
    mm <- body_add_par(mm, "")
  }

for (t in .REG$tables)
  if (t$section == "main") {
    mm <- body_add_fpar(mm, fpar(
      ftext(paste0(tbl(t$id), ". "), prop = fp_bold),
      ftext(trimws(t$caption), prop = fp_norm)))
    mm <- body_add_par(mm, "")
  }
mm <- body_add_break(mm, pos = "after")

# Main figures — legend text FIRST, then one panel per page at 300 DPI
missing_figs <- character(0)
embed_panel_pdf <- function(doc, pdf_path) {
  img <- if (tolower(tools::file_ext(pdf_path)) == "pdf")
           magick::image_read_pdf(pdf_path, density = 300, pages = 1)
         else
           magick::image_read(pdf_path)
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp))
  magick::image_write(img, tmp, format = "png")
  info <- magick::image_info(img)
  # For raster sources image_info returns nominal pixel dims; assume 300 DPI
  # for layout scaling (conservative — always downscaled to fit the page).
  dpi_assumed <- 300
  native_w <- info$width  / dpi_assumed
  native_h <- info$height / dpi_assumed
  scale <- min(1, 9 / native_w, 6.5 / native_h)
  w <- native_w * scale
  h <- native_h * scale
  body_add_img(doc, tmp, width = w, height = h)
}

for (f in .REG$figures) {
  if (f$section != "main") next
  panels <- fig_panels[[f$id]]
  if (is.null(panels)) {
    .dbg("main_materials: ", f$id, " MISSING panels \u2014 inserting placeholder")
    missing_figs <- c(missing_figs, basename(here(f$output$pdf)))
    mm <- body_add_fpar(mm, fpar(ftext(
      paste0("[MISSING: ", basename(here(f$output$pdf)), "]"), prop = fp_norm)))
    mm <- body_add_break(mm, pos = "after")
    next
  }
  .dbg("main_materials \u2192 embedding ", f$id)
  t0 <- Sys.time()
  # Full legend text before panels
  mm <- add_legend_fpar(mm, paste0(fig(f$id), ". "),
                        full_legend_text(f), fp_bold, fp_norm)
  mm <- body_add_par(mm, "")
  mm <- body_add_break(mm, pos = "after")
  # Embed composite image if available, otherwise individual panels
  composite_png <- fig_composites[[f$id]]
  if (!is.null(composite_png) && file.exists(composite_png)) {
    .dbg("  composite path \u2192 ", basename(composite_png))
    cinfo <- magick::image_info(magick::image_read(composite_png))
    native_w <- cinfo$width  / 300
    native_h <- cinfo$height / 300
    scale <- min(1, 9 / native_w, 6.5 / native_h)
    mm <- body_add_img(mm, composite_png,
                       width = native_w * scale, height = native_h * scale)
    mm <- body_add_break(mm, pos = "after")
  } else {
    .dbg("  no composite PNG \u2014 falling back to ", length(panels), " per-panel PDF embeds (slow)")
    for (j in seq_along(panels)) {
      panel_pdf <- panels[[j]]
      .dbg("    embed panel ", j, "/", length(panels), " \u2192 ", basename(panel_pdf))
      tp <- Sys.time()
      mm <- embed_panel_pdf(mm, panel_pdf)
      mm <- body_add_break(mm, pos = "after")
      .dbg("    embed \u2713 in ",
           sprintf("%.2fs", as.numeric(difftime(Sys.time(), tp, units = "secs"))))
    }
  }
  .dbg("main_materials \u2713 ", f$id, " in ",
       sprintf("%.2fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# Main tables — each must have a CSV or build fails
missing_tbls <- character(0)
for (t in .REG$tables) {
  if (t$section != "main") next
  mm <- body_add_fpar(mm, fpar(
    ftext(paste0(tbl(t$id), ". "), prop = fp_bold),
    ftext(trimws(t$caption), prop = fp_norm)))
  mm <- body_add_par(mm, "")
  csv_p <- here(t$output$csv)
  if (!file.exists(csv_p)) {
    missing_tbls <- c(missing_tbls, basename(csv_p))
    mm <- body_add_fpar(mm, fpar(ftext(
      paste0("[MISSING CSV: ", basename(csv_p), "]"), prop = fp_norm)))
  } else {
    d  <- read.csv(csv_p, check.names = FALSE)
    ft <- flextable::flextable(d) |>
            style_flextable() |>
            flextable::align(j = 2:ncol(d), align = "center", part = "body") |>
            flextable::align(j = 2:ncol(d), align = "center", part = "header")
    mm <- flextable::body_add_flextable(mm, ft)
  }
  mm <- body_add_break(mm, pos = "after")
}

mm_path <- file.path(DIR_SUBMISSION, "main_materials.docx")
print(mm, target = mm_path)
if (!file.exists(mm_path)) stop("main_materials.docx was not created at: ", mm_path)
cat("  →", mm_path, sprintf("(%.0f KB)\n", file.info(mm_path)$size / 1024))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n══ build complete ══\n")
cat("Output root:", DIR_SUBMISSION, "\n")

outputs <- list(
  "main_text/manuscript.docx"                = file.path(DIR_MAIN_TEXT, "manuscript.docx"),
  "main_materials.docx"                      = mm_path,
  "supplementary/Supplementary_Material.docx"= file.path(DIR_SUPP, "Supplementary_Material.docx")
)
for (nm in names(outputs)) {
  p <- outputs[[nm]]
  if (file.exists(p)) {
    cat(sprintf("  ✓  %-50s  %.0f KB\n", nm, file.info(p)$size / 1024))
  } else {
    cat(sprintf("  ✗  %-50s  NOT FOUND\n", nm))
  }
}

# Warn (not fail) if figures or tables were missing at bundle time
if (length(missing_figs) > 0)
  warning("main_materials.docx built with missing figures: ",
          paste(missing_figs, collapse = ", "),
          "\nRun the figure analysis scripts first.")
if (length(missing_tbls) > 0)
  warning("main_materials.docx built with missing table CSVs: ",
          paste(missing_tbls, collapse = ", "),
          "\nRun the table analysis scripts first.")

cat("\n")
