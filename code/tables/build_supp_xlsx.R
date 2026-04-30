# ============================================================================
# build_supp_xlsx.R — Pack all supplementary tables into one .xlsx
# ============================================================================
# Iterates supp tables in the registry, reads their staged CSVs, and writes
# one sheet per table into supplementary/Supplementary_Tables.xlsx, with a
# Contents sheet listing sheet name -> table number -> caption.
#
# Sheet names use the already-assigned `tbl()` number (so they match what the
# manuscript text and the Supplementary_Material.docx contents listing show),
# truncated to Excel's 31-character limit.
#
# Requires `stage_table_csvs()` to have run first so output$csv paths exist.
# ============================================================================

build_supp_xlsx <- function(out_path = file.path(DIR_SUPP,
                                                 "Supplementary_Tables.xlsx")) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("openxlsx required for build_supp_xlsx()")

  wb <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(textDecoration = "bold",
                                     fgFill = "#F2F2F2",
                                     border = "Bottom")

  openxlsx::addWorksheet(wb, "Contents")

  contents <- data.frame(Sheet = character(0), Table = character(0),
                         Caption = character(0), stringsAsFactors = FALSE)
  missing  <- character(0)
  n <- 0L
  for (t in .REG$tables) {
    if (!identical(t$section, "supplement")) next
    csv_p <- here(t$output$csv)
    if (!file.exists(csv_p)) {
      missing <- c(missing, sprintf("%s (%s)", t$id, basename(csv_p)))
      next
    }
    d <- read.csv(csv_p, check.names = FALSE)

    label <- tbl(t$id)                              # "Supplementary Table N"
    num   <- sub(".*Table\\s+", "", label)          # "N"
    short <- sub("^supp(lementary)?_?table_?", "",  # strip redundant prefix
                 t$slug, ignore.case = TRUE)
    sheet <- substr(paste0("S", num, "_", short), 1L, 31L)
    sheet <- gsub("[\\\\/?*\\[\\]:]", "_", sheet)

    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, d, headerStyle = hdr_style)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(d)), widths = "auto")

    contents[nrow(contents) + 1L, ] <- list(sheet, label, trimws(t$caption))
    n <- n + 1L
  }

  openxlsx::writeData(wb, "Contents", contents, headerStyle = hdr_style)
  openxlsx::setColWidths(wb, "Contents", cols = 1:3,
                         widths = c(28, 26, 100))
  openxlsx::freezePane(wb, "Contents", firstRow = TRUE)

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  cat(sprintf("  wrote %d sheets → %s\n", n, basename(out_path)))
  if (length(missing) > 0)
    warning("build_supp_xlsx: ", length(missing),
            " supp table(s) missing CSV: ",
            paste(missing, collapse = "; "), call. = FALSE)
  invisible(out_path)
}
