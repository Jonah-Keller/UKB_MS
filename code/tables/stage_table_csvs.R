# ============================================================================
# stage_table_csvs.R — Copy registered tables' source CSVs into _build_intermediates
# ============================================================================
# For each table in the registry that defines `source_csv:`, copy (or row-bind
# if a list of paths) the source file(s) into the registry-defined output$csv
# path. Loose copy: no column reshaping. Missing sources are warned, not fatal.
#
# Call once in build_submission.R before build_tables_docx.R / build_supplementary.R.
# ============================================================================

stage_table_csvs <- function() {
  staged   <- 0L
  missing  <- character(0)
  for (t in .REG$tables) {
    src <- t$source_csv
    if (is.null(src)) next
    out_csv <- here(t$output$csv)
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

    src_paths <- vapply(as.list(src), function(p) here(as.character(p)),
                        character(1))
    exist <- file.exists(src_paths)
    if (!any(exist)) {
      missing <- c(missing, sprintf("%s ← %s", t$id,
                                     paste(src_paths[!exist], collapse = ", ")))
      next
    }
    if (length(src_paths) == 1L) {
      file.copy(src_paths[1], out_csv, overwrite = TRUE)
    } else {
      # rbind multiple sources — fill missing columns with NA across files
      dfs <- lapply(src_paths[exist], read.csv, check.names = FALSE)
      all_cols <- unique(unlist(lapply(dfs, names)))
      dfs <- lapply(dfs, function(d) {
        for (cc in setdiff(all_cols, names(d))) d[[cc]] <- NA
        d[, all_cols, drop = FALSE]
      })
      combined <- do.call(rbind, dfs)
      write.csv(combined, out_csv, row.names = FALSE)
    }
    staged <- staged + 1L
  }
  cat(sprintf("  staged %d table CSV(s) into %s\n",
              staged, sub(here(), "", DIR_BUILD_INT, fixed = TRUE)))
  if (length(missing) > 0) {
    warning("stage_table_csvs: missing source CSVs for ", length(missing),
            " table(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
  invisible(staged)
}
