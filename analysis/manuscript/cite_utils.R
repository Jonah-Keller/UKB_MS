# ============================================================================
# cite_utils.R — Tokenized reference system for manuscript.Rmd
#
# Usage in Rmd text:
#   `r cite("goris2022")`               → ^11^
#   `r cite("walton2020", "gbd2019")`   → ^1,2^
#   `r render_references()`             → numbered reference list (at end of doc)
#
# Reference numbers are assigned in order of FIRST appearance in the document.
# Moving, deleting, or adding cite() calls automatically renumbers on next knit.
# All keys must exist in REFERENCES (defined in references.R).
# ============================================================================

# ── Internal state (environment survives across Rmd chunks) ─────────────────
.cite_state <- new.env(parent = emptyenv())
.cite_state$counter  <- 0L
.cite_state$assigned <- list()   # key → integer number

# ── Assign a number to a key (idempotent) ────────────────────────────────────
.cite_assign <- function(key) {
  if (!is.null(.cite_state$assigned[[key]])) {
    return(.cite_state$assigned[[key]])
  }
  if (!key %in% names(REFERENCES)) {
    stop("cite(): unknown reference key '", key,
         "'. Add it to analysis/manuscript/references.R")
  }
  .cite_state$counter <- .cite_state$counter + 1L
  n <- .cite_state$counter
  .cite_state$assigned[[key]] <- n
  n
}

# ── Main cite() function ─────────────────────────────────────────────────────
# Returns a superscript markdown string, e.g. "^1^" or "^1,2^"
cite <- function(...) {
  keys <- c(...)
  nums <- vapply(keys, .cite_assign, integer(1L))
  paste0("^", paste(nums, collapse = ","), "^")
}

# ── Render numbered reference list ───────────────────────────────────────────
# Call as `r render_references()` in the References section of the Rmd.
render_references <- function() {
  if (length(.cite_state$assigned) == 0L) return("")
  nums  <- unlist(as.list(.cite_state$assigned))
  keys  <- names(nums)[order(nums)]
  lines <- vapply(seq_along(keys), function(i) {
    paste0(i, ". ", REFERENCES[[ keys[[i]] ]])
  }, character(1L))
  paste(lines, collapse = "\n\n")
}
