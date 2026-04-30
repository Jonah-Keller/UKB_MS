# ============================================================================
# figure_registry.R — Load figures.yaml; assign numbers by first citation order
# ============================================================================
# API:
#   registry_load("figures.yaml")   # parses metadata; does NOT assign numbers
#   registry_reset()                # clear numbering state (call in Rmd setup)
#   fig("fig1_id")                  # -> "Figure 1"  (number assigned on first call)
#   fig("fig1_id", "b")             # -> "Figure 1b"
#   fig_caption("fig1_id")          # -> "**Figure 1.** <caption>"
#   tbl("tbl1_id")                  # -> "Table 1"   (number assigned on first call)
#   out_path("fig1_id", "pdf")      # resolved output path (static, from YAML)
#   registry_lint(md_path)          # hard-fail on panel-interleave or order violations
#   registry_summary()              # data.frame of assigned numbers (post-knit)
#
# NUMBERING: Numbers are assigned lazily, in the order fig()/tbl() are first
# called during knit.  Main figures and supplement figures are numbered
# independently.  registry_reset() must be called at the top of every Rmd
# setup chunk to guarantee a clean count on each knit.
#
# LINT RULE: Two checks, both hard-fail by default:
#   1. First-introduction order — the first citation of Figure N must precede
#      the first citation of Figure N+1.
#   2. Panel-interleave — once any panel of Figure N+1 appears in the text,
#      no panel of Figure N may follow it.  Back-references to a bare
#      "Figure N" (without a panel letter) are always permitted.
#      Example violation: "Figure 2a … Figure 1d"
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(here)
})

# ── Registry store (metadata only, no numbers) ──────────────────────────────
.REG <- new.env(parent = emptyenv())
.REG$figures <- list()
.REG$tables  <- list()
.REG$loaded  <- FALSE

# ── Dynamic numbering state (reset between knits) ────────────────────────────
.FIG_NUM <- new.env(parent = emptyenv())
.FIG_NUM$fig_main_n <- 0L
.FIG_NUM$fig_supp_n <- 0L
.FIG_NUM$tbl_main_n <- 0L
.FIG_NUM$tbl_supp_n <- 0L
.FIG_NUM$figs <- list()   # id -> integer
.FIG_NUM$tbls <- list()   # id -> integer

#' Reset numbering state.  Call once in every Rmd setup chunk.
registry_reset <- function() {
  .FIG_NUM$fig_main_n <- 0L
  .FIG_NUM$fig_supp_n <- 0L
  .FIG_NUM$tbl_main_n <- 0L
  .FIG_NUM$tbl_supp_n <- 0L
  .FIG_NUM$figs <- list()
  .FIG_NUM$tbls <- list()
  invisible(NULL)
}

# ── Internal helpers ─────────────────────────────────────────────────────────

.slugify <- function(s) {
  s <- tolower(s); s <- gsub("[^a-z0-9]+", "_", s); s <- gsub("^_|_$", "", s); s
}

.resolve_paths <- function(entry, slug) {
  if (is.null(entry$output)) return(entry)
  entry$output <- lapply(entry$output, function(p) gsub("\\{slug\\}", slug, p))
  entry
}

.assign_panels <- function(fig) {
  if (is.null(fig$panels) || length(fig$panels) == 0) return(fig)
  auto_letters <- letters[seq_along(fig$panels)]
  for (i in seq_along(fig$panels)) {
    if (is.null(fig$panels[[i]]$letter))
      fig$panels[[i]]$letter <- auto_letters[i]
    fig$panels[[i]]$index <- i
  }
  fig
}

# ── Public API ───────────────────────────────────────────────────────────────

#' Load registry metadata from figures.yaml.
#' Does NOT assign numbers — numbers come from first fig()/tbl() call order.
registry_load <- function(yaml_path = here("figures.yaml")) {
  raw <- yaml::read_yaml(yaml_path)

  for (kind in c("figures", "tables")) {
    entries <- raw[[kind]] %||% list()
    for (i in seq_along(entries)) {
      e <- entries[[i]]
      if (is.null(e$id))
        stop("Registry entry missing `id` in ", kind, "[", i, "]")
      if (is.null(e$section))
        stop("Registry entry missing `section` in ", e$id)
      if (!e$section %in% c("main", "supplement"))
        stop("section must be 'main' or 'supplement', got '", e$section,
             "' in ", e$id)
      e$slug <- .slugify(e$id)
      e <- .resolve_paths(e, e$slug)
      if (kind == "figures") e <- .assign_panels(e)
      entries[[i]] <- e
    }
    names(entries) <- vapply(entries, `[[`, character(1), "id")
    .REG[[kind]] <- entries
  }
  .REG$loaded <- TRUE
  invisible(.REG)
}

.ensure_loaded <- function() {
  if (!.REG$loaded) registry_load()
}

.lookup <- function(id, kind) {
  .ensure_loaded()
  e <- .REG[[kind]][[id]]
  if (is.null(e)) stop(kind, " id '", id, "' not in figures.yaml")
  e
}

# Lazy number assignment for figures
.fig_num <- function(id) {
  if (!is.null(.FIG_NUM$figs[[id]])) return(.FIG_NUM$figs[[id]])
  e <- .lookup(id, "figures")
  if (e$section == "main") {
    .FIG_NUM$fig_main_n <- .FIG_NUM$fig_main_n + 1L
    n <- .FIG_NUM$fig_main_n
  } else {
    .FIG_NUM$fig_supp_n <- .FIG_NUM$fig_supp_n + 1L
    n <- .FIG_NUM$fig_supp_n
  }
  .FIG_NUM$figs[[id]] <- n
  n
}

# Lazy number assignment for tables
.tbl_num <- function(id) {
  if (!is.null(.FIG_NUM$tbls[[id]])) return(.FIG_NUM$tbls[[id]])
  e <- .lookup(id, "tables")
  if (e$section == "main") {
    .FIG_NUM$tbl_main_n <- .FIG_NUM$tbl_main_n + 1L
    n <- .FIG_NUM$tbl_main_n
  } else {
    .FIG_NUM$tbl_supp_n <- .FIG_NUM$tbl_supp_n + 1L
    n <- .FIG_NUM$tbl_supp_n
  }
  .FIG_NUM$tbls[[id]] <- n
  n
}

#' Figure reference.  Number assigned on first call (text order = figure order).
#' fig("id")        -> "Figure 1"
#' fig("id", "b")   -> "Figure 1b"   (letter directly)
#' fig("id", "p1b") -> "Figure 1b"   (panel id lookup)
fig <- function(id, panel = NULL) {
  e <- .lookup(id, "figures")
  n  <- .fig_num(id)
  prefix <- if (e$section == "supplement") "Extended Data Fig. " else "Figure "
  out <- paste0(prefix, n)
  if (!is.null(panel)) {
    if (panel %in% c(letters, LETTERS)) {
      out <- paste0(out, panel)
    } else {
      p <- Filter(function(pn) pn$id == panel, e$panels)
      if (length(p) == 0)
        stop("panel id '", panel, "' not found in figure '", id, "'")
      out <- paste0(out, p[[1]]$letter)
    }
  }
  out
}

#' Full caption line: "**Figure 1.** <caption text>"
fig_caption <- function(id) {
  e <- .lookup(id, "figures")
  paste0("**", fig(id), ".** ", trimws(e$caption))
}

#' Table reference.  Number assigned on first call (text order = table order).
tbl <- function(id) {
  e <- .lookup(id, "tables")
  n  <- .tbl_num(id)
  prefix <- if (e$section == "supplement") "Supplementary Table " else "Table "
  paste0(prefix, n)
}

tbl_caption <- function(id) {
  e <- .lookup(id, "tables")
  paste0("**", tbl(id), ".** ", trimws(e$caption))
}

#' Resolved static output path for a registered figure or table.
out_path <- function(id, format = c("pdf", "png", "svg", "csv", "docx")) {
  format <- match.arg(format)
  e <- tryCatch(.lookup(id, "figures"), error = function(e) .lookup(id, "tables"))
  p <- e$output[[format]]
  if (is.null(p)) stop("format '", format, "' not registered for '", id, "'")
  here(p)
}

#' Summary data.frame of all assigned numbers (meaningful only after knit).
registry_summary <- function() {
  .ensure_loaded()
  rows <- list()
  for (e in .REG$figures) {
    n   <- .FIG_NUM$figs[[e$id]]
    ref <- if (!is.null(n)) fig(e$id) else "(unassigned)"
    rows[[length(rows) + 1]] <- data.frame(
      kind = "figure", id = e$id, section = e$section,
      number = n %||% NA_integer_, ref = ref,
      n_panels = length(e$panels %||% list()),
      stringsAsFactors = FALSE)
  }
  for (e in .REG$tables) {
    n   <- .FIG_NUM$tbls[[e$id]]
    ref <- if (!is.null(n)) tbl(e$id) else "(unassigned)"
    rows[[length(rows) + 1]] <- data.frame(
      kind = "table", id = e$id, section = e$section,
      number = n %||% NA_integer_, ref = ref,
      n_panels = NA_integer_,
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# ── Citation-order linter ────────────────────────────────────────────────────

#' Lint a rendered-markdown file for citation-order violations.
#'
#' Three rules, all hard-fail by default:
#'
#'   Rule 1 — First-introduction order:
#'     The first citation of Figure N must precede the first citation of Figure N+1.
#'     (Figures 1..N-1 must all be introduced before Figure N is introduced.)
#'
#'   Rule 2 — Panel-interleave (across figures):
#'     Once any PANEL citation of Figure N+1 appears, no panel citation of
#'     Figure N may follow it.  Bare figure refs (no panel) are always allowed.
#'     Example violation: "Figure 2a" then later "Figure 1d".
#'
#'   Rule 3 — Within-figure panel order + completeness:
#'     (3a) Within a single figure, the first citation of panel X must precede
#'     the first citation of panel Y whenever X < Y in the alphabet.
#'     (3b) Every panel letter defined in figures.yaml for a figure must be
#'     cited at least once in the body (outside the Legends block).  Panels
#'     defined but never cited are reported as missing.
#'     Range citations (e.g., "Figure 2a\u2013c") introduce every letter in
#'     the range at their position.  Later re-citations of an already-
#'     introduced panel are always allowed.
#'     Example violation: "Figure 3b" first cited, later "Figure 3a".
#'
#'   Rule 4 — Cite \u2192 registry (reverse check):
#'     Every panel letter cited in the body must correspond to a panel
#'     actually defined in figures.yaml for that figure.  Catches typos and
#'     stale citations that point to panels that no longer exist.
#'     Example violation: manuscript cites "Figure 1z" but figures.yaml
#'     defines only panels a\u2013i for Figure 1.
#'
#'   Rule 5 — Registry \u2192 composite (source resolvability):
#'     Every panel defined in figures.yaml must resolve to an existing
#'     source PDF using the same logic as compose_figure.R:
#'       \u2022 if `source` contains "/", it is a project-relative path
#'       \u2022 otherwise `source` (or `letter` if unset) is used as a
#'         case-insensitive pattern "^panel_[xX]([._]|$)" inside `panels_dir`
#'     Catches broken composites where a yaml panel points to a missing
#'     upstream PDF that would be silently dropped by compose_figure.R.
#'
#'   (Note: re-citing the same panel multiple times in the body is
#'   permitted.  Each (figure, letter) pair resolves to exactly one panel
#'   by yaml construction, so duplicate citations always refer to the
#'   same panel and do not need to be flagged.)
#'
#' The Table/Figure Legends block at the end of the document is excluded
#' from Rule 3 because legends list panels in alphabetical order by design.
#'
#' Numbers in the markdown reflect text-order assignment (from fig()/tbl()),
#' so Rule 1 catches genuine out-of-order introductions, and Rule 2 catches
#' interleaved panel citations across adjacent figures.
#'
#' @param md_path  path to the post-transform markdown (figure refs are bolded)
#' @param strict   if TRUE (default) stop() on violation; if FALSE warn() only
#' @return invisible list of violations
registry_lint <- function(md_path, strict = TRUE) {
  .ensure_loaded()
  md <- paste(readLines(md_path, warn = FALSE), collapse = "\n")

  pat <- "\\*\\*((?:Supplementary Table|Extended Data Fig\\.|e?Figure|e?Table))\\s+(\\d+)([a-zA-Z])?(?:[\u2013\\-]([a-zA-Z]))?\\*\\*"
  gr  <- gregexpr(pat, md, perl = TRUE, ignore.case = TRUE)[[1]]
  if (gr[1] == -1L) {
    message("[registry_lint] no citations found in ", md_path)
    return(invisible(list()))
  }
  m         <- regmatches(md, list(gr))[[1]]
  positions <- as.integer(gr)

  # Identify the "Legends" block so Rule 3 can exclude it (legends list panels
  # alphabetically by design and would otherwise always match ascending order)
  legend_m <- regexpr("\\*\\*(?:Table|Figure) Legends\\*\\*", md, perl = TRUE)
  legend_start <- if (legend_m[1] > 0) as.integer(legend_m[1]) else .Machine$integer.max

  parse_one <- function(s, pos) {
    mm       <- regmatches(s, regexec(pat, s, perl = TRUE, ignore.case = TRUE))[[1]]
    full_pfx <- mm[2]
    kind     <- if (grepl("fig", full_pfx, ignore.case = TRUE)) "figure" else "table"
    section  <- if (grepl("^(?:Supplementary|Extended|e)", full_pfx, perl = TRUE)) "supplement" else "main"
    num      <- as.integer(mm[3])
    p_start  <- if (nzchar(mm[4])) tolower(mm[4]) else NA_character_
    p_end    <- if (nzchar(mm[5])) tolower(mm[5]) else p_start
    list(raw = s, pos = pos, kind = kind, section = section,
         number = num, panel = p_start, panel_to = p_end)
  }
  cites <- Map(parse_one, m, positions)

  violations <- list()

  # Rule 1 — first-introduction order (unchanged from original logic)
  check_first_mention_order <- function(subset, label) {
    if (length(subset) < 2) return()
    first_seen    <- integer(0)
    max_introduced <- 0L
    for (cur in subset) {
      n <- cur$number
      if (!(n %in% first_seen)) {
        if (n < max_introduced)
          violations[[length(violations) + 1]] <<- list(
            reason = sprintf(
              "Rule 1 [%s]: %s introduced out of order — Figure %d was already introduced",
              label, cur$raw, max_introduced))
        first_seen     <- c(first_seen, n)
        max_introduced <- max(max_introduced, n)
      }
    }
  }

  # Rule 2 — no panel of Figure N after any panel of Figure N+1
  check_panel_interleave <- function(subset, label) {
    panel_cites <- Filter(function(x) !is.na(x$panel), subset)
    if (length(panel_cites) < 2) return()
    max_panel_fig <- 0L
    for (cur in panel_cites) {
      n <- cur$number
      if (n < max_panel_fig)
        violations[[length(violations) + 1]] <<- list(
          reason = sprintf(
            "Rule 2 [%s]: panel citation %s appears after panels of Figure %d — interleaving not allowed",
            label, cur$raw, max_panel_fig))
      max_panel_fig <- max(max_panel_fig, n)
    }
  }

  # Rule 3 — within a single figure, first citations of panels must appear in
  # alphabetical order (3a), and every panel letter defined in figures.yaml
  # must be cited at least once in the body (3b).
  # Range citations introduce every letter in the range.
  # Legends (everything at/after the "Figure Legends" heading) are excluded.
  check_within_figure_panel_order <- function(subset, label, section) {
    body_cites <- Filter(function(x) x$pos < legend_start, subset)
    # Build reverse map: rendered number -> yaml id (within this section)
    num_to_id <- list()
    for (id in names(.FIG_NUM$figs)) {
      e <- .REG$figures[[id]]
      if (!is.null(e) && e$section == section)
        num_to_id[[as.character(.FIG_NUM$figs[[id]])]] <- id
    }
    # Every figure we might want to check, based on BOTH cited numbers and
    # numbers assigned in the registry (so a figure with zero body citations
    # is still flagged if it has yaml panels).
    cited_nums <- unique(vapply(body_cites, `[[`, integer(1), "number"))
    all_nums   <- unique(c(cited_nums, as.integer(names(num_to_id))))
    for (nn in all_nums) {
      fig_cites <- Filter(function(x) x$number == nn && !is.na(x$panel),
                          body_cites)
      # --- (3a) ascending-order check -----------------------------------
      seen           <- character(0)
      max_introduced <- ""
      for (cur in fig_cites) {
        a <- utf8ToInt(cur$panel)
        b <- utf8ToInt(cur$panel_to)
        for (code in seq(a, b)) {
          L <- intToUtf8(code)
          if (L %in% seen) next           # already introduced → fine
          if (nzchar(max_introduced) && L < max_introduced)
            violations[[length(violations) + 1]] <<- list(
              reason = sprintf(
                "Rule 3a [%s]: panel %s of %s first cited after panel %s (panels must be introduced in alphabetical order)",
                label, L, sub("[a-z]$", "", cur$raw), max_introduced))
          seen <- c(seen, L)
          if (L > max_introduced) max_introduced <- L
        }
      }
      # --- (3b) completeness check: every yaml panel letter is cited ----
      id <- num_to_id[[as.character(nn)]]
      if (is.null(id)) next
      yaml_letters <- tolower(vapply(.REG$figures[[id]]$panels,
                                     function(p) as.character(p$letter %||% ""),
                                     character(1)))
      yaml_letters <- yaml_letters[nzchar(yaml_letters)]
      missing <- setdiff(yaml_letters, seen)
      if (length(missing) > 0) {
        prefix <- if (section == "supplement") "Extended Data Fig. " else "Figure "
        violations[[length(violations) + 1]] <<- list(
          reason = sprintf(
            "Rule 3b [%s]: %s%d missing panel citation(s) for: %s (defined in figures.yaml but never cited in body)",
            label, prefix, nn, paste(missing, collapse = ", ")))
      }
      # --- (4) reverse check: every cited letter is defined in yaml -----
      unknown <- setdiff(seen, yaml_letters)
      if (length(unknown) > 0) {
        prefix <- if (section == "supplement") "Extended Data Fig. " else "Figure "
        violations[[length(violations) + 1]] <<- list(
          reason = sprintf(
            "Rule 4 [%s]: %s%d cites panel(s) %s but figures.yaml defines only: %s",
            label, prefix, nn, paste(sort(unknown), collapse = ", "),
            paste(yaml_letters, collapse = ", ")))
      }
    }
  }

  # Rule 5 — every yaml panel must resolve to an existing source PDF,
  # mirroring compose_figure.R's resolution logic so the lint hard-fails
  # in the same cases the composite silently drops panels.
  resolve_panel_source <- function(entry, panel) {
    panels_dir <- entry$panels_dir
    src <- if (!is.null(panel$source)) panel$source else panel$letter
    if (is.null(src) || !nzchar(as.character(src))) return(list(ok = FALSE, tried = "(no source or letter)"))
    src <- as.character(src)
    if (grepl("/", src, fixed = TRUE)) {
      full <- if (requireNamespace("here", quietly = TRUE)) here::here(src) else src
      return(list(ok = file.exists(full), tried = full))
    }
    if (is.null(panels_dir) || !nzchar(panels_dir))
      return(list(ok = FALSE, tried = sprintf("(source='%s' but no panels_dir)", src)))
    pd <- if (requireNamespace("here", quietly = TRUE)) here::here(panels_dir) else panels_dir
    if (!dir.exists(pd))
      return(list(ok = FALSE, tried = sprintf("%s (directory missing)", pd)))
    pat <- paste0("^panel_[", tolower(src), toupper(src), "]([._]|$)")
    fs  <- list.files(pd, pattern = pat, full.names = FALSE)
    list(ok = length(fs) > 0,
         tried = sprintf("%s/panel_%s.*", pd, src))
  }

  check_composite_sources <- function() {
    for (id in names(.REG$figures)) {
      e <- .REG$figures[[id]]
      panels <- e$panels %||% list()
      if (length(panels) == 0) next
      n <- .FIG_NUM$figs[[id]]
      prefix <- if (identical(e$section, "supplement")) "Extended Data Fig. " else "Figure "
      label <- if (!is.null(n)) paste0(prefix, n) else id
      for (p in panels) {
        r <- resolve_panel_source(e, p)
        if (!isTRUE(r$ok)) {
          violations[[length(violations) + 1]] <<- list(
            reason = sprintf(
              "Rule 5 [composite]: %s panel %s (id=%s) not resolvable \u2192 %s",
              label, p$letter %||% "?", p$id %||% "?", r$tried))
        }
      }
    }
  }

  # Bounds check — cited number must not exceed total assigned
  n_main_fig <- .FIG_NUM$fig_main_n
  n_supp_fig <- .FIG_NUM$fig_supp_n
  n_main_tbl <- .FIG_NUM$tbl_main_n
  n_supp_tbl <- .FIG_NUM$tbl_supp_n

  for (sect in c("main", "supplement")) {
    figs <- Filter(function(x) x$kind == "figure" && x$section == sect, cites)
    tbls <- Filter(function(x) x$kind == "table"  && x$section == sect, cites)
    check_first_mention_order(figs, paste0(sect, " figures"))
    check_first_mention_order(tbls, paste0(sect, " tables"))
    check_panel_interleave(figs, paste0(sect, " figures"))
    check_panel_interleave(tbls, paste0(sect, " tables"))
    check_within_figure_panel_order(figs, paste0(sect, " figures"), sect)
  }

  # Rule 5 — run once across the whole registry (not section-dependent)
  check_composite_sources()

  for (c in cites) {
    cap <- if      (c$kind == "figure" && c$section == "main")       n_main_fig
           else if (c$kind == "figure" && c$section == "supplement") n_supp_fig
           else if (c$kind == "table"  && c$section == "main")       n_main_tbl
           else                                                        n_supp_tbl
    if (c$number > cap)
      violations[[length(violations) + 1]] <- list(
        reason = sprintf("Bounds: %s cited but only %d %s %s assigned",
                         c$raw, cap, c$section,
                         if (c$kind == "figure") "figures" else "tables"))
  }

  if (length(violations) == 0) {
    message("[registry_lint] \u2713 citation order OK (", length(cites), " refs checked)")
    return(invisible(list()))
  }

  msg <- paste(
    c("[registry_lint] \u2717 violations:",
      vapply(violations, function(v) paste0("  \u2022 ", v$reason), character(1))),
    collapse = "\n")
  if (strict) stop(msg, call. = FALSE)
  warning(msg, call. = FALSE)
  invisible(violations)
}

# ── Null coalesce ─────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (is.null(a)) b else a
