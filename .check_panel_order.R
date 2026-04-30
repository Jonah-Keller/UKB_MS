#!/usr/bin/env Rscript
# Verify panel-level ascending citation order in the rendered manuscript.md.
#
# Rules enforced:
#   1. Within a single figure, the first citation of panel X must precede the
#      first citation of panel Y whenever X < Y in the alphabet.
#      (e.g., Figure 2a must first appear before Figure 2b's first appearance.)
#   2. Across figures, the first citation of Figure N must precede the first
#      citation of Figure N+1. (delegated to existing registry_lint Rule 1)
#   3. Panels of Figure N+1 must not interleave before Figure N's last panel.
#      (delegated to existing registry_lint Rule 2)
#
# This script complements registry_lint by adding the within-figure panel
# ordering check, which the existing linter does not enforce.

suppressPackageStartupMessages({
  library(here)
})

md_path <- here("manuscript", "submission", "_build_intermediates", "manuscript.md")
stopifnot(file.exists(md_path))

# Only inspect the body text, not the Table/Figure Legends section at the end,
# because the legend block repeats every panel in alphabetical order by design.
md <- readLines(md_path, warn = FALSE)
legend_start <- grep("^\\*\\*(?:Table|Figure) Legends\\*\\*", md)
if (length(legend_start) > 0) md <- md[seq_len(legend_start[1] - 1L)]
txt <- paste(md, collapse = "\n")

# Citation regex: bolded figure/table refs with panel letters.
# Matches "**Figure 2b**", "**Figure 2f\u2013h**", "**Extended Data Fig. 5a**", etc.
pat <- "\\*\\*((?:Supplementary Table|Extended Data Fig\\.|Figure|Table))\\s+(\\d+)([a-zA-Z])?(?:[\u2013\\-]([a-zA-Z]))?\\*\\*"

matches <- gregexpr(pat, txt, perl = TRUE)[[1]]
if (matches[1] == -1) {
  cat("[check_panel_order] no bolded citations found\n")
  quit(status = 0)
}

cites <- list()
for (i in seq_along(matches)) {
  s  <- matches[i]
  ln <- attr(matches, "match.length")[i]
  raw <- substring(txt, s, s + ln - 1L)
  mm <- regmatches(raw, regexec(pat, raw, perl = TRUE))[[1]]
  prefix     <- mm[2]
  number     <- as.integer(mm[3])
  panel_from <- if (nzchar(mm[4])) tolower(mm[4]) else NA_character_
  panel_to   <- if (nzchar(mm[5])) tolower(mm[5]) else panel_from

  # Figure-like vs table-like grouping
  is_fig <- grepl("Fig|Figure", prefix)
  section <- if (grepl("^(?:Supplementary|Extended)", prefix)) "supplement" else "main"
  fig_key <- paste(prefix, number, sep = "|")

  cites[[length(cites) + 1L]] <- list(
    pos = s, raw = raw, prefix = prefix, number = number,
    section = section, is_fig = is_fig, fig_key = fig_key,
    panel_from = panel_from, panel_to = panel_to)
}

# ── Rule 3: within-figure panel ordering ─────────────────────────────────────
violations <- list()

fig_keys <- unique(vapply(cites, `[[`, character(1), "fig_key"))
for (fk in fig_keys) {
  sub <- Filter(function(c) c$fig_key == fk, cites)
  panel_cites <- Filter(function(c) !is.na(c$panel_from), sub)
  if (length(panel_cites) < 2) next

  # First-seen position of each panel letter (based on panel_from; ranges are
  # expanded to every contained letter since each member is 'introduced' by
  # the range citation).
  first_seen <- list()
  for (pc in panel_cites) {
    letters_in_cite <- {
      a <- utf8ToInt(pc$panel_from)
      b <- utf8ToInt(pc$panel_to)
      intToUtf8(seq(a, b), multiple = TRUE)
    }
    for (L in letters_in_cite) {
      if (is.null(first_seen[[L]])) first_seen[[L]] <- pc$pos
    }
  }

  # Check ordering: sort panels alphabetically; their first-seen positions
  # must be strictly increasing.
  ordered_letters <- sort(names(first_seen))
  positions <- vapply(ordered_letters, function(L) first_seen[[L]], numeric(1))
  # Find inversions: if a later (alphabetically) panel is first seen before an
  # earlier one.
  for (i in 2:length(ordered_letters)) {
    for (j in 1:(i - 1)) {
      if (positions[i] < positions[j]) {
        violations[[length(violations) + 1L]] <- list(
          reason = sprintf(
            "Panel order [%s]: panel %s first cited before panel %s (pos %d vs %d)",
            fk, ordered_letters[i], ordered_letters[j],
            positions[i], positions[j]))
      }
    }
  }
}

# ── Summary ──────────────────────────────────────────────────────────────────
# List of all citations in text order, grouped by figure, for manual inspection
cat("══ Citation audit (body text only) ══\n\n")
for (fk in fig_keys) {
  sub <- Filter(function(c) c$fig_key == fk, cites)
  panel_cites <- Filter(function(c) !is.na(c$panel_from), sub)
  if (length(panel_cites) == 0) next
  panel_seq <- vapply(panel_cites, function(c) {
    if (c$panel_from == c$panel_to) c$panel_from
    else paste0(c$panel_from, "\u2013", c$panel_to)
  }, character(1))
  cat(sprintf("  %-28s  %s\n", fk, paste(panel_seq, collapse = ", ")))
}

cat("\n")
if (length(violations) == 0) {
  cat("\u2713 All panel citations appear in ascending alphabetical order within each figure.\n")
  quit(status = 0)
} else {
  cat("\u2717 ", length(violations), " panel-ordering violation(s):\n", sep = "")
  for (v in violations) cat("  \u2022 ", v$reason, "\n", sep = "")
  quit(status = 1)
}
