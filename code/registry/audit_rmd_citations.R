# Source-level audit: parse manuscript.Rmd and figures.yaml directly to
# reconstruct the text-order citation sequence WITHOUT a full knit.
# Useful for quickly re-checking panel order after .Rmd edits.

suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

rmd_path  <- here("analysis", "manuscript", "manuscript.Rmd")
yaml_path <- here("configs", "figures.yaml")
rmd <- readLines(rmd_path, warn = FALSE)

# Strip YAML front matter (first two "---" lines)
yaml_ends <- which(rmd == "---")
if (length(yaml_ends) >= 2) rmd <- rmd[-(1:yaml_ends[2])]

# Strip R code fences
in_code <- FALSE
keep <- logical(length(rmd))
for (i in seq_along(rmd)) {
  if (grepl("^```", rmd[i])) { in_code <- !in_code; keep[i] <- FALSE }
  else keep[i] <- !in_code
}
rmd <- rmd[keep]

# Cut off at the Figure / Table Legends heading — legends list panels
# alphabetically by design and would pollute the ordering audit.
legend_hdr <- grep("^##\\s+(Figure|Table)\\s+Legends", rmd, perl = TRUE)
if (length(legend_hdr) > 0) rmd <- rmd[seq_len(legend_hdr[1] - 1L)]

txt <- paste(rmd, collapse = "\n")

# Load registry for panel lookup
reg <- yaml::read_yaml(yaml_path)
fig_entries <- reg$figures
names(fig_entries) <- vapply(fig_entries, `[[`, character(1), "id")

# Resolve a (fig_id, panel_arg) to a panel letter
resolve_letter <- function(fig_id, panel_arg) {
  if (is.na(panel_arg) || !nzchar(panel_arg)) return(NA_character_)
  if (panel_arg %in% c(letters, LETTERS)) return(tolower(panel_arg))
  fe <- fig_entries[[fig_id]]
  if (is.null(fe)) return(NA_character_)
  hit <- Filter(function(p) identical(p$id, panel_arg), fe$panels)
  if (length(hit) == 0) return(NA_character_)
  tolower(hit[[1]]$letter)
}

# Extract all fig() calls in text order.
# Patterns: fig("ID") or fig("ID", "PANEL")
pat <- "fig\\(\\s*\"([^\"]+)\"(?:\\s*,\\s*\"([^\"]+)\")?\\s*\\)"
m <- gregexpr(pat, txt, perl = TRUE)[[1]]
starts <- as.integer(m)
lens   <- attr(m, "match.length")

# Assign figure numbers in order of first citation (mimics fig_registry)
main_n <- 0L; supp_n <- 0L
num_for_id <- list()
cites <- list()

for (i in seq_along(starts)) {
  raw <- substring(txt, starts[i], starts[i] + lens[i] - 1L)
  gp  <- regmatches(raw, regexec(pat, raw, perl = TRUE))[[1]]
  fig_id   <- gp[2]
  panel_arg <- if (nzchar(gp[3])) gp[3] else NA_character_
  fe <- fig_entries[[fig_id]]
  if (is.null(fe)) next
  section <- fe$section
  if (is.null(num_for_id[[fig_id]])) {
    if (section == "main") { main_n <- main_n + 1L; num_for_id[[fig_id]] <- main_n }
    else                   { supp_n <- supp_n + 1L; num_for_id[[fig_id]] <- supp_n }
  }
  n <- num_for_id[[fig_id]]
  letter <- resolve_letter(fig_id, panel_arg)
  prefix <- if (section == "main") "Figure" else "Extended Data Fig."
  cites[[length(cites) + 1L]] <- list(
    pos = starts[i], fig_id = fig_id, section = section,
    number = n, letter = letter,
    rendered = paste0(prefix, " ", n, if (!is.na(letter)) letter else ""))
}

# Per-figure audit: ordering (3a) + completeness (3b)
cat("══ Source-level citation audit (Rmd body, pre-knit) ══\n\n")
fig_keys <- unique(vapply(cites, function(x) paste(x$section, x$number, sep = ":"), character(1)))

# Reverse map: (section, number) -> fig_id for yaml-panel lookup
key_to_fig_id <- list()
for (nm in names(num_for_id)) {
  fe <- fig_entries[[nm]]
  n  <- num_for_id[[nm]]
  key_to_fig_id[[paste(fe$section, n, sep = ":")]] <- nm
}

violations <- list()
for (fk in fig_keys) {
  sub <- Filter(function(c) paste(c$section, c$number, sep = ":") == fk, cites)
  panel_cites <- Filter(function(c) !is.na(c$letter), sub)
  if (length(panel_cites) == 0) next

  label <- panel_cites[[1]]$rendered
  label <- sub("[a-z]$", "", label)
  letter_seq <- vapply(panel_cites, `[[`, character(1), "letter")
  cat(sprintf("  %-30s  %s\n", label, paste(letter_seq, collapse = ", ")))

  # Rule 3a — ascending-order
  seen <- character()
  max_letter <- ""
  for (pc in panel_cites) {
    L <- pc$letter
    if (L %in% seen) next
    if (nzchar(max_letter) && L < max_letter)
      violations[[length(violations) + 1L]] <- list(
        rule = "3a", figure = label,
        msg  = sprintf("panel %s first cited after panel %s", L, max_letter))
    seen <- c(seen, L)
    if (L > max_letter) max_letter <- L
  }

  # Rule 3b — completeness: every yaml panel must be cited
  fig_id <- key_to_fig_id[[fk]]
  if (!is.null(fig_id)) {
    yaml_letters <- tolower(vapply(fig_entries[[fig_id]]$panels,
                                   function(p) as.character(p$letter %||% ""),
                                   character(1)))
    yaml_letters <- yaml_letters[nzchar(yaml_letters)]
    missing <- setdiff(yaml_letters, seen)
    if (length(missing) > 0)
      violations[[length(violations) + 1L]] <- list(
        rule = "3b", figure = label,
        msg  = sprintf("missing citation(s) for panel(s): %s (defined in figures.yaml but never cited)",
                       paste(missing, collapse = ", ")))
  }
}

cat("\n")
if (length(violations) == 0) {
  cat("\u2713 All panels cited in ascending order and completely.\n")
  quit(status = 0)
}
cat("\u2717 ", length(violations), " violation(s):\n", sep = "")
for (v in violations)
  cat(sprintf("  \u2022 Rule %-3s %s: %s\n", v$rule, v$figure, v$msg))
quit(status = 1)
