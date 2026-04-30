# ============================================================================
# render_figures_yaml.R — Jinja-style placeholder substitution for figures.yaml
# ============================================================================
# Reads configs/figures.yaml as text, substitutes {{var}} placeholders using
# scalar values from configs/disease.yaml, and returns the rendered YAML.
# Called by code/registry/figure_registry.R::registry_load() before parsing.
#
# Supported placeholders are any scalar (length-1 character/numeric) field in
# disease.yaml. Lists and nested objects are skipped (can't be inlined as
# strings into YAML structure).
#
# Example: a caption "{{disease_short_caps}} vs HC" with cohort_short="ms"
#          and disease_short_caps="MS" renders to "MS vs HC".
#
# Re-running this function is idempotent — substitutions only fire when the
# placeholder text is present.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
})

source(here("analysis", "helpers", "disease_config.R"))

render_figures_yaml <- function(yaml_path = here("configs", "figures.yaml"),
                                cfg = NULL) {
  if (is.null(cfg)) cfg <- load_disease_config()
  text <- paste(readLines(yaml_path, warn = FALSE), collapse = "\n")

  # Substitute every scalar field from disease.yaml. Skip lists/nested objects.
  for (key in names(cfg)) {
    val <- cfg[[key]]
    if (length(val) != 1L) next
    if (!is.character(val) && !is.numeric(val) && !is.logical(val)) next
    placeholder <- paste0("\\{\\{\\s*", key, "\\s*\\}\\}")
    text <- gsub(placeholder, as.character(val), text, perl = TRUE)
  }

  # Warn on any unresolved placeholders so typos surface immediately
  unresolved <- regmatches(text,
                           gregexpr("\\{\\{[^}]+\\}\\}", text, perl = TRUE))[[1]]
  if (length(unresolved) > 0) {
    warning("render_figures_yaml: ", length(unresolved),
            " unresolved placeholder(s): ",
            paste(unique(unresolved), collapse = ", "), call. = FALSE)
  }

  text
}

if (sys.nframe() == 0L) {
  rendered <- render_figures_yaml()
  cat("rendered length:", nchar(rendered), "chars\n")
  cat("first 400 chars:\n")
  cat(substr(rendered, 1L, 400L), "\n")
}
