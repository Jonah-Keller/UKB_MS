# ============================================================================
# disease_config.R — R-side loader for configs/disease.yaml
# ============================================================================
# To replicate this study on a different cohort, edit configs/disease.yaml
# (do not modify any analysis script). Every R script that needs a disease-,
# HLA-, or PRS-specific value sources this file and reads from `cfg`.
#
# Usage:
#   source(here::here("analysis", "helpers", "disease_config.R"))
#   cfg <- load_disease_config()
#   cohort_short <- cfg$cohort_short
#   icd_codes    <- cfg$icd_codes
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

.disease_cfg_cache <- NULL

load_disease_config <- function(path = here::here("configs", "disease.yaml"),
                                refresh = FALSE) {
  if (!refresh && !is.null(.disease_cfg_cache)) return(.disease_cfg_cache)
  cfg <- yaml::read_yaml(path)

  cfg$all_exclusion_codes <- unique(unlist(cfg$control_exclusion_codes,
                                           use.names = FALSE))

  assign(".disease_cfg_cache", cfg, envir = topenv())
  cfg
}

if (sys.nframe() == 0L) {
  cfg <- load_disease_config()
  cat("cohort_short        =", cfg$cohort_short, "\n")
  cat("disease_long        =", cfg$disease_long, "\n")
  cat("icd_codes           =", paste(cfg$icd_codes, collapse = ","), "\n")
  cat("hla_allele          =", cfg$hla_allele, "\n")
  cat("prs_pgs_ids         =", paste(cfg$prs_pgs_ids, collapse = ","), "\n")
  cat("all_exclusion_codes =", length(cfg$all_exclusion_codes), "codes\n")
}
