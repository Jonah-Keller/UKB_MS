#!/usr/bin/env Rscript
# ============================================================================
# 00_config.R — UKB MS proteomics project configuration
# ============================================================================

suppressPackageStartupMessages({ library(here) })

PROJECT_TITLE   <- "Pre-diagnostic comorbidity endophenotypes reveal proteomic and genetic heterogeneity in multiple sclerosis"
PROJECT_SHORT   <- "ukb_ms_proteomics"
PROJECT_AUTHORS <- c(
  "Jonah N. Keller1,2",
  "Jolie A. Kantor3",
  "Fanny M. Elahi2*"
)
PROJECT_AFFIL <- c(
  "1 Perelman School of Medicine at the University of Pennsylvania, Philadelphia, Pennsylvania, USA.",
  "2 Departments of Neurology, Neuroscience, and AI and Human Health, Icahn School of Medicine at Mount Sinai, New York, New York, USA.",
  "3 New York University School of Medicine, New York, New York, USA."
)
PROJECT_CORRESPONDING <- list(
  name  = "Fanny M. Elahi, MD, PhD",
  email = "fanny.elahi@mssm.edu"
)
TARGET_JOURNAL <- "Nature Medicine"

# ── Paths ────────────────────────────────────────────────────────────────────
DIR_DATA_RAW        <- here("data", "raw")
DIR_DATA_PROCESSED  <- here("data", "ukb")
DIR_RESULTS         <- here("results")
DIR_RESULTS_FIGURES <- here("results", "figures")

DIR_SUBMISSION      <- here("manuscript", "submission")
DIR_MAIN_TEXT       <- file.path(DIR_SUBMISSION, "main_text")
DIR_MAIN_FIGURES    <- file.path(DIR_SUBMISSION, "main_figures")
DIR_MAIN_TABLES     <- file.path(DIR_SUBMISSION, "main_tables")
DIR_SUPP            <- file.path(DIR_SUBMISSION, "supplementary")
DIR_SUPP_EXTDATA   <- file.path(DIR_SUPP, "extended_data")
DIR_SUPP_TABLES    <- file.path(DIR_SUPP, "supplementary_tables")
DIR_BUILD_INT       <- file.path(DIR_SUBMISSION, "_build_intermediates")

DIR_TEMPLATES  <- here("manuscript", "templates")
REFERENCE_DOCX <- file.path(DIR_TEMPLATES, "reference.docx")
FIGURES_YAML   <- here("figures.yaml")

for (.d in c(DIR_MAIN_TEXT, DIR_MAIN_FIGURES, DIR_MAIN_TABLES,
             DIR_SUPP_EXTDATA, DIR_SUPP_TABLES, DIR_BUILD_INT))
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)

# ── Analysis parameters ──────────────────────────────────────────────────────
ALPHA      <- 0.05
FDR_METHOD <- "BH"
SEED       <- 42
set.seed(SEED)

# ── Style ────────────────────────────────────────────────────────────────────
source(here("analysis", "helpers", "ukb_theme.R"))

cat("[00_config] project =", PROJECT_SHORT, "| target =", TARGET_JOURNAL, "\n")
