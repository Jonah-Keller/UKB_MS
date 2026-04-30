#!/usr/bin/env Rscript
# 00_merge_prs_scores.R
# Merge per-chromosome plink1 .profile files from RAP into a single PRS file
#
# Run AFTER dx download of /data/prs/results/*.profile files
#
# Usage:
#   Rscript 00_merge_prs_scores.R
#
# Input:  data/ukb/genetics/prs_profiles/prs_{PGS_ID}_c{CHR}.profile
# Output: data/ukb/genetics/ms_prs_scores.csv  (eid, prs_PGS000809, prs_PGS004699)

suppressPackageStartupMessages(library(data.table))

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

PROFILE_DIR <- file.path(PROJ_DIR, "data", "ukb", "genetics", "prs_profiles")
OUT_FILE    <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")

if (!dir.exists(PROFILE_DIR)) {
    message("═══════════════════════════════════════════════════════════════")
    message("Profile directory not found: ", PROFILE_DIR)
    message("")
    message("First download RAP results:")
    message("  mkdir -p data/ukb/genetics/prs_profiles")
    message("  dx download '/data/prs/results/*.profile' -o data/ukb/genetics/prs_profiles/")
    message("")
    message("Check RAP job status:")
    message("  dx find jobs --state done --name 'Swiss Army Knife' | head -50")
    message("═══════════════════════════════════════════════════════════════")
    quit(save = "no", status = 0)
}

merge_prs <- function(pgs_id) {
    pattern <- file.path(PROFILE_DIR, sprintf("prs_%s_c*.profile", pgs_id))
    files   <- Sys.glob(pattern)
    if (length(files) == 0) {
        message(sprintf("No .profile files found for %s (pattern: %s)", pgs_id, pattern))
        return(NULL)
    }
    cat(sprintf("%s: merging %d chromosomes...\n", pgs_id, length(files)))

    # plink1 .profile columns: FID IID PHENO CNT CNT2 SCORE
    # SCORE = sum(weight * allele_count) / (2 * n_scored)
    # To get total PRS, need to re-weight by n_scored per chromosome
    # Simpler: sum SCORE * CNT2 across chromosomes, divide by total CNT2
    chrom_list <- lapply(files, function(f) {
        dt <- fread(f, showProgress = FALSE)
        setnames(dt, tolower(names(dt)))
        # Columns: fid iid pheno cnt cnt2 score
        if (!"iid" %in% names(dt)) setnames(dt, 2, "iid")
        dt[, eid := as.integer(iid)]
        dt[, weighted_sum := score * cnt2]
        dt[, .(eid, cnt2, weighted_sum)]
    })
    combined <- rbindlist(chrom_list)
    agg <- combined[, .(
        total_cnt2 = sum(cnt2, na.rm = TRUE),
        total_wsum = sum(weighted_sum, na.rm = TRUE)
    ), by = eid]
    agg[, prs_raw := total_wsum / total_cnt2]
    # Standardize to mean=0, SD=1
    agg[, prs_z := (prs_raw - mean(prs_raw, na.rm=TRUE)) / sd(prs_raw, na.rm=TRUE)]
    setnames(agg, "prs_z", paste0("prs_", pgs_id))
    agg[, .(eid, get(paste0("prs_", pgs_id)))]
}

cat("Merging PRS scores from per-chromosome plink1 profiles...\n")

prs809 <- merge_prs("PGS000809")
prs699 <- merge_prs("PGS004699")

if (is.null(prs809) && is.null(prs699)) {
    stop("No profile files found for either PRS. Run dx download first.")
}

# Combine
if (!is.null(prs809) && !is.null(prs699)) {
    setnames(prs809, 2, "prs_PGS000809")
    setnames(prs699, 2, "prs_PGS004699")
    prs_all <- merge(prs809, prs699, by = "eid", all = TRUE)
} else {
    prs_all <- if (!is.null(prs809)) {
        setnames(prs809, 2, "prs_PGS000809"); prs809
    } else {
        setnames(prs699, 2, "prs_PGS004699"); prs699
    }
}

# Also create a combined column: prs_score = mean of available z-scores
score_cols <- intersect(c("prs_PGS000809","prs_PGS004699"), names(prs_all))
prs_all[, prs_score := rowMeans(.SD, na.rm=TRUE), .SDcols = score_cols]

fwrite(prs_all, OUT_FILE)
cat(sprintf("Saved: %d participants, columns: %s\n", nrow(prs_all), paste(names(prs_all), collapse=", ")))

# Quick QC
cat(sprintf("PRS distributions (z-scored):\n"))
for (col in score_cols) {
    cat(sprintf("  %s: mean=%.3f, SD=%.3f, range=[%.2f, %.2f]\n",
                col, mean(prs_all[[col]], na.rm=TRUE), sd(prs_all[[col]], na.rm=TRUE),
                min(prs_all[[col]], na.rm=TRUE), max(prs_all[[col]], na.rm=TRUE)))
}
cat("\n00_merge_prs_scores.R complete.\n")
