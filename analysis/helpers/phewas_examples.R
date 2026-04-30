#!/usr/bin/env Rscript
# phewas_examples.R — Examples: variant / protein / cluster → UKB PheWAS
#
# Run any of these blocks independently after sourcing phewas_runner.R.
# Each example saves results CSV + Manhattan PDF to results/phewas/<label>/

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

# Resolve project root relative to this file
args     <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

source(file.path(PROJ_DIR, "analysis", "helpers", "phewas_runner.R"))

QC_FILE  <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
OUT_BASE <- file.path(PROJ_DIR, "results", "phewas")

# ── Shared cohort ─────────────────────────────────────────────────────────────
cat("Loading cohort...\n")
qc     <- fread(QC_FILE, showProgress = FALSE)
cohort <- qc[qc_outlier == FALSE & !is.na(age_at_sampling) & !is.na(sex)]
cohort_cov <- cohort[, .(eid, age = age_at_sampling, sex)]

# ═══════════════════════════════════════════════════════════════════════════════
# Example 1: Case-control PheWAS (MS status)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n── Example 1: MS case-control PheWAS ──────────────────────────────────\n")

trait_ms <- cohort[, .(eid, trait = as.integer(ms_status %in% c("pre_onset","post_onset")))]

res_ms <- run_phewas(
    trait       = trait_ms,
    cohort      = cohort_cov,
    min_cases   = 10L,
    propagate   = TRUE,
    timing_col  = NULL       # set to "age_at_diagnosis" to compute pre/post timing
)

save_phewas(res_ms, file.path(OUT_BASE, "ms_casecontrol"), "ms_casecontrol")

# ═══════════════════════════════════════════════════════════════════════════════
# Example 2: Protein PheWAS (any Olink protein as continuous trait)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n── Example 2: Protein PheWAS (GFAP) ───────────────────────────────────\n")

# Olink matrix is wide; NPX values are in columns named by protein
olink_file <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed",
                         "ms_olink_qc.csv")   # reuse QC file which has NPX columns
protein_cols <- setdiff(names(qc),
                         c("eid","ms_status","age_at_sampling","sex",
                           "age_at_diagnosis","years_to_diagnosis","qc_outlier",
                           "plate_id","sample_id","n_prot_missing","n_prot_outlier"))
# pick GFAP as example — swap to any protein name
EXAMPLE_PROTEIN <- "GFAP"
if (EXAMPLE_PROTEIN %in% names(qc)) {
    trait_prot <- cohort[, .(eid, trait = get(EXAMPLE_PROTEIN))]
    trait_prot <- trait_prot[!is.na(trait)]

    res_prot <- run_phewas(
        trait     = trait_prot,
        cohort    = cohort_cov,
        min_cases = 10L,
        propagate = TRUE
    )

    save_phewas(res_prot, file.path(OUT_BASE, paste0("protein_", EXAMPLE_PROTEIN)),
                paste0("protein_", EXAMPLE_PROTEIN))
} else {
    cat("  GFAP not found in QC file — skipping protein example\n")
    cat("  Available protein columns:", paste(head(protein_cols, 10), collapse=", "), "...\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Example 3: Cluster membership PheWAS (which ICD-10 phenotypes define cluster 3?)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n── Example 3: Cluster-3 membership PheWAS ─────────────────────────────\n")

clust_file <- file.path(PROJ_DIR, "results", "endophenotype",
                         "ms_prems_cluster_assignments.csv")
if (file.exists(clust_file)) {
    clust <- fread(clust_file, showProgress = FALSE)
    # Restrict to MS cases only, compare cluster 3 vs all other MS clusters
    ms_eids <- cohort[ms_status %in% c("pre_onset","post_onset"), eid]
    clust_ms <- clust[eid %in% ms_eids]

    for (clust_id in 0:4) {
        trait_clust <- clust_ms[, .(eid, trait = as.integer(cluster == clust_id))]
        n_clust <- sum(trait_clust$trait)
        if (n_clust < 20L) {
            cat(sprintf("  Cluster %d: n=%d — skipping (< 20)\n", clust_id, n_clust))
            next
        }
        cat(sprintf("  Running PheWAS for cluster %d (n=%d vs rest)\n", clust_id, n_clust))
        res_clust <- run_phewas(
            trait     = trait_clust,
            cohort    = cohort_cov,
            min_cases = 5L,
            propagate = TRUE,
            verbose   = FALSE
        )
        label <- sprintf("cluster_%d", clust_id)
        save_phewas(res_clust, file.path(OUT_BASE, label), label)
        cat(sprintf("  Cluster %d: %d tested, %d FDR<0.05\n",
                    clust_id, nrow(res_clust), sum(res_clust$fdr < 0.05)))
    }
} else {
    cat("  Cluster assignments file not found — run 15_endophenotype first\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Example 4: PRS PheWAS (what phenotypes co-associate with MS genetic risk?)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n── Example 4: MS PRS PheWAS ────────────────────────────────────────────\n")

prs_file <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
if (file.exists(prs_file)) {
    prs <- fread(prs_file, showProgress = FALSE)
    prs[, prs_z := scale(prs_score)[,1]]   # standardise
    trait_prs <- prs[, .(eid, trait = prs_z)]

    res_prs <- run_phewas(
        trait     = trait_prs,
        cohort    = cohort_cov,
        min_cases = 10L,
        propagate = TRUE
    )

    save_phewas(res_prs, file.path(OUT_BASE, "ms_prs"), "ms_prs")
} else {
    cat("  PRS file not found — skipping PRS example\n")
}

cat("\n\nAll examples complete. Results in:", OUT_BASE, "\n")
