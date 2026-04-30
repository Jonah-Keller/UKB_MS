# manuscript_stats.R
# ----------------------------------------------------------------------------
# Precompute numeric statistics that appear inline in Results prose.
# Exposes two entry points:
#   manuscript_stats_load(proj_dir)  — call once in Rmd setup chunk
#   st(key)                          — inline retrieval, errors on unknown key
#
# Convention: values are stored as pre-formatted character strings
# (thousands separators applied) so that `r st("n_hc")` renders as "47,801".
# Raw numeric variants are stored under the same key prefixed with "raw_".
# ----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(here)
    library(data.table)
})

.STATS <- new.env(parent = emptyenv())

.fmt_int <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
    format(as.integer(x), big.mark = ",", trim = TRUE, scientific = FALSE)
}

st <- function(key) {
    if (!exists(key, envir = .STATS))
        stop("manuscript stat '", key, "' not defined; ",
             "add it to manuscript_stats.R")
    .STATS[[key]]
}

.set <- function(key, value, raw = NULL) {
    .STATS[[key]] <- value
    if (!is.null(raw)) .STATS[[paste0("raw_", key)]] <- raw
}

manuscript_stats_load <- function(proj_dir = here()) {
    rassign <- function(path) fread(file.path(proj_dir, path), showProgress = FALSE)

    # ── Cohort counts (QC-filtered) ──────────────────────────────────────────
    qc <- rassign("data/ukb/olink/processed/ms_olink_qc.csv")
    meta_cols <- c("eid","ms_status","age_at_sampling","age_at_diagnosis",
                   "years_to_diagnosis","sex","olink_instance","qc_outlier",
                   "UMAP1","UMAP2","mean_npx")
    protein_cols <- setdiff(names(qc), meta_cols)
    qc <- qc[qc_outlier == FALSE]

    .set("n_pre_onset",  .fmt_int(sum(qc$ms_status == "pre_onset")))
    .set("n_post_onset", .fmt_int(sum(qc$ms_status == "post_onset")))
    .set("n_hc",         .fmt_int(sum(qc$ms_status == "control")))
    .set("n_ms_total",   .fmt_int(sum(qc$ms_status != "control")))
    .set("n_proteins_retained",   .fmt_int(length(protein_cols)))
    .set("n_proteins_full_panel", .fmt_int(2923))  # Olink Explore 3072 panel

    # Sex in MS cases (UKB: 0=Female, 1=Male)
    ms <- qc[ms_status != "control" & !is.na(sex)]
    .set("n_female_ms", .fmt_int(sum(ms$sex == 0)))
    .set("n_male_ms",   .fmt_int(sum(ms$sex == 1)))

    # ── DEP counts ───────────────────────────────────────────────────────────
    pre_de  <- rassign("results/differential/ms_pre_vs_hc.csv")
    post_de <- rassign("results/differential/ms_post_vs_hc.csv")
    comb_de <- tryCatch(rassign("results/differential/ms_combined_vs_hc.csv"),
                        error = function(e) NULL)
    .set("n_pre_dep_fdr05",  .fmt_int(sum(pre_de$adj.P.Val  < 0.05, na.rm = TRUE)))
    .set("n_post_dep_fdr05", .fmt_int(sum(post_de$adj.P.Val < 0.05, na.rm = TRUE)))
    if (!is.null(comb_de))
        .set("n_combined_dep_fdr05",
             .fmt_int(sum(comb_de$adj.P.Val < 0.05, na.rm = TRUE)))

    # ── Cox survival (pre-diagnostic) ────────────────────────────────────────
    cox <- rassign("results/survival/ms_protein_cox_results.csv")
    fdr_col <- if ("fdr" %in% names(cox)) "fdr" else "adj.P.Val"
    .set("n_cox_fdr05", .fmt_int(sum(cox[[fdr_col]] < 0.05, na.rm = TRUE)))

    # ── Sex-stratified concordance ───────────────────────────────────────────
    fem <- rassign("results/sex_stratified/ms_female_vs_hc.csv")
    mal <- rassign("results/sex_stratified/ms_male_vs_hc.csv")
    fem_sig <- fem[adj.P.Val < 0.05, protein]
    mal_sig <- mal[adj.P.Val < 0.05, protein]
    .set("n_sex_concordant", .fmt_int(length(intersect(fem_sig, mal_sig))))
    .set("n_female_dep_fdr05", .fmt_int(length(fem_sig)))
    .set("n_male_dep_fdr05",   .fmt_int(length(mal_sig)))

    # ── PheWAS of pre-MS ICD-10 conditions ───────────────────────────────────
    phew <- rassign("results/comorbidity/ms_prems_phewas.csv")
    sig  <- phew[fdr < 0.05]
    .set("n_prems_phewas_fdr05",    .fmt_int(nrow(sig)))
    .set("n_prems_phewas_enriched", .fmt_int(sum(sig$OR > 1, na.rm = TRUE)))
    .set("n_prems_phewas_depleted", .fmt_int(sum(sig$OR < 1, na.rm = TRUE)))

    # ── Endophenotype clusters ───────────────────────────────────────────────
    clus <- rassign("results/endophenotype/ms_prems_cluster_assignments.csv")
    clus_valid <- clus[!is.na(cluster) & as.character(cluster) != ""]
    .set("n_clusters", .fmt_int(length(unique(clus_valid$cluster))))
    .set("n_clustered_preonset", .fmt_int(nrow(clus_valid)))
    for (k in sort(unique(as.integer(clus_valid$cluster)))) {
        .set(paste0("n_c", k, "_cases"),
             .fmt_int(sum(as.integer(clus_valid$cluster) == k)))
    }

    # ── GMM co-presenter rate (fraction spanning multiple clusters) ──────────
    copres <- tryCatch(
        rassign("results/endophenotype/gmm_copresentation/gmm_copresenters.csv"),
        error = function(e) NULL)
    if (!is.null(copres) && "co_presenter" %in% names(copres)) {
        pct_copres <- if (nrow(copres) > 0)
            100 * mean(copres$co_presenter, na.rm = TRUE) else NA_real_
        .set("pct_copresenters",
             if (is.na(pct_copres))         "NA"
             else if (pct_copres == 0)      "0"
             else if (pct_copres < 1)       "<1"
             else format(round(pct_copres, 1), trim = TRUE, nsmall = 1))
    }

    # ── HLA-DRB1*15:01 carrier vs non-carrier DEPs (within MS cases) ─────────
    hla <- rassign("results/hla/ms_hla_subgroup_results.csv")
    .set("n_hla_dep_fdr05", .fmt_int(sum(hla$fdr < 0.05, na.rm = TRUE)))

    # ── HLA-stratified MS-vs-HC DEP counts ───────────────────────────────────
    hla_strat <- rassign("results/endophenotype/ms_hla_stratified_logfc.csv")
    .set("n_hla_carrier_msdep_fdr05",
         .fmt_int(sum(hla_strat$fdr_hlapos < 0.05, na.rm = TRUE)))
    .set("n_hla_noncarrier_msdep_fdr05",
         .fmt_int(sum(hla_strat$fdr_hlaneg < 0.05, na.rm = TRUE)))

    # ── PRS-correlated proteins (population scale) ───────────────────────────
    prs <- rassign("results/prs/ms_prs_results.csv")
    .set("n_prs_proteome_fdr05", .fmt_int(sum(prs$fdr < 0.05, na.rm = TRUE)))

    # ── Spline divergence timepoints ─────────────────────────────────────────
    # Source: results/differential/cns_trajectories_divergence.csv written by
    # analysis/04_differential/03_cns_trajectories.R. diverge_yr is the first
    # grid point where the MS spline 95% CI separates from the HC mean ±1.96 SE
    # band; recovery_yr is the first post-onset grid point where the divergence
    # condition no longer holds (suppression lifts or elevation resolves).
    fmt_yr <- function(x) {
        if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
        if (isTRUE(all.equal(x, round(x)))) format(round(x), trim = TRUE)
        else format(round(x, 2), trim = TRUE, nsmall = 2)
    }
    div <- rassign("results/differential/cns_trajectories_divergence.csv")
    nefl  <- div[toupper(protein) == "NEFL"]
    gfap  <- div[toupper(protein) == "GFAP"]
    erbb2 <- div[toupper(protein) == "ERBB2"]
    .set("nefl_divergence_years",
         fmt_yr(if (nrow(nefl)  > 0) abs(nefl$diverge_yr)    else NA_real_))
    .set("gfap_emergence_years",
         fmt_yr(if (nrow(gfap)  > 0) gfap$diverge_yr         else NA_real_))
    .set("erbb2_recovery_years",
         fmt_yr(if (nrow(erbb2) > 0) erbb2$recovery_yr       else NA_real_))

    invisible(.STATS)
}
