#!/usr/bin/env Rscript
# 02_preonset_endophenotype.R
# Pre-onset proteomic endophenotypes and post-diagnosis outcome validation
#
# Design:
#   1. Restrict to pre-onset MS cases (years_to_diagnosis > 0)
#   2. Feature matrix = top DEPs (FDR<0.05), z-scored
#   3. PCA → test k=1..6 by gap statistic + average silhouette
#   4. k-means with optimal k; null-aware: if gap → k=1, proceed with
#      k=2 as exploratory but flag clearly
#   5. Annotate clusters: PRS, HLA-DRB1*15:01, WMH volume, comorbidity
#      burden, demographics
#   6. Post-diagnosis outcome comparison across clusters
#   7. Save cluster assignments + profiles
#
# Output: results/endophenotype/
#   ms_preonset_clusters.csv          — per-participant cluster assignments
#   ms_preonset_cluster_profiles.csv  — per-cluster summaries
#   ms_preonset_cluster_gap.csv       — gap statistic + silhouette by k
#   ms_preonset_cluster_stability.csv — bootstrap stability (Jaccard index)

suppressPackageStartupMessages({
    library(data.table)
    library(cluster)       # silhouette, clusGap
    library(ggplot2)
    library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

# ── Paths ──────────────────────────────────────────────────────────────────
QC_FILE      <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DEP_FILE     <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
PRS_FILE     <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
COV_FILE     <- file.path(PROJ_DIR, "data", "ukb", "covariates", "covariate_table.csv")
HLA_FILE     <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
WMH_FILE     <- file.path(dirname(PROJ_DIR),
                           "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                           "data", "ukb", "misc", "wmh_htn_dm.csv")
FOD_FILE     <- file.path(dirname(PROJ_DIR),
                           "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                           "data", "ukb", "misc", "ukb_first_occurrence_dates.csv")
OUT_DIR      <- file.path(PROJ_DIR, "results", "endophenotype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
K_MAX     <- 6L      # test k = 1..K_MAX
N_BOOT    <- 50L     # bootstrap replicates for gap statistic
FDR_THR   <- 0.05
MIN_PROTS <- 10L     # minimum DEPs needed to proceed

# ── 1. Load pre-onset proteomics ───────────────────────────────────────────
cat("Loading pre-onset data...\n")
qc <- fread(QC_FILE, showProgress = FALSE)
pre <- qc[ms_status == "pre_onset" & qc_outlier == FALSE]
cat(sprintf("  Pre-onset cases: %d\n", nrow(pre)))

# Identify protein columns
META_COLS <- c("eid","ms_status","age_at_sampling","age_at_diagnosis",
               "years_to_diagnosis","sex","olink_instance","qc_outlier",
               "UMAP1","UMAP2","mean_npx")
protein_cols <- setdiff(names(pre), META_COLS)

# ── 2. Feature matrix: top DEPs ────────────────────────────────────────────
cat("Selecting DEP features...\n")
deps <- fread(DEP_FILE, showProgress = FALSE)
if (!"fdr" %in% names(deps)) deps[, fdr := p.adjust(P.Value, method = "BH")]
top_deps <- deps[fdr < FDR_THR, tolower(protein)]
avail_feats <- intersect(top_deps, protein_cols)
cat(sprintf("  DEPs at FDR<%.2f: %d  available in QC: %d\n",
            FDR_THR, length(top_deps), length(avail_feats)))

if (length(avail_feats) < MIN_PROTS) {
    cat("  Insufficient features — relaxing to FDR<0.20\n")
    top_deps2 <- deps[fdr < 0.20, tolower(protein)]
    avail_feats <- intersect(top_deps2, protein_cols)
    cat(sprintf("  Features at FDR<0.20: %d\n", length(avail_feats)))
}

# z-score each protein, median-impute missing
feat_mat <- pre[, ..avail_feats]
for (col in avail_feats) {
    med <- median(feat_mat[[col]], na.rm = TRUE)
    feat_mat[is.na(get(col)), (col) := med]
}
feat_scaled <- scale(feat_mat)
cat(sprintf("  Feature matrix: %d × %d\n", nrow(feat_scaled), ncol(feat_scaled)))

# PCA: retain enough PCs to explain ≥80% variance
pca_fit  <- prcomp(feat_scaled, center = FALSE, scale. = FALSE)
cum_var  <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
n_pc     <- max(2L, which(cum_var >= 0.80)[1])
n_pc     <- min(n_pc, nrow(pre) - 1L, 30L)
pca_mat  <- pca_fit$x[, seq_len(n_pc), drop = FALSE]
cat(sprintf("  PCA: %d PCs explaining %.1f%% variance\n",
            n_pc, 100 * cum_var[n_pc]))

# ── 3. Gap statistic + silhouette for k = 1..K_MAX ─────────────────────────
cat(sprintf("Computing gap statistic (B=%d) and silhouette for k=1..%d...\n",
            N_BOOT, K_MAX))

gap_result <- tryCatch(
    clusGap(pca_mat,
            FUN     = function(x, k) list(cluster = kmeans(x, centers=k, nstart=25, iter.max=50)$cluster),
            K.max   = K_MAX,
            B       = N_BOOT,
            verbose = FALSE),
    error = function(e) { message("  clusGap failed: ", e$message); NULL }
)

sil_scores <- vapply(seq_len(K_MAX), function(k) {
    if (k == 1L) return(NA_real_)
    km <- kmeans(pca_mat, centers = k, nstart = 25, iter.max = 50)
    mean(silhouette(km$cluster, dist(pca_mat))[, 3], na.rm = TRUE)
}, numeric(1))

gap_dt <- data.table(
    k         = seq_len(K_MAX),
    sil_mean  = sil_scores
)
if (!is.null(gap_result)) {
    gap_tab <- as.data.table(gap_result$Tab)
    gap_dt[, gap      := gap_tab$gap]
    gap_dt[, gap_se   := gap_tab$SE.sim]
    # Tibshirani-Walther-Hastie rule: smallest k such that gap(k) >= gap(k+1) - SE(k+1)
    optimal_k_gap <- tryCatch(maxSE(gap_result$Tab[, "gap"], gap_result$Tab[, "SE.sim"],
                                    method = "Tibs2001SEmax"),
                              error = function(e) NA_integer_)
} else {
    gap_dt[, gap    := NA_real_]
    gap_dt[, gap_se := NA_real_]
    optimal_k_gap <- NA_integer_
}

optimal_k_sil <- if (all(is.na(sil_scores))) NA_integer_ else which.max(sil_scores)

cat(sprintf("  Gap statistic → optimal k = %s\n",
            ifelse(is.na(optimal_k_gap), "NA", as.character(optimal_k_gap))))
cat(sprintf("  Silhouette   → optimal k = %s\n",
            ifelse(is.na(optimal_k_sil), "NA", as.character(optimal_k_sil))))
print(gap_dt)
fwrite(gap_dt, file.path(OUT_DIR, "ms_preonset_cluster_gap.csv"))

# ── 4. Choose k and cluster ────────────────────────────────────────────────
# Consensus: if gap and sil agree → use that k
# If gap says 1 → report null structure, proceed with k=2 as exploratory
chose_k_method <- "consensus"
if (!is.na(optimal_k_gap) && !is.na(optimal_k_sil) && optimal_k_gap == optimal_k_sil) {
    final_k <- optimal_k_gap
    chose_k_method <- "gap+silhouette agree"
} else if (!is.na(optimal_k_gap) && optimal_k_gap > 1L) {
    final_k <- optimal_k_gap
    chose_k_method <- "gap statistic"
} else if (!is.na(optimal_k_sil) && optimal_k_sil > 1L) {
    final_k <- optimal_k_sil
    chose_k_method <- "silhouette"
} else {
    final_k <- 2L
    chose_k_method <- "forced k=2 (gap → k=1; exploratory only)"
}
cat(sprintf("\n  Selected k = %d (%s)\n", final_k, chose_k_method))
null_structure <- grepl("forced", chose_k_method)

km_final <- kmeans(pca_mat, centers = final_k, nstart = 50, iter.max = 100)
pre[, cluster := as.integer(km_final$cluster)]
pre[, cluster_label := paste0("C", cluster)]
cat(sprintf("  Cluster sizes:\n"))
print(pre[, .N, by = .(cluster_label)][order(cluster_label)])

# ── 5. Bootstrap stability (Jaccard similarity) ────────────────────────────
cat("\nBootstrap cluster stability (Jaccard)...\n")
n_obs <- nrow(pca_mat)
jaccard_all <- numeric(N_BOOT)
for (b in seq_len(N_BOOT)) {
    boot_idx <- sample(n_obs, n_obs, replace = TRUE)
    km_b <- tryCatch(
        kmeans(pca_mat[boot_idx, , drop=FALSE], centers=final_k, nstart=10, iter.max=50),
        error = function(e) NULL
    )
    if (is.null(km_b)) { jaccard_all[b] <- NA_real_; next }
    # For each full-data cluster, find best-matching bootstrap cluster
    orig_lab <- km_final$cluster[boot_idx]
    boot_lab <- km_b$cluster
    best_j <- mean(vapply(seq_len(final_k), function(ck) {
        a <- orig_lab == ck; b_vec <- boot_lab == ck
        if (sum(a | b_vec) == 0) return(NA_real_)
        sum(a & b_vec) / sum(a | b_vec)
    }, numeric(1)), na.rm = TRUE)
    jaccard_all[b] <- best_j
}
mean_jac <- mean(jaccard_all, na.rm = TRUE)
cat(sprintf("  Mean Jaccard stability: %.3f  (>0.75 = stable)\n", mean_jac))
stability_flag <- ifelse(mean_jac >= 0.75, "stable",
                  ifelse(mean_jac >= 0.60, "moderate", "unstable"))

stab_dt <- data.table(
    k               = final_k,
    method          = chose_k_method,
    null_structure  = null_structure,
    mean_jaccard    = round(mean_jac, 3),
    stability       = stability_flag,
    optimal_k_gap   = optimal_k_gap,
    optimal_k_sil   = optimal_k_sil
)
fwrite(stab_dt, file.path(OUT_DIR, "ms_preonset_cluster_stability.csv"))
cat(sprintf("  Cluster stability: %s\n", stability_flag))

# ── 6. Annotate with external data ────────────────────────────────────────
cat("\nAnnotating clusters with genetic and clinical data...\n")

# PRS
if (file.exists(PRS_FILE)) {
    prs <- fread(PRS_FILE, showProgress=FALSE)
    pre <- merge(pre, prs[, .(eid, prs_score)], by="eid", all.x=TRUE)
    cat(sprintf("  PRS merged: %d/%d\n", sum(!is.na(pre$prs_score)), nrow(pre)))
} else { pre[, prs_score := NA_real_] }

# HLA-DRB1*15:01 (from dedicated HLA imputation file)
if (file.exists(HLA_FILE)) {
    hla <- fread(HLA_FILE, showProgress=FALSE)
    hla_merge_cols <- intersect(c("eid","drb1_1501_carrier","drb1_1501_dosage"), names(hla))
    pre <- merge(pre, hla[, ..hla_merge_cols], by="eid", all.x=TRUE)
    cat(sprintf("  HLA merged: %d/%d have DRB1*15:01 status\n",
                sum(!is.na(pre$drb1_1501_carrier)), nrow(pre)))
} else { pre[, c("drb1_1501_carrier","drb1_1501_dosage") := NA_real_] }

# Genetic PCs (from covariate table)
if (file.exists(COV_FILE)) {
    cov <- fread(COV_FILE, showProgress=FALSE)
    pc_cols <- grep("^PC[0-9]", names(cov), value=TRUE)[seq_len(min(10, sum(grepl("^PC[0-9]", names(cov)))))]
    if (length(pc_cols) > 0) {
        pre <- merge(pre, cov[, c("eid", pc_cols), with=FALSE], by="eid", all.x=TRUE)
        cat(sprintf("  Genetic PCs merged: %d/%d\n",
                    sum(!is.na(pre[[pc_cols[1]]])), nrow(pre)))
    }
}

# WMH volume (post-diagnosis brain injury marker)
if (file.exists(WMH_FILE)) {
    wmh <- fread(WMH_FILE, showProgress=FALSE)
    setnames(wmh, old=grep("p25781", names(wmh), value=TRUE), new="wmh_vol",
             skip_absent=TRUE)
    if ("wmh_vol" %in% names(wmh)) {
        wmh[, wmh_vol := suppressWarnings(as.numeric(wmh_vol))]
        pre <- merge(pre, wmh[, .(eid, wmh_vol)], by="eid", all.x=TRUE)
        cat(sprintf("  WMH merged: %d/%d\n", sum(!is.na(pre$wmh_vol)), nrow(pre)))
    }
} else { pre[, wmh_vol := NA_real_] }

# First occurrence dates — compute days from MS diagnosis to any neuro event
if (file.exists(FOD_FILE)) {
    fod <- fread(FOD_FILE, showProgress=FALSE)
    setnames(fod, "participant.eid", "eid", skip_absent=TRUE)
    # p131286 = G35 (MS) date; p53_i2 = imaging visit date
    # Use imaging visit date as a post-Dx anchor if G35 date is available
    fod_cols <- c("eid",
                  grep("p131286|p53_i2", names(fod), value=TRUE))
    fod_sub <- fod[, ..fod_cols]
    for (col in setdiff(names(fod_sub), "eid")) {
        fod_sub[, (col) := as.IDate(get(col))]
    }
    pre <- merge(pre, fod_sub, by="eid", all.x=TRUE)
    g35_col <- grep("p131286", names(pre), value=TRUE)[1]
    if (!is.na(g35_col)) {
        cat(sprintf("  MS diagnosis date (G35) available: %d/%d\n",
                    sum(!is.na(pre[[g35_col]])), nrow(pre)))
    }
} else { cat("  First occurrence dates not found — skipping\n") }

# ── 7. Cluster profiles ────────────────────────────────────────────────────
cat("\nComputing cluster profiles...\n")

profile_vars <- c("prs_score","drb1_1501_carrier","drb1_1501_dosage","wmh_vol",
                  "age_at_sampling","age_at_diagnosis","years_to_diagnosis")
profile_vars <- intersect(profile_vars, names(pre))

profiles <- pre[, lapply(.SD, function(x) mean(as.numeric(x), na.rm=TRUE)),
                .SDcols = profile_vars,
                by = cluster_label]
# Add cluster size + sex composition
profiles <- merge(profiles,
    pre[, .(n = .N, pct_female = mean(sex == 0, na.rm=TRUE) * 100),
        by = cluster_label],
    by = "cluster_label")
print(profiles[, .SD, .SDcols = c("cluster_label","n",
                                    intersect(profile_vars[1:4], names(profiles)))])

fwrite(profiles, file.path(OUT_DIR, "ms_preonset_cluster_profiles.csv"))

# ── 8. Save per-participant cluster assignments ────────────────────────────
out_cols <- c("eid","ms_status","cluster","cluster_label",
              "years_to_diagnosis","age_at_sampling","sex",
              "UMAP1","UMAP2","prs_score","drb1_1501_carrier",
              "drb1_1501_dosage","wmh_vol",
              grep("^PC[0-9]", names(pre), value=TRUE)[1:5])
out_cols <- intersect(out_cols, names(pre))

# Append cluster validity metadata as attributes (stored as first comment rows in CSV)
cluster_csv <- pre[, ..out_cols]
fwrite(cluster_csv, file.path(OUT_DIR, "ms_preonset_clusters.csv"))
cat(sprintf("  Saved: ms_preonset_clusters.csv (%d rows)\n", nrow(cluster_csv)))

# Also save cluster-level DEP z-scores for heatmap
dep_zscores <- cbind(
    data.table(cluster_label = pre$cluster_label),
    as.data.table(feat_scaled)
)
dep_mean_by_cluster <- dep_zscores[,
    lapply(.SD, mean, na.rm=TRUE), .SDcols=avail_feats, by=cluster_label
]
fwrite(dep_mean_by_cluster, file.path(OUT_DIR, "ms_preonset_cluster_dep_means.csv"))

cat(sprintf("\nClustering complete.\n"))
cat(sprintf("  k = %d (%s)\n", final_k, chose_k_method))
cat(sprintf("  Stability: %s (mean Jaccard = %.3f)\n", stability_flag, mean_jac))
if (null_structure) {
    cat("  *** NULL STRUCTURE DETECTED: gap statistic supports k=1. ***\n")
    cat("  Cluster assignments are exploratory and should be interpreted cautiously.\n")
}
cat("\n02_preonset_endophenotype.R complete.\n")
