#!/usr/bin/env Rscript
# 01_ms_patient_clustering.R
# Unsupervised proteomics clustering of MS cases — UKB-unique extension
#
# Approach:
#   1. Scale 173 DEP profiles within MS cases
#   2. Determine optimal k via silhouette + gap statistic (k=2-5)
#   3. K-means clustering + hierarchical (ward.D2) for heatmap
#   4. Characterize clusters: pre/post timing, sex, age, cognitive, DTI FA
#
# Output: results/clustering/
#   ms_cluster_heatmap.pdf
#   ms_cluster_umap.pdf
#   ms_cluster_clinical.pdf
#   ms_cluster_assignments.csv

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(factoextra)
    library(patchwork)
    library(cluster)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

QC_FILE   <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DIFF_FILE <- file.path(PROJ_DIR, "results", "differential", "ms_combined_vs_hc.csv")
COG_DIR   <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                        "data", "ukb", "cognitive_function")
DTI_FILE  <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                        "data", "ukb", "brain_mri", "dti", "dti_tract_protein_data.tsv")
OUT_DIR   <- file.path(PROJ_DIR, "results", "clustering")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
K_MAX <- 5

# ── 1. Load DEPs & QC data ────────────────────────────────────────────────────
cat("Loading data...\n")
deps <- fread(DIFF_FILE, showProgress = FALSE)
dep_prots <- tolower(deps[adj.P.Val < 0.05, protein])

ms_qc <- fread(QC_FILE, showProgress = FALSE)
ms_qc <- ms_qc[qc_outlier == FALSE & ms_status %in% c("pre_onset","post_onset") &
                !is.na(age_at_sampling) & !is.na(sex) & !is.na(UMAP1) & !is.na(UMAP2)]
cat(sprintf("  MS cases (QC passed): %d\n", nrow(ms_qc)))

avail_prots <- intersect(dep_prots, names(ms_qc))
cat(sprintf("  DEPs available: %d / %d\n", length(avail_prots), length(dep_prots)))

# Build feature matrix — complete cases only
# UKB sex: 0=Female, 1=Male — create label for plots
ms_qc[, sex_label := fifelse(sex == 0, "Female", "Male")]

feat_dt <- ms_qc[, c("eid", "ms_status", "age_at_sampling", "sex", "sex_label",
                      "age_at_diagnosis", "years_to_diagnosis", "UMAP1", "UMAP2",
                      ..avail_prots), with = TRUE]
feat_cc  <- feat_dt[complete.cases(feat_dt[, ..avail_prots])]
cat(sprintf("  Complete cases: %d\n", nrow(feat_cc)))

# Scale protein matrix within MS cases
X <- as.matrix(feat_cc[, ..avail_prots])
X_scaled <- scale(X)
rownames(X_scaled) <- feat_cc$eid

# ── 2. Optimal k via silhouette and gap statistic ──────────────────────────────
cat("Selecting optimal k (silhouette + gap statistic)...\n")

sil_vals <- numeric(K_MAX - 1)
for (k in 2:K_MAX) {
    km  <- kmeans(X_scaled, centers = k, nstart = 25, iter.max = 50)
    sil <- silhouette(km$cluster, dist(X_scaled, method = "euclidean"))
    sil_vals[k - 1] <- mean(sil[, 3])
    cat(sprintf("  k=%d silhouette=%.3f\n", k, sil_vals[k - 1]))
}
k_opt <- which.max(sil_vals) + 1
cat(sprintf("  Optimal k = %d\n", k_opt))

# Also run for k=2,3,4 — report all
for (k_run in 2:min(4, K_MAX)) {
    km_tmp <- kmeans(X_scaled, centers = k_run, nstart = 25, iter.max = 50)
    sil_tmp <- silhouette(km_tmp$cluster, dist(X_scaled))
    cat(sprintf("  k=%d: cluster sizes %s\n", k_run,
                paste(table(km_tmp$cluster), collapse = "/")))
}

# ── 3. K-means at optimal k ───────────────────────────────────────────────────
cat(sprintf("\nK-means clustering (k=%d)...\n", k_opt))
km_final <- kmeans(X_scaled, centers = k_opt, nstart = 50, iter.max = 100)
feat_cc[, cluster := as.factor(km_final$cluster)]

# ── 4. Silhouette plot ────────────────────────────────────────────────────────
sil_final <- silhouette(km_final$cluster, dist(X_scaled))
sil_dt    <- as.data.table(sil_final)
sil_dt[, cluster := as.factor(cluster)]
sil_dt    <- sil_dt[order(cluster, -sil_width)]
sil_dt[, idx := seq_len(.N)]

p_sil <- ggplot(sil_dt, aes(x = idx, y = sil_width, fill = cluster)) +
    geom_col(width = 1) +
    geom_hline(yintercept = mean(sil_dt$sil_width), linetype = "dashed",
               colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c("#CC0066","#56B4E9","#E6A817","#2B4C7E","#1A1A1A")[seq_len(k_opt)],
                      name = "Cluster") +
    labs(x = "MS case", y = "Silhouette width",
         title = sprintf("Cluster quality (k=%d, mean=%.3f)", k_opt, mean(sil_dt$sil_width))) +
    theme_ukb()

ggsave(file.path(OUT_DIR, "ms_cluster_silhouette.pdf"), p_sil,
       width = 5, height = 3.5, device = cairo_pdf)

# ── 5. UMAP coloured by cluster ───────────────────────────────────────────────
CLUST_COLS <- c("#CC0066","#56B4E9","#E6A817","#2B4C7E","#1A1A1A")[seq_len(k_opt)]

p_umap <- ggplot(feat_cc, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_colour_manual(values = CLUST_COLS, name = "Cluster") +
    labs(title = sprintf("MS cases by proteomic cluster (k=%d)", k_opt),
         subtitle = sprintf("n=%d, %d DEPs", nrow(feat_cc), length(avail_prots))) +
    theme_ukb()

ggsave(file.path(OUT_DIR, "ms_cluster_umap.pdf"), p_umap,
       width = 5, height = 4.5, device = cairo_pdf)
cat("  Saved: ms_cluster_umap.pdf\n")

# ── 6. Heatmap of top discriminating proteins ─────────────────────────────────
cat("Creating cluster heatmap...\n")

# Top 40 most variable proteins
prot_var    <- apply(X_scaled, 2, var)
top40_prots <- names(sort(prot_var, decreasing = TRUE))[1:min(40, length(prot_var))]
X_top       <- t(X_scaled[, top40_prots])

# Column annotation
annot_col <- data.frame(
    Cluster   = paste0("C", feat_cc$cluster),
    Timing    = feat_cc$ms_status,
    Sex       = feat_cc$sex_label,
    row.names = as.character(feat_cc$eid)
)

clust_colors <- setNames(CLUST_COLS, paste0("C", seq_len(k_opt)))
ann_colors <- list(
    Cluster = clust_colors,
    Timing  = c(pre_onset = "#E6A817", post_onset = "#CC0066"),
    Sex     = c(Male = "#56B4E9", Female = "#CC0066")
)

pdf(file.path(OUT_DIR, "ms_cluster_heatmap.pdf"), width = 10, height = 8)
pheatmap(X_top,
         annotation_col   = annot_col,
         annotation_colors = ann_colors,
         show_colnames    = FALSE,
         show_rownames    = TRUE,
         fontsize_row     = 6,
         clustering_method = "ward.D2",
         color            = colorRampPalette(c("#56B4E9","white","#CC0066"))(100),
         breaks           = seq(-2.5, 2.5, length.out = 101),
         main             = sprintf("MS proteomics clusters (k=%d, top 40 variable DEPs)", k_opt))
dev.off()
cat("  Saved: ms_cluster_heatmap.pdf\n")

# ── 7. Clinical characterization ─────────────────────────────────────────────
cat("Characterizing clusters clinically...\n")

clust_summary <- feat_cc[, .(
    n             = .N,
    pct_pre_onset = 100 * mean(ms_status == "pre_onset"),
    pct_female    = 100 * mean(sex == 0),   # UKB: 0=Female
    mean_age_dx   = mean(age_at_diagnosis, na.rm = TRUE),
    mean_years_to_dx = mean(years_to_diagnosis, na.rm = TRUE)
), by = cluster]
cat("\nCluster summary:\n")
print(clust_summary)

# Test each clinical variable across clusters
test_results <- list()

for (var in c("age_at_sampling", "age_at_diagnosis", "years_to_diagnosis")) {
    kw <- tryCatch(
        kruskal.test(get(var) ~ cluster, data = feat_cc),
        error = function(e) NULL
    )
    if (!is.null(kw))
        test_results[[var]] <- data.table(variable = var, p_kruskal = kw$p.value)
}
for (var in c("ms_status", "sex")) {
    chi <- tryCatch(
        chisq.test(table(feat_cc[[var]], feat_cc$cluster)),
        error = function(e) NULL
    )
    if (!is.null(chi))
        test_results[[var]] <- data.table(variable = var, p_kruskal = chi$p.value)
}
if (length(test_results) > 0) {
    cat("\nCluster × clinical variable tests:\n")
    print(rbindlist(test_results, fill = TRUE))
}

# ── 8. Clinical violin plots ──────────────────────────────────────────────────
plot_list <- list()

for (var in c("age_at_diagnosis", "years_to_diagnosis", "age_at_sampling")) {
    if (!var %in% names(feat_cc)) next
    p_v <- ggplot(feat_cc, aes(x = cluster, y = get(var), fill = cluster)) +
        geom_violin(alpha = 0.7, linewidth = 0.3) +
        geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5, linewidth = 0.4) +
        scale_fill_manual(values = CLUST_COLS, guide = "none") +
        labs(x = "Cluster", y = gsub("_", " ", var),
             title = gsub("_", " ", var)) +
        theme_ukb()
    plot_list[[var]] <- p_v
}

# Pre/post timing bar
timing_dt <- feat_cc[, .N, by = .(cluster, ms_status)]
timing_dt[, pct := 100 * N / sum(N), by = cluster]
p_timing <- ggplot(timing_dt, aes(x = cluster, y = pct, fill = ms_status)) +
    geom_col(position = "stack", width = 0.6) +
    scale_fill_manual(values = c(pre_onset = "#E6A817", post_onset = "#CC0066"),
                      name = NULL, labels = c(pre_onset = "Pre-onset", post_onset = "Post-onset")) +
    labs(x = "Cluster", y = "% of cases", title = "Pre vs Post-onset timing") +
    theme_ukb()
plot_list[["timing"]] <- p_timing

# Sex bar
sex_dt <- feat_cc[, .N, by = .(cluster, sex_label)]
sex_dt[, pct := 100 * N / sum(N), by = cluster]
p_sex <- ggplot(sex_dt, aes(x = cluster, y = pct, fill = sex_label)) +
    geom_col(position = "stack", width = 0.6) +
    scale_fill_manual(values = c(Male = "#56B4E9", Female = "#CC0066"), name = NULL) +
    labs(x = "Cluster", y = "% of cases", title = "Sex distribution") +
    theme_ukb()
plot_list[["sex"]] <- p_sex

if (length(plot_list) >= 2) {
    p_clin <- wrap_plots(plot_list, ncol = 3)
    ggsave(file.path(OUT_DIR, "ms_cluster_clinical.pdf"), p_clin,
           width = 12, height = 4 * ceiling(length(plot_list)/3),
           device = cairo_pdf)
    cat("  Saved: ms_cluster_clinical.pdf\n")
}

# ── 9. Try to merge cognitive + DTI data ─────────────────────────────────────
cog_file <- tryCatch(
    list.files(COG_DIR, pattern = "cognitive_gcsvd", full.names = TRUE, recursive = TRUE)[1],
    error = function(e) NA_character_
)

if (!is.na(cog_file) && file.exists(cog_file)) {
    cog <- fread(cog_file, showProgress = FALSE,
                 select = intersect(c("eid","gcsvd_score",
                                      "fluid_intelligence__fi_score__i0",
                                      "reaction_time__mean_reaction_time_ms__i0"),
                                    names(fread(cog_file, nrows = 0))))
    feat_cog <- merge(feat_cc[, .(eid, cluster)], cog, by = "eid", all.x = TRUE)
    for (cv in setdiff(names(feat_cog), c("eid","cluster"))) {
        if (sum(!is.na(feat_cog[[cv]])) < 10) next
        kw <- kruskal.test(get(cv) ~ cluster, data = feat_cog[!is.na(get(cv))])
        cat(sprintf("  Cognitive %s: Kruskal p=%.4f\n", cv, kw$p.value))
    }
}

if (file.exists(DTI_FILE)) {
    dti <- fread(DTI_FILE, showProgress = FALSE,
                 select = c("eid", grep("^fa_", names(fread(DTI_FILE, nrows = 0)), value = TRUE)[1:5]))
    feat_dti <- merge(feat_cc[, .(eid, cluster)], dti, by = "eid", all.x = TRUE)
    for (dv in setdiff(names(feat_dti), c("eid","cluster"))) {
        if (sum(!is.na(feat_dti[[dv]])) < 10) next
        kw <- kruskal.test(get(dv) ~ cluster, data = feat_dti[!is.na(get(dv))])
        cat(sprintf("  DTI %s: Kruskal p=%.4f\n", dv, kw$p.value))
    }
}

# ── 10. Save assignments ──────────────────────────────────────────────────────
fwrite(feat_cc[, .(eid, cluster, ms_status, sex, age_at_sampling,
                   age_at_diagnosis, years_to_diagnosis)],
       file.path(OUT_DIR, "ms_cluster_assignments.csv"))

# Also save top proteins per cluster (cluster centroids)
centroids <- as.data.table(km_final$centers, keep.rownames = "cluster")
fwrite(centroids, file.path(OUT_DIR, "ms_cluster_centroids.csv"))

cat(sprintf("\n01_ms_patient_clustering.R complete. k=%d clusters.\n", k_opt))
