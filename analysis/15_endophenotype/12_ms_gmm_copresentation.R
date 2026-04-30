#!/usr/bin/env Rscript
# 12_ms_gmm_copresentation.R — Multimodal GMM co-presentation model
#
# Replace the hard k-means cluster assignments (phenotype-only, binary) with a
# Gaussian Mixture Model fit to the full multimodal feature matrix
# (ProtPC1-10 + ICD binary codes from {cohort_short}_endophenotype_features.csv).
#
# Output:
#   1. Soft cluster membership probabilities P(C0|pt), P(C1|pt), P(C2|pt)
#   2. Co-presenter flag: patients with bimodal/mixed membership
#   3. Comparison of hard (k-means) vs soft (GMM) assignments
#   4. Co-presenter enrichment: are co-presenters enriched for specific HLA/PRS?
#
# GMM approach:
#   - mclust package (Bayesian GMM, covariance structure selected by BIC)
#   - G = 3 components to align with k-means C0/C1/C2
#   - Also fit G = 2..5 and compare BIC (panel o)
#   - Co-presenter threshold: max(posterior) < CO_PRES_THRESH (default 0.65)
#
# Outputs (results/endophenotype/gmm_copresentation/):
#   gmm_posterior_probs.csv     — eid + P(C0) P(C1) P(C2) + hard GMM label
#   gmm_copresenters.csv        — eid + co-presenter flag + dominant cluster
#   gmm_vs_kmeans_confusion.csv — confusion matrix hard GMM vs original k-means
#   panel_w.pdf                 — BIC selection curve
#   panel_x.pdf                 — ternary / simplex membership probability plot
#   panel_y.pdf                 — co-presenter HLA/PRS enrichment

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(here)
    library(glue)
    library(mclust)
    library(patchwork)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))

cfg <- load_disease_config()

FEAT_FILE   <- here::here("results", "endophenotype",
                          glue("{cfg$cohort_short}_endophenotype_features.csv"))
CLUST_FILE  <- here::here("results", "endophenotype",
                          glue("{cfg$cohort_short}_prems_cluster_assignments.csv"))
QC_FILE     <- here::here("data", "ukb", "olink", "processed",
                          glue("{cfg$cohort_short}_olink_qc.csv"))
OUT_DIR     <- here::here("results", "endophenotype", "gmm_copresentation")
FIG_DIR     <- here::here("results", "figures", "5")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
G_TARGET      <- 3L          # align with k-means C0/C1/C2
G_RANGE       <- 2:5         # BIC sweep range
if (!(G_TARGET %in% G_RANGE)) stop("G_TARGET must be in G_RANGE")
CO_PRES_THRESH <- 0.90       # max posterior < threshold → co-presenter (GMM is highly confident; use 0.90)
CLUST_COLS["Co-presenter"]    <- "#7B2FBE"
CLUST_COLS["None/unassigned"] <- "grey65"

# ── 1. Load feature matrix ────────────────────────────────────────────────────
cat("Loading multimodal feature matrix...\n")
feat  <- fread(FEAT_FILE)
clust <- fread(CLUST_FILE)
qc    <- fread(QC_FILE, showProgress = FALSE)

# Feature columns: ProtPC1-10 + ICD binary codes (everything except eid)
feat_cols <- setdiff(names(feat), "eid")
prot_cols <- grep("^ProtPC", feat_cols, value = TRUE)
icd_cols  <- setdiff(feat_cols, prot_cols)
cat(sprintf("  Features: %d ProtPCs + %d ICD blocks = %d total\n",
            length(prot_cols), length(icd_cols), length(feat_cols)))

# Scale ProtPCs (already normalised but re-scale for joint modelling with binary ICD)
feat_mat       <- as.matrix(feat[, ..feat_cols])
feat_mat_scale <- feat_mat
feat_mat_scale[, prot_cols] <- scale(feat_mat[, prot_cols])
rownames(feat_mat_scale)    <- feat$eid

cat(sprintf("  Feature matrix: %d patients × %d features\n",
            nrow(feat_mat_scale), ncol(feat_mat_scale)))

# ── 2. BIC sweep: G = 2..5 ───────────────────────────────────────────────────
cat("Fitting GMM across G = 2..5 (BIC sweep)...\n")
bic_results <- list()
mc_objects  <- list()

for (g in G_RANGE) {
    mc <- tryCatch(
        Mclust(feat_mat_scale, G = g, verbose = FALSE),
        error = function(e) { cat(sprintf("  G=%d failed: %s\n", g, e$message)); NULL }
    )
    if (!is.null(mc)) {
        mc_objects[[as.character(g)]] <- mc
        bic_val <- mc$BIC[1L, mc$modelName]
        bic_results[[as.character(g)]] <- data.table(G = g, BIC = bic_val)
        cat(sprintf("  G=%d  BIC=%.1f  model=%s\n", g, bic_val, mc$modelName))
    }
}

bic_dt <- rbindlist(bic_results)
fwrite(bic_dt, file.path(OUT_DIR, "gmm_bic_sweep.csv"))

# Panel o: BIC curve
p_bic <- ggplot(bic_dt, aes(x = G, y = BIC)) +
    geom_line(linewidth = 0.7, colour = "#CC0066") +
    geom_point(size = 2.5, colour = "#CC0066") +
    geom_vline(xintercept = G_TARGET, linetype = "dashed",
               linewidth = 0.35, colour = "grey55") +
    annotate("text", x = G_TARGET + 0.1, y = max(bic_dt$BIC, na.rm=TRUE) * 0.99,
             label = sprintf("G=%d (target)", G_TARGET),
             hjust = 0, size = 2.5, colour = "grey45") +
    scale_x_continuous(breaks = G_RANGE) +
    labs(x = "Number of components (G)",
         y = "BIC",
         title    = "o  GMM component selection",
         subtitle = "mclust | ProtPC1-10 + ICD binary features | higher BIC = better") +
    theme_ukb(base_size = 9)

ggsave(file.path(FIG_DIR, "panel_w.pdf"), p_bic,
       width = 3.5, height = 3.0, device = cairo_pdf)
cat("  panel_w.pdf\n")

# ── 3. Extract target GMM (G = 3) ────────────────────────────────────────────
cat(sprintf("\nExtracting GMM G=%d from sweep...\n", G_TARGET))
mc3 <- mc_objects[[as.character(G_TARGET)]]
if (is.null(mc3)) stop("GMM fit failed for G=", G_TARGET)
cat(sprintf("  Model: %s  BIC: %.1f\n", mc3$modelName, mc3$BIC[1L, mc3$modelName]))

# Posterior probabilities (n × G matrix)
post <- mc3$z  # rows = patients, columns = GMM components
colnames(post) <- paste0("P_comp", seq_len(G_TARGET))

# Hard GMM assignment (MAP)
hard_gmm <- mc3$classification

# ── 4. Align GMM components to k-means clusters ──────────────────────────────
# Match GMM component → k-means cluster by majority vote
cat("Aligning GMM components to k-means clusters...\n")

feat_with_clust <- merge(feat, clust[, .(eid, cluster)], by = "eid", all.x = TRUE)
feat_with_clust[, kmeans_label := ifelse(is.na(cluster), "None", paste0("C", cluster))]

# Confusion matrix: GMM component (rows) vs k-means label (cols)
gmm_dt <- data.table(
    eid      = feat$eid,
    gmm_comp = hard_gmm
)
gmm_dt <- merge(gmm_dt, feat_with_clust[, .(eid, kmeans_label)], by = "eid")

conf <- gmm_dt[kmeans_label != "None",
               .N, by = .(gmm_comp, kmeans_label)]
conf_wide <- dcast(conf, gmm_comp ~ kmeans_label, value.var = "N", fill = 0L)
cat("  GMM component vs k-means confusion (non-None patients):\n")
print(conf_wide)
fwrite(conf_wide, file.path(OUT_DIR, "gmm_vs_kmeans_confusion.csv"))

# Map each GMM component to the plurality k-means cluster
kmeans_clusters <- c("C0","C1","C2")
comp_to_cluster <- vapply(seq_len(G_TARGET), function(g) {
    row <- conf[gmm_comp == g]
    if (nrow(row) == 0L) return(NA_character_)
    row[kmeans_label %in% kmeans_clusters][which.max(N), kmeans_label]
}, character(1))
names(comp_to_cluster) <- paste0("comp", seq_len(G_TARGET))
cat("  Component → cluster mapping:", paste(names(comp_to_cluster), comp_to_cluster, sep="→", collapse=", "), "\n")

# Rename posterior columns to matched cluster labels
post_dt <- as.data.table(post)
colnames(post_dt) <- paste0("P_", comp_to_cluster)
post_dt[, eid := feat$eid]

# Handle duplicate column names if two components map to same cluster (shouldn't happen)
if (anyDuplicated(names(post_dt)) > 0) {
    colnames(post_dt) <- make.unique(colnames(post_dt))
    cat("  WARNING: duplicate component-to-cluster mappings — check GMM alignment\n")
}

# ── 5. Co-presenter classification ───────────────────────────────────────────
cat("\nClassifying co-presenters...\n")

p_cols     <- grep("^P_C", names(post_dt), value = TRUE)
max_post   <- apply(as.matrix(post_dt[, ..p_cols]), 1, max)
dom_comp   <- apply(as.matrix(post_dt[, ..p_cols]), 1, which.max)
dom_label  <- gsub("P_", "", p_cols[dom_comp])

post_dt[, max_posterior   := max_post]
post_dt[, dominant_cluster := dom_label]
post_dt[, co_presenter    := max_posterior < CO_PRES_THRESH]
post_dt[, hard_gmm_label  := dom_label]

n_copres <- sum(post_dt$co_presenter)
cat(sprintf("  Co-presenters (max posterior < %.2f): %d / %d (%.1f%%)\n",
            CO_PRES_THRESH, n_copres, nrow(post_dt),
            100 * n_copres / nrow(post_dt)))

fwrite(post_dt[order(eid)], file.path(OUT_DIR, "gmm_posterior_probs.csv"))

copres_dt <- post_dt[, .(eid, dominant_cluster, co_presenter, max_posterior)]
copres_dt <- merge(copres_dt, clust[, .(eid, cluster)], by = "eid", all.x = TRUE)
copres_dt[, kmeans_label := ifelse(is.na(cluster), "None", paste0("C", cluster))]
fwrite(copres_dt, file.path(OUT_DIR, "gmm_copresenters.csv"))

# ── 6. Panel p: ternary / 2D simplex probability plot ────────────────────────
cat("\nBuilding membership probability plot (panel p)...\n")

if (length(p_cols) >= 2L) {
    # Project 3D simplex onto 2D via first two principal membership axes
    pca_post    <- prcomp(as.matrix(post_dt[, ..p_cols]), scale. = FALSE)
    plot_dt     <- data.table(
        PC1       = pca_post$x[, 1],
        PC2       = pca_post$x[, 2],
        co_pres   = post_dt$co_presenter,
        dom_clust = post_dt$dominant_cluster,
        max_post  = post_dt$max_posterior
    )
    plot_dt[, point_label := ifelse(co_pres, "Co-presenter", dom_clust)]
    plot_dt[, point_label := factor(point_label,
                                     levels = c("C0","C1","C2","Co-presenter"))]
    point_cols <- CLUST_COLS[c("C0", "C1", "C2", "Co-presenter")]

    p_p <- ggplot(plot_dt, aes(x = PC1, y = PC2,
                                colour = point_label,
                                alpha  = max_post,
                                size   = ifelse(co_pres, 1.8, 1.0))) +
        geom_point(stroke = 0) +
        scale_colour_manual(values = point_cols, name = "GMM label") +
        scale_alpha_continuous(range = c(0.4, 0.95), guide = "none") +
        scale_size_identity() +
        labs(x        = "PC1 (membership)",
             y        = "PC2 (membership)",
             title    = "p  GMM soft membership probability space",
             subtitle = sprintf("mclust G=3 | G=%d used | co-presenters: max posterior < %.2f | n=%d",
                                G_TARGET, CO_PRES_THRESH, n_copres)) +
        theme_ukb(base_size = 9) +
        theme(legend.position = "right")

    ggsave(file.path(FIG_DIR, "panel_x.pdf"), p_p,
           width = 4.5, height = 3.8, device = cairo_pdf)
    cat("  panel_x.pdf\n")
}

# ── 7. Panel q: Co-presenter enrichment for HLA / PRS ────────────────────────
cat("\nBuilding co-presenter HLA/PRS enrichment (panel q)...\n")

if ("hla_carrier" %in% names(clust) && "prs_score" %in% names(clust)) {
    enr_dt <- merge(copres_dt[, .(eid, co_presenter, dominant_cluster)],
                    clust[, .(eid, hla_carrier, prs_score)], by = "eid", all.x = TRUE)
    enr_dt <- enr_dt[!is.na(hla_carrier) | !is.na(prs_score)]

    if (nrow(enr_dt) > 0L) {
        # HLA carrier rate by co-presenter status
        hla_dt <- enr_dt[!is.na(hla_carrier), .(
            hla_rate = mean(hla_carrier == "carrier", na.rm = TRUE),
            n        = .N
        ), by = co_presenter]
        hla_dt[, group := ifelse(co_presenter, "Co-presenter", "Single-cluster")]

        # PRS score by co-presenter status
        prs_dt <- enr_dt[!is.na(prs_score), .(
            prs_mean = mean(prs_score, na.rm = TRUE),
            prs_se   = sd(prs_score, na.rm = TRUE) / sqrt(.N),
            n        = .N
        ), by = co_presenter]
        prs_dt[, group := ifelse(co_presenter, "Co-presenter", "Single-cluster")]

        p_hla <- ggplot(hla_dt, aes(x = group, y = hla_rate, fill = group)) +
            geom_col(width = 0.55, alpha = 0.9) +
            geom_text(aes(label = sprintf("n=%d", n)),
                      vjust = -0.4, size = 2.5, colour = "grey30") +
            scale_fill_manual(values = c("Co-presenter" = "#7B2FBE",
                                         "Single-cluster" = "grey65"),
                              guide = "none") +
            scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
            labs(x = NULL, y = glue("HLA-{cfg$hla_allele} carrier rate"),
                 title = "q  Co-presenter HLA + PRS",
                 subtitle = glue("HLA carrier rate (left) | mean {cfg$prs_label} (right)")) +
            theme_ukb(base_size = 9) +
            theme(panel.grid.major.x = element_blank())

        p_prs <- ggplot(prs_dt, aes(x = group, y = prs_mean, fill = group)) +
            geom_col(width = 0.55, alpha = 0.9) +
            geom_errorbar(aes(ymin = prs_mean - prs_se, ymax = prs_mean + prs_se),
                          width = 0.2, linewidth = 0.5, colour = "grey30") +
            scale_fill_manual(values = c("Co-presenter" = "#7B2FBE",
                                         "Single-cluster" = "grey65"),
                              guide = "none") +
            labs(x = NULL, y = glue("Mean {cfg$prs_label} (±SE)")) +
            theme_ukb(base_size = 9) +
            theme(panel.grid.major.x = element_blank())

        p_q <- p_hla + p_prs + patchwork::plot_layout(nrow = 1)
        ggsave(file.path(FIG_DIR, "panel_y.pdf"), p_q,
               width = 5.5, height = 3.5, device = cairo_pdf)
        cat("  panel_y.pdf\n")
    }
} else {
    cat("  hla_carrier or prs_score missing from cluster assignments — skipping panel q\n")
}

cat("\n12_ms_gmm_copresentation.R complete.\n")
cat("Figure 5 multimodal extension complete.\n")
cat("Outputs:\n")
cat("  results/endophenotype/gmm_copresentation/\n")
cat("  results/figures/5/panel_w.pdf  (BIC)\n")
cat("  results/figures/5/panel_x.pdf  (membership space)\n")
cat("  results/figures/5/panel_y.pdf  (HLA/PRS enrichment)\n")
