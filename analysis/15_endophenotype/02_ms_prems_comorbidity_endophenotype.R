#!/usr/bin/env Rscript
# 02_ms_prems_comorbidity_endophenotype.R — Figure 5, pre-MS clustering
#
# FEATURES
#   Enriched blocks (OR>1, FDR<0.05, ≥5 pre-MS patients): ICD 2-char block
#   composites used as primary clustering features. Block = 1 iff patient has
#   ≥1 OR>1 FDR<0.05 code in that block >1yr before MS diagnosis.
#   Depleted blocks (OR<1, FDR<0.05) retained as post-hoc cluster annotations.
#   Full 3-char code matrix also saved for specificity analysis.
#
# DIMENSIONALITY REDUCTION
#   Precomputed Jaccard distance on patients with ≥1 enriched block.
#   Parameter sweep: n_neighbors × min_dist × k, selected by mean silhouette.
#
# OUTPUTS (results/endophenotype/)
#   ms_prems_block_features.csv       block × patient matrix + cluster labels
#   ms_prems_cluster_assignments.csv  cluster + outcomes per patient
#   ms_prems_cluster_summary.csv      per-cluster summary stats
#   ms_prems_umap_config.csv          winning hyperparameters
#   ms_prems_block_overrep.csv        Fisher enrichment per block × cluster
#   ms_prems_code_specificity.csv     3-char code Fisher per cluster (sensitivity)
#   fig5_row1_prems_umap.pdf
#   ms_prems_overrep_heatmap.pdf
#   ms_prems_cluster_outcomes.pdf

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(cluster)
    library(uwot)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

ICD_DIR     <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                         "data", "ukb", "diagnoses", "icd_codes")
QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
HLA_FILE    <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE    <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
PHEWAS_FILE <- file.path(PROJ_DIR, "results", "comorbidity", "ms_prems_phewas.csv")
OUT_DIR     <- file.path(PROJ_DIR, "results", "endophenotype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
MS_EXCL   <- c("G35", "G36", "G37")
MIN_PRE_N <- 5L   # minimum patients with an enriched block pre-MS

CLUST_COLS <- c("0"="grey65","1"="#E6A817","2"="#CC0066",
                "3"="#2B4C7E","4"="#56B4E9","5"="#1A1A1A")

# ── 1. MS cohort ──────────────────────────────────────────────────────────────
cat("Loading MS cohort...\n")
qc <- fread(QC_FILE, showProgress = FALSE)
ms <- qc[qc_outlier == FALSE &
         ms_status %in% c("pre_onset", "post_onset") &
         !is.na(age_at_diagnosis)]
cat(sprintf("  n=%d MS cases (pre-onset: %d, post-onset: %d)\n",
            nrow(ms), sum(ms$ms_status == "pre_onset"), sum(ms$ms_status == "post_onset")))
ms_dx <- ms[, .(eid, age_at_diagnosis)]

# ── 2. Split PheWAS results into enriched / depleted ──────────────────────────
cat("\nLoading pre-MS PheWAS results...\n")
phewas <- fread(PHEWAS_FILE)
phewas_sig <- phewas[fdr < 0.05 & !icd10 %in% MS_EXCL]

enriched_codes <- phewas_sig[OR > 1, icd10]
depleted_codes <- phewas_sig[OR < 1, icd10]
cat(sprintf("  FDR<0.05: %d enriched (OR>1), %d depleted (OR<1)\n",
            length(enriched_codes), length(depleted_codes)))

# ── 3. Load pre-MS events ─────────────────────────────────────────────────────
load_code_prems <- function(code) {
    ch  <- substr(code, 1, 1)
    blk <- substr(code, 1, 2)
    fp  <- file.path(ICD_DIR, ch, blk, code, paste0(code, ".csv"))
    if (!file.exists(fp)) return(NULL)
    dt  <- fread(fp, showProgress = FALSE)
    names(dt) <- c("eid", "ever_present", "age_at_event")
    dt  <- dt[ever_present == 1L & eid %in% ms$eid,
              .(eid           = as.integer(eid),
                age_at_event  = as.numeric(age_at_event),
                icd10         = code,
                block         = blk)]
    if (nrow(dt) == 0L) return(NULL)
    dt
}

cat("Loading enriched code events...\n")
raw_enriched <- rbindlist(Filter(Negate(is.null),
                                  lapply(enriched_codes, load_code_prems)), fill = TRUE)
prems_enriched <- merge(raw_enriched, ms_dx, by = "eid")
prems_enriched <- prems_enriched[age_at_event < age_at_diagnosis - 1]

cat("Loading depleted code events...\n")
raw_depleted <- rbindlist(Filter(Negate(is.null),
                                  lapply(depleted_codes, load_code_prems)), fill = TRUE)
prems_depleted <- merge(raw_depleted, ms_dx, by = "eid")
prems_depleted <- prems_depleted[age_at_event < age_at_diagnosis - 1]

cat(sprintf("  Enriched pre-MS: %d events | %d patients | %d codes\n",
            nrow(prems_enriched), uniqueN(prems_enriched$eid),
            uniqueN(prems_enriched$icd10)))
cat(sprintf("  Depleted pre-MS: %d events | %d patients | %d codes\n",
            nrow(prems_depleted), uniqueN(prems_depleted$eid),
            uniqueN(prems_depleted$icd10)))

# ── 4. Build enriched block composites (primary clustering features) ───────────
cat("\nBuilding enriched block composites...\n")
enr_blocks_all <- unique(substr(enriched_codes, 1, 2))

block_dt <- data.table(eid = ms$eid)
for (blk in enr_blocks_all) {
    codes_in_blk <- enriched_codes[substr(enriched_codes, 1, 2) == blk]
    eids_with    <- prems_enriched[icd10 %in% codes_in_blk, unique(eid)]
    block_dt[, (blk) := as.integer(eid %in% eids_with)]
}

# Filter to blocks with ≥ MIN_PRE_N patients
block_n <- sapply(enr_blocks_all, function(b) sum(block_dt[[b]]))
enr_blocks <- enr_blocks_all[block_n >= MIN_PRE_N]
cat(sprintf("  %d enriched blocks pass min_N=%d:\n", length(enr_blocks), MIN_PRE_N))
blk_tbl <- data.table(block = enr_blocks, n = block_n[enr_blocks],
                       codes = sapply(enr_blocks, function(b)
                           paste(enriched_codes[substr(enriched_codes,1,2)==b],
                                 collapse=",")))
print(blk_tbl[order(-n)])

# ── 5. Build depleted block annotations ──────────────────────────────────────
dep_blocks_all <- unique(substr(depleted_codes, 1, 2))
for (blk in dep_blocks_all) {
    codes_in_blk <- depleted_codes[substr(depleted_codes, 1, 2) == blk]
    eids_with    <- prems_depleted[icd10 %in% codes_in_blk, unique(eid)]
    block_dt[, paste0("dep_", blk) := as.integer(eid %in% eids_with)]
}
dep_cols <- paste0("dep_", dep_blocks_all)

# ── 6. Code-level matrix for sensitivity analysis ────────────────────────────
cat("\nBuilding code-level matrix (sensitivity layer)...\n")
all_sig_codes <- c(enriched_codes, depleted_codes)
loaded_codes  <- unique(c(prems_enriched$icd10, prems_depleted$icd10))
code_dt <- data.table(eid = ms$eid)
for (cc in loaded_codes) {
    src <- if (cc %in% enriched_codes) prems_enriched else prems_depleted
    code_dt[, (cc) := as.integer(eid %in% src[icd10 == cc, unique(eid)])]
}

# ── 7. UMAP parameter sweep on enriched blocks ───────────────────────────────
cat("\n=== UMAP parameter sweep (enriched blocks, Jaccard distance) ===\n")
bmat_full  <- as.matrix(block_dt[, ..enr_blocks])
rownames(bmat_full) <- block_dt$eid
has_feat   <- rowSums(bmat_full) > 0
bmat       <- bmat_full[has_feat, , drop = FALSE]
eids_sub   <- block_dt$eid[has_feat]
cat(sprintf("  Input: %d patients with ≥1 enriched block | %d blocks\n",
            nrow(bmat), length(enr_blocks)))

int_mat  <- tcrossprod(bmat)
rs       <- rowSums(bmat)
union_m  <- outer(rs, rs, "+") - int_mat
jac_mat  <- 1 - int_mat / pmax(union_m, 1e-10)
diag(jac_mat) <- 0
jac_dist <- as.dist(jac_mat)
cat(sprintf("  Jaccard dist: mean=%.3f\n",
            mean(jac_mat[upper.tri(jac_mat)])))

nn_grid <- c(10L, 15L, 20L, 30L)
md_grid <- c(0.05, 0.10, 0.20)
k_grid  <- 2:5

best_sil    <- -Inf
best_config <- NULL
best_umap   <- NULL
best_km     <- NULL
best_k      <- 2L

max_nn <- nrow(bmat) - 1L
for (nn in pmin(nn_grid, max_nn)) {
    for (md in md_grid) {
        set.seed(42)
        um <- tryCatch(
            umap(jac_dist, n_neighbors = nn, min_dist = md,
                 n_threads = 2, verbose = FALSE),
            error = function(e) NULL
        )
        if (is.null(um)) next
        for (k in k_grid) {
            if (k >= nrow(bmat)) next
            km  <- kmeans(um, centers = k, nstart = 50, iter.max = 200)
            sil <- mean(silhouette(km$cluster, dist(um))[, 3])
            cat(sprintf("  nn=%2d md=%.2f k=%d sil=%.3f\n", nn, md, k, sil))
            if (sil > best_sil) {
                best_sil    <- sil
                best_config <- list(nn = nn, md = md, k = k)
                best_umap   <- um
                best_km     <- km
                best_k      <- k
            }
        }
    }
}

if (is.null(best_umap)) stop("All UMAP configurations failed.")
cat(sprintf("\n  WINNER: nn=%d  md=%.2f  k=%d  sil=%.3f\n",
            best_config$nn, best_config$md, best_config$k, best_sil))

# ── 8. Assemble annotation table ──────────────────────────────────────────────
cat("\nAssembling annotation table...\n")
umap_sub <- data.table(eid     = eids_sub,
                        U1      = best_umap[, 1],
                        U2      = best_umap[, 2],
                        cluster = best_km$cluster - 1L)
umap_all <- merge(data.table(eid = ms$eid), umap_sub, by = "eid", all.x = TRUE)

hla <- tryCatch(fread(HLA_FILE, showProgress = FALSE), error = function(e) NULL)
prs <- tryCatch(fread(PRS_FILE, showProgress = FALSE), error = function(e) NULL)

ann <- merge(block_dt[, c("eid", enr_blocks, dep_cols), with = FALSE],
             umap_all, by = "eid")
ann <- merge(ann, ms[, .(eid, ms_status, age_at_diagnosis, sex, years_to_diagnosis)],
             by = "eid")
if (!is.null(hla)) ann <- merge(ann, hla[, .(eid, hla_carrier = drb1_1501_carrier)],
                                 by = "eid", all.x = TRUE)
if (!is.null(prs)) ann <- merge(ann, prs[, .(eid, prs_score)], by = "eid", all.x = TRUE)

ann[, n_enr_blocks := rowSums(.SD), .SDcols = enr_blocks]
ann[, burden := fcase(n_enr_blocks == 0, "None",
                      n_enr_blocks == 1, "1 block",
                      n_enr_blocks >= 2, "\u22652 blocks")]
ann[, burden    := factor(burden, levels = c("None", "1 block", "\u22652 blocks"))]
ann[, cluster_f := factor(ifelse(is.na(cluster), "None", paste0("C", cluster)))]
ann[, sex_label := ifelse(sex == 0, "Female", "Male")]
ann[, age_dx_q  := cut(age_at_diagnosis, breaks = c(-Inf, 40, 50, Inf),
                        labels = c("Early (<40)", "Mid (40\u201350)", "Late (>50)"))]
if ("hla_carrier" %in% names(ann))
    ann[, hla_label := factor(ifelse(hla_carrier == 1, "HLA+", "HLA-"))]

cluster_levels <- c(paste0("C", sort(unique(umap_sub$cluster))), "None")
ann[, cluster_f := factor(cluster_f, levels = cluster_levels)]
clust_col_map  <- c(
    setNames(CLUST_COLS[as.character(sort(unique(umap_sub$cluster)))],
             paste0("C", sort(unique(umap_sub$cluster)))),
    "None" = "grey88"
)

# ── 9. Panel A: burden overlay ────────────────────────────────────────────────
burden_cols <- c("None" = "grey88", "1 block" = "#56B4E9", "\u22652 blocks" = "#CC0066")

pA <- ggplot(ann, aes(U1, U2, colour = burden)) +
    geom_point(size = 0.55, alpha = 0.65, stroke = 0) +
    scale_colour_manual(values = burden_cols, name = "Pre-MS enriched blocks",
                        na.value = "grey88") +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = "a  Pre-MS comorbidity burden",
         subtitle = sprintf("n=%d MS cases | %.0f%% with \u22651 enriched block pre-Dx",
                            nrow(ann), 100 * mean(ann$n_enr_blocks > 0))) +
    theme_ukb(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), aspect.ratio = 1) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# ── 10. Panel B: cluster UMAP ─────────────────────────────────────────────────
clust_centroids <- ann[!is.na(U1) & cluster_f != "None",
                        .(U1 = median(U1), U2 = median(U2)), by = cluster_f]

pB <- ggplot(ann, aes(U1, U2, colour = cluster_f)) +
    geom_point(size = 0.55, alpha = 0.65, stroke = 0) +
    geom_label(data = clust_centroids, aes(U1, U2, label = cluster_f),
               colour = "grey15", fill = alpha("white", 0.8),
               size = 2.5, linewidth = 0, inherit.aes = FALSE) +
    scale_colour_manual(values = clust_col_map, name = "Cluster", na.value = "grey88") +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = "b  Pre-MS cluster assignments",
         subtitle = sprintf("k=%d | %d blocks | nn=%d md=%.2f | sil=%.2f",
                            best_config$k, length(enr_blocks),
                            best_config$nn, best_config$md, best_sil)) +
    theme_ukb(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), aspect.ratio = 1) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# ── 11. Panel C: enriched block small multiples ───────────────────────────────
blk_long <- melt(ann[, c("eid", "U1", "U2", enr_blocks), with = FALSE],
                 id.vars = c("eid", "U1", "U2"),
                 variable.name = "block", value.name = "present")

pC <- ggplot(blk_long, aes(U1, U2, colour = factor(present))) +
    geom_point(size = 0.35, alpha = 0.45, stroke = 0) +
    facet_wrap(~block, ncol = ceiling(sqrt(length(enr_blocks)))) +
    scale_colour_manual(values = c("0" = "grey78", "1" = "#CC0066"),
                        labels = c("Absent", "Present"), name = NULL) +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = "c  Enriched block distribution (pre-MS, OR>1)") +
    theme_ukb(base_size = 7.5) +
    theme(strip.text = element_text(size = 5.5, face = "bold"),
          axis.text = element_blank(), axis.ticks = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# ── 12. Panel D: biological overlays ─────────────────────────────────────────
make_overlay <- function(data, col_var, col_vals, title_str) {
    data_sub <- data[!is.na(get(col_var))]
    data_sub[, .col := as.character(get(col_var))]
    ggplot(data_sub, aes(U1, U2, colour = .col)) +
        geom_point(size = 0.55, alpha = 0.55, stroke = 0) +
        scale_colour_manual(values = col_vals, name = NULL, na.value = "grey88") +
        labs(x = "UMAP 1", y = "UMAP 2", title = title_str) +
        theme_ukb(base_size = 8) +
        theme(axis.text = element_blank(), axis.ticks = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
}

pD1 <- make_overlay(ann, "sex_label",
                    c("Female" = "#CC0066", "Male" = "#2B4C7E"), "Sex")
pD2 <- if ("hla_label" %in% names(ann)) {
    make_overlay(ann, "hla_label",
                 c("HLA+" = "#E6A817", "HLA-" = "grey60"), "HLA-DRB1*15:01")
} else plot_spacer()
pD3 <- make_overlay(ann, "age_dx_q",
                    c("Early (<40)" = "#56B4E9", "Mid (40\u201350)" = "#E6A817",
                      "Late (>50)" = "#CC0066"), "Age at Dx")
pD4 <- make_overlay(ann, "ms_status",
                    c("pre_onset" = "#CC0066", "post_onset" = "grey60"),
                    "Sampling vs onset")
pD <- (pD1 | pD2 | pD3 | pD4) +
    plot_annotation(title = "d  Biological annotations",
                    theme = theme(plot.title = element_text(size = 9)))

# ── 13. Block overrepresentation per cluster ──────────────────────────────────
cat("\nBlock overrepresentation (Fisher, clustered patients only)...\n")
ann_c <- ann[!is.na(cluster)]
fish_blk <- list()
for (blk in enr_blocks) {
    for (cl in sort(unique(ann_c$cluster))) {
        in_cl  <- ann_c$cluster == cl
        has_b  <- ann_c[[blk]] == 1L
        a <- sum( in_cl &  has_b); b <- sum(!in_cl &  has_b)
        cv <- sum( in_cl & !has_b); d <- sum(!in_cl & !has_b)
        if (a < 2L) next
        ft <- fisher.test(matrix(c(a, b, cv, d), 2, 2), alternative = "greater")
        fish_blk[[length(fish_blk) + 1]] <- data.table(
            block = blk, cluster = cl,
            n_in = a + cv, n_with = a,
            pct_in = round(100 * a / (a + cv), 1),
            pct_out = round(100 * b / (b + d), 1),
            OR = ft$estimate, pval = ft$p.value
        )
    }
}
fish_blk_dt <- rbindlist(fish_blk)
fish_blk_dt[, fdr := p.adjust(pval, method = "BH")]
fish_blk_dt <- fish_blk_dt[order(pval)]
cat(sprintf("  Significant (FDR<0.05): %d block\u00d7cluster pairs\n",
            sum(fish_blk_dt$fdr < 0.05, na.rm = TRUE)))
cat("\nTop block enrichments:\n")
print(fish_blk_dt[fdr < 0.05, .(block, cluster = paste0("C", cluster),
                                  n_with, pct_in, pct_out,
                                  OR = round(OR, 2), fdr = round(fdr, 3))])
fwrite(fish_blk_dt, file.path(OUT_DIR, "ms_prems_block_overrep.csv"))

over_heat <- fish_blk_dt[, .(block, cluster_f = paste0("C", cluster),
                               OR_cap = pmin(OR, 10), sig = fdr < 0.05)]
over_heat[, blk_f   := factor(block, levels = rev(enr_blocks))]
over_heat[, clust_f := factor(cluster_f,
                               levels = paste0("C", sort(unique(ann_c$cluster))))]

pOverrep <- ggplot(over_heat, aes(x = clust_f, y = blk_f, fill = OR_cap)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(data = over_heat[sig == TRUE],
              aes(label = sprintf("%.1f\u00d7", OR_cap)),
              size = 2.2, colour = "grey10") +
    scale_fill_gradient(low = "white", high = "#CC0066",
                        limits = c(1, 10), oob = scales::squish,
                        na.value = "grey95", name = "OR (cap 10)") +
    labs(x = NULL, y = NULL, title = "e  Enriched block \u00d7 cluster",
         subtitle = "labels: FDR<0.05") +
    theme_ukb(base_size = 9) +
    theme(axis.text.x = element_text(face = "bold", size = 8),
          axis.text.y = element_text(size = 7))

# ── 14. Outcome panels ────────────────────────────────────────────────────────
pOut1 <- ggplot(ann[!is.na(age_at_diagnosis)],
                aes(x = cluster_f, y = age_at_diagnosis, fill = cluster_f)) +
    geom_boxplot(width = 0.55, outlier.size = 0.5, linewidth = 0.35, alpha = 0.85) +
    scale_fill_manual(values = clust_col_map, guide = "none") +
    labs(x = NULL, y = "Age at MS Dx (yrs)", title = "f  Age at diagnosis") +
    theme_ukb(base_size = 9)

pOut2 <- if ("hla_carrier" %in% names(ann)) {
    hd <- ann[!is.na(hla_carrier),
               .(pct = 100 * mean(hla_carrier == 1), n = .N), by = cluster_f]
    hd[, cluster_f := factor(cluster_f, levels = cluster_levels)]
    ggplot(hd, aes(x = cluster_f, y = pct, fill = cluster_f)) +
        geom_col(width = 0.6, alpha = 0.9) +
        geom_text(aes(label = sprintf("%.0f%%", pct)),
                  vjust = -0.3, size = 2.5, colour = "grey30") +
        geom_hline(yintercept = 100 * mean(ann$hla_carrier == 1, na.rm = TRUE),
                   linetype = "dashed", linewidth = 0.3, colour = "grey45") +
        scale_fill_manual(values = clust_col_map, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = NULL, y = "HLA-DRB1*15:01 (%)", title = "g  HLA carrier rate") +
        theme_ukb(base_size = 9)
} else plot_spacer()

pOut3 <- if ("prs_score" %in% names(ann)) {
    ggplot(ann[!is.na(prs_score)],
           aes(x = cluster_f, y = prs_score, fill = cluster_f)) +
    geom_boxplot(width = 0.55, outlier.size = 0.5, linewidth = 0.35, alpha = 0.85) +
    scale_fill_manual(values = clust_col_map, guide = "none") +
    labs(x = NULL, y = "MS PRS", title = "h  Polygenic risk score") +
    theme_ukb(base_size = 9)
} else plot_spacer()

pOut4 <- {
    sd <- ann[!is.na(sex), .(pct = 100 * mean(sex == 0), n = .N), by = cluster_f]
    sd[, cluster_f := factor(cluster_f, levels = cluster_levels)]
    ggplot(sd, aes(x = cluster_f, y = pct, fill = cluster_f)) +
        geom_col(width = 0.6, alpha = 0.9) +
        geom_text(aes(label = sprintf("%.0f%%", pct)),
                  vjust = -0.3, size = 2.5, colour = "grey30") +
        geom_hline(yintercept = 100 * mean(ann$sex == 0, na.rm = TRUE),
                   linetype = "dashed", linewidth = 0.3, colour = "grey45") +
        scale_fill_manual(values = clust_col_map, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = NULL, y = "Female (%)", title = "i  Sex") +
        theme_ukb(base_size = 9)
}

pOut5 <- ggplot(ann[!is.na(years_to_diagnosis)],
                aes(x = cluster_f, y = years_to_diagnosis, fill = cluster_f)) +
    geom_boxplot(width = 0.55, outlier.size = 0.5, linewidth = 0.35, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
    scale_fill_manual(values = clust_col_map, guide = "none") +
    labs(x = NULL, y = "Years: sampling \u2212 Dx", title = "j  Sampling timing") +
    theme_ukb(base_size = 9)

# ── 15. Assemble and save figures ─────────────────────────────────────────────
cat("\nAssembling figures...\n")

fig_umap <- (pA | pB) / pC / pD + plot_layout(heights = c(1, 2, 0.85))
ggsave(file.path(OUT_DIR, "fig5_row1_prems_umap.pdf"),
       fig_umap, width = 14, height = 18, device = cairo_pdf)
cat("  Saved: fig5_row1_prems_umap.pdf\n")

ggsave(file.path(OUT_DIR, "ms_prems_overrep_heatmap.pdf"),
       pOverrep,
       width  = max(4, length(cluster_levels) * 1.2 + 2),
       height = max(4, length(enr_blocks) * 0.38 + 2),
       device = cairo_pdf)
cat("  Saved: ms_prems_overrep_heatmap.pdf\n")

outcomes_fig <- (pOut1 | pOut2 | pOut3) / (pOut4 | pOut5 | plot_spacer()) +
    plot_annotation(
        title    = "Pre-MS cluster outcomes",
        subtitle = sprintf("n=%d | %d clusters | post-hoc", nrow(ann), best_config$k),
        theme    = theme(plot.title    = element_text(size = 10),
                         plot.subtitle = element_text(size = 8, colour = "grey40"))
    )
ggsave(file.path(OUT_DIR, "ms_prems_cluster_outcomes.pdf"),
       outcomes_fig, width = 12, height = 8, device = cairo_pdf)
cat("  Saved: ms_prems_cluster_outcomes.pdf\n")

# ── 16. Code-level specificity within clusters ───────────────────────────────
cat("\nCode-level cluster specificity (sensitivity layer)...\n")
code_ann   <- merge(code_dt, ann[, .(eid, cluster)], by = "eid")
code_ann_c <- code_ann[!is.na(cluster)]

fish_code <- list()
for (cc in loaded_codes) {
    if (!cc %in% names(code_ann_c)) next
    for (cl in sort(unique(code_ann_c$cluster))) {
        in_cl  <- code_ann_c$cluster == cl
        has_cc <- code_ann_c[[cc]] == 1L
        a <- sum( in_cl &  has_cc); b <- sum(!in_cl &  has_cc)
        cv <- sum( in_cl & !has_cc); d <- sum(!in_cl & !has_cc)
        if (a < 2L) next
        ft <- fisher.test(matrix(c(a, b, cv, d), 2, 2), alternative = "greater")
        fish_code[[length(fish_code) + 1]] <- data.table(
            icd10   = cc,
            block   = substr(cc, 1, 2),
            direction = if (cc %in% enriched_codes) "enriched" else "depleted",
            cluster = cl,
            n_in = a + cv, n_with = a,
            pct_in  = round(100 * a / (a + cv), 1),
            pct_out = round(100 * b / (b + d), 1),
            OR = ft$estimate, pval = ft$p.value
        )
    }
}
fish_code_dt <- rbindlist(fish_code)
fish_code_dt[, fdr := p.adjust(pval, method = "BH")]
fish_code_dt <- fish_code_dt[order(pval)]
cat(sprintf("  Significant code-level (FDR<0.05): %d\n",
            sum(fish_code_dt$fdr < 0.05, na.rm = TRUE)))
fwrite(fish_code_dt, file.path(OUT_DIR, "ms_prems_code_specificity.csv"))

# ── 17. Save data ─────────────────────────────────────────────────────────────
fwrite(ann, file.path(OUT_DIR, "ms_prems_block_features.csv"))

clust_assign <- ann[, .(eid, cluster, n_enr_blocks,
                          ms_status, age_at_diagnosis, sex, years_to_diagnosis)]
if ("hla_carrier" %in% names(ann)) clust_assign[, hla_carrier := ann$hla_carrier]
if ("prs_score"   %in% names(ann)) clust_assign[, prs_score   := ann$prs_score]
fwrite(clust_assign, file.path(OUT_DIR, "ms_prems_cluster_assignments.csv"))

clust_summ <- ann[, .(
    n             = .N,
    pct_any_prems = round(100 * mean(n_enr_blocks > 0), 1),
    pct_female    = round(100 * mean(sex == 0, na.rm = TRUE), 1),
    pct_pre_onset = round(100 * mean(ms_status == "pre_onset"), 1),
    mean_age_dx   = round(mean(age_at_diagnosis, na.rm = TRUE), 1),
    mean_yrs_dx   = round(mean(years_to_diagnosis, na.rm = TRUE), 1),
    pct_hla       = if ("hla_carrier" %in% names(.SD))
                        round(100 * mean(hla_carrier == 1, na.rm = TRUE), 1) else NA_real_,
    mean_prs      = if ("prs_score" %in% names(.SD))
                        round(mean(prs_score, na.rm = TRUE), 3) else NA_real_
), by = cluster_f][order(cluster_f)]
cat("\nCluster summary:\n"); print(clust_summ)
fwrite(clust_summ, file.path(OUT_DIR, "ms_prems_cluster_summary.csv"))

config_dt <- data.table(
    n_enr_blocks = length(enr_blocks),
    blocks       = paste(enr_blocks, collapse = ","),
    min_pre_n    = MIN_PRE_N,
    metric       = "jaccard_precomputed",
    n_neighbors  = best_config$nn,
    min_dist     = best_config$md,
    k            = best_config$k,
    silhouette   = round(best_sil, 4),
    n_patients   = nrow(ann),
    n_clustered  = sum(!is.na(ann$cluster))
)
fwrite(config_dt, file.path(OUT_DIR, "ms_prems_umap_config.csv"))
cat("\nUMAP config:\n"); print(config_dt)
cat("\n02_ms_prems_comorbidity_endophenotype.R complete.\n")
