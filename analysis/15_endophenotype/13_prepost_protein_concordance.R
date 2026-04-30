#!/usr/bin/env Rscript
# 13_prepost_protein_concordance.R — Pre/post protein signature concordance
#
# Addresses the question: do cluster-specific protein signatures found across
# all MS patients (script 04) replicate in post-diagnosis MS alone (script 08)?
#
# Since UK Biobank has no within-person longitudinal proteomics (each
# participant has one blood draw), concordance is tested cross-sectionally by
# correlating per-protein logFC estimates from the two analyses:
#   x-axis: logFC from 04_ms_cluster_proteomics.R   (pre + post-onset MS)
#   y-axis: logFC from 08_postms_cluster_proteomics.R (post-onset MS only)
#
# High Pearson r → pre-MS endophenotype signature persists in established disease
# Low r           → signature is largely driven by pre-diagnostic biology
#
# Outputs (results/endophenotype/prepost_concordance/):
#   concordance_stats.csv             — r, rho, n_proteins per contrast
#   panel_l/m/n.pdf                   → Figure 5S

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(here)
    library(glue)
})

source(here::here("analysis", "helpers", "disease_config.R"))
source(here::here("analysis", "helpers", "ukb_theme.R"))

cfg <- load_disease_config()

ALL_DEP_FILE  <- here::here("results", "endophenotype",
                            "cluster_proteomics",       "cluster_deps_all_contrasts.csv")
POST_DEP_FILE <- here::here("results", "endophenotype",
                            "postms_cluster_proteomics", "postms_deps_all_contrasts.csv")
OUT_DIR <- here::here("results", "endophenotype", "prepost_concordance")
FIG_DIR <- here::here("results", "figures", "5S")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

N_LABEL <- 8L

all_deps  <- fread(ALL_DEP_FILE)
post_deps <- fread(POST_DEP_FILE)

contrasts <- c("C0_vs_None", "C1_vs_None", "C2_vs_None")
panel_labels <- c("C0_vs_None" = "l", "C1_vs_None" = "m", "C2_vs_None" = "n")
cluster_names <- c("C0_vs_None" = "C0", "C1_vs_None" = "C1", "C2_vs_None" = "C2")

stats_list <- list()

for (cname in contrasts) {
    a <- all_deps[contrast  == cname, .(protein, logFC_all  = logFC, P_all  = P.Value,
                                         fdr_all  = adj.P.Val)]
    p <- post_deps[contrast == cname, .(protein, logFC_post = logFC, P_post = P.Value,
                                         fdr_post = adj.P.Val)]

    m <- merge(a, p, by = "protein")
    m <- m[!is.na(logFC_all) & !is.na(logFC_post)]

    r   <- cor(m$logFC_all, m$logFC_post, method = "pearson")
    rho <- cor(m$logFC_all, m$logFC_post, method = "spearman")

    stats_list[[cname]] <- data.table(
        contrast = cname, cluster = cluster_names[cname],
        n_proteins = nrow(m), pearson_r = round(r, 3), spearman_rho = round(rho, 3)
    )
    cat(sprintf("  %s: n=%d  r=%.3f  rho=%.3f\n", cname, nrow(m), r, rho))

    sig_all_label <- glue("FDR<0.05 (all {cfg$disease_short_caps})")
    m[, sig := fcase(
        fdr_all  < 0.05 & fdr_post < 0.05, "FDR<0.05 (both)",
        fdr_all  < 0.05,                   sig_all_label,
        fdr_post < 0.05,                   "FDR<0.05 (post-onset)",
        P_all    < 0.05 | P_post < 0.05,   "Nominal (either)",
        default = "NS"
    )]
    m[, sig := factor(sig, levels = c("FDR<0.05 (both)", sig_all_label,
                                       "FDR<0.05 (post-onset)", "Nominal (either)", "NS"))]

    sig_cols <- setNames(
        c("#7B2FBE",
          unname(CLUST_COLS[cluster_names[cname]]),
          "#E6740A",
          "grey60",
          "grey88"),
        c("FDR<0.05 (both)", sig_all_label,
          "FDR<0.05 (post-onset)", "Nominal (either)", "NS")
    )
    m[, alpha_val := fcase(sig == "NS", 0.25, default = 0.75)]

    label_set <- rbind(
        m[fdr_all  < 0.05][order(abs(logFC_all),  decreasing = TRUE)][seq_len(min(N_LABEL, .N))],
        m[fdr_post < 0.05][order(abs(logFC_post), decreasing = TRUE)][seq_len(min(N_LABEL, .N))]
    )
    label_set <- unique(label_set, by = "protein")

    xy_lim <- max(abs(c(m$logFC_all, m$logFC_post)), na.rm = TRUE) * 1.05

    p_plot <- ggplot(m, aes(x = logFC_all, y = logFC_post, colour = sig)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    linewidth = 0.35, colour = "grey50") +
        geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
        geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
        geom_point(aes(alpha = alpha_val), size = 1.2) +
        geom_text_repel(
            data = label_set, aes(label = protein),
            size = 2.3, max.overlaps = 18, seed = 42,
            segment.size = 0.25, segment.colour = "grey50"
        ) +
        scale_colour_manual(values = sig_cols, name = NULL,
                            guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        scale_alpha_identity(guide = "none") +
        annotate("text",
                 x = -xy_lim * 0.95, y = xy_lim * 0.95,
                 label = sprintf("r = %.3f\n\u03c1 = %.3f\nn = %d", r, rho, nrow(m)),
                 hjust = 0, vjust = 1, size = 2.8, colour = "grey30") +
        coord_fixed(xlim = c(-xy_lim, xy_lim), ylim = c(-xy_lim, xy_lim)) +
        labs(
            x        = glue("logFC  (all {cfg$disease_short_caps}: pre + post-onset)"),
            y        = glue("logFC  (post-onset {cfg$disease_short_caps} only)"),
            title    = sprintf("%s  %s — protein signature concordance",
                               panel_labels[cname], cluster_names[cname]),
            subtitle = "Each point = one protein; dashed line = perfect concordance"
        ) +
        theme_ukb(base_size = 9) +
        theme(
            aspect.ratio    = 1,
            legend.position = "bottom",
            legend.text     = element_text(size = 7),
            plot.subtitle   = element_text(size = 7, colour = "grey40")
        )

    out <- file.path(FIG_DIR, sprintf("panel_%s.pdf", panel_labels[cname]))
    ggsave(out, p_plot, width = 4.2, height = 4.8, device = cairo_pdf)
    cat(sprintf("  panel_%s.pdf\n", panel_labels[cname]))
}

stats_dt <- rbindlist(stats_list)
fwrite(stats_dt, file.path(OUT_DIR, "concordance_stats.csv"))
cat("\nConcordance stats:\n"); print(stats_dt)

cat("\n13_prepost_protein_concordance.R complete.\n")
cat("Panels l/m/n.pdf → results/figures/5S/\n")
