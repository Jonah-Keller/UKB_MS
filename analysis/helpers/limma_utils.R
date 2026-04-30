#!/usr/bin/env Rscript
# limma_utils.R — safe wrappers for limma::topTable
#
# topTable_safe() prevents the protein-to-statistic misalignment that occurs
# when topTable is called twice with different sort.by values and rownames from
# one call are paired with statistics from the other.
#
# Usage:
#   source(file.path(PROJ_DIR, "analysis", "helpers", "limma_utils.R"))
#   res <- topTable_safe(fit, coef = "ms - control")
#   res <- topTable_safe(fit, coef = "is_ms:sex_fFemale")

# ── compute_control_pcs() ─────────────────────────────────────────────────────
# Fit PCA on healthy-control samples only (center=TRUE, scale.=FALSE), then
# project all samples into that space.  Returns a data.table with columns
# (eid, PC1, PC2).  Mean-imputes missing values per protein before PCA.
#
# Rationale (matches 04_differential/01_limma_ms_vs_hc.R): control-derived PCs
# capture technical/batch variance without disease-signal contamination and
# avoid collinearity with the outcome variable in downstream DE/Cox models.
compute_control_pcs <- function(dt_filt, protein_cols, control_value = "control",
                                 n_pc = 2L, use_irlba = TRUE) {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table not available")
    prot_mat <- t(as.matrix(dt_filt[, ..protein_cols]))
    colnames(prot_mat) <- as.character(dt_filt$eid)
    row_means <- rowMeans(prot_mat, na.rm = TRUE)
    for (i in seq_len(nrow(prot_mat))) {
        nas <- is.na(prot_mat[i, ])
        if (any(nas)) prot_mat[i, nas] <- row_means[i]
    }
    ctrl_eids <- as.character(dt_filt[ms_status == control_value, eid])
    ctrl_x    <- t(prot_mat[, ctrl_eids])                 # samples × proteins
    ctrl_ctr  <- colMeans(ctrl_x)
    ctrl_xc   <- sweep(ctrl_x, 2L, ctrl_ctr, "-")         # HC-centred

    # Truncated SVD when irlba is available — ~100× faster than full prcomp on
    # 48K × 2911 matrices and yields identical top-k PCs up to sign flips.
    # Falls back to base prcomp otherwise.
    if (use_irlba && requireNamespace("irlba", quietly = TRUE) &&
        n_pc < min(dim(ctrl_xc)) - 1L) {
        sv <- irlba::irlba(ctrl_xc, nv = n_pc)
        rotation <- sv$v                                  # proteins × n_pc
        sdev_pc  <- sv$d / sqrt(max(1L, nrow(ctrl_xc) - 1L))
        # Align signs to prcomp convention (largest |loading| is positive).
        for (k in seq_len(n_pc)) {
            i_max <- which.max(abs(rotation[, k]))
            if (rotation[i_max, k] < 0) rotation[, k] <- -rotation[, k]
        }
        total_var <- sum(apply(ctrl_xc, 2L, var))
        pc_var    <- 100 * (sdev_pc ^ 2) / total_var
    } else {
        pc_fit   <- prcomp(ctrl_x, center = TRUE, scale. = FALSE)
        rotation <- pc_fit$rotation[, seq_len(n_pc), drop = FALSE]
        pc_var   <- 100 * pc_fit$sdev[seq_len(n_pc)] ^ 2 /
                           sum(pc_fit$sdev ^ 2)
    }

    # Project every sample (including non-controls) into the HC PC space.
    all_x   <- sweep(t(prot_mat), 2L, ctrl_ctr, "-")
    all_pcs <- all_x %*% rotation

    message(sprintf("  Control PC1: %.1f%%  PC2: %.1f%% variance explained",
                    pc_var[1], if (n_pc >= 2) pc_var[2] else NA_real_))

    out <- data.table::data.table(eid = as.integer(colnames(prot_mat)))
    for (k in seq_len(n_pc))
        data.table::set(out, j = paste0("PC", k), value = all_pcs[, k])
    out
}

topTable_safe <- function(fit, coef, number = Inf, ...) {
    if (!requireNamespace("limma", quietly = TRUE)) stop("limma not available")
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table not available")
    res_raw <- limma::topTable(fit, coef = coef, number = number, sort.by = "none", ...)
    dt <- data.table::as.data.table(res_raw)
    dt[, protein := rownames(res_raw)]
    dt[order(P.Value)]
}

# ── make_volcano_plot() ───────────────────────────────────────────────────────
# Standard cluster-vs-None volcano used across Figures 5 and 5S.
# Caller must have sourced ukb_theme.R and loaded ggplot2 + ggrepel.
make_volcano_plot <- function(
    tt,               # data.table from topTable_safe with sig column
    title_str,
    panel_label,
    clust_col,        # colour for up-regulated proteins
    bonf,             # -log10(Bonferroni threshold), precomputed
    n_label         = 10L,
    subtitle_suffix = "ref = MS-None (no pre-MS comorbidities)"
) {
    tt <- data.table::copy(tt)
    tt[, log10p := -log10(P.Value)]
    tt[, colour := data.table::fcase(
        sig & logFC  > 0, clust_col,
        sig & logFC  < 0, COL_DOWN,
        !sig & P.Value < 0.05 & logFC > 0, scales::alpha(clust_col, 0.65),
        !sig & P.Value < 0.05 & logFC < 0, scales::alpha(COL_DOWN, 0.65),
        default = "grey72"
    )]
    n_up    <- sum(tt$sig & tt$logFC > 0, na.rm = TRUE)
    n_down  <- sum(tt$sig & tt$logFC < 0, na.rm = TRUE)
    # Label FDR-sig hits; if none, fall back to top nominal hits
    top_lbl <- if (sum(tt$sig, na.rm = TRUE) > 0)
        tt[sig == TRUE][order(-abs(logFC))][seq_len(min(n_label, .N))]
    else
        tt[P.Value < 0.05][order(P.Value)][seq_len(min(n_label, .N))]

    ggplot2::ggplot(tt, ggplot2::aes(logFC, log10p, colour = colour)) +
        ggplot2::geom_hline(yintercept = bonf, linetype = "dashed",
                            linewidth = 0.3, colour = "grey45") +
        ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
        ggplot2::geom_point(size = 0.55, alpha = 0.55, stroke = 0) +
        ggrepel::geom_text_repel(
            data = top_lbl,
            ggplot2::aes(label = toupper(protein)),
            colour = "grey20",
            size = 2.0, fontface = "bold", max.overlaps = 20,
            segment.size = 0.22, segment.colour = "grey55",
            box.padding = 0.28, show.legend = FALSE
        ) +
        ggplot2::annotate("text",
                          x = max(tt$logFC, na.rm = TRUE) * 0.9,
                          y = bonf + 0.35, label = "Bonferroni",
                          size = 2.0, hjust = 1, colour = "grey45") +
        ggplot2::scale_colour_identity() +
        ggplot2::labs(
            x        = expression(log[2]~"fold change  (cluster vs MS-None)"),
            y        = expression(-log[10](italic(P))),
            title    = sprintf("%s  %s", panel_label, title_str),
            subtitle = sprintf("Up: %d  |  Down: %d  (FDR<0.05, n=%d proteins)  |  %s",
                               n_up, n_down, nrow(tt), subtitle_suffix)
        ) +
        theme_ukb(base_size = 9) +
        ggplot2::theme(legend.position = "none")
}

# ── run_cluster_psm_limma() ───────────────────────────────────────────────────
# PSM (1:psm_ratio NN) + limma for C0/C1/C2 vs None within an MS cohort.
# time_var: column name of the temporal covariate ("years_to_diagnosis" for
# pre-onset, "years_since_diagnosis" for post-onset). The same column is used
# in both the PSM formula and the limma design matrix.
#
# Returns: list(all_dep, bonf, psm_data)
#   all_dep  — named list of data.tables (one per contrast) from topTable_safe
#   bonf     — Bonferroni -log10 threshold for volcanos
#   psm_data — the PSM-matched analysis data.table
run_cluster_psm_limma <- function(
    ms_data,
    prot_cols,
    time_var,
    psm_ratio = 10L,
    psm_cal   = 0.2,
    fdr_thr   = 0.05,
    bonf_thr  = 0.05
) {
    ms_data <- data.table::copy(ms_data)
    if (!"sex_num" %in% names(ms_data))
        ms_data[, sex_num := as.integer(sex == "Female")]

    psm_ids <- list()
    for (cl in c("C0", "C1", "C2")) {
        sub_all  <- ms_data[cluster_f %in% c(cl, "None")]
        sub_all  <- sub_all[!is.na(get(time_var))]
        psm_df   <- data.frame(
            eid             = sub_all$eid,
            treat           = as.integer(sub_all$cluster_f == cl),
            age_at_sampling = sub_all$age_at_sampling,
            sex_num         = sub_all$sex_num
        )
        psm_df[[time_var]] <- sub_all[[time_var]]

        if (sum(psm_df$treat) < 3L) {
            cat(sprintf("  %s: too few treated, skipping PSM\n", cl)); next
        }
        psm_form <- as.formula(
            paste("treat ~ age_at_sampling + sex_num +", time_var)
        )
        m <- tryCatch(
            MatchIt::matchit(psm_form, data = psm_df,
                             method = "nearest", ratio = psm_ratio,
                             caliper = psm_cal, std.caliper = TRUE),
            error = function(e) {
                cat(sprintf("  %s PSM failed: %s\n", cl, conditionMessage(e))); NULL
            }
        )
        if (is.null(m)) next
        psm_ids[[cl]] <- sub_all$eid[m$weights > 0]
        n_t <- sum(psm_df$treat[m$weights > 0])
        n_c <- sum((1 - psm_df$treat)[m$weights > 0])
        cat(sprintf("  %s vs None: %d treated, %d controls after PSM\n", cl, n_t, n_c))
    }
    if (length(psm_ids) == 0L) stop("PSM failed for all clusters")

    psm_data <- ms_data[eid %in% union(ms_data[cluster_f != "None", eid],
                                        unique(unlist(psm_ids)))]
    cat(sprintf("  PSM set: %d (%d cluster, %d None)\n",
                nrow(psm_data), sum(psm_data$cluster_f != "None"),
                sum(psm_data$cluster_f == "None")))

    mat    <- t(as.matrix(psm_data[, ..prot_cols]))
    colnames(mat) <- psm_data$eid
    design <- model.matrix(
        as.formula(paste("~ 0 + cluster_f + age_at_sampling + sex +", time_var)),
        data = psm_data
    )
    colnames(design) <- gsub("cluster_f", "", colnames(design))

    cmat <- limma::makeContrasts(
        C0_vs_None = C0 - None, C1_vs_None = C1 - None, C2_vs_None = C2 - None,
        C1_vs_C0   = C1 - C0,  C2_vs_C0   = C2 - C0,  C2_vs_C1   = C2 - C1,
        levels     = design
    )
    fit2 <- limma::eBayes(
        limma::contrasts.fit(limma::lmFit(mat, design), cmat),
        trend = TRUE, robust = TRUE
    )

    bonf    <- -log10(bonf_thr / nrow(mat))
    all_dep <- lapply(colnames(cmat), function(cname) {
        tt <- topTable_safe(fit2, coef = cname, adjust.method = "BH")
        tt[, contrast := cname][, sig := adj.P.Val < fdr_thr]
    })
    names(all_dep) <- colnames(cmat)

    list(all_dep = all_dep, bonf = bonf, psm_data = psm_data)
}
