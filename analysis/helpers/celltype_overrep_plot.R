#!/usr/bin/env Rscript
# celltype_overrep_plot.R — cell-type enrichment helpers + plots.
#
# ─── WHAT CHANGED (2025 refactor) ────────────────────────────────────────────
# The previous pipeline filtered markers with a raw abundance threshold
# (HPA `nTPM >= 10`) and tested them with Fisher's exact on up/down DEP lists.
# Marker sets covered 35-50% of the Olink universe and overlapped massively,
# which made the Fisher test measure "is the DEP list enriched for tissue-
# expressed proteins in general" rather than "is it enriched for cell-type-
# *specific* markers." Significance was driven by marker-set size, not biology.
#
# The current pipeline fixes both issues:
#   (1) MARKERS: specificity-filtered using a tau + per-cell-type z-score on
#       log1p(nTPM). Genes qualify as markers of cell type ct only if they are
#       tissue-specific (tau high) AND highly expressed in ct relative to the
#       gene's own across-cell-type mean (z high). Typical set sizes drop from
#       ~1100 to ~30-400 per cell type, with minimal cross-cell-type overlap.
#   (2) TEST: signed, rank-based preranked GSEA (`fgsea::fgsea`) on the full
#       moderated-t statistic from limma. This uses direction natively (no
#       up/down double-testing), preserves sub-threshold signal, and returns
#       a single signed NES per (contrast x cell type) with BH-adjusted p.
#       An optional `cameraPR` pass (accounts for inter-gene correlation) is
#       available for confirmation.
#
# Exports (primary):
#   CT_LABELS, CT_ORDER, CT_VASC, CT_CNS, CT_IMMUNE
#   load_celltype_markers()  — specificity-filtered marker sets (default)
#   run_celltype_gsea()       — preranked fgsea on moderated-t
#   run_celltype_camerapr()   — optional cameraPR companion (inter-gene corr)
#   make_celltype_gsea_plot() — signed-NES dot plot (1 col per contrast)
#   make_celltype_heatmap()   — multi-contrast signed-NES heatmap
#
# INTERPRETATION (important; use this wording in figures and prose):
#   The ranking statistic is the signed limma moderated-t from the stated
#   contrast (typically case vs matched healthy control). Positive t means
#   the protein is HIGHER in cases than in HC; negative t means LOWER.
#   fgsea asks, per cell-type marker set, whether the set is concentrated
#   at one end of that ranking.
#     NES > 0 :: markers of this cell type are coordinately HIGHER in cases
#                than in HC (i.e., the cell type's plasma footprint is
#                elevated in cases).
#     NES < 0 :: markers are coordinately LOWER in cases than in HC (the
#                cell type's plasma footprint is reduced in cases).
#   This is a plasma-level readout. It does NOT imply the cells themselves
#   are dying, expanding, or functionally activated — plasma abundance is
#   downstream of secretion, clearance, barrier permeability, and turnover.
#   "Depleted markers" language from the old Fisher-on-up/down pipeline is
#   obsolete and should not be used; use "lower/higher in cases vs HC".
#
# Legacy (retained for back-compatibility only; not recommended for new work):
#   run_ora(), make_celltype_overrep_plot()
#
# Required columns in gsea_dt passed to make_celltype_gsea_plot():
#   cell_type, contrast (or analysis), NES, pval, padj
#
# Required columns in heat_dt passed to make_celltype_heatmap():
#   cluster, cell_type, NES, padj   (legacy fallback: direction, odds_ratio)
#
# Caller must have sourced ukb_theme.R and loaded ggplot2 + data.table.

suppressPackageStartupMessages({
    if (!requireNamespace("fgsea",    quietly = TRUE)) stop("fgsea not installed")
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table not installed")
})

# ── Canonical cell-type labels and order ─────────────────────────────────────
# Key = raw marker-source label (HPA long-format cell_type / Walchli col tag)
# Value = display label used in figures.
CT_LABELS <- c(
    "Endothelial"                     = "Endothelial",
    "Pericytes"                       = "Pericytes",
    "Smooth muscle"                   = "Smooth muscle",
    "Fibroblasts"                     = "Fibroblasts",
    "Microglial cells"                = "Microglia",
    "Astrocytes"                      = "Astrocytes",
    "Oligodendrocytes"                = "Oligodendrocytes",
    "Oligodendrocyte precursor cells" = "OPCs",
    "Excitatory neurons"              = "Exc. neurons",
    "Inhibitory neurons"              = "Inh. neurons",
    "T-cells"                         = "T cells",
    "B-cells"                         = "B cells",
    "NK-cells"                        = "NK cells",
    "monocytes"                       = "Monocytes",
    "dendritic cells"                 = "Dendritic cells",
    "Macrophages"                     = "Macrophages",
    "Plasma cells"                    = "Plasma cells",
    "granulocytes"                    = "Granulocytes"
)

CT_ORDER <- c(
    "Endothelial", "Pericytes", "Smooth muscle", "Fibroblasts",
    "Microglia", "Astrocytes", "Oligodendrocytes", "OPCs",
    "Exc. neurons", "Inh. neurons",
    "T cells", "B cells", "NK cells", "Monocytes",
    "Dendritic cells", "Macrophages", "Plasma cells", "Granulocytes"
)

CT_VASC   <- c("Endothelial", "Pericytes", "Smooth muscle", "Fibroblasts")
CT_CNS    <- c("Microglia", "Astrocytes", "Oligodendrocytes", "OPCs",
               "Exc. neurons", "Inh. neurons")
CT_IMMUNE <- c("T cells", "B cells", "NK cells", "Monocytes",
               "Dendritic cells", "Macrophages", "Plasma cells", "Granulocytes")

# Mapping between Walchli wide-format column names and helper keys.
.WALCHLI_COL2KEY <- c(
    "Endothelial.cells"   = "Endothelial",
    "Pericytes"           = "Pericytes",
    "Smooth.muscle.cells" = "Smooth muscle",
    "Fibroblasts"         = "Fibroblasts"
)
.WALCHLI_ALL_COLS <- c(
    "Endothelial.cells", "Pericytes", "Smooth.muscle.cells", "Fibroblasts",
    "Oligodendrocytes", "Microglia.and.Macrophages",
    "Astrocytes", "Neurons", "Stem.cells", "Neuron.progenitor", "T.cells"
)

# HPA cell types we export markers for (raw labels, match CT_LABELS keys).
.HPA_CT_CNS    <- c("Microglial cells", "Astrocytes", "Oligodendrocytes",
                    "Oligodendrocyte precursor cells",
                    "Excitatory neurons", "Inhibitory neurons")
.HPA_CT_IMMUNE <- c("T-cells", "B-cells", "NK-cells", "monocytes",
                    "dendritic cells", "Macrophages", "Plasma cells", "granulocytes")

# ── Internal: tau (Yanai 2005) + per-ct z-score from a gene x cell type matrix
# Both computed on log1p(nTPM). Rows with max expression < `min_max_log` are
# returned with tau = NA (effectively "not expressed anywhere").
.compute_specificity <- function(mat, min_max_log = log1p(1)) {
    stopifnot(is.matrix(mat))
    row_max  <- apply(mat, 1, max, na.rm = TRUE)
    row_mean <- rowMeans(mat, na.rm = TRUE)
    row_sd   <- apply(mat, 1, sd, na.rm = TRUE)
    n_ct     <- ncol(mat)
    xhat     <- mat / pmax(row_max, 1e-9)
    tau      <- rowSums(1 - xhat, na.rm = TRUE) / (n_ct - 1)
    tau[row_max < min_max_log] <- NA_real_
    z <- (mat - row_mean) / pmax(row_sd, 1e-9)
    list(tau = tau, z = z, row_max = row_max)
}

# ── Internal: HPA long-format TSV → specificity + marker sets ────────────────
.load_hpa_markers <- function(hpa_file, cts_keep,
                               tau_thr   = 0.60,
                               z_thr     = 1.50,
                               min_ntpm  = 1.0,
                               cache_tsv = NULL,
                               verbose   = TRUE) {
    # Decompress once if given a zip.
    hpa_tsv <- sub("\\.zip$", ".tsv", hpa_file)
    if (file.exists(hpa_tsv)) {
        hpa <- data.table::fread(hpa_tsv, showProgress = FALSE)
    } else if (file.exists(hpa_file)) {
        if (verbose) cat("  Decompressing HPA zip (once)...\n")
        hpa <- data.table::fread(cmd = paste("unzip -p", shQuote(hpa_file)),
                                  showProgress = FALSE)
        data.table::fwrite(hpa, hpa_tsv, sep = "\t")
    } else {
        stop("HPA file not found: ", hpa_file)
    }
    data.table::setnames(hpa,
        c("Gene", "Gene name", "Cell type", "nTPM"),
        c("ensembl", "gene_name", "cell_type", "nTPM"),
        skip_absent = TRUE
    )
    hpa[, gene_name := toupper(gene_name)]

    # Pivot: genes x cell types, log1p(nTPM). Any cell type absent for a gene
    # is filled with 0 (HPA is typically fully rectangular, but guard anyway).
    wide <- data.table::dcast(hpa, gene_name ~ cell_type,
                               value.var = "nTPM", fun.aggregate = mean,
                               fill = 0)
    genes  <- wide$gene_name
    ct_all <- setdiff(names(wide), "gene_name")
    mat    <- log1p(as.matrix(wide[, ..ct_all]))
    rownames(mat) <- genes

    spec <- .compute_specificity(mat, min_max_log = log1p(min_ntpm))

    if (!is.null(cache_tsv)) {
        cache_dt <- data.table::data.table(
            gene_name = genes,
            tau       = spec$tau,
            max_log1p_ntpm = spec$row_max
        )
        data.table::fwrite(cache_dt, cache_tsv, sep = "\t")
    }

    markers <- list()
    for (ct in cts_keep) {
        if (!ct %in% ct_all) next
        z_ct    <- spec$z[, ct]
        expr_ct <- mat[, ct]
        keep    <- !is.na(spec$tau) &
                   spec$tau >= tau_thr &
                   z_ct    >= z_thr &
                   expr_ct >= log1p(min_ntpm)
        markers[[ct]] <- sort(unique(genes[keep]))
        if (verbose) {
            cat(sprintf("  HPA %-35s tau>=%.2f z>=%.2f -> %d markers\n",
                        ct, tau_thr, z_thr, length(markers[[ct]])))
        }
    }
    markers
}

# ── Internal: Walchli wide CSV → specificity + marker sets ───────────────────
.load_walchli_markers <- function(walchli_csv,
                                   tau_thr  = 0.45,
                                   z_thr    = 1.00,
                                   min_ntpm = 0.5,
                                   verbose  = TRUE) {
    if (!file.exists(walchli_csv)) {
        if (verbose) cat("  Walchli CSV missing - skipping vascular\n")
        return(list())
    }
    w <- data.table::fread(walchli_csv, showProgress = FALSE)
    data.table::setnames(w, "Gene", "gene_name")
    w[, gene_name := toupper(gene_name)]
    avail_all <- intersect(.WALCHLI_ALL_COLS, names(w))
    if (length(avail_all) < 3) {
        if (verbose) cat("  Walchli columns not as expected - skipping\n")
        return(list())
    }
    mat <- log1p(as.matrix(w[, ..avail_all]))
    rownames(mat) <- w$gene_name

    spec <- .compute_specificity(mat, min_max_log = log1p(min_ntpm))

    markers <- list()
    for (col in names(.WALCHLI_COL2KEY)) {
        if (!col %in% avail_all) next
        key     <- .WALCHLI_COL2KEY[[col]]
        z_ct    <- spec$z[, col]
        expr_ct <- mat[, col]
        keep    <- !is.na(spec$tau) &
                   spec$tau >= tau_thr &
                   z_ct    >= z_thr &
                   expr_ct >= log1p(min_ntpm)
        markers[[key]] <- sort(unique(w$gene_name[keep]))
        if (verbose) {
            cat(sprintf("  Walchli %-20s tau>=%.2f z>=%.2f -> %d markers\n",
                        key, tau_thr, z_thr, length(markers[[key]])))
        }
    }
    markers
}

# ── load_celltype_markers() ───────────────────────────────────────────────────
# Build specificity-filtered marker sets for Walchli vascular + HPA CNS/immune.
#
# Args:
#   walchli_csv   Path to data/reference/walchli_avg_expr_by_celltype.csv
#   hpa_file       Path to HPA v23 single-cell TSV (zip or tsv).
#   hpa_tau_thr    HPA tau cutoff (default 0.60; increase -> fewer, more specific)
#   hpa_z_thr      HPA per-cell-type z cutoff (default 1.50)
#   walchli_tau_thr / walchli_z_thr  Walchli equivalents (looser because only
#                  11 cell types in Walchli matrix).
#   method         "specificity" (default, new) or "legacy_ntpm" (old behaviour
#                  kept for backwards reproducibility: nTPM >= 10 on HPA +
#                  top-percentile on Walchli).
#   hpa_ntpm_thr   Used only when method = "legacy_ntpm".
#   walchli_pctile Used only when method = "legacy_ntpm".
#   walchli_min_lp Used only when method = "legacy_ntpm".
#   verbose        Print set sizes to stdout.
#
# Returns a named list of marker gene-symbol vectors (ALL UPPERCASE).
load_celltype_markers <- function(walchli_csv, hpa_file,
                                   method          = c("specificity", "legacy_ntpm"),
                                   hpa_tau_thr     = 0.60,
                                   hpa_z_thr       = 1.50,
                                   hpa_min_ntpm    = 1.0,
                                   walchli_tau_thr = 0.45,
                                   walchli_z_thr   = 1.00,
                                   walchli_min_ntpm = 0.5,
                                   hpa_ntpm_thr    = 10,
                                   walchli_pctile  = 0.65,
                                   walchli_min_lp  = 1.0,
                                   verbose         = TRUE) {
    method <- match.arg(method)

    if (method == "legacy_ntpm") {
        return(.load_celltype_markers_legacy(
            walchli_csv, hpa_file,
            walchli_pctile = walchli_pctile,
            walchli_min_lp = walchli_min_lp,
            hpa_ntpm_thr   = hpa_ntpm_thr,
            verbose        = verbose
        ))
    }

    # Specificity-based (primary pathway).
    cts_hpa <- c(.HPA_CT_CNS, .HPA_CT_IMMUNE)
    markers_hpa <- .load_hpa_markers(
        hpa_file, cts_keep = cts_hpa,
        tau_thr  = hpa_tau_thr,
        z_thr    = hpa_z_thr,
        min_ntpm = hpa_min_ntpm,
        verbose  = verbose
    )
    markers_walchli <- .load_walchli_markers(
        walchli_csv,
        tau_thr  = walchli_tau_thr,
        z_thr    = walchli_z_thr,
        min_ntpm = walchli_min_ntpm,
        verbose  = verbose
    )
    # Walchli vascular sets take precedence for vascular labels; HPA for the
    # rest. If HPA also has "Fibroblasts" (it does), keep the Walchli version
    # for brain-vascular specificity.
    out <- c(markers_walchli, markers_hpa[setdiff(names(markers_hpa), names(markers_walchli))])
    out <- out[!vapply(out, is.null, logical(1L))]
    out <- out[vapply(out, length, integer(1L)) > 0L]
    if (verbose) {
        total <- sum(vapply(out, length, integer(1L)))
        cat(sprintf("  load_celltype_markers: %d cell types, %d total markers (method=%s)\n",
                    length(out), total, method))
    }
    out
}

# Legacy pathway kept only for side-by-side reproducibility of pre-refactor runs.
.load_celltype_markers_legacy <- function(walchli_csv, hpa_file,
                                           walchli_pctile = 0.65,
                                           walchli_min_lp = 1.0,
                                           hpa_ntpm_thr   = 10,
                                           verbose        = TRUE) {
    WALCHLI_VASC_COLS <- .WALCHLI_COL2KEY
    WALCHLI_ALL_COLS  <- c("Endothelial.cells", "Pericytes", "Smooth.muscle.cells",
                           "Fibroblasts", "Oligodendrocytes", "Microglia.and.Macrophages",
                           "Astrocytes", "Neurons")
    marker_sets <- list()
    if (file.exists(walchli_csv)) {
        w <- data.table::fread(walchli_csv, showProgress = FALSE)
        data.table::setnames(w, "Gene", "gene_name")
        w[, gene_name := toupper(gene_name)]
        avail_all <- intersect(WALCHLI_ALL_COLS, names(w))
        log_mat   <- log1p(as.matrix(w[, ..avail_all]))
        row_med   <- apply(log_mat, 1, median, na.rm = TRUE)
        for (col in names(WALCHLI_VASC_COLS)) {
            if (!col %in% names(w)) next
            ct_lbl  <- WALCHLI_VASC_COLS[[col]]
            lp_col  <- log1p(w[[col]])
            pct_thr <- quantile(lp_col, walchli_pctile, na.rm = TRUE)
            keep    <- lp_col >= max(pct_thr, walchli_min_lp) & lp_col >= row_med
            marker_sets[[ct_lbl]] <- w$gene_name[keep]
        }
    }
    hpa_tsv <- sub("\\.zip$", ".tsv", hpa_file)
    if (file.exists(hpa_tsv)) {
        hpa <- data.table::fread(hpa_tsv, showProgress = FALSE)
    } else {
        hpa <- data.table::fread(cmd = paste("unzip -p", shQuote(hpa_file)),
                                  showProgress = FALSE)
    }
    data.table::setnames(hpa,
        c("Gene", "Gene name", "Cell type", "nTPM"),
        c("ensembl", "gene_name", "cell_type", "nTPM"),
        skip_absent = TRUE
    )
    cts <- intersect(c(.HPA_CT_CNS, .HPA_CT_IMMUNE), unique(hpa$cell_type))
    for (ct in cts) {
        marker_sets[[ct]] <- toupper(hpa[cell_type == ct & nTPM >= hpa_ntpm_thr,
                                          unique(gene_name)])
    }
    if (verbose) {
        cat(sprintf("  load_celltype_markers (legacy): %d cell types\n",
                    length(marker_sets)))
    }
    marker_sets
}

# ── run_celltype_gsea() ──────────────────────────────────────────────────────
# Preranked fgsea of cell-type marker sets against a signed ranking statistic
# (typically limma's moderated t). Returns a long data.table.
#
# Args:
#   stats        Named numeric vector: names = UPPERCASE gene symbols; values
#                = signed ranking statistic. Duplicate names are averaged with
#                a warning; NAs are dropped.
#   marker_sets  Named list from load_celltype_markers() (UPPERCASE symbols).
#   min_size     Drop cell-type sets smaller than this (after universe intersect).
#   max_size     Drop cell-type sets larger than this.
#   n_perm       fgsea eps=0 permutation budget per test (default 10000 via eps).
#   seed         RNG seed for reproducibility.
#
# Returns data.table with columns:
#   cell_type, size, NES, ES, pval, padj, leading_edge (comma-joined symbols)
run_celltype_gsea <- function(stats, marker_sets,
                               min_size = 10, max_size = 500,
                               eps       = 0,
                               seed      = 42) {
    if (!is.numeric(stats)) stop("`stats` must be a named numeric vector")
    if (is.null(names(stats))) stop("`stats` must be named with gene symbols")
    if (length(marker_sets) == 0L) return(NULL)
    names(stats) <- toupper(names(stats))

    # Handle duplicate gene symbols by averaging (fgsea requires unique names).
    if (anyDuplicated(names(stats))) {
        dt <- data.table::data.table(g = names(stats), s = stats)
        dt <- dt[!is.na(s), .(s = mean(s)), by = g]
        stats <- setNames(dt$s, dt$g)
    } else {
        stats <- stats[!is.na(stats)]
    }

    # Intersect marker sets with universe, drop sets outside size bounds.
    marker_sets_u <- lapply(marker_sets, function(m) intersect(toupper(m), names(stats)))
    set_sizes     <- vapply(marker_sets_u, length, integer(1L))
    keep          <- set_sizes >= min_size & set_sizes <= max_size
    if (!any(keep)) {
        warning(sprintf("run_celltype_gsea: no marker set within size bounds [%d, %d]",
                        min_size, max_size))
        return(NULL)
    }
    marker_sets_u <- marker_sets_u[keep]

    set.seed(seed)
    res <- fgsea::fgsea(
        pathways = marker_sets_u,
        stats    = stats,
        minSize  = min_size,
        maxSize  = max_size,
        eps      = eps
    )
    if (is.null(res) || nrow(res) == 0L) return(NULL)
    res <- as.data.table(res)
    data.table::setnames(res, "pathway", "cell_type")
    # Fold leadingEdge list-col into a comma-joined string for CSV portability.
    res[, leading_edge := vapply(leadingEdge, function(x) paste(x, collapse = ","),
                                  character(1L))]
    res[, leadingEdge := NULL]
    data.table::setcolorder(res, c("cell_type", "size", "NES", "ES",
                                     "pval", "padj", "leading_edge"))
    res[]
}

# ── run_celltype_camerapr() ─────────────────────────────────────────────────
# Optional confirmatory test using limma::cameraPR, which accounts for
# inter-gene correlation within each marker set. Returns a data.table
# with columns: cell_type, NGenes, direction, pval, padj.
run_celltype_camerapr <- function(stats, marker_sets,
                                    inter_gene_cor = 0.01,
                                    min_size = 10, max_size = 500) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        warning("limma not installed - skipping cameraPR")
        return(NULL)
    }
    if (is.null(names(stats))) stop("`stats` must be named with gene symbols")
    names(stats) <- toupper(names(stats))
    stats <- stats[!is.na(stats)]
    marker_sets_u <- lapply(marker_sets, function(m) intersect(toupper(m), names(stats)))
    set_sizes <- vapply(marker_sets_u, length, integer(1L))
    keep <- set_sizes >= min_size & set_sizes <= max_size
    if (!any(keep)) return(NULL)
    marker_sets_u <- marker_sets_u[keep]
    idx <- lapply(marker_sets_u, function(gs) match(gs, names(stats)))
    cam <- limma::cameraPR(statistic = stats, index = idx,
                            inter.gene.cor = inter_gene_cor, use.ranks = FALSE)
    out <- data.table::as.data.table(cam, keep.rownames = "cell_type")
    data.table::setnames(out, c("PValue", "FDR", "Direction"),
                          c("pval", "padj", "direction"), skip_absent = TRUE)
    out[]
}

# ── build_rank_stat() ───────────────────────────────────────────────────────
# Convenience: convert a limma topTable-style data.table into the named-vector
# stats required by run_celltype_gsea(). Uses moderated t by default (which
# combines effect size and precision); falls back to signed -log10(P.Value) if
# no t column is present.
build_rank_stat <- function(dt, gene_col = "protein",
                              stat_col = c("t", "signed_logp")) {
    stat_col <- match.arg(stat_col)
    if (!gene_col %in% names(dt)) stop("`gene_col` not found: ", gene_col)
    genes <- toupper(dt[[gene_col]])
    if (stat_col == "t" && "t" %in% names(dt)) {
        s <- as.numeric(dt$t)
    } else if ("logFC" %in% names(dt) && "P.Value" %in% names(dt)) {
        s <- sign(as.numeric(dt$logFC)) * -log10(pmax(as.numeric(dt$P.Value), 1e-300))
    } else {
        stop("Need either a `t` column or logFC+P.Value columns to build rank stat")
    }
    setNames(s, genes)
}

# ── run_ora() ─────────────────────────────────────────────────────────────────
# LEGACY: one-sided Fisher ORA of module_genes against each set in gs_list.
# Kept unchanged for backwards compatibility and side-by-side comparison with
# the new GSEA pipeline. New code should use run_celltype_gsea() instead.
run_ora <- function(module_genes, bg_genes, gs_list) {
    mod_u <- intersect(toupper(module_genes), toupper(bg_genes))
    bg_n  <- length(unique(toupper(bg_genes)))
    m_n   <- length(mod_u)
    if (m_n == 0L) return(NULL)

    res <- data.table::rbindlist(lapply(names(gs_list), function(ct) {
        markers <- intersect(toupper(gs_list[[ct]]), toupper(bg_genes))
        k <- length(markers)
        if (k < 3L) return(NULL)
        x   <- length(intersect(mod_u, markers))
        tab <- matrix(c(x, m_n - x, k - x, bg_n - m_n - k + x), nrow = 2)
        ft1 <- fisher.test(tab, alternative = "greater")
        ft2 <- fisher.test(tab, alternative = "two.sided")
        data.table::data.table(
            cell_type     = ct,
            n_universe    = bg_n,
            n_markers     = k,
            n_DEPs        = m_n,
            n_DEP_markers = x,
            odds_ratio    = ft1$estimate,
            p_overrep     = ft1$p.value,
            p_twosided    = ft2$p.value,
            overlap_genes = paste(sort(intersect(mod_u, markers)), collapse = ",")
        )
    }))
    if (is.null(res) || nrow(res) == 0L) return(NULL)
    res[, p_adj := p.adjust(p_twosided, method = "BH")]
    res[order(p_overrep)]
}

# ── make_celltype_gsea_plot() ────────────────────────────────────────────────
# Signed-NES dot plot. One column per contrast; dot colour = NES (red = up,
# blue = down in that contrast); dot size = -log10(padj); asterisks mark FDR.
#
# Args:
#   gsea_dt        Long data.table; cols required: cell_type, contrast (or
#                  analysis), NES, pval, padj. Optional: size (leading edge size).
#   contrast_col   Column to use for x-axis. Defaults to "contrast"; will
#                  fall back to "analysis" if "contrast" absent.
#   contrast_order Optional character vector giving x-axis order.
#   contrast_labels Optional named character vector mapping raw contrast values
#                  to display labels.
#   high_col / low_col   NES extreme colours.
#   title_str / subtitle_str   Plot labels.
make_celltype_gsea_plot <- function(
    gsea_dt,
    contrast_col     = NULL,
    contrast_order   = NULL,
    contrast_labels  = NULL,
    title_str        = "",
    subtitle_str     = paste0(
        "Preranked fgsea on limma moderated-t | ",
        "NES>0: markers higher in cases vs HC | ",
        "NES<0: lower | * FDR<0.05"
    ),
    high_col         = "#CC0066",
    low_col          = "#56B4E9",
    nes_clip         = 3,
    size_max         = 4
) {
    dt <- data.table::copy(gsea_dt)
    if (is.null(contrast_col)) {
        contrast_col <- if ("contrast" %in% names(dt)) "contrast" else "analysis"
    }
    if (!contrast_col %in% names(dt))
        stop("gsea_dt missing contrast column: ", contrast_col)
    req <- c("cell_type", "NES", "pval", "padj", contrast_col)
    missing <- setdiff(req, names(dt))
    if (length(missing))
        stop("gsea_dt missing columns: ", paste(missing, collapse = ", "))

    # Rename contrast column to the constant `ctx` to simplify aes().
    data.table::setnames(dt, contrast_col, "ctx")

    # Resolve display labels.
    dt[, ct_label := ifelse(cell_type %in% names(CT_LABELS),
                             CT_LABELS[cell_type], cell_type)]
    ct_order_actual <- CT_ORDER[CT_ORDER %in% dt$ct_label]
    if (length(ct_order_actual) == 0L) {
        message("make_celltype_gsea_plot: no recognised cell types - returning NULL")
        return(NULL)
    }
    dt[, ct_f := factor(ct_label, levels = rev(ct_order_actual))]

    if (!is.null(contrast_labels)) {
        dt[, ctx_label := ifelse(ctx %in% names(contrast_labels),
                                  contrast_labels[ctx], ctx)]
    } else {
        dt[, ctx_label := ctx]
    }
    if (is.null(contrast_order)) {
        contrast_order <- unique(dt$ctx)
    }
    lbl_order <- if (!is.null(contrast_labels))
                    contrast_labels[contrast_order] else contrast_order
    lbl_order <- unique(lbl_order)
    dt[, ctx_f := factor(ctx_label, levels = lbl_order)]

    dt[, NES_clip := pmax(pmin(NES, nes_clip), -nes_clip)]
    dt[, neglog_padj := pmin(-log10(pmax(padj, 1e-12)), size_max)]
    dt[, sig_label := data.table::fcase(
        padj < 0.001, "***",
        padj < 0.01,  "**",
        padj < 0.05,  "*",
        default = ""
    )]
    dt[, compartment := data.table::fcase(
        ct_label %in% CT_VASC,   "Vascular",
        ct_label %in% CT_CNS,    "CNS",
        ct_label %in% CT_IMMUNE, "Immune",
        default = "Other"
    )]
    dt[, compartment := factor(compartment, levels = c("Vascular", "CNS", "Immune"))]

    ggplot2::ggplot(dt, ggplot2::aes(x = ctx_f, y = ct_f,
                                      size = neglog_padj, colour = NES_clip)) +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::geom_text(ggplot2::aes(label = sig_label),
                            size = 3, colour = "black",
                            vjust = 0.5, hjust = 0.5, show.legend = FALSE) +
        ggplot2::facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_colour_gradient2(
            low = low_col, mid = "grey92", high = high_col,
            midpoint = 0,
            name   = "NES",
            limits = c(-nes_clip, nes_clip), oob = scales::squish,
            guide  = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 3.5)
        ) +
        ggplot2::scale_size_continuous(
            name   = expression(-log[10] ~ italic(FDR)),
            range  = c(2, 8),
            limits = c(0, size_max)
        ) +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.8)) +
        ggplot2::labs(title = title_str, subtitle = subtitle_str,
                       x = NULL, y = NULL) +
        theme_ukb(base_size = 9) +
        ggplot2::theme(
            axis.text.x      = ggplot2::element_text(size = 8.5, face = "bold"),
            axis.text.y      = ggplot2::element_text(size = 8),
            strip.text       = ggplot2::element_text(size = 8, face = "bold"),
            panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
            legend.position  = "right"
        )
}

# ── make_celltype_overrep_plot() ─────────────────────────────────────────────
# LEGACY: old-style up/down-split Fisher dot plot. Kept so supplementary
# "legacy pipeline" panels can still be regenerated. For new panels use
# make_celltype_gsea_plot() above.
make_celltype_overrep_plot <- function(
    overrep_dt,
    title_str    = "",
    subtitle_str = "Fisher's exact (legacy) | Walchli + HPA markers | * p_adj<0.05",
    high_col,
    low_col      = COL_DOWN,
    up_label     = "Up-regulated DEPs",
    down_label   = "Down-regulated DEPs"
) {
    dt <- data.table::copy(overrep_dt)

    dt[, ct_label := ifelse(cell_type %in% names(CT_LABELS),
                             CT_LABELS[cell_type], cell_type)]
    ct_order_actual <- CT_ORDER[CT_ORDER %in% dt$ct_label]
    if (length(ct_order_actual) == 0L) {
        message("make_celltype_overrep_plot: no recognised cell types - returning NULL")
        return(NULL)
    }

    dt[, ct_f      := factor(ct_label, levels = rev(ct_order_actual))]
    dt[, dir_label := ifelse(direction == "up", up_label, down_label)]
    dt[, dir_f     := factor(dir_label, levels = c(up_label, down_label))]
    dt[, log2_OR_clamp := pmax(pmin(log2(pmax(odds_ratio, 0.01)), 3), -3)]
    dt[, compartment := data.table::fcase(
        ct_label %in% CT_VASC,   "Vascular",
        ct_label %in% CT_CNS,    "CNS",
        ct_label %in% CT_IMMUNE, "Immune",
        default = "Other"
    )]
    dt[, compartment := factor(compartment, levels = c("Vascular", "CNS", "Immune"))]

    ggplot2::ggplot(dt, ggplot2::aes(x = dir_f, y = ct_f,
                                      size   = pmin(-log10(p_twosided), 4),
                                      colour = log2_OR_clamp)) +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_colour_gradient2(
            low = low_col, mid = "grey92", high = high_col,
            midpoint = 0,
            name   = expression(log[2] ~ OR),
            limits = c(-3, 3), oob = scales::squish,
            guide  = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 3.5)
        ) +
        ggplot2::scale_size_continuous(
            name   = expression(-log[10] ~ italic(p)),
            range  = c(2, 8), breaks = c(1, 2, 3, 4),
            labels = c("1", "2", "3", "\u22654")
        ) +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.8)) +
        ggplot2::labs(title = title_str, subtitle = subtitle_str,
                       x = NULL, y = NULL) +
        theme_ukb(base_size = 9) +
        ggplot2::theme(
            axis.text.x      = ggplot2::element_text(size = 8.5, face = "bold"),
            axis.text.y      = ggplot2::element_text(size = 8),
            strip.text       = ggplot2::element_text(size = 8, face = "bold"),
            panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
            legend.position  = "right"
        )
}

# ── make_celltype_heatmap() ──────────────────────────────────────────────────
# Multi-contrast cell-type heatmap. Supports both the new GSEA format
# (NES + padj columns; one row per cluster x cell_type) and the legacy Fisher
# format (direction == "up" rows with odds_ratio + p_adj). Back-compatible.
make_celltype_heatmap <- function(
    ora_all,
    fdr_thr      = 0.20,
    title_str    = "",
    subtitle_str = NULL,
    cluster_col  = "cluster",
    value_limits = c(-3, 3)
) {
    dt <- data.table::copy(ora_all)
    if (!cluster_col %in% names(dt))
        stop("make_celltype_heatmap: missing column ", cluster_col)

    dt[, ct_label := ifelse(cell_type %in% names(CT_LABELS),
                             CT_LABELS[cell_type], cell_type)]
    dt <- dt[ct_label %in% CT_ORDER]
    if (nrow(dt) == 0L) return(NULL)

    # Detect format: new (NES) vs legacy (odds_ratio + direction).
    is_new <- "NES" %in% names(dt) && "padj" %in% names(dt)
    if (is_new) {
        heat_dt <- dt[, .(cluster = get(cluster_col), ct_label,
                           value = NES,
                           sig   = padj < fdr_thr)]
        value_lab <- "NES"
    } else {
        dt <- dt[direction == "up"]
        heat_dt <- dt[, .(cluster = get(cluster_col), ct_label,
                           value = log2(pmax(odds_ratio, 1e-3)),
                           sig   = p_adj < fdr_thr)]
        value_lab <- expression(log[2](OR))
    }
    heat_dt[, ct_f := factor(ct_label, levels = rev(CT_ORDER))]
    heat_dt[, compartment := data.table::fcase(
        ct_label %in% CT_VASC,   "Vascular",
        ct_label %in% CT_CNS,    "CNS",
        ct_label %in% CT_IMMUNE, "Immune",
        default = "Other"
    )]
    heat_dt[, compartment := factor(compartment, levels = c("Vascular", "CNS", "Immune"))]

    high_col <- if (exists("CLUST_COLS")) CLUST_COLS["C1"] else "#CC0066"
    low_col  <- if (exists("CLUST_COLS")) CLUST_COLS["C2"] else "#56B4E9"

    ggplot2::ggplot(heat_dt[!is.na(value)],
                     ggplot2::aes(x = cluster, y = ct_f, fill = value)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::geom_text(data = heat_dt[sig == TRUE], ggplot2::aes(label = "*"),
                            colour = "black", size = 3, vjust = 0.65) +
        ggplot2::facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
        ggplot2::scale_fill_gradient2(
            low      = low_col, mid = "grey96", high = high_col,
            midpoint = 0, limits = value_limits, oob = scales::squish,
            na.value = "grey92", name = value_lab
        ) +
        ggplot2::labs(x = NULL, y = NULL, title = title_str, subtitle = subtitle_str) +
        theme_ukb(base_size = 9) +
        ggplot2::theme(
            axis.text.x       = ggplot2::element_text(face = "bold", size = 8, lineheight = 1.1),
            axis.text.y       = ggplot2::element_text(size = 8),
            strip.text        = ggplot2::element_text(size = 8, face = "bold"),
            legend.key.height = ggplot2::unit(0.5, "cm")
        )
}
