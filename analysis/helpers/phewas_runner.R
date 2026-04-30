#!/usr/bin/env Rscript
# phewas_runner.R — Reusable UKB PheWAS engine (PEACOK ICD-10 tree)
#
# Adapts PEACOK's ICD-10 hierarchy propagation (datacode-19.tsv) for use
# with our per-code icd_codes/ directory.  Any trait can be tested against
# all available ICD-10 phenotypes in one call.
#
# USAGE ──────────────────────────────────────────────────────────────────────
#
#   source(file.path(PROJ_DIR, "analysis", "helpers", "phewas_runner.R"))
#
#   # Protein PheWAS (continuous trait → binary ICD outcomes)
#   res <- run_phewas(
#       trait   = olink[, .(eid, trait = GFAP)],
#       cohort  = cohort[, .(eid, age, sex)],
#       icd_dir = ICD_DIR
#   )
#
#   # Cluster membership PheWAS (binary trait)
#   res <- run_phewas(
#       trait       = clust[, .(eid, trait = as.integer(cluster == 3))],
#       cohort      = cohort[, .(eid, age, sex)],
#       trait_label = "Cluster-3 membership"
#   )
#
#   # PRS PheWAS with genetic PCs as covariates
#   res <- run_phewas(
#       trait   = prs[, .(eid, trait = prs_score)],
#       cohort  = cohort[, .(eid, age, sex, PC1, PC2, PC3, PC4, PC5)],
#       covariates = c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")
#   )
#
#   # Visualise
#   p <- plot_phewas(res, title = "MS PRS PheWAS")
#   ggsave("phewas_prs.pdf", p, width = 12, height = 5)
#
# ARGUMENTS (run_phewas) ─────────────────────────────────────────────────────
#
#   trait         data.table(eid, trait)  — any numeric/binary vector with eids
#   cohort        data.table with eid + covariate columns (at minimum age, sex)
#   icd_dir       path to icd_codes/ root directory  [auto-detected]
#   covariates    character vector of covariate column names in cohort
#                 [default: c("age","sex")]
#   min_cases     integer  — minimum cases per ICD phenotype tested  [10]
#   fdr_method    "BH" | "bonferroni"  [default "BH"]
#   propagate     logical  — propagate diagnoses up ICD-10 tree via PEACOK
#                 datacode-19.tsv (adds parent-level phenotypes)  [TRUE]
#   timing_col    name of cohort column holding reference date (e.g. age at
#                 MS diagnosis) for pre/post timing.  Set NULL to skip.  [NULL]
#   verbose       logical  [TRUE]
#
# RETURNS ────────────────────────────────────────────────────────────────────
#
#   data.table with columns:
#     icd10, meaning, chapter, level (leaf/parent), n_cases, pct_cases,
#     OR (or beta for continuous), OR_lo, OR_hi, pval, fdr,
#     pct_pre_timing (if timing_col supplied)
#
# REQUIRES:  data.table, ggplot2, ggrepel

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
})

# ── Internal: default paths (can be overridden by caller before sourcing) ────
if (!exists(".phewas_env") || !is.environment(.phewas_env)) {
    .phewas_env <- new.env(parent = emptyenv())
}
if (!exists("PEACOK_CODES", envir = .phewas_env, inherits = FALSE)) {
    # Resolve relative to this file when sourced via source("path/to/phewas_runner.R")
    .phewas_this_file <- tryCatch(
        normalizePath(sys.frames()[[length(sys.frames())]]$ofile),
        error = function(e) NULL
    )
    if (is.null(.phewas_this_file) || !nzchar(.phewas_this_file)) {
        # Fallback: assume CWD is project root
        .phewas_helpers <- "analysis/helpers"
    } else {
        .phewas_helpers <- dirname(normalizePath(.phewas_this_file))
    }
    .phewas_proj    <- normalizePath(file.path(.phewas_helpers, "..", ".."))
    .phewas_cadasil <- normalizePath(
        file.path(.phewas_proj, "..", "CADASIL_Proteome_ML_Keller_2024_Rebuttal"),
        mustWork = FALSE
    )
    .phewas_env$PEACOK_CODES      <- file.path(.phewas_helpers, "PEACOK",
                                                "variable-info", "data-codes")
    .phewas_env$ICD_DIR_DEFAULT   <- file.path(.phewas_cadasil, "data", "ukb",
                                                "diagnoses", "icd_codes")
    rm(.phewas_this_file, .phewas_helpers, .phewas_proj, .phewas_cadasil)
}

# ── ICD-10 tree (PEACOK datacode-19.tsv) ─────────────────────────────────────

.build_icd_tree <- function(peacok_codes_dir = .phewas_env$PEACOK_CODES) {
    f <- file.path(peacok_codes_dir, "datacode-19.tsv")
    if (!file.exists(f)) stop("PEACOK datacode-19.tsv not found at: ", f)

    tree <- fread(f, sep = "\t", header = TRUE, na.strings = c("","NA"),
                  showProgress = FALSE)
    # Columns: coding, meaning, node_id, parent_id, selectable

    # Build level (depth) via BFS from roots (parent_id == 0)
    tree[, level    := NA_integer_]
    tree[, root_nid := NA_integer_]
    tree[parent_id == 0L, level    := 0L]
    tree[parent_id == 0L, root_nid := node_id]

    # node_id → level/root lookup table (indexed by node_id)
    for (iter in seq_len(20L)) {
        # Find nodes whose parent is already assigned
        assigned <- tree[!is.na(level), .(node_id, level, root_nid)]
        # Join: for each unassigned node, look up its parent
        unassigned <- tree[is.na(level)]
        if (nrow(unassigned) == 0L) break
        merged <- unassigned[assigned, on = .(parent_id = node_id), nomatch = 0L,
                              .(row_node = x.node_id,
                                new_level = i.level + 1L,
                                new_root  = i.root_nid)]
        if (nrow(merged) == 0L) break
        # Assign via row_node (suppressed: harmless recycling when node has multiple children)
        suppressWarnings(
            tree[merged, on = .(node_id = row_node), `:=`(level = i.new_level, root_nid = i.new_root)]
        )
    }

    # icd_clean: remove dots (A00.1 → A001) for matching against our file names
    tree[, icd_clean := gsub("\\.", "", coding)]
    tree[]
}

# Build once and cache
.get_icd_tree <- local({
    cache <- NULL
    function(peacok_codes_dir = .phewas_env$PEACOK_CODES) {
        if (is.null(cache)) cache <<- .build_icd_tree(peacok_codes_dir)
        cache
    }
})

# ── Internal: load all per-code files and propagate up the tree ───────────────

.load_icd_phenotypes <- function(icd_dir, eids, tree, propagate = TRUE, verbose = TRUE) {
    all_files <- list.files(icd_dir, pattern = "^[A-Z][0-9]+\\.csv$",
                             recursive = TRUE, full.names = TRUE)
    if (verbose) cat(sprintf("  Loading %d ICD code files...\n", length(all_files)))

    # Load all present cases  →  data.table(eid, icd10, age_at_event)
    load_one <- function(fp) {
        code <- tools::file_path_sans_ext(basename(fp))
        dt <- tryCatch(fread(fp, showProgress = FALSE), error = function(e) NULL)
        if (is.null(dt) || nrow(dt) == 0L) return(NULL)
        if (ncol(dt) < 2L) return(NULL)
        names(dt)[1:2] <- c("eid","present")
        if (ncol(dt) >= 3L) names(dt)[3] <- "age_at_event"
        else dt[, age_at_event := NA_real_]
        dt[present == 1L & eid %in% eids,
           .(eid = as.integer(eid),
             icd10 = code,
             age_at_event = as.numeric(age_at_event))]
    }

    events <- rbindlist(Filter(Negate(is.null), lapply(all_files, load_one)),
                        fill = TRUE)
    if (verbose) cat(sprintf("  Loaded %d events across %d codes\n",
                              nrow(events), uniqueN(events$icd10)))

    if (!propagate || nrow(events) == 0L) return(events)

    # Propagate up: for each leaf code, add all ancestor codes for the same eid
    if (verbose) cat("  Propagating diagnoses up ICD-10 tree...\n")

    # Build child→ancestors map from tree
    parent_map <- tree[!is.na(parent_id) & parent_id != 0L,
                       .(icd_clean, parent_icd = tree[match(parent_id, node_id), icd_clean])]
    parent_map <- parent_map[!is.na(parent_icd)]

    # Walk up the tree for each observed leaf code
    leaf_codes   <- unique(events$icd10)
    ancestor_map <- list()  # leaf_code → vector of ancestor codes

    for (lc in leaf_codes) {
        current <- lc
        ancestors <- character(0)
        for (depth in seq_len(8L)) {
            par <- parent_map[icd_clean == current, parent_icd]
            if (length(par) == 0L || is.na(par)) break
            ancestors <- c(ancestors, par)
            current   <- par
        }
        if (length(ancestors) > 0L) ancestor_map[[lc]] <- ancestors
    }

    # Expand events table with ancestor entries
    extras <- rbindlist(lapply(names(ancestor_map), function(lc) {
        ancs <- ancestor_map[[lc]]
        base <- events[icd10 == lc, .(eid, age_at_event)]
        rbindlist(lapply(ancs, function(a) {
            base[, .(eid, icd10 = a, age_at_event)]
        }))
    }), fill = TRUE)

    out <- rbindlist(list(events, extras), fill = TRUE)
    out <- out[!duplicated(out[, .(eid, icd10)])]  # one row per eid-code
    if (verbose) {
        n_codes_after <- uniqueN(out$icd10)
        cat(sprintf("  After propagation: %d codes (%d leaf + %d parent)\n",
                    n_codes_after, uniqueN(events$icd10),
                    n_codes_after - uniqueN(events$icd10)))
    }
    out[]
}

# ── Main: run_phewas ──────────────────────────────────────────────────────────

run_phewas <- function(
    trait,
    cohort,
    icd_dir      = .phewas_env$ICD_DIR_DEFAULT,
    covariates   = c("age", "sex"),
    min_cases    = 10L,
    fdr_method   = "BH",
    propagate    = TRUE,
    timing_col   = NULL,
    peacok_dir   = .phewas_env$PEACOK_CODES,
    verbose      = TRUE
) {
    stopifnot(is.data.table(trait) || is.data.frame(trait))
    stopifnot(is.data.table(cohort) || is.data.frame(cohort))
    setDT(trait); setDT(cohort)

    if (!all(c("eid","trait") %in% names(trait)))
        stop("trait must be a data.table with columns eid and trait")
    if (!"eid" %in% names(cohort))
        stop("cohort must contain column eid")

    miss_cov <- setdiff(covariates, names(cohort))
    if (length(miss_cov) > 0L)
        stop("covariates missing from cohort: ", paste(miss_cov, collapse=", "))

    # Merge trait + cohort
    keep_cols <- unique(c("eid", covariates, timing_col))
    ana <- merge(trait, cohort[, ..keep_cols], by = "eid")
    ana <- ana[!is.na(trait)]
    for (v in covariates) ana <- ana[!is.na(get(v))]
    if (verbose) cat(sprintf("Analytic sample: %d individuals\n", nrow(ana)))

    # Determine trait type
    is_binary_trait <- length(unique(na.omit(ana$trait))) == 2L
    if (verbose) cat(sprintf("Trait type: %s\n",
                              if (is_binary_trait) "binary" else "continuous"))

    # Load ICD-10 tree
    tree <- .get_icd_tree(peacok_dir)

    # Load + propagate ICD phenotypes
    events <- .load_icd_phenotypes(icd_dir, ana$eid, tree,
                                    propagate = propagate, verbose = verbose)
    if (nrow(events) == 0L) stop("No ICD-10 events found for the analytic sample")

    all_codes <- unique(events$icd10)
    N_ANA     <- nrow(ana)
    n_total   <- 0L

    if (verbose) cat(sprintf("Testing %d ICD-10 phenotypes (min cases = %d)...\n",
                              length(all_codes), min_cases))

    test_one <- function(code) {
        ev_sub <- events[icd10 == code]
        n_cases <- sum(ana$eid %in% ev_sub$eid)
        if (n_cases < min_cases || (N_ANA - n_cases) < min_cases) return(NULL)

        dt <- copy(ana)
        dt[, present := as.integer(eid %in% ev_sub$eid)]

        formula_str <- paste("present ~ trait +",
                              paste(covariates, collapse = " + "))
        fit <- tryCatch(
            suppressWarnings(glm(as.formula(formula_str), data = dt,
                                 family = binomial)),
            error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
        if (is.null(cf) || !"trait" %in% rownames(cf)) return(NULL)

        beta <- cf["trait","Estimate"]
        se   <- cf["trait","Std. Error"]
        pval <- cf["trait","Pr(>|z|)"]

        # Timing: % with event before reference date (timing_col)
        pct_pre <- NA_real_
        if (!is.null(timing_col) && timing_col %in% names(ana)) {
            ms <- ev_sub[eid %in% ana[!is.na(get(timing_col)), eid]]
            if (nrow(ms) >= 5L && any(!is.na(ms$age_at_event))) {
                tref <- ana[, .(eid, t_ref = get(timing_col))]
                tdt  <- merge(ms, tref, by = "eid")
                tdt[, delta := age_at_event - t_ref]
                pct_pre <- 100 * mean(tdt$delta < -1, na.rm = TRUE)
            }
        }

        # ICD meaning + level from tree
        tree_row <- tree[icd_clean == code]
        meaning  <- if (nrow(tree_row) > 0L) tree_row$meaning[1L] else code
        # Strip ICD code prefix from meaning ("G35 Multiple sclerosis" → "Multiple sclerosis")
        meaning <- sub("^[A-Z][0-9]+(?:\\.[0-9]+)?\\s+", "", meaning)
        lvl      <- if (nrow(tree_row) > 0L) tree_row$level[1L] else NA_integer_

        data.table(
            icd10        = code,
            meaning      = meaning,
            chapter      = substr(code, 1, 1),
            tree_level   = lvl,
            is_leaf      = code %in% all_codes[nchar(all_codes) >= 3L],
            n_cases      = n_cases,
            pct_cases    = 100 * n_cases / N_ANA,
            beta         = beta,
            se           = se,
            OR           = exp(beta),
            OR_lo        = exp(beta - 1.96 * se),
            OR_hi        = exp(beta + 1.96 * se),
            pval         = pval,
            pct_pre      = pct_pre
        )
    }

    results <- rbindlist(
        Filter(Negate(is.null), lapply(all_codes, test_one)),
        fill = TRUE
    )

    if (nrow(results) == 0L) {
        warning("No phenotypes passed min_cases filter")
        return(results)
    }

    results[, fdr := p.adjust(pval, method = fdr_method)]
    setorder(results, pval)

    if (verbose) {
        cat(sprintf("\nResults: %d phenotypes tested\n", nrow(results)))
        cat(sprintf("FDR < 0.05: %d phenotypes\n", sum(results$fdr < 0.05)))
        cat(sprintf("Nominal p < 0.05: %d phenotypes\n", sum(results$pval < 0.05)))
        cat("\nTop 10 associations:\n")
        print(results[1:min(10L, .N),
                      .(icd10, meaning = substr(meaning, 1, 35), OR = round(OR, 2),
                        pval = signif(pval, 3), fdr = round(fdr, 3))])
    }
    results[]
}

# ── Plot: plot_phewas ─────────────────────────────────────────────────────────

plot_phewas <- function(
    results,
    title       = "PheWAS",
    fdr_thr     = 0.05,
    n_label     = 20L,
    leaf_only   = FALSE,
    chap_order  = c("G","H","I","J","K","L","M","N","O","E","F","C","D","P","Q","A","B"),
    peacok_dir  = .phewas_env$PEACOK_CODES
) {
    stopifnot(is.data.table(results))

    dt <- copy(results)
    if (leaf_only) dt <- dt[is_leaf == TRUE]
    if (nrow(dt) == 0L) stop("No rows to plot")

    # ── Chapter colours ──
    CHAP_COLS <- c(
        A="#9E2A2B", B="#E05C5C",
        C="grey70",  D="#FF9933",
        E="#006D77", F="#83C5BE",
        G="#2B4C7E", H="#5B9BD5",
        I="#E6A817", J="#4CAF50",
        K="#A0522D", L="#CC0066",
        M="#9370DB", N="#20B2AA",
        O="#DC143C", P="#CDDC39", Q="#808080"
    )

    chap_present  <- chap_order[chap_order %in% dt$chapter]
    dt[, chapter  := factor(chapter, levels = chap_present)]
    dt            <- dt[!is.na(chapter)]
    setorder(dt, chapter, pval)

    # ── X positions with chapter gap ──
    dt[, x_pos := seq_len(.N), by = chapter]
    GAP <- 4L
    offset <- 0L
    for (ch in chap_present) {
        idx <- dt$chapter == ch
        dt[idx, x_pos := x_pos + offset]
        offset <- offset + max(dt[idx, x_pos]) + GAP
    }
    dt[, x_pos := as.numeric(x_pos)]

    # Chapter midpoints
    chap_mids <- dt[, .(x_mid = mean(x_pos)), by = chapter]
    CHAP_NAMES <- c(
        G="Neurological", H="Eye/Ear",    I="Cardiovascular", J="Respiratory",
        K="Digestive",    L="Skin",       M="Musculoskeletal", N="Genitourinary",
        O="Obstetric",    E="Metabolic",  F="Psychiatric",
        C="Neoplasms",    D="Blood",      P="Perinatal", Q="Congenital",
        A="Infections",   B="Viral"
    )
    chap_mids[, chap_name := CHAP_NAMES[as.character(chapter)]]

    # FDR threshold horizontal line
    fdr_y <- {
        ord <- dt[order(fdr)]
        first_ns <- which(ord$fdr >= fdr_thr)[1L]
        if (is.na(first_ns)) NA_real_ else -log10(ord$pval[first_ns])
    }

    dt[, log10p   := -log10(pval)]
    dt[, sig      := fdr < fdr_thr]
    dt[, fill_col := CHAP_COLS[as.character(chapter)]]
    dt[sig == FALSE, fill_col := "grey82"]
    dt[, pt_size  := pmin(abs(log2(pmax(OR, 0.01))), 4) * 0.6 + 0.3]

    # Labels: top n_label FDR-sig
    label_dt <- dt[sig == TRUE][order(pval)][seq_len(min(n_label, sum(dt$sig)))]
    label_dt[, lbl := paste0(icd10, " ", substr(meaning, 1, 28))]

    ggplot(dt, aes(x = x_pos, y = log10p)) +
        { if (!is.na(fdr_y))
            geom_hline(yintercept = fdr_y,
                       linewidth = 0.3, linetype = "dashed", colour = "grey50") } +
        geom_point(aes(size = pt_size, colour = fill_col),
                   alpha = ifelse(dt$sig, 0.9, 0.4), stroke = 0) +
        geom_text_repel(
            data      = label_dt,
            aes(label = lbl),
            size      = 2.1, segment.size = 0.25, max.overlaps = 30,
            box.padding = 0.3, colour = "grey10"
        ) +
        scale_colour_identity() +
        scale_size_identity() +
        scale_x_continuous(
            breaks = chap_mids$x_mid,
            labels = chap_mids$chap_name,
            expand = expansion(add = c(2, 2))
        ) +
        annotate("text",
                 x     = min(dt$x_pos),
                 y     = max(dt$log10p) * 0.97,
                 label = sprintf("n=%d tested | %d FDR<%.2f",
                                  nrow(dt), sum(dt$sig), fdr_thr),
                 size  = 2.5, hjust = 0, vjust = 1, colour = "grey30") +
        labs(
            x = "ICD-10 chapter",
            y = expression(-log[10](P)),
            title = title
        ) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            axis.text.x        = element_text(size = 8, angle = 35, hjust = 1),
            plot.title         = element_text(size = 11, face = "bold")
        )
}

# ── Convenience wrapper: save results + plot to a directory ──────────────────

save_phewas <- function(results, out_dir, label, w = 12, h = 5, fdr_thr = 0.05) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    fwrite(results, file.path(out_dir, paste0(label, "_phewas_results.csv")))
    p <- plot_phewas(results, title = label, fdr_thr = fdr_thr)
    ggsave(file.path(out_dir, paste0(label, "_phewas_manhattan.pdf")),
           p, width = w, height = h, device = cairo_pdf)
    invisible(results)
}
