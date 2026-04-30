# ============================================================================
# top_n_selector.R — programmatic top-N protein selection with FDR fallback
# ============================================================================
# Replaces hardcoded protein name lists in figure scripts.  Returns top-N up
# and top-N down DEPs from a limma/Cox/PheWAS results table, ranked by
# adj.P.Val when there are enough FDR-significant hits, else by nominal P.
#
# Usage:
#   source(here::here("analysis", "helpers", "top_n_selector.R"))
#   sel <- top_n_by_direction(dt, n_each = 20)
#   # sel$up    — character vector, top N up-regulated proteins (lowercase)
#   # sel$down  — character vector, top N down-regulated proteins
#   # sel$all   — c(up, down)
#   # attr(sel, "ordered_by") — "adj.P.Val" or "P.Value" depending on fallback
#
# Single-direction helper (e.g. Cox HR or PheWAS where logFC has no meaning):
#   sel <- top_n_one_direction(dt, n = 20, p_col = "P.Value")
#
# Behaviour:
#   - When >= min_sig_for_fdr proteins meet adj.P.Val < adj_threshold in a
#     given direction, the helper picks top-N from that significant subset
#     ordered by adj.P.Val (P.Value as tiebreaker).
#   - When the FDR-significant subset is too small (e.g. pre-onset MS with
#     only 4 FDR-sig down-regulated DEPs), it falls back to the full
#     direction subset ordered by nominal P.Value.  The fallback ensures
#     that under-powered contrasts still show their best-evidence proteins
#     instead of producing an empty plot.
# ============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

# Default-when-NULL operator (R 4.4+ ships ?: but base R does not).
# Idempotent so loading this helper multiple times is safe.
if (!exists("%||%", mode = "function")) {
    `%||%` <- function(x, y) if (is.null(x)) y else x
}

#' Pick top-N proteins per direction (up- and down-regulated separately).
#'
#' @param dt data.table or data.frame with columns:
#'        protein, logFC (or effect-size column), P.Value, adj.P.Val
#' @param n_each integer, number of proteins to pick per direction
#' @param adj_threshold FDR threshold for "significant" subset (default 0.05)
#' @param min_sig_for_fdr integer.  If fewer than this number of proteins meet
#'        adj.P.Val < adj_threshold in a given direction, the function falls
#'        back to nominal P.Value ordering on the full direction subset.
#'        Default 5.
#' @param effect_col column to use for direction sign (default "logFC")
#' @param protein_col column with protein identifiers (default "protein")
#' @param p_col nominal p-value column (default "P.Value")
#' @param adj_col adjusted p-value column (default "adj.P.Val")
#' @param case_fn function applied to protein names before return
#'        (default tolower so values match Olink lowercase column convention)
#'
#' @return list with elements `up`, `down`, `all`.  Attribute "ordered_by"
#'         records the actual ordering criterion per direction.
top_n_by_direction <- function(dt,
                                n_each = 20L,
                                adj_threshold = 0.05,
                                min_sig_for_fdr = 5L,
                                effect_col   = "logFC",
                                protein_col  = "protein",
                                p_col        = "P.Value",
                                adj_col      = "adj.P.Val",
                                case_fn      = tolower) {
    dt <- as.data.table(dt)
    needed <- c(protein_col, effect_col, p_col, adj_col)
    missing_cols <- setdiff(needed, names(dt))
    if (length(missing_cols) > 0) {
        stop(sprintf("top_n_by_direction: missing column(s) %s in input",
                     paste(missing_cols, collapse = ", ")))
    }

    pick_dir <- function(d, dir_name, sign_fn) {
        d_dir <- d[sign_fn(get(effect_col))]
        if (nrow(d_dir) == 0L) {
            return(list(picked = character(0), ordered_by = NA_character_))
        }
        sig_subset <- d_dir[get(adj_col) < adj_threshold]
        if (nrow(sig_subset) >= min_sig_for_fdr) {
            picked_dt <- sig_subset[order(get(adj_col), get(p_col))][
                seq_len(min(n_each, .N))]
            return(list(picked     = case_fn(picked_dt[[protein_col]]),
                        ordered_by = adj_col))
        }
        picked_dt <- d_dir[order(get(p_col))][seq_len(min(n_each, .N))]
        list(picked     = case_fn(picked_dt[[protein_col]]),
             ordered_by = p_col)
    }

    up_res   <- pick_dir(dt, "up",   function(x) x > 0)
    down_res <- pick_dir(dt, "down", function(x) x < 0)

    out <- list(up   = up_res$picked,
                down = down_res$picked,
                all  = c(up_res$picked, down_res$picked))
    attr(out, "ordered_by") <- c(up = up_res$ordered_by,
                                  down = down_res$ordered_by)
    attr(out, "n_up")   <- length(out$up)
    attr(out, "n_down") <- length(out$down)
    out
}

#' Pick top-N proteins by p-value with FDR-aware fallback (no direction split).
#' Useful for Cox PH (HR has no sign convention here), PheWAS by ICD code, etc.
top_n_one_direction <- function(dt,
                                 n = 20L,
                                 adj_threshold = 0.05,
                                 min_sig_for_fdr = 5L,
                                 protein_col = "protein",
                                 p_col       = "P.Value",
                                 adj_col     = "adj.P.Val",
                                 case_fn     = tolower) {
    dt <- as.data.table(dt)
    needed <- c(protein_col, p_col, adj_col)
    missing_cols <- setdiff(needed, names(dt))
    if (length(missing_cols) > 0) {
        stop(sprintf("top_n_one_direction: missing column(s) %s",
                     paste(missing_cols, collapse = ", ")))
    }
    sig_subset <- dt[get(adj_col) < adj_threshold]
    if (nrow(sig_subset) >= min_sig_for_fdr) {
        picked <- sig_subset[order(get(adj_col), get(p_col))][
            seq_len(min(n, .N))]
        attr_ord <- adj_col
    } else {
        picked <- dt[order(get(p_col))][seq_len(min(n, .N))]
        attr_ord <- p_col
    }
    out <- case_fn(picked[[protein_col]])
    attr(out, "ordered_by") <- attr_ord
    attr(out, "n") <- length(out)
    out
}

# Self-test --------------------------------------------------------------------
if (sys.nframe() == 0L) {
    set.seed(1)
    test_dt <- data.table::data.table(
        protein  = paste0("P", 1:30),
        logFC    = c(runif(15, 0, 2), runif(15, -2, 0)),
        P.Value  = sort(runif(30, 0, 0.5)),
        adj.P.Val = sort(runif(30, 0, 0.6))
    )
    sel <- top_n_by_direction(test_dt, n_each = 5L)
    cat("top_n_by_direction self-test:\n")
    cat("  up:   ",  paste(sel$up, collapse = ", "), "\n")
    cat("  down: ",  paste(sel$down, collapse = ", "), "\n")
    cat("  ordered_by:", paste(attr(sel, "ordered_by"), collapse = " / "), "\n")
}
