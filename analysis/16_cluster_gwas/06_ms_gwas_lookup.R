#!/usr/bin/env Rscript
# 06_ms_gwas_lookup.R
# Post-process cluster GWAS results and compare hits to MS GWAS catalog.
#
# Steps:
#   1. Load compiled REGENIE output for C1 and C2
#   2. Identify genome-wide significant loci (p < 5e-8), clump by 500kb windows
#   3. Download IMSGC 2023 MS GWAS summary stats from GWAS Catalog
#      (accession GCST90012877 = Brum et al. 2023 Nature Genetics)
#   4. Check overlap: do cluster GWAS hits co-localise with MS risk loci?
#   5. Generate Manhattan + QQ plots
#   6. Output tables: cluster_gwas_hits.csv, ms_overlap_C1.csv, ms_overlap_C2.csv
#
# Output: results/gwas_cluster/

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

GWAS_DIR  <- file.path(PROJ_DIR, "results", "gwas_cluster")
OUT_DIR   <- GWAS_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GW_THR    <- 5e-8
SUGGESTIVE <- 1e-6
WINDOW_KB  <- 500

# IMSGC 2023 MS GWAS Catalog accession (Brum et al. Nat Genet 2023)
MS_GWAS_URL <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012877/harmonised/GCST90012877_buildGRCh37.tsv.gz"
MS_GWAS_LOCAL <- file.path(PROJ_DIR, "data", "external", "ms_imsgc2023_gwas.tsv.gz")
dir.create(dirname(MS_GWAS_LOCAL), showWarnings=FALSE, recursive=TRUE)

# ── 1. Load REGENIE results ────────────────────────────────────────────────────
load_regenie <- function(gz_file, phenotype) {
    if (!file.exists(gz_file)) {
        warning(sprintf("GWAS results not found: %s", gz_file))
        return(NULL)
    }
    dt <- fread(gz_file, showProgress=FALSE)
    # REGENIE output columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P
    # Convert LOG10P to p-value
    dt[, p_value := 10^(-LOG10P)]
    dt[, phenotype := phenotype]
    cat(sprintf("  %s: %d variants loaded\n", phenotype, nrow(dt)))
    dt
}

cat("Loading GWAS results...\n")
gwas_c1 <- load_regenie(file.path(GWAS_DIR, "C1_gwas_all_chr.tsv.gz"), "C1")
gwas_c2 <- load_regenie(file.path(GWAS_DIR, "C2_gwas_all_chr.tsv.gz"), "C2")

if (is.null(gwas_c1) && is.null(gwas_c2)) {
    cat("No GWAS results found. Run 05_collect_results.sh first.\n")
    quit(status=0)
}

# ── 2. Identify significant loci ──────────────────────────────────────────────
clump_loci <- function(dt, thr = GW_THR, window_kb = WINDOW_KB) {
    if (is.null(dt)) return(NULL)
    sig <- dt[p_value < thr][order(CHROM, GENPOS)]
    if (nrow(sig) == 0) return(data.table())
    # Greedy window clumping
    sig[, locus := 0L]
    locus_id <- 1L
    prev_chr  <- -1L
    prev_pos  <- -Inf
    for (i in seq_len(nrow(sig))) {
        if (sig$CHROM[i] != prev_chr || (sig$GENPOS[i] - prev_pos) > window_kb * 1000) {
            locus_id <- locus_id + 1L
            prev_chr <- sig$CHROM[i]
            prev_pos <- sig$GENPOS[i]
        }
        sig$locus[i] <- locus_id
    }
    # Keep index SNP (lowest p) per locus
    sig[, .SD[which.min(p_value)], by=locus]
}

cat("\nIdentifying genome-wide significant loci...\n")
hits_c1 <- clump_loci(gwas_c1)
hits_c2 <- clump_loci(gwas_c2)
cat(sprintf("  C1 GW hits (p<5e-8): %d loci\n",
            if (is.null(hits_c1)) 0L else nrow(hits_c1)))
cat(sprintf("  C2 GW hits (p<5e-8): %d loci\n",
            if (is.null(hits_c2)) 0L else nrow(hits_c2)))

# ── 3. Download MS GWAS catalog ───────────────────────────────────────────────
cat("\nLoading IMSGC 2023 MS GWAS...\n")
if (!file.exists(MS_GWAS_LOCAL)) {
    cat("  Downloading IMSGC 2023 GWAS summary stats...\n")
    tryCatch(
        download.file(MS_GWAS_URL, MS_GWAS_LOCAL, mode="wb", quiet=TRUE),
        error = function(e) {
            cat(sprintf("  WARNING: Could not download MS GWAS: %s\n", e$message))
        }
    )
}

ms_gwas <- NULL
if (file.exists(MS_GWAS_LOCAL)) {
    ms_gwas <- fread(MS_GWAS_LOCAL, showProgress=FALSE)
    cat(sprintf("  MS GWAS loaded: %d variants\n", nrow(ms_gwas)))
    # Harmonise column names (GWAS Catalog format)
    if ("hm_chrom" %in% names(ms_gwas)) {
        setnames(ms_gwas, c("hm_chrom","hm_pos","p_value"),
                          c("CHR","POS","ms_pvalue"), skip_absent=TRUE)
    }
    ms_sig <- ms_gwas[ms_pvalue < 5e-8]
    cat(sprintf("  MS GW loci: %d variants\n", nrow(ms_sig)))
}

# ── 4. Overlap analysis ───────────────────────────────────────────────────────
check_ms_overlap <- function(hits, ms_sig, window_kb = 500) {
    if (is.null(hits) || nrow(hits) == 0 || is.null(ms_sig)) return(data.table())
    results <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
        chr <- hits$CHROM[i]
        pos <- hits$GENPOS[i]
        nearby_ms <- ms_sig[CHR == chr &
                             abs(POS - pos) <= window_kb * 1000]
        if (nrow(nearby_ms) == 0) {
            return(data.table(
                cluster_snp=hits$ID[i], CHROM=chr, GENPOS=pos,
                cluster_p=hits$p_value[i], ms_overlap=FALSE,
                nearest_ms_snp=NA, nearest_ms_dist=NA, nearest_ms_p=NA
            ))
        }
        best <- nearby_ms[which.min(ms_pvalue)]
        data.table(
            cluster_snp=hits$ID[i], CHROM=chr, GENPOS=pos,
            cluster_p=hits$p_value[i], ms_overlap=TRUE,
            nearest_ms_snp=best$hm_rsid[1], nearest_ms_dist=abs(best$POS[1]-pos),
            nearest_ms_p=best$ms_pvalue[1]
        )
    }))
    results
}

if (!is.null(ms_gwas)) {
    overlap_c1 <- check_ms_overlap(hits_c1, ms_sig)
    overlap_c2 <- check_ms_overlap(hits_c2, ms_sig)
    cat(sprintf("\n  C1 hits overlapping MS loci (±%dkb): %d / %d\n",
                WINDOW_KB,
                if (nrow(overlap_c1) > 0) sum(overlap_c1$ms_overlap) else 0,
                if (is.null(hits_c1)) 0 else nrow(hits_c1)))
    cat(sprintf("  C2 hits overlapping MS loci (±%dkb): %d / %d\n",
                WINDOW_KB,
                if (nrow(overlap_c2) > 0) sum(overlap_c2$ms_overlap) else 0,
                if (is.null(hits_c2)) 0 else nrow(hits_c2)))

    fwrite(overlap_c1, file.path(OUT_DIR, "ms_overlap_C1.csv"))
    fwrite(overlap_c2, file.path(OUT_DIR, "ms_overlap_C2.csv"))
}

# ── 5. Manhattan + QQ plots ───────────────────────────────────────────────────
make_manhattan <- function(dt, phenotype, hits) {
    if (is.null(dt)) return(NULL)
    # Downsample non-significant variants for plotting speed
    set.seed(42)
    ns  <- dt[p_value >= 0.01]
    sig <- dt[p_value < 0.01]
    ns_sample <- ns[sample(.N, min(200000L, .N))]
    plot_dt <- rbind(ns_sample, sig)

    # Compute x-axis cumulative position
    chr_sizes <- plot_dt[, .(max_pos = max(GENPOS)), by=CHROM][order(CHROM)]
    chr_sizes[, cum_start := cumsum(shift(max_pos, fill=0L))]
    plot_dt <- merge(plot_dt, chr_sizes[, .(CHROM, cum_start)], by="CHROM")
    plot_dt[, plot_pos := GENPOS + cum_start]
    chr_ticks <- chr_sizes[, .(CHROM, tick = cum_start + max_pos / 2)]

    plot_dt[, chr_col := fifelse(CHROM %% 2 == 0, "grey60", "grey30")]
    if (!is.null(hits) && nrow(hits) > 0) {
        hits_plot <- merge(hits[, .(ID, CHROM, GENPOS, p_value)],
                           chr_sizes[, .(CHROM, cum_start)], by="CHROM")
        hits_plot[, plot_pos := GENPOS + cum_start]
    }

    p <- ggplot(plot_dt, aes(plot_pos, -log10(p_value))) +
        geom_point(aes(colour=chr_col), size=0.4, alpha=0.5, show.legend=FALSE) +
        scale_colour_identity() +
        geom_hline(yintercept=-log10(GW_THR), linetype="dashed", colour="red", linewidth=0.5) +
        geom_hline(yintercept=-log10(SUGGESTIVE), linetype="dotted", colour="orange", linewidth=0.4) +
        scale_x_continuous(breaks=chr_ticks$tick,
                           labels=chr_ticks$CHROM, expand=c(0.01,0)) +
        labs(title=sprintf("Manhattan plot: %s vs None", phenotype),
             x="Chromosome", y="-log10(p)") +
        theme_ukb(base_size=10) +
        theme(panel.grid.major.x = element_blank(),
              axis.text.x = element_text(size=7))

    if (!is.null(hits) && nrow(hits) > 0) {
        p <- p + geom_point(data=hits_plot,
                            aes(plot_pos, -log10(p_value)),
                            colour="red", size=1.5)
    }
    p
}

make_qq <- function(dt, phenotype) {
    if (is.null(dt)) return(NULL)
    set.seed(42)
    obs  <- sort(-log10(dt[sample(.N, min(100000L, .N)), p_value]))
    exp  <- sort(-log10(seq(1, length(obs)) / (length(obs) + 1)))
    qq_dt <- data.table(exp=exp, obs=obs)
    lambda <- median(dt$p_value < 1, na.rm=TRUE)  # rough lambda placeholder
    ggplot(qq_dt, aes(exp, obs)) +
        geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
        geom_point(size=0.5, alpha=0.5, colour=COL_PRE_UP) +
        labs(title=sprintf("QQ plot: %s", phenotype),
             x="Expected -log10(p)", y="Observed -log10(p)") +
        theme_ukb(base_size=10)
}

cat("\nGenerating plots...\n")
plots <- list()
for (nm in c("C1","C2")) {
    dt   <- get(sprintf("gwas_%s", tolower(nm)))
    hits <- get(sprintf("hits_%s", tolower(nm)))
    if (is.null(dt)) next
    plots[[paste0(nm,"_manhattan")]] <- make_manhattan(dt, nm, hits)
    plots[[paste0(nm,"_qq")]]        <- make_qq(dt, nm)
}

if (length(plots) > 0) {
    pdf(file.path(OUT_DIR, "cluster_gwas_manhattan.pdf"), width=14, height=10)
    non_null <- Filter(Negate(is.null), plots)
    if (length(non_null) >= 2) print(wrap_plots(non_null, ncol=2))
    dev.off()
    cat("  Saved: cluster_gwas_manhattan.pdf\n")
}

# ── 6. Save hit tables ────────────────────────────────────────────────────────
all_hits <- rbind(
    if (!is.null(hits_c1) && nrow(hits_c1) > 0) cbind(hits_c1, phenotype="C1") else NULL,
    if (!is.null(hits_c2) && nrow(hits_c2) > 0) cbind(hits_c2, phenotype="C2") else NULL
)
if (!is.null(all_hits)) {
    fwrite(all_hits, file.path(OUT_DIR, "cluster_gwas_hits.csv"))
    cat(sprintf("  Saved: cluster_gwas_hits.csv (%d total loci)\n", nrow(all_hits)))
}

cat("\nDone. Results in", OUT_DIR, "\n")
