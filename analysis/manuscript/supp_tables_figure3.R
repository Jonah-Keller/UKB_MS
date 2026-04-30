#!/usr/bin/env Rscript
# supp_tables_figure3.R
# Generate supplementary tables for Figure 3 (HLA/PRS genetic heterogeneity)
# Outputs: results/supp_tables/supp_figure3_tables.xlsx (+ console summary)

suppressPackageStartupMessages({
    library(data.table)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

ENDO_DIR  <- file.path(PROJ_DIR, "results", "endophenotype")
GEN_DIR   <- file.path(PROJ_DIR, "data", "ukb", "genetics")
PRS_DIR   <- file.path(PROJ_DIR, "results", "prs")
OUT_DIR   <- file.path(PROJ_DIR, "results", "supp_tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

HAS_XLSX <- requireNamespace("openxlsx", quietly = TRUE)
if (HAS_XLSX) suppressPackageStartupMessages(library(openxlsx))

# =============================================================================
# Load all figure 3 result files
# =============================================================================
cat("Loading result files...\n")

load_or_null <- function(path, label) {
    if (file.exists(path)) {
        dt <- fread(path, showProgress = FALSE)
        cat(sprintf("  %-40s %d rows\n", label, nrow(dt)))
        return(dt)
    } else {
        cat(sprintf("  %-40s NOT FOUND\n", label))
        return(NULL)
    }
}

hla_dep     <- load_or_null(file.path(ENDO_DIR, "ms_hla_dep_results.csv"),        "within-MS HLA DEPs")
hla_int_all <- load_or_null(file.path(ENDO_DIR, "ms_hla_interaction_all.csv"),    "HLA×MS interaction (all)")
hla_int_pre <- load_or_null(file.path(ENDO_DIR, "ms_hla_interaction_pre.csv"),    "HLA×MS interaction (pre)")
hla_int_pst <- load_or_null(file.path(ENDO_DIR, "ms_hla_interaction_post.csv"),   "HLA×MS interaction (post)")
hla_strat   <- load_or_null(file.path(ENDO_DIR, "ms_hla_stratified_logfc.csv"),   "HLA-stratified MS logFC")
sex_strat   <- load_or_null(file.path(ENDO_DIR, "ms_sex_stratified_genetic.csv"), "Sex-stratified genetic")
prs_int     <- load_or_null(file.path(ENDO_DIR, "ms_prs_interaction_genome.csv"), "PRS×MS interaction")
comor       <- load_or_null(file.path(ENDO_DIR, "ms_hla_prems_comorbidity.csv"),  "Pre-MS comorbidity by HLA")
hla_enr     <- load_or_null(file.path(GEN_DIR,  "hla_allele_enrichment.csv"),     "HLA allele enrichment")
hla_enr_p   <- load_or_null(file.path(GEN_DIR,  "hla_allele_enrichment_pruned.csv"), "HLA allele enrichment (pruned)")
prs_res     <- load_or_null(file.path(PRS_DIR,   "ms_prs_results.csv"),           "PRS-proteome association")

# =============================================================================
# Console summary: story-relevant statistics
# =============================================================================
cat("\n=== FIGURE 3 NARRATIVE SUMMARY ===\n\n")

# --- HLA allele landscape ---
if (!is.null(hla_enr_p)) {
    cat("HLA ALLELE LANDSCAPE (pruned, LD-independent):\n")
    cat(sprintf("  Total alleles tested:       %d\n", nrow(hla_enr_p)))
    cat(sprintf("  FDR<0.05:                   %d\n", sum(hla_enr_p$fdr < 0.05, na.rm=TRUE)))
    risk <- hla_enr_p[fdr < 0.05 & enrichment > 0][order(-enrichment)]
    prot <- hla_enr_p[fdr < 0.05 & enrichment < 0][order(enrichment)]
    cat(sprintf("  Risk alleles (FDR<0.05):    %d\n", nrow(risk)))
    cat(sprintf("  Protective alleles:          %d\n", nrow(prot)))
    if (nrow(risk) > 0) {
        cat("  Top risk alleles:\n")
        print(risk[, .(allele_name = if("allele_name" %in% names(risk)) allele_name else paste0("idx",idx),
                        enrichment = round(enrichment,1), fdr = signif(fdr,3))][seq_len(min(5,.N))], row.names=FALSE)
    }
    cat("\n")
}

# --- Within-MS HLA DEPs ---
if (!is.null(hla_dep)) {
    if (!"fdr" %in% names(hla_dep)) hla_dep[, fdr := p.adjust(P.Value, method="BH")]
    cat("WITHIN-MS HLA+ vs HLA- DIFFERENTIAL PROTEINS:\n")
    cat(sprintf("  Total proteins tested:      %d\n", nrow(hla_dep)))
    cat(sprintf("  FDR<0.05:                   %d\n", sum(hla_dep$fdr < 0.05)))
    cat(sprintf("  Nominal p<0.05:             %d\n", sum(hla_dep$P.Value < 0.05)))
    cat(sprintf("  Nominal p<0.01:             %d\n", sum(hla_dep$P.Value < 0.01)))
    top_up   <- hla_dep[fdr < 0.05 & logFC > 0][order(-logFC)]
    top_down <- hla_dep[fdr < 0.05 & logFC < 0][order(logFC)]
    cat(sprintf("  FDR-sig up in HLA+:         %d\n", nrow(top_up)))
    cat(sprintf("  FDR-sig down in HLA+:       %d\n", nrow(top_down)))
    cat("  Top FDR-significant proteins (by logFC):\n")
    top_show <- hla_dep[fdr < 0.05][order(-abs(logFC))]
    print(top_show[, .(protein = toupper(protein),
                        logFC = round(logFC,3),
                        P.Value = signif(P.Value,3),
                        fdr = signif(fdr,3))][seq_len(min(10,.N))], row.names=FALSE)
    cat("\n")
}

# --- Sex-stratified within-MS HLA ---
if (!is.null(sex_strat) && "hla_within_ms" %in% sex_strat$analysis) {
    cat("SEX-STRATIFIED within-MS HLA DEPs:\n")
    hw <- sex_strat[analysis == "hla_within_ms"]
    hw_f <- hw[sex_group == "Female"]
    hw_m <- hw[sex_group == "Male"]
    cat(sprintf("  Female: FDR<0.05=%d  nominal p<0.05=%d\n",
        sum(hw_f$fdr < 0.05), sum(hw_f$P.Value < 0.05)))
    cat(sprintf("  Male:   FDR<0.05=%d  nominal p<0.05=%d\n",
        sum(hw_m$fdr < 0.05), sum(hw_m$P.Value < 0.05)))
    # Overlap
    f_nom <- hw_f[P.Value < 0.05, toupper(protein)]
    m_nom <- hw_m[P.Value < 0.05, toupper(protein)]
    cat(sprintf("  Overlap nominal p<0.05:     %d proteins in both sexes\n",
        length(intersect(f_nom, m_nom))))
    cat("  Top proteins significant in BOTH sexes (FDR):\n")
    both_fdr <- intersect(hw_f[fdr < 0.05, toupper(protein)],
                           hw_m[fdr < 0.05, toupper(protein)])
    cat(sprintf("    %s\n", if(length(both_fdr)>0) paste(both_fdr, collapse=", ") else "(none)"))
    cat("\n")
}

# --- HLA×MS interaction (formal) ---
if (!is.null(hla_int_all)) {
    if (!"fdr" %in% names(hla_int_all)) hla_int_all[, fdr := p.adjust(P.Value, method="BH")]
    cat("HLA×MS FORMAL INTERACTION (all MS vs HC):\n")
    cat(sprintf("  FDR<0.05:       %d\n", sum(hla_int_all$fdr < 0.05)))
    cat(sprintf("  Nominal p<0.05: %d\n", sum(hla_int_all$P.Value < 0.05)))
    top_int <- hla_int_all[P.Value < 0.05][order(P.Value)][seq_len(min(8,.N))]
    if (nrow(top_int) > 0) {
        cat("  Top nominal interaction hits:\n")
        print(top_int[, .(protein=toupper(protein), logFC=round(logFC,3),
                           P.Value=signif(P.Value,3), fdr=signif(fdr,3))], row.names=FALSE)
    }
    cat("\n")
}

# --- HLA-stratified MS logFC concordance ---
if (!is.null(hla_strat)) {
    cat("HLA-STRATIFIED MS-EFFECT CONCORDANCE:\n")
    r_val <- cor(hla_strat$logFC_hlaneg, hla_strat$logFC_hlapos, use="complete.obs")
    cat(sprintf("  Overall r (HLA- vs HLA+ MS effect):  %.3f\n", r_val))
    cat(sprintf("  DEPs in HLA- background (FDR<0.05):  %d\n",
        sum(hla_strat$fdr_hlaneg < 0.05, na.rm=TRUE)))
    cat(sprintf("  DEPs in HLA+ background (FDR<0.05):  %d\n",
        sum(hla_strat$fdr_hlapos < 0.05, na.rm=TRUE)))
    cat(sprintf("  DEPs in BOTH backgrounds (FDR<0.05): %d\n",
        sum(hla_strat$fdr_hlaneg < 0.05 & hla_strat$fdr_hlapos < 0.05, na.rm=TRUE)))
    # Proteins amplified specifically in HLA+ background
    if ("sig_int" %in% names(hla_strat)) {
        cat(sprintf("  Interaction-significant proteins:     %d\n",
            sum(hla_strat$sig_int, na.rm=TRUE)))
    }
    cat("\n")
}

# --- PRS×MS interaction ---
if (!is.null(prs_int)) {
    if (!"fdr" %in% names(prs_int)) prs_int[, fdr := p.adjust(P.Value, method="BH")]
    cat("PRS×MS INTERACTION (genome-wide, all proteins):\n")
    cat(sprintf("  FDR<0.05:       %d\n", sum(prs_int$fdr < 0.05)))
    cat(sprintf("  Nominal p<0.05: %d\n", sum(prs_int$P.Value < 0.05)))
    top_prs <- prs_int[order(P.Value)][seq_len(min(8,.N))]
    cat("  Top nominal PRS×MS interaction proteins:\n")
    print(top_prs[, .(protein=toupper(protein), logFC=round(logFC,3),
                       P.Value=signif(P.Value,3), fdr=signif(fdr,3))], row.names=FALSE)
    # Overlap between PRS interaction and HLA DEPs
    if (!is.null(hla_dep)) {
        hla_nom  <- hla_dep[P.Value < 0.05, tolower(protein)]
        prs_nom  <- prs_int[P.Value < 0.05, tolower(protein)]
        overlap  <- intersect(hla_nom, prs_nom)
        cat(sprintf("  Overlap HLA DEP (p<0.05) & PRS int (p<0.05): %d proteins\n",
            length(overlap)))
        if (length(overlap) > 0)
            cat(sprintf("    %s\n", paste(toupper(overlap[seq_len(min(8,length(overlap)))]),
                                          collapse=", ")))
    }
    cat("\n")
}

# --- Pre-MS comorbidity by HLA ---
if (!is.null(comor)) {
    cat("PRE-MS COMORBIDITY ENRICHMENT BY HLA STATUS:\n")
    cat(sprintf("  Comorbidities tested:  %d\n", nrow(comor)))
    sig_c <- comor[fdr < 0.2]
    cat(sprintf("  FDR<0.20 (suggestive): %d\n", nrow(sig_c)))
    cat(sprintf("  FDR<0.05:              %d\n", sum(comor$fdr < 0.05, na.rm=TRUE)))
    comor_show <- comor[order(fdr)]
    print(comor_show[, .(label, OR=round(OR,2), OR_lo=round(OR_lo,2), OR_hi=round(OR_hi,2),
                          fdr=signif(fdr,3))][seq_len(min(10,.N))], row.names=FALSE)
    cat("\n")
}

# --- PRS-proteome association ---
if (!is.null(prs_res)) {
    cat("PRS-PROTEOME ASSOCIATION (HLA-adjusted, within MS cases):\n")
    cat(sprintf("  Proteins tested:    %d\n", nrow(prs_res)))
    cat(sprintf("  FDR<0.05:           %d\n", sum(prs_res$fdr < 0.05, na.rm=TRUE)))
    cat(sprintf("  Nominal p<0.05:     %d\n", sum(prs_res$pval < 0.05, na.rm=TRUE)))
    top_prs_r <- prs_res[order(pval)][seq_len(min(8,.N))]
    print(top_prs_r[, .(protein=toupper(protein), rho=round(rho,3),
                          pval=signif(pval,3), fdr=signif(fdr,3))], row.names=FALSE)
    cat("\n")
}

# =============================================================================
# Write supplementary tables (Excel if openxlsx available, else CSVs)
# =============================================================================
cat("Writing supplementary tables...\n")

clean_hla_dep <- function(dt) {
    if (is.null(dt)) return(NULL)
    dt <- copy(dt)
    if (!"fdr" %in% names(dt)) dt[, fdr := p.adjust(P.Value, method="BH")]
    dt[, protein := toupper(protein)]
    setorder(dt, P.Value)
    dt[, .(protein, logFC = round(logFC,4), AveExpr = round(AveExpr,3),
           t = round(t,3), P.Value = signif(P.Value,4), FDR = signif(fdr,4))]
}

tables <- list(
    "S3a_HLA_DEPs"         = clean_hla_dep(hla_dep),
    "S3b_HLA_int_all"      = clean_hla_dep(hla_int_all),
    "S3c_HLA_int_pre"      = clean_hla_dep(hla_int_pre),
    "S3d_HLA_int_post"     = clean_hla_dep(hla_int_pst),
    "S3e_HLA_strat_logFC"  = if (!is.null(hla_strat)) {
        hs <- copy(hla_strat); hs[, protein := toupper(protein)]
        hs[, abs_lfc := abs(logFC_hlaneg)]; setorder(hs, -abs_lfc); hs[, abs_lfc := NULL]
        hs
    } else NULL,
    "S3f_Sex_stratified"   = if (!is.null(sex_strat)) {
        ss <- copy(sex_strat)[analysis == "hla_within_ms"]
        ss[, protein := toupper(protein)]
        setorder(ss, P.Value)
        ss
    } else NULL,
    "S3g_PRS_MS_int"       = clean_hla_dep(prs_int),
    "S3h_PRS_proteome"     = if (!is.null(prs_res)) {
        pr <- copy(prs_res); pr[, protein := toupper(protein)]
        setorder(pr, pval); pr
    } else NULL,
    "S3i_Comorbidity_HLA"  = comor,
    "S3j_HLA_allele_enr"   = if (!is.null(hla_enr_p)) {
        ha <- copy(hla_enr_p)
        setorder(ha, fdr)
        ha
    } else NULL
)
tables <- Filter(Negate(is.null), tables)

if (HAS_XLSX) {
    wb <- createWorkbook()
    for (nm in names(tables)) {
        addWorksheet(wb, nm)
        writeData(wb, nm, tables[[nm]])
        setColWidths(wb, nm, cols = seq_len(ncol(tables[[nm]])), widths = "auto")
    }
    out_xlsx <- file.path(OUT_DIR, "supp_figure3_tables.xlsx")
    saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    cat(sprintf("  Wrote: %s\n", out_xlsx))
} else {
    for (nm in names(tables)) {
        out_csv <- file.path(OUT_DIR, paste0(nm, ".csv"))
        fwrite(tables[[nm]], out_csv)
        cat(sprintf("  Wrote: %s\n", basename(out_csv)))
    }
}

cat("\nsupp_tables_figure3.R complete.\n")
