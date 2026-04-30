#!/usr/bin/env Rscript
# 01_build_phenotype.R
# Build REGENIE phenotype + covariate files for cluster GWAS
#
# Phenotype definitions (ICD-10 blocks enriched in each comorbidity cluster,
# derived from ms_prems_block_overrep.csv):
#
#   C1 (cranial nerve/demyelinating):
#     H4 = H40-H49 (glaucoma & other eye disorders)
#     G5 = G50-G59 (nerve, nerve root & plexus disorders)
#     M4 = M40-M49 (deforming dorsopathies / spondylosis)
#     CASE  = has ≥1 code in any of H4, G5, M4
#
#   C2 (neurological symptoms + EBV):
#     R2 = R20-R29 (symptoms & signs: nervous & musculoskeletal)
#     B2 = B20-B29 (viral infections; B27 = infectious mononucleosis/EBV)
#     CASE  = has ≥1 code in R2 AND ≥1 code in B2
#
#   CONTROL for both GWASs = no code in ANY of H4, G5, M4, R2, B2
#   (excludes the M5/Z8/H2/N1/K7 C0-cluster codes too to keep controls clean)
#   Participants with only C0-defining codes are dropped from both GWASs.
#
# Covariates: age (instance 0), sex, PC1-PC10, genotype array (UKBB vs UKBL)
#
# Output: data/ukb/genetics/gwas_cluster/
#   pheno_C1_vs_none.tsv    REGENIE-format (FID IID C1)
#   pheno_C2_vs_none.tsv    REGENIE-format (FID IID C2)
#   covariates.tsv          REGENIE-format (FID IID AGE SEX PC1..PC10 ARRAY)
#   sample_counts.txt       sanity check

suppressPackageStartupMessages({
    library(data.table)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))

ICD_DIR   <- file.path(dirname(PROJ_DIR),
                        "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                        "data", "ukb", "diagnoses", "icd_codes")
COV_FILE  <- file.path(PROJ_DIR, "data", "ukb", "covariates", "covariate_table.csv")
RAP_COV   <- file.path(PROJ_DIR, "data", "ukb", "rap_extraction", "covariates_pcs.csv")
OUT_DIR   <- file.path(PROJ_DIR, "data", "ukb", "genetics", "gwas_cluster")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ICD blocks defining each cluster
C1_BLOCKS  <- c("H4", "G5", "M4")
C2_BLOCKS  <- c("R2", "B2")
C0_BLOCKS  <- c("M5", "Z8", "H2")
EXCL_BLOCKS <- c(C0_BLOCKS, "N1", "K7")   # C0/partial C2 — drop from controls

# ── 1. Load ICD block events ───────────────────────────────────────────────────
load_block_eids <- function(block) {
    ch     <- substr(block, 1, 1)
    codes  <- list.dirs(file.path(ICD_DIR, ch, block), recursive=FALSE, full.names=FALSE)
    eids   <- integer(0)
    for (code in codes) {
        fp <- file.path(ICD_DIR, ch, block, code, paste0(code, ".csv"))
        if (!file.exists(fp)) next
        dt <- fread(fp, showProgress=FALSE, select=1:2)
        setnames(dt, c("eid","ever_present"))
        eids <- c(eids, dt[ever_present == 1L, as.integer(eid)])
    }
    unique(eids)
}

cat("Loading ICD block data for all UKB participants...\n")
block_eids <- list()
for (blk in c(C1_BLOCKS, C2_BLOCKS, EXCL_BLOCKS)) {
    block_eids[[blk]] <- load_block_eids(blk)
    cat(sprintf("  Block %-3s: %d participants\n", blk, length(block_eids[[blk]])))
}

c1_eids   <- Reduce(union, block_eids[C1_BLOCKS])
c2_r2_ids <- block_eids[["R2"]]
c2_b2_ids <- block_eids[["B2"]]
c2_eids   <- intersect(c2_r2_ids, c2_b2_ids)   # must have BOTH R2 and B2

excl_eids <- Reduce(union, block_eids[EXCL_BLOCKS])

cat(sprintf("\nC1 cases (H4 OR G5 OR M4): %d\n", length(c1_eids)))
cat(sprintf("C2 cases (R2 AND B2):       %d\n", length(c2_eids)))
cat(sprintf("C2 partial (R2 only):       %d\n", length(setdiff(c2_r2_ids, c2_b2_ids))))
cat(sprintf("Excluded (C0/N1/K7 only):   %d\n",
            length(setdiff(excl_eids, c(c1_eids, c2_eids)))))

# "None" controls: no code in any cluster-defining block
all_case_eids <- union(c1_eids, c2_eids)
all_excl_eids <- union(all_case_eids, excl_eids)

# ── 2. Load covariates ─────────────────────────────────────────────────────────
cat("\nLoading covariates...\n")
# RAP extraction has sex (p31) and genetic PCs (p22009_a1-10)
rap <- fread(RAP_COV, showProgress=FALSE)
setnames(rap, c("participant.eid","participant.p31"),
              c("eid","sex"), skip_absent=TRUE)
pc_cols_old <- paste0("participant.p22009_a", 1:10)
pc_cols_new <- paste0("PC", 1:10)
setnames(rap, pc_cols_old, pc_cols_new, skip_absent=TRUE)
rap <- rap[, c("eid","sex", pc_cols_new), with=FALSE]
setnames(rap, pc_cols_new, pc_cols_new)   # no-op, just for clarity

# Age at instance 0 from covariate_table
age_dt <- fread(COV_FILE, showProgress=FALSE, select=c("eid","age_i0"))

covar_dt <- merge(rap, age_dt, by="eid", all=FALSE)
setnames(covar_dt, "age_i0", "age")

# Genotype array: UKBB=1 (ukb22418 starts with chr; chip 1=Axiom),
# distinguish by checking if the EID was typed on the UKBB or UKBL array.
# UKB field 22006: genetic ethnic grouping; 22019: sex chromosome aneuploidy;
# 22021: in kinship inference exclusions. We use PC loading as a proxy.
# For now, set ARRAY as constant placeholder — replace if array FAM is available.
covar_dt[, ARRAY := 0L]   # TODO: fill from ukb22000 (chip ID) if needed

# Restrict to complete cases
covar_dt <- covar_dt[complete.cases(covar_dt[, .(age, sex, PC1)])]
cat(sprintf("  Participants with complete covariates: %d\n", nrow(covar_dt)))

# ── 3. Build REGENIE phenotype files ──────────────────────────────────────────
cat("\nBuilding phenotype files...\n")
all_eids <- covar_dt$eid
none_eids <- setdiff(all_eids, all_excl_eids)

cat(sprintf("  None (clean controls): %d\n", length(none_eids)))
cat(sprintf("  C1 cases with genotype: %d\n",
            sum(c1_eids %in% all_eids)))
cat(sprintf("  C2 cases with genotype: %d\n",
            sum(c2_eids %in% all_eids)))

make_pheno_file <- function(case_eids, none_eids, pheno_name, out_file) {
    case_set <- intersect(case_eids, all_eids)
    ctrl_set <- none_eids
    eids     <- c(case_set, ctrl_set)
    dt <- data.table(FID=eids, IID=eids)
    dt[, (pheno_name) := c(rep(1L, length(case_set)), rep(0L, length(ctrl_set)))]  # REGENIE: 1=case, 0=ctrl, NA=missing
    fwrite(dt, out_file, sep="\t", na="NA")
    cat(sprintf("  %s: %d cases / %d controls → %s\n",
                pheno_name, length(case_set), length(ctrl_set), basename(out_file)))
}

make_pheno_file(c1_eids, none_eids, "C1",
                file.path(OUT_DIR, "pheno_C1_vs_none.tsv"))
make_pheno_file(c2_eids, none_eids, "C2",
                file.path(OUT_DIR, "pheno_C2_vs_none.tsv"))

# ── 4. Write covariate file (REGENIE format) ───────────────────────────────────
cat("\nWriting covariate file...\n")
covar_out <- covar_dt[, .(FID=eid, IID=eid, AGE=age, SEX=sex,
                           PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                           ARRAY)]
fwrite(covar_out, file.path(OUT_DIR, "covariates.tsv"), sep="\t")
cat(sprintf("  Covariates: %d participants\n", nrow(covar_out)))

# ── 5. Sanity check summary ────────────────────────────────────────────────────
sink(file.path(OUT_DIR, "sample_counts.txt"))
cat("=== GWAS Phenotype Summary ===\n")
cat(sprintf("Date: %s\n\n", Sys.Date()))
cat(sprintf("Total UKB participants with covariates: %d\n", nrow(covar_dt)))
cat(sprintf("\nC1 (H4 or G5 or M4): %d cases\n", sum(c1_eids %in% all_eids)))
cat(sprintf("C2 (R2 and B2):       %d cases\n", sum(c2_eids %in% all_eids)))
cat(sprintf("None (clean ctrl):    %d\n", length(none_eids)))
cat(sprintf("\nNote: C1 and C2 phenotypes share the same control group.\n"))
cat(sprintf("REGENIE: 1=control, 2=case, -9=NA/excluded\n"))
sink()

cat("\nDone. Files in", OUT_DIR, "\n")
