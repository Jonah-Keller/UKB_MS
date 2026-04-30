#!/usr/bin/env Rscript
# 01_ms_comorbidity_endophenotype.R
# Comorbidity-based MS endophenotype analysis
#
# Design:
#   1. Load UKB first-occurrence ICD10 data for a curated panel of comorbidities
#      spanning autoimmune, neurological, cardiovascular, respiratory,
#      digestive, and musculoskeletal domains
#   2. PheWAS: logistic regression per condition ~ ms_status + age + sex
#   3. Build comorbidity matrix for MS cases (presence × timing)
#   4. Hierarchical/k-means clustering on significant comorbidities
#   5. Cross-tabulate comorbidity clusters with:
#        - Proteomic cluster (k=2 from 12_clustering)
#        - HLA-DRB1*15:01 carrier status
#        - MS PRS quartile
#   6. SENSITIVITY: re-run clustering restricted to:
#        - Pre-MS comorbidities (>1 yr before diagnosis)
#        - Post-MS comorbidities (>1 yr after diagnosis)
#      Tests whether endophenotype is etiological vs consequential.
#
# Timing:
#   Comorbidity timing = age_at_comorbidity - age_at_ms_diagnosis (years)
#   Negative = before MS diagnosis, positive = after
#   Note: ages used (not calendar dates), so timing is relative, not absolute
#
# Available chapters: A–Q, R–Z (full ICD-10 coverage as of 2026-04-18)
#
# Output: results/comorbidity/

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(cluster)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

ICD_DIR   <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                        "data", "ukb", "diagnoses", "icd_codes")
DIAG_FILE <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                        "data", "ukb", "diagnoses",
                        "ukb_first_occurrence_FINAL_20251031_144424.csv")
QC_FILE   <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
CLUST_FILE <- file.path(PROJ_DIR, "results", "clustering", "ms_cluster_assignments.csv")
HLA_FILE  <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE  <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
OUT_DIR   <- file.path(PROJ_DIR, "results", "comorbidity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ── 0. Curated comorbidity panel ──────────────────────────────────────────────
# ICD10 code → (label, domain, csv_path)
# Chapters available: G, H, I, J, K, L, M (E and F not parsed)

comorbidities <- list(
    # Autoimmune / inflammatory
    M05 = list(label="Seropositive RA",            domain="Autoimmune"),
    M06 = list(label="Other RA",                   domain="Autoimmune"),
    M32 = list(label="SLE",                        domain="Autoimmune"),
    M35 = list(label="Connective tissue disease",  domain="Autoimmune"),
    K50 = list(label="Crohn's disease",            domain="Autoimmune"),
    K51 = list(label="Ulcerative colitis",         domain="Autoimmune"),
    L40 = list(label="Psoriasis",                  domain="Autoimmune"),
    H20 = list(label="Uveitis/iridocyclitis",      domain="Autoimmune"),
    G61 = list(label="Inflammatory polyneuropathy",domain="Autoimmune"),
    G70 = list(label="Myasthenia gravis",          domain="Autoimmune"),
    E05 = list(label="Thyrotoxicosis (Graves')",   domain="Autoimmune"),
    E06 = list(label="Thyroiditis (Hashimoto's)",  domain="Autoimmune"),

    # Neurological (non-MS)
    G43 = list(label="Migraine",                   domain="Neurological"),
    G40 = list(label="Epilepsy",                   domain="Neurological"),
    G45 = list(label="TIA",                        domain="Neurological"),
    H46 = list(label="Optic neuritis",             domain="Neurological"),
    H81 = list(label="Vestibular disorder",        domain="Neurological"),
    G62 = list(label="Other polyneuropathy",       domain="Neurological"),
    G50 = list(label="Trigeminal neuralgia",       domain="Neurological"),

    # Psychiatric
    F32 = list(label="Depression",                 domain="Psychiatric"),
    F41 = list(label="Anxiety disorders",          domain="Psychiatric"),
    F45 = list(label="Somatoform disorders",       domain="Psychiatric"),

    # Cardiovascular
    I10 = list(label="Hypertension",               domain="Cardiovascular"),
    I25 = list(label="Chronic IHD",                domain="Cardiovascular"),
    I21 = list(label="Myocardial infarction",      domain="Cardiovascular"),
    I48 = list(label="Atrial fibrillation",        domain="Cardiovascular"),
    I50 = list(label="Heart failure",              domain="Cardiovascular"),
    I63 = list(label="Cerebral infarction",        domain="Cardiovascular"),

    # Respiratory
    J45 = list(label="Asthma",                     domain="Respiratory"),
    J44 = list(label="COPD",                       domain="Respiratory"),

    # Digestive
    K21 = list(label="GERD",                       domain="Digestive"),
    K58 = list(label="IBS",                        domain="Digestive"),
    K57 = list(label="Diverticular disease",       domain="Digestive"),

    # Musculoskeletal (non-autoimmune)
    M54 = list(label="Back pain",                  domain="Musculoskeletal"),
    M25 = list(label="Joint disorders NEC",        domain="Musculoskeletal"),
    M79 = list(label="Soft tissue disorders",      domain="Musculoskeletal"),
    M47 = list(label="Spondylosis",                domain="Musculoskeletal"),

    # Eye
    H04 = list(label="Lacrimal system disorder",   domain="Other")
)

# Helper: load per-code CSV
load_icd <- function(code) {
    ch  <- substr(code, 1, 1)
    gr  <- substr(code, 1, 2)
    fp  <- file.path(ICD_DIR, ch, gr, code, paste0(code, ".csv"))
    if (!file.exists(fp)) return(NULL)
    dt <- fread(fp, showProgress=FALSE)
    # Cols: eid, diagnosis (1=present), age_at_diagnosis
    setnames(dt, c("eid","diagnosis","age_at_event"), skip_absent=TRUE)
    if (!"age_at_event" %in% names(dt))
        setnames(dt, names(dt)[3], "age_at_event")
    dt[, icd10 := code]
    dt[, .(eid=as.integer(eid), present=as.integer(diagnosis),
           age_at_event=as.numeric(age_at_event), icd10)]
}

# ── 1. Load MS cohort ─────────────────────────────────────────────────────────
cat("Loading MS cohort...\n")
ms_qc <- fread(QC_FILE, showProgress=FALSE)
ms_qc <- ms_qc[qc_outlier==FALSE & !is.na(age_at_sampling) & !is.na(sex)]
ms_qc[, is_ms := as.integer(ms_status %in% c("pre_onset","post_onset"))]

# age_at_diagnosis: from ms_qc (age_at_diagnosis column)
# timing reference: age at MS diagnosis
ms_cases <- ms_qc[is_ms==1 & !is.na(age_at_diagnosis)]
cat(sprintf("  MS cases with diagnosis age: %d\n", nrow(ms_cases)))
cat(sprintf("  Controls: %d\n", sum(ms_qc$is_ms==0)))

# ── 2. Load all comorbidity data ───────────────────────────────────────────────
cat("Loading comorbidity data...\n")
codes_to_load <- names(comorbidities)
icd_list <- lapply(codes_to_load, function(code) {
    dt <- load_icd(code)
    if (is.null(dt)) {
        cat(sprintf("  MISSING: %s (%s)\n", code,
                    comorbidities[[code]]$label))
        return(NULL)
    }
    cat(sprintf("  %s (%s): %d participants\n", code,
                comorbidities[[code]]$label, sum(dt$present)))
    dt
})
icd_dt <- rbindlist(Filter(Negate(is.null), icd_list), fill=TRUE)
cat(sprintf("  Loaded: %d ICD10 codes\n", length(unique(icd_dt$icd10))))

# Create wide matrix: eid × icd10 → present (0/1) + age_at_event
# For each eid × icd10 keep one row (some codes may have multiple entries — take first)
icd_wide_pres <- dcast(icd_dt[, .(eid, icd10, present)],
                        eid ~ icd10, value.var="present",
                        fun.aggregate=max, fill=0L)
icd_wide_age  <- dcast(icd_dt[, .(eid, icd10, age_at_event)],
                        eid ~ icd10, value.var="age_at_event",
                        fun.aggregate=min, fill=NA_real_)

# Available codes after loading
avail_codes <- setdiff(names(icd_wide_pres), "eid")
cat(sprintf("  Available codes: %d\n", length(avail_codes)))

# Merge with cohort
cohort_wide <- merge(ms_qc[, .(eid, ms_status, is_ms, age_at_sampling,
                                age_at_diagnosis, sex)],
                     icd_wide_pres, by="eid", all.x=TRUE)
# Fill NA → 0 (not diagnosed = absent)
for (cc in avail_codes) {
    set(cohort_wide, which(is.na(cohort_wide[[cc]])), cc, 0L)
}

cat(sprintf("  Total cohort with comorbidity data: %d\n", nrow(cohort_wide)))

# ── 3. PheWAS: comorbidity enrichment in MS vs HC ─────────────────────────────
cat("\nPheWAS: testing comorbidity enrichment in MS vs HC...\n")

phewas_list <- lapply(avail_codes, function(cc) {
    sub <- cohort_wide[!is.na(get(cc))]
    n_ms_with <- sum(sub[is_ms==1, get(cc)])
    n_hc_with <- sum(sub[is_ms==0, get(cc)])
    if (n_ms_with < 3) return(NULL)

    fit <- tryCatch(
        glm(get(cc) ~ is_ms + age_at_sampling + sex,
            data=sub, family=binomial),
        error=function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    cf <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
    if (is.null(cf) || !"is_ms" %in% rownames(cf)) return(NULL)

    or  <- exp(cf["is_ms","Estimate"])
    pval <- cf["is_ms","Pr(>|z|)"]
    ci   <- tryCatch(confint(fit)["is_ms",], error=function(e) c(NA,NA))

    data.table(
        icd10    = cc,
        label    = comorbidities[[cc]]$label,
        domain   = comorbidities[[cc]]$domain,
        OR       = or,
        OR_lo    = exp(ci[1]),
        OR_hi    = exp(ci[2]),
        pval     = pval,
        n_ms     = n_ms_with,
        n_hc     = n_hc_with,
        pct_ms   = 100*n_ms_with/sum(sub$is_ms),
        pct_hc   = 100*n_hc_with/sum(!sub$is_ms)
    )
})
phewas_dt <- rbindlist(Filter(Negate(is.null), phewas_list))
phewas_dt[, fdr := p.adjust(pval, method="BH")]
phewas_dt <- phewas_dt[order(pval)]

cat("\nPheWAS results (all conditions):\n")
print(phewas_dt[, .(icd10, label, OR=round(OR,2), pct_ms=round(pct_ms,1),
                     pct_hc=round(pct_hc,1), pval=signif(pval,3), fdr=round(fdr,3))])

fwrite(phewas_dt, file.path(OUT_DIR, "ms_phewas_comorbidity.csv"))

# ── 4. PheWAS forest plot ──────────────────────────────────────────────────────
cat("\nPlotting PheWAS forest plot...\n")

phewas_dt[, sig := fcase(
    fdr  < 0.05, "FDR<0.05",
    pval < 0.05, "p<0.05",
    default     = "NS"
)]
# Domain colours — all derived from UKB_PALETTE + UKB_NAVY + neutrals
DOMAIN_COLS <- c(
    Autoimmune      = "#CC0066",   # magenta  (UKB_PALETTE[4])
    Neurological    = "#2B4C7E",   # navy     (UKB_NAVY)
    Psychiatric     = "#7B2D8B",   # purple
    Cardiovascular  = "#E6A817",   # amber    (UKB_PALETTE[2])
    Respiratory     = "#56B4E9",   # sky-blue (UKB_PALETTE[3])
    Digestive       = "#1A1A1A",   # near-black (UKB_PALETTE[1])
    Musculoskeletal = "grey55",    # neutral mid-grey
    Other           = "grey80"     # light grey
)
phewas_dt[, cond_f := factor(label, levels=rev(label))]
phewas_dt[, domain_f := factor(domain, levels=names(DOMAIN_COLS))]

phewas_dt[, log2or := log2(pmax(OR, 0.01))]
phewas_dt[, cond_f := factor(label, levels=label[order(log2or)])]
phewas_dt[, domain_f := factor(domain, levels=names(DOMAIN_COLS))]

p_phewas <- ggplot(phewas_dt, aes(x=log2or, y=cond_f, colour=domain_f)) +
    geom_vline(xintercept=0, linewidth=0.35, colour="grey50") +
    geom_segment(aes(x=0, xend=log2or, y=cond_f, yend=cond_f,
                     colour=domain_f, alpha=sig),
                 linewidth=0.6) +
    geom_point(aes(size=sig, alpha=sig)) +
    geom_text(data=phewas_dt[pval<0.05],
              aes(x=log2or + sign(log2or)*0.08,
                  label=sprintf("%.0f%%", pct_ms)),
              size=2.2, hjust=ifelse(phewas_dt[pval<0.05, log2or]>0,0,1),
              colour="grey30") +
    scale_colour_manual(values=DOMAIN_COLS, name="Domain") +
    scale_alpha_manual(values=c("FDR<0.05"=1,"p<0.05"=0.75,"NS"=0.35),
                       guide="none") +
    scale_size_manual(values=c("FDR<0.05"=3,"p<0.05"=2,"NS"=1.2),
                      guide="none") +
    labs(x="log\u2082 OR (MS vs HC, adjusted for age + sex)",
         y=NULL,
         title="Comorbidity enrichment in MS",
         subtitle=sprintf("MS n=%d, HC n=%d | BH FDR | %% = MS prevalence",
                          sum(cohort_wide$is_ms), sum(!cohort_wide$is_ms))) +
    theme_ukb(base_size=9) +
    theme(axis.text.y=element_text(size=8))

ggsave(file.path(OUT_DIR, "ms_phewas_forest.pdf"), p_phewas,
       width=7.5, height=max(5, nrow(phewas_dt)*0.28+2), device=cairo_pdf)
cat("  Saved: ms_phewas_forest.pdf\n")

# ── 5. Comorbidity matrix for MS cases ────────────────────────────────────────
cat("\nBuilding comorbidity matrix for MS cases...\n")

# Use all conditions with p<0.10 OR pre-specified clinical interest
sig_codes <- phewas_dt[pval < 0.10, icd10]
# Always include these clinically important ones even if NS
always_include <- c("G43","G40","H46","M05","K50","K51","L40","M32","I10","J45","K58")
use_codes <- unique(c(sig_codes, intersect(always_include, avail_codes)))
cat(sprintf("  Conditions for clustering: %d\n", length(use_codes)))

ms_only <- cohort_wide[is_ms==1 & !is.na(age_at_diagnosis)]
X_comorbid <- as.matrix(ms_only[, ..use_codes])
rownames(X_comorbid) <- ms_only$eid

# Remove conditions with <3 cases in MS
col_sums <- colSums(X_comorbid, na.rm=TRUE)
use_codes2 <- use_codes[col_sums[use_codes] >= 3]
X_comorbid <- X_comorbid[, use_codes2, drop=FALSE]
cat(sprintf("  After prevalence filter (≥3 MS cases): %d conditions\n", ncol(X_comorbid)))

# ── 6. Timing matrix: years relative to MS diagnosis ─────────────────────────
cat("Computing comorbidity timing relative to MS diagnosis...\n")

# Merge age at event with MS cases
ms_ages <- icd_wide_age[eid %in% ms_only$eid]
ms_timing <- merge(ms_only[, .(eid, age_at_diagnosis)],
                   ms_ages[, c("eid", use_codes2), with=FALSE],
                   by="eid", all.x=TRUE)

# timing = age_at_comorbidity - age_at_ms_diagnosis (negative = before MS)
timing_mat <- ms_timing[, ..use_codes2] - ms_timing$age_at_diagnosis
timing_mat <- as.matrix(timing_mat)
rownames(timing_mat) <- ms_timing$eid

# Classify each comorbidity × MS case as: pre (< -1yr), peri (-1 to +1yr), post (>+1yr)
pre_mat  <- (!is.na(timing_mat) & timing_mat < -1) * X_comorbid[rownames(timing_mat),]
post_mat <- (!is.na(timing_mat) & timing_mat >  1) * X_comorbid[rownames(timing_mat),]

# Timing summary per condition (among MS cases who HAVE the condition)
timing_summary <- lapply(use_codes2, function(cc) {
    has_it <- which(!is.na(timing_mat[,cc]) & X_comorbid[,cc]==1)
    if (length(has_it) < 2) return(NULL)
    t_vals <- timing_mat[has_it, cc]
    data.table(
        icd10   = cc,
        label   = comorbidities[[cc]]$label,
        n_with  = length(has_it),
        pct_pre = 100*mean(t_vals < -1, na.rm=TRUE),
        pct_peri= 100*mean(abs(t_vals) <= 1, na.rm=TRUE),
        pct_post= 100*mean(t_vals > 1, na.rm=TRUE),
        median_timing = median(t_vals, na.rm=TRUE)
    )
})
timing_dt <- rbindlist(Filter(Negate(is.null), timing_summary))
timing_dt <- timing_dt[order(median_timing)]
cat("\nComorbidity timing relative to MS diagnosis:\n")
print(timing_dt[, .(icd10, label, n_with, pct_pre=round(pct_pre,0),
                     pct_post=round(pct_post,0), median_timing=round(median_timing,1))])
fwrite(timing_dt, file.path(OUT_DIR, "ms_comorbidity_timing.csv"))

# ── 7. Timing violin plot ─────────────────────────────────────────────────────
cat("\nPlotting comorbidity timing...\n")

timing_long <- melt(as.data.table(timing_mat)[, eid := ms_timing$eid],
                    id.vars="eid", variable.name="icd10", value.name="years_rel")
timing_long[, present := {
    ri <- match(eid, rownames(X_comorbid))
    ci <- match(as.character(icd10), colnames(X_comorbid))
    X_comorbid[cbind(ri, ci)]
}]
timing_long <- timing_long[present==1 & !is.na(years_rel)]
timing_long <- merge(timing_long,
                     phewas_dt[, .(icd10, label, domain, fdr)],
                     by="icd10", all.x=TRUE)
timing_long[, label_f := factor(label, levels=timing_dt$label)]

p_timing <- ggplot(timing_long[!is.na(label_f)],
                   aes(x=years_rel, y=label_f, fill=domain)) +
    geom_vline(xintercept=0, linewidth=0.6, colour="#CC0066",
               linetype="dashed") +
    geom_vline(xintercept=c(-1,1), linewidth=0.3, colour="grey60",
               linetype="dotted") +
    geom_boxplot(width=0.5, outlier.size=0.5, linewidth=0.35, alpha=0.7) +
    scale_fill_manual(values=DOMAIN_COLS, name="Domain") +
    labs(x="Years relative to MS diagnosis (negative = before)",
         y=NULL,
         title="Comorbidity timing relative to MS diagnosis",
         caption="Red dashed line = MS diagnosis; ±1yr grey dotted = peri-MS window") +
    theme_ukb() +
    theme(axis.text.y=element_text(size=8))

ggsave(file.path(OUT_DIR, "ms_comorbidity_timing.pdf"), p_timing,
       width=9, height=max(5, nrow(timing_dt)*0.35+2), device=cairo_pdf)
cat("  Saved: ms_comorbidity_timing.pdf\n")

# ── 8. Clustering function ────────────────────────────────────────────────────
run_comorbidity_clustering <- function(mat, label_suffix="", k_max=4) {
    if (nrow(mat) < 20 || ncol(mat) < 3) {
        cat(sprintf("  Skipping clustering [%s]: too few cases/conditions\n", label_suffix))
        return(NULL)
    }

    # Replace any residual NAs with 0 (absent) before clustering
    mat[is.na(mat)] <- 0L

    # Optimal k by silhouette on Jaccard distances
    jac_dist <- dist(mat, method="binary")  # Jaccard distance for binary data

    sil_vals <- numeric(k_max - 1)
    for (k in 2:k_max) {
        km <- tryCatch(
            kmeans(mat, centers=k, nstart=50, iter.max=200),
            error=function(e) NULL
        )
        if (is.null(km)) next
        sil <- silhouette(km$cluster, jac_dist)
        sil_vals[k-1] <- mean(sil[,3])
        cat(sprintf("  [%s] k=%d: silhouette=%.3f, sizes=%s\n",
                    label_suffix, k,
                    sil_vals[k-1], paste(table(km$cluster), collapse="/")))
    }
    k_opt <- which.max(sil_vals) + 1
    cat(sprintf("  [%s] Optimal k=%d\n", label_suffix, k_opt))

    km_final <- kmeans(mat, centers=k_opt, nstart=100, iter.max=200)
    data.table(eid=as.integer(rownames(mat)),
               comorbidity_cluster=km_final$cluster,
               cluster_label=label_suffix)
}

# ── 9. Primary clustering: all comorbidities ──────────────────────────────────
cat("\nClustering MS cases by comorbidity profile...\n")
clust_all <- run_comorbidity_clustering(X_comorbid, "all_comorbidities")

# ── 10. Sensitivity 1: pre-MS comorbidities only ──────────────────────────────
cat("\nSensitivity: pre-MS comorbidities only (>1yr before diagnosis)...\n")
pre_mat2 <- pre_mat[, use_codes2, drop=FALSE]
pre_col_sums <- colSums(pre_mat2, na.rm=TRUE)
pre_codes <- use_codes2[pre_col_sums >= 3]
pre_mat2 <- pre_mat2[, pre_codes, drop=FALSE]
cat(sprintf("  Pre-MS conditions with ≥3 MS cases: %d\n", ncol(pre_mat2)))
clust_pre <- run_comorbidity_clustering(pre_mat2, "pre_MS_only")

# ── 11. Sensitivity 2: post-MS comorbidities only ─────────────────────────────
cat("\nSensitivity: post-MS comorbidities only (>1yr after diagnosis)...\n")
post_mat2 <- post_mat[, use_codes2, drop=FALSE]
post_col_sums <- colSums(post_mat2, na.rm=TRUE)
post_codes <- use_codes2[post_col_sums >= 3]
post_mat2 <- post_mat2[, post_codes, drop=FALSE]
cat(sprintf("  Post-MS conditions with ≥3 MS cases: %d\n", ncol(post_mat2)))
clust_post <- run_comorbidity_clustering(post_mat2, "post_MS_only")

# ── 12. Cross-tabulation with proteomic cluster + HLA + PRS ──────────────────
cat("\nCross-tabulating comorbidity clusters with molecular data...\n")

# Merge all cluster assignments
merge_clusters <- function(comorbidity_clust) {
    if (is.null(comorbidity_clust)) return(NULL)
    dt <- copy(comorbidity_clust)

    # Proteomic cluster
    if (file.exists(CLUST_FILE)) {
        prot_c <- fread(CLUST_FILE, select=c("eid","cluster"))
        setnames(prot_c, "cluster", "prot_cluster")
        dt <- merge(dt, prot_c, by="eid", all.x=TRUE)
    }
    # HLA
    if (file.exists(HLA_FILE)) {
        hla <- fread(HLA_FILE, select=c("eid","drb1_1501_carrier"))
        dt <- merge(dt, hla, by="eid", all.x=TRUE)
    }
    # PRS
    if (file.exists(PRS_FILE)) {
        prs <- fread(PRS_FILE, select=c("eid","prs_score"))
        prs[, prs_quartile := cut(prs_score, quantile(prs_score, 0:4/4),
                                   labels=c("Q1","Q2","Q3","Q4"),
                                   include.lowest=TRUE)]
        dt <- merge(dt, prs[, .(eid, prs_quartile)], by="eid", all.x=TRUE)
    }
    # MS status + timing
    dt <- merge(dt, ms_only[, .(eid, ms_status, age_at_diagnosis,
                                 age_at_sampling, sex)],
                by="eid", all.x=TRUE)
    dt
}

ct_all <- merge_clusters(clust_all)

if (!is.null(ct_all) && "prot_cluster" %in% names(ct_all)) {
    cat("\nCross-tabulation: comorbidity cluster × proteomic cluster\n")
    xtab <- table(ct_all$comorbidity_cluster, ct_all$prot_cluster,
                  dnn=c("Comorbidity", "Proteomic"))
    print(xtab)
    chi_prot <- chisq.test(xtab)
    cat(sprintf("  Chi-sq p=%.3g\n", chi_prot$p.value))
}

if (!is.null(ct_all) && "drb1_1501_carrier" %in% names(ct_all)) {
    cat("\nCross-tabulation: comorbidity cluster × HLA-DRB1*15:01\n")
    xtab_hla <- table(ct_all$comorbidity_cluster, ct_all$drb1_1501_carrier,
                      dnn=c("Comorbidity", "HLA"))
    print(xtab_hla)
    chi_hla <- chisq.test(xtab_hla)
    cat(sprintf("  Chi-sq p=%.3g\n", chi_hla$p.value))
}

# ── 13. Cross-tabulation heatmap ──────────────────────────────────────────────
if (!is.null(ct_all) && ncol(ct_all) > 3) {
    # Characterize each comorbidity cluster
    clust_char <- ct_all[, .(
        n             = .N,
        pct_pre_onset = 100*mean(ms_status=="pre_onset", na.rm=TRUE),
        pct_female    = 100*mean(sex==0, na.rm=TRUE),
        mean_age_dx   = mean(age_at_diagnosis, na.rm=TRUE),
        pct_hla_carrier = if("drb1_1501_carrier" %in% names(.SD))
                            100*mean(drb1_1501_carrier, na.rm=TRUE) else NA_real_,
        pct_prot_c1   = if("prot_cluster" %in% names(.SD))
                            100*mean(prot_cluster==1, na.rm=TRUE) else NA_real_,
        mean_prs      = if("prs_quartile" %in% names(.SD))
                            mean(as.integer(prs_quartile), na.rm=TRUE) else NA_real_
    ), by=comorbidity_cluster]
    cat("\nComorbidity cluster characteristics:\n")
    print(clust_char)
    fwrite(clust_char, file.path(OUT_DIR, "ms_comorbidity_cluster_chars.csv"))
}

# ── 14. Comorbidity heatmap (geom_tile) ──────────────────────────────────────
cat("\nCreating comorbidity heatmap...\n")
if (!is.null(clust_all) && nrow(X_comorbid) > 10) {

    # Build long-form data for geom_tile
    heat_dt <- as.data.table(X_comorbid, keep.rownames="eid")
    heat_dt[, eid := as.integer(eid)]
    heat_long <- melt(heat_dt, id.vars="eid",
                      variable.name="icd10", value.name="present")
    # Attach cluster and label
    heat_long <- merge(heat_long,
                       clust_all[, .(eid, comorbidity_cluster)],
                       by="eid", all.x=TRUE)
    heat_long[, cond_label := sapply(as.character(icd10), function(cc)
        if (cc %in% names(comorbidities)) comorbidities[[cc]]$label else cc)]
    heat_long[, domain := sapply(as.character(icd10), function(cc)
        if (cc %in% names(comorbidities)) comorbidities[[cc]]$domain else "Other")]
    # Order patients by cluster, conditions by domain then prevalence
    cond_order <- heat_long[, .(pct=mean(present)), by=.(icd10, domain)][order(domain, -pct)]
    heat_long[, cond_f := factor(cond_label,
                                  levels=sapply(cond_order$icd10, function(cc)
                                      if (cc %in% names(comorbidities))
                                          comorbidities[[cc]]$label else cc))]
    heat_long[, clust_f := factor(paste0("C", comorbidity_cluster))]

    # Compute prevalence per (cluster × condition) for heatmap
    prev_dt <- heat_long[!is.na(comorbidity_cluster),
                          .(pct=100*mean(present)), by=.(clust_f, cond_label, domain)]

    p_heatmap <- ggplot(prev_dt,
                        aes(x=clust_f, y=cond_label, fill=pct)) +
        geom_tile(colour="white", linewidth=0.4) +
        geom_text(aes(label=ifelse(pct>=5, sprintf("%.0f%%",pct), "")),
                  size=2.3, colour="grey10") +
        facet_grid(domain ~ ., scales="free_y", space="free_y") +
        scale_fill_gradient2(low="white", mid="#56B4E9", high="#CC0066",
                             midpoint=10, limits=c(0, NA),
                             name="Prevalence (%)") +
        labs(x="Comorbidity cluster", y=NULL,
             title=sprintf("MS comorbidity prevalence by cluster (n=%d)", nrow(X_comorbid))) +
        theme_ukb(base_size=8) +
        theme(strip.text.y=element_text(angle=0, size=7, face="bold"),
              axis.text.x=element_text(face="bold"))

    ggsave(file.path(OUT_DIR, "ms_comorbidity_heatmap.pdf"), p_heatmap,
           width=12, height=6, device=cairo_pdf)
    cat("  Saved: ms_comorbidity_heatmap.pdf\n")
}

# ── 15. Sensitivity comparison: pre vs post vs all clustering ─────────────────
cat("\nSensitivity analysis summary:\n")
sens_compare <- list()
for (cc_name in c("all","pre","post")) {
    clust_obj <- get(paste0("clust_",cc_name))
    if (!is.null(clust_obj)) {
        k_used  <- max(clust_obj$comorbidity_cluster)
        sizes   <- paste(as.integer(table(clust_obj$comorbidity_cluster)), collapse="/")
        sens_compare[[cc_name]] <- data.table(
            analysis=cc_name, k=k_used, cluster_sizes=sizes
        )
    }
}
if (length(sens_compare) > 0) {
    sens_dt <- rbindlist(sens_compare)
    cat("\n"); print(sens_dt)
    fwrite(sens_dt, file.path(OUT_DIR, "ms_sensitivity_summary.csv"))
}

# Compare: do pre-MS and all-comorbidity clusters agree?
if (!is.null(clust_all) && !is.null(clust_pre)) {
    shared_eids <- intersect(clust_all$eid, clust_pre$eid)
    comp_dt <- merge(clust_all[eid %in% shared_eids, .(eid, comorbidity_cluster)],
                     clust_pre[eid %in% shared_eids, .(eid, comorbidity_cluster)],
                     by="eid", suffixes=c("_all","_pre"))
    xtab_sens <- table(comp_dt$comorbidity_cluster_all,
                       comp_dt$comorbidity_cluster_pre)
    cat("\nAll-comorbidity vs pre-MS-only cluster agreement:\n")
    print(xtab_sens)
    if (nrow(xtab_sens) > 1 && ncol(xtab_sens) > 1) {
        chi_s <- chisq.test(xtab_sens)
        cat(sprintf("  Agreement chi-sq p=%.3g (concordant clusters?)\n", chi_s$p.value))
    }
}

# ── 16. Save merged cluster assignment ───────────────────────────────────────
if (!is.null(ct_all)) {
    fwrite(ct_all, file.path(OUT_DIR, "ms_endophenotype_assignments.csv"))
    cat(sprintf("\nSaved endophenotype assignments: %d MS cases\n", nrow(ct_all)))
}

cat("\n01_ms_comorbidity_endophenotype.R complete.\n")
