#!/usr/bin/env Rscript
# 01_ms_multimodal_endophenotype.R  —  Figure 5
# Multi-modal MS endophenotype discovery
#
# DESIGN PRINCIPLE: strict input/outcome separation
#
#   INPUT features (what we cluster on — independent, phenotypic):
#     - PROTEOMICS  top 10 PCs from the 200 most variable proteins in MS cases
#     - COMORBIDITY 24-condition binary matrix (all timing windows)
#   These are clinical and molecular *observations* collected before/at sampling.
#
#   OUTCOMES (pre-specified, tested AFTER clustering — never seen by the algorithm):
#     - HLA-DRB1*15:01 carrier status
#     - MS PRS (polygenic risk score)
#     - Pre/post-MS onset timing (prodromal vs. prevalent)
#     - Age at diagnosis, sex, years to diagnosis
#     - Independent proteomic cluster (k=2, from analysis 12)
#
#   Rationale: including HLA or PRS as inputs causes the UMAP to cluster on
#   genotype rather than discovering novel multi-modal structure.
#
# Figure 5 panels:
#   A — UMAP coloured by multi-modal endophenotype cluster
#   B — UMAP coloured by proteomic cluster (outcome overlay)
#   C — UMAP coloured by HLA-DRB1*15:01 (outcome overlay)
#   D — UMAP coloured by MS PRS quartile (outcome overlay)
#   E — UMAP coloured by pre/post-onset timing (outcome overlay)
#   F — UMAP coloured by comorbidity cluster (outcome overlay)
#   G — Feature profile heatmap per endophenotype
#   H — Clinical validation: age at dx, sex ratio, prodromal window
#
# Sensitivity: re-run with pre-MS comorbidities only (>1yr before diagnosis)

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

QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
DEP_FILE    <- file.path(PROJ_DIR, "results", "differential", "ms_all_vs_hc.csv")
CLUST_FILE  <- file.path(PROJ_DIR, "results", "clustering", "ms_cluster_assignments.csv")
COMORBID_FILE <- file.path(PROJ_DIR, "results", "comorbidity", "ms_endophenotype_assignments.csv")
ICD_DIR     <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal",
                          "data", "ukb", "diagnoses", "icd_codes")
HLA_FILE    <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE    <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
OUT_DIR     <- file.path(PROJ_DIR, "results", "endophenotype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

ENDO_COLS <- c("1"="#CC0066","2"="#2B4C7E","3"="#E6A817","4"="#56B4E9","5"="#1A1A1A")

# ── 1. INPUT A: Proteomic features ────────────────────────────────────────────
cat("=== INPUT A: Proteomics ===\n")
qc       <- fread(QC_FILE, showProgress=FALSE)
ms_cases <- qc[qc_outlier==FALSE &
               ms_status %in% c("pre_onset","post_onset") &
               !is.na(age_at_diagnosis)]

meta_cols <- c("eid","ms_status","age_at_sampling","age_at_diagnosis",
               "years_to_diagnosis","sex","olink_instance","qc_outlier",
               "UMAP1","UMAP2","mean_npx","sd_npx",
               "lower_iqr","upper_iqr","iqr_npx")
prot_cols <- setdiff(names(ms_cases), meta_cols)
cat(sprintf("  MS cases: %d | Proteins available: %d\n",
            nrow(ms_cases), length(prot_cols)))

# Select top 200 most variable proteins WITHIN MS cases
# (not DEPs selected for MS vs HC difference — captures intra-MS heterogeneity)
prot_mat <- as.matrix(ms_cases[, ..prot_cols])
prot_sds <- apply(prot_mat, 2, sd, na.rm=TRUE)
top200   <- names(sort(prot_sds, decreasing=TRUE))[seq_len(min(200, length(prot_sds)))]
cat(sprintf("  Top 200 most variable proteins (SD range: %.2f – %.2f)\n",
            min(prot_sds[top200]), max(prot_sds[top200])))

prot_mat200 <- as.matrix(ms_cases[, ..top200])
prot_sc     <- scale(prot_mat200)
prot_sc[is.na(prot_sc)] <- 0   # impute residual missing with column mean (=0 after scale)

pca_res  <- prcomp(prot_sc, center=FALSE, scale.=FALSE)
n_pcs    <- 10
prot_pcs <- scale(pca_res$x[, seq_len(n_pcs)])   # z-score each PC
pct_var  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
colnames(prot_pcs) <- paste0("ProtPC", seq_len(n_pcs))
rownames(prot_pcs) <- ms_cases$eid
cat(sprintf("  PC1-10 variance explained: %.1f%%\n", sum(pct_var[seq_len(n_pcs)])))

# ── 2. INPUT B: Comorbidity features ─────────────────────────────────────────
cat("\n=== INPUT B: Comorbidity ===\n")

# Reload the 24-condition panel (same as analysis 14)
comorbidities <- list(
    M05=list(label="Seropositive RA",           domain="Autoimmune"),
    M06=list(label="Other RA",                  domain="Autoimmune"),
    M32=list(label="SLE",                       domain="Autoimmune"),
    M35=list(label="Connective tissue disease", domain="Autoimmune"),
    K50=list(label="Crohn's disease",           domain="Autoimmune"),
    K51=list(label="Ulcerative colitis",        domain="Autoimmune"),
    L40=list(label="Psoriasis",                 domain="Autoimmune"),
    H20=list(label="Uveitis/iridocyclitis",     domain="Autoimmune"),
    G61=list(label="Inflammatory polyneuropathy",domain="Autoimmune"),
    G70=list(label="Myasthenia gravis",         domain="Autoimmune"),
    G43=list(label="Migraine",                  domain="Neurological"),
    G40=list(label="Epilepsy",                  domain="Neurological"),
    G45=list(label="TIA",                       domain="Neurological"),
    H46=list(label="Optic neuritis",            domain="Neurological"),
    H81=list(label="Vestibular disorder",       domain="Neurological"),
    G62=list(label="Other polyneuropathy",      domain="Neurological"),
    G50=list(label="Trigeminal neuralgia",      domain="Neurological"),
    I10=list(label="Hypertension",              domain="Cardiovascular"),
    I25=list(label="Chronic IHD",               domain="Cardiovascular"),
    I21=list(label="Myocardial infarction",     domain="Cardiovascular"),
    I48=list(label="Atrial fibrillation",       domain="Cardiovascular"),
    I50=list(label="Heart failure",             domain="Cardiovascular"),
    I63=list(label="Cerebral infarction",       domain="Cardiovascular"),
    J45=list(label="Asthma",                    domain="Respiratory"),
    J44=list(label="COPD",                      domain="Respiratory"),
    K21=list(label="GERD",                      domain="Digestive"),
    K58=list(label="IBS",                       domain="Digestive"),
    K57=list(label="Diverticular disease",      domain="Digestive"),
    M54=list(label="Back pain",                 domain="Musculoskeletal"),
    M25=list(label="Joint disorders NEC",       domain="Musculoskeletal"),
    M79=list(label="Soft tissue disorders",     domain="Musculoskeletal"),
    M47=list(label="Spondylosis",               domain="Musculoskeletal"),
    H04=list(label="Lacrimal system disorder",  domain="Other")
)

load_icd_simple <- function(code) {
    ch <- substr(code, 1, 1); gr <- substr(code, 1, 2)
    fp <- file.path(ICD_DIR, ch, gr, code, paste0(code, ".csv"))
    if (!file.exists(fp)) return(NULL)
    dt <- fread(fp, showProgress=FALSE)
    data.table(eid=as.integer(dt[[1]]),
               present=as.integer(dt[[2]]),
               age_at_event=as.numeric(dt[[3]]),
               icd10=code)
}

all_codes <- names(comorbidities)
icd_list  <- lapply(all_codes, load_icd_simple)
icd_dt    <- rbindlist(Filter(Negate(is.null), icd_list), fill=TRUE)

# Wide presence matrix
icd_wide_pres <- dcast(icd_dt[, .(eid, icd10, present)],
                        eid ~ icd10, value.var="present",
                        fun.aggregate=max, fill=0L)
# Wide age matrix (for pre-MS sensitivity)
icd_wide_age  <- dcast(icd_dt[, .(eid, icd10, age_at_event)],
                        eid ~ icd10, value.var="age_at_event",
                        fun.aggregate=min, fill=NA_real_)

ms_eids <- ms_cases$eid
avail   <- setdiff(names(icd_wide_pres), "eid")

# All-timing comorbidity matrix (input B)
comorbid_base <- merge(data.table(eid=ms_eids),
                        icd_wide_pres, by="eid", all.x=TRUE)
for (cc in avail) set(comorbid_base, which(is.na(comorbid_base[[cc]])), cc, 0L)
col_sums <- colSums(comorbid_base[, ..avail])
use_codes <- avail[col_sums >= 3]
comorbid_mat <- as.matrix(comorbid_base[, ..use_codes])
rownames(comorbid_mat) <- comorbid_base$eid
cat(sprintf("  Comorbidity conditions (≥3 MS cases): %d\n", ncol(comorbid_mat)))

# Pre-MS comorbidity matrix (>1yr before diagnosis) — used in sensitivity
ms_ages   <- merge(data.table(eid=ms_eids, age_at_dx=ms_cases$age_at_diagnosis),
                    icd_wide_age[, c("eid", use_codes), with=FALSE],
                    by="eid", all.x=TRUE)
timing_mat <- as.matrix(ms_ages[, ..use_codes]) - ms_ages$age_at_dx
pre_mat    <- ((!is.na(timing_mat)) & (timing_mat < -1)) * comorbid_mat
pre_mat[is.na(pre_mat)] <- 0L

pre_col_sums <- colSums(pre_mat)
pre_codes    <- use_codes[pre_col_sums >= 3]
pre_mat_filt <- pre_mat[, pre_codes, drop=FALSE]
cat(sprintf("  Pre-MS conditions (>1yr before dx, ≥3 cases): %d\n", ncol(pre_mat_filt)))

# ── 3. OUTCOMES (genetics + clinical — never used in clustering) ───────────────
cat("\n=== OUTCOMES (genetics + clinical) ===\n")
hla <- fread(HLA_FILE, showProgress=FALSE)
prs <- fread(PRS_FILE, showProgress=FALSE)
prot_clust <- fread(CLUST_FILE, showProgress=FALSE)
comorbid_dt <- fread(COMORBID_FILE, showProgress=FALSE)

outcomes <- data.table(eid=ms_eids)
outcomes <- merge(outcomes,
    hla[, .(eid, hla_carrier=drb1_1501_carrier, hla_dosage=drb1_1501_dosage)],
    by="eid", all.x=TRUE)
outcomes <- merge(outcomes,
    prs[, .(eid, prs_score)],
    by="eid", all.x=TRUE)
outcomes[, prs_quartile := cut(prs_score,
    quantile(prs_score, 0:4/4, na.rm=TRUE),
    labels=c("Q1","Q2","Q3","Q4"), include.lowest=TRUE)]
outcomes <- merge(outcomes,
    prot_clust[, .(eid, prot_cluster=paste0("PC",cluster))],
    by="eid", all.x=TRUE)
outcomes <- merge(outcomes,
    comorbid_dt[, .(eid, comorbidity_cluster=paste0("CC",comorbidity_cluster))],
    by="eid", all.x=TRUE)
outcomes <- merge(outcomes,
    ms_cases[, .(eid, ms_status, age_at_diagnosis, sex, years_to_diagnosis)],
    by="eid", all.x=TRUE)
outcomes[, hla_label  := ifelse(hla_carrier==1, "HLA+","HLA-")]
outcomes[, onset_label := ifelse(ms_status=="pre_onset","Pre-onset","Post-onset")]

cat(sprintf("  HLA data: %d | PRS data: %d | Proteomic cluster: %d\n",
            sum(!is.na(outcomes$hla_carrier)),
            sum(!is.na(outcomes$prs_score)),
            sum(!is.na(outcomes$prot_cluster))))
cat(sprintf("  HLA carrier rate in MS: %.1f%%\n",
            100*mean(outcomes$hla_carrier, na.rm=TRUE)))

# ── 4. Joint feature matrix (proteomics + comorbidity ONLY) ──────────────────
cat("\nBuilding joint feature matrix (inputs only)...\n")
feat_dt <- data.table(eid=ms_eids)
feat_dt <- merge(feat_dt,
    as.data.table(prot_pcs, keep.rownames="eid")[, eid := as.integer(eid)],
    by="eid", all.x=TRUE)
feat_dt <- merge(feat_dt,
    as.data.table(comorbid_mat, keep.rownames="eid")[, eid := as.integer(eid)],
    by="eid", all.x=TRUE)

feat_mat <- as.matrix(feat_dt[, -"eid"])
feat_mat[is.na(feat_mat)] <- 0
cat(sprintf("  Joint matrix: %d cases × %d features (%d proteomic PCs + %d comorbidity)\n",
            nrow(feat_mat), ncol(feat_mat), n_pcs, ncol(comorbid_mat)))

fwrite(as.data.table(feat_mat)[, eid := ms_eids],
       file.path(OUT_DIR, "ms_endophenotype_features.csv"))

# ── 5. Primary UMAP + clustering ─────────────────────────────────────────────
cat("\nRunning multi-modal UMAP (proteomics + comorbidity inputs)...\n")
umap_res <- umap(feat_mat, n_neighbors=15, min_dist=0.10,
                  metric="euclidean", n_components=2, verbose=FALSE)

umap_dt <- data.table(eid=ms_eids, UMAP1=umap_res[,1], UMAP2=umap_res[,2])

cat("Selecting optimal k...\n")
umap_mat  <- as.matrix(umap_dt[, .(UMAP1, UMAP2)])
umap_dist <- dist(umap_mat)
k_range   <- 2:5
sil_vals  <- numeric(length(k_range))
for (i in seq_along(k_range)) {
    k  <- k_range[i]
    km <- tryCatch(kmeans(umap_mat, centers=k, nstart=100, iter.max=300),
                   error=function(e) NULL)
    if (is.null(km)) next
    sil_vals[i] <- mean(silhouette(km$cluster, umap_dist)[,3])
    cat(sprintf("  k=%d: silhouette=%.3f, sizes=%s\n", k, sil_vals[i],
                paste(table(km$cluster), collapse="/")))
}
k_opt    <- k_range[which.max(sil_vals)]
km_final <- kmeans(umap_mat, centers=k_opt, nstart=200, iter.max=300)
umap_dt[, endophenotype := as.character(km_final$cluster)]
cat(sprintf("  Optimal k = %d\n", k_opt))

# Attach all outcome variables
umap_dt <- merge(umap_dt, outcomes, by="eid", all.x=TRUE)

fwrite(umap_dt, file.path(OUT_DIR, "ms_endophenotype_clusters.csv"))

# ── 6. Sensitivity: re-run with pre-MS comorbidities only ────────────────────
cat("\nSensitivity UMAP: pre-MS comorbidities only...\n")
feat_dt_pre <- data.table(eid=ms_eids)
feat_dt_pre <- merge(feat_dt_pre,
    as.data.table(prot_pcs, keep.rownames="eid")[, eid := as.integer(eid)],
    by="eid", all.x=TRUE)
pre_mat_dt  <- as.data.table(pre_mat_filt)
pre_mat_dt[, eid := as.integer(rownames(pre_mat_filt))]
feat_dt_pre <- merge(feat_dt_pre, pre_mat_dt, by="eid", all.x=TRUE)
feat_mat_pre <- as.matrix(feat_dt_pre[, -"eid"])
feat_mat_pre[is.na(feat_mat_pre)] <- 0
cat(sprintf("  Pre-MS matrix: %d cases × %d features\n",
            nrow(feat_mat_pre), ncol(feat_mat_pre)))

umap_pre <- umap(feat_mat_pre, n_neighbors=15, min_dist=0.10,
                  metric="euclidean", n_components=2, verbose=FALSE)
umap_dt_pre <- data.table(eid=ms_eids, UMAP1=umap_pre[,1], UMAP2=umap_pre[,2])
umap_mat_pre <- as.matrix(umap_dt_pre[, .(UMAP1, UMAP2)])
sil_pre <- sapply(2:4, function(k) {
    km <- tryCatch(kmeans(umap_mat_pre, centers=k, nstart=100, iter.max=300),
                   error=function(e) NULL)
    if (is.null(km)) return(NA)
    mean(silhouette(km$cluster, dist(umap_mat_pre))[,3])
})
k_pre <- (2:4)[which.max(sil_pre)]
km_pre <- kmeans(umap_mat_pre, centers=k_pre, nstart=200, iter.max=300)
umap_dt_pre[, endophenotype_pre := as.character(km_pre$cluster)]
cat(sprintf("  Pre-MS sensitivity: k=%d, silhouette=%.3f\n",
            k_pre, max(sil_pre, na.rm=TRUE)))

# ── 7. Panels A-F: UMAP with outcome overlays ─────────────────────────────────
cat("\nPlotting UMAP panels A-F...\n")

base_umap <- function(dt, col_var, col_map, title_str) {
    ggplot(dt[!is.na(get(col_var))],
           aes(x=UMAP1, y=UMAP2, colour=.data[[col_var]])) +
        geom_point(size=0.9, alpha=0.65, stroke=0) +
        scale_colour_manual(values=col_map, na.value="grey85", name=NULL) +
        labs(x="UMAP 1", y="UMAP 2", title=title_str) +
        theme_ukb(base_size=9) +
        guides(colour=guide_legend(override.aes=list(size=2.5, alpha=1)))
}

pA <- base_umap(umap_dt, "endophenotype",
                ENDO_COLS[seq_len(k_opt)],
                "A  Multi-modal endophenotype")

pB <- base_umap(umap_dt[!is.na(prot_cluster)], "prot_cluster",
                c(PC1=UKB_NAVY, PC2=UKB_PALETTE[4]),
                "B  Proteomic cluster")

pC <- base_umap(umap_dt[!is.na(hla_label)], "hla_label",
                c("HLA+"=UKB_PALETTE[4], "HLA-"="grey65"),
                "C  HLA-DRB1*15:01 status")

pD <- base_umap(umap_dt[!is.na(prs_quartile)], "prs_quartile",
                c(Q1="#56B4E9", Q2="grey65", Q3="#E6A817", Q4="#CC0066"),
                "D  MS PRS quartile")

pE <- base_umap(umap_dt[!is.na(onset_label)], "onset_label",
                c("Pre-onset"=UKB_PALETTE[2], "Post-onset"=UKB_NAVY),
                "E  Pre/post-onset")

pF <- base_umap(umap_dt[!is.na(comorbidity_cluster)], "comorbidity_cluster",
                setNames(ENDO_COLS[1:4], paste0("CC",1:4)),
                "F  Comorbidity cluster")

p_umap_grid <- (pA | pB | pC) / (pD | pE | pF) +
    plot_layout(guides="collect") &
    theme(legend.position="bottom")

ggsave(file.path(OUT_DIR, "ms_endophenotype_umap.pdf"),
       p_umap_grid, width=13, height=9, device=cairo_pdf)
cat("  Saved: ms_endophenotype_umap.pdf\n")

# ── 8. Panel G: feature profile per endophenotype ────────────────────────────
cat("Building characterisation heatmap...\n")

# Key features: top 5 proteomic PCs + top comorbidities by variance
top_comorb <- use_codes[order(-apply(comorbid_mat, 2, var))][seq_len(min(12, length(use_codes)))]
key_feats  <- c(paste0("ProtPC", 1:5), top_comorb)
key_feats  <- intersect(key_feats, names(feat_dt))

char_long <- melt(
    merge(feat_dt[, c("eid", key_feats), with=FALSE],
          umap_dt[, .(eid, endophenotype)], by="eid"),
    id.vars=c("eid","endophenotype"),
    variable.name="feature", value.name="value"
)

cond_label_map <- sapply(top_comorb, function(cc)
    if (cc %in% names(comorbidities)) comorbidities[[cc]]$label else cc)
pc_label_map   <- setNames(paste0("Prot-PC", 1:5), paste0("ProtPC", 1:5))
full_label_map <- c(pc_label_map, cond_label_map)

char_mean <- char_long[, .(mean_val=mean(value, na.rm=TRUE)),
                        by=.(endophenotype, feature)]
char_mean[, feat_label := full_label_map[as.character(feature)]]
char_mean[, feat_label := ifelse(is.na(feat_label), as.character(feature), feat_label)]
char_mean[, feat_type  := ifelse(grepl("^ProtPC", feature), "Proteomic PCs", "Comorbidity")]
char_mean[, feat_type_f := factor(feat_type, levels=c("Proteomic PCs","Comorbidity"))]

p_char <- ggplot(char_mean, aes(x=endophenotype, y=feat_label, fill=mean_val)) +
    geom_tile(colour="white", linewidth=0.3) +
    facet_grid(feat_type_f ~ ., scales="free_y", space="free_y") +
    scale_fill_gradient2(low="#56B4E9", mid="white", high="#CC0066",
                         midpoint=0, name="Mean\n(z-score)") +
    labs(x="Endophenotype cluster", y=NULL,
         title="G  Molecular + comorbidity profile per endophenotype") +
    theme_ukb(base_size=9) +
    theme(axis.text.y=element_text(size=7),
          strip.text.y=element_text(angle=0, size=7, face="bold"))

# ── 9. Panel H: clinical validation ──────────────────────────────────────────
cat("Building clinical validation...\n")

make_violin <- function(y_var, y_lab, title_letter) {
    ggplot(umap_dt[!is.na(endophenotype) & !is.na(get(y_var))],
           aes(x=endophenotype, y=get(y_var), fill=endophenotype)) +
        geom_violin(alpha=0.5, linewidth=0.3) +
        geom_boxplot(width=0.15, fill="white",
                     outlier.size=0.4, linewidth=0.35, colour="grey30") +
        scale_fill_manual(values=ENDO_COLS[seq_len(k_opt)], guide="none") +
        labs(x="Endophenotype", y=y_lab, title=title_letter) +
        theme_ukb(base_size=9)
}
pH1 <- make_violin("age_at_diagnosis", "Age at diagnosis (yr)", "H1  Age at diagnosis")
pH2 <- make_violin("years_to_diagnosis", "Years to diagnosis", "H2  Prodromal window")
pH3 <- ggplot(umap_dt[!is.na(endophenotype)],
              aes(x=endophenotype, fill=onset_label)) +
    geom_bar(position="fill", width=0.7) +
    scale_fill_manual(values=c("Pre-onset"=UKB_PALETTE[2],"Post-onset"=UKB_NAVY),
                       name=NULL) +
    scale_y_continuous(labels=scales::percent_format()) +
    labs(x="Endophenotype", y="% of cases", title="H3  Pre/post onset") +
    theme_ukb(base_size=9)

p_charact <- (p_char | (pH1 / pH2 / pH3)) + plot_layout(widths=c(1.4, 1))
ggsave(file.path(OUT_DIR, "ms_endophenotype_characterisation.pdf"),
       p_charact, width=14, height=8, device=cairo_pdf)
cat("  Saved: ms_endophenotype_characterisation.pdf\n")

# ── 10. Outcome association tests ─────────────────────────────────────────────
cat("\n=== Outcome association tests ===\n")

clust_summary <- umap_dt[!is.na(endophenotype) & !is.na(age_at_diagnosis), .(
    n              = .N,
    pct_pre_onset  = round(100*mean(ms_status=="pre_onset", na.rm=TRUE), 1),
    pct_female     = round(100*mean(sex==0, na.rm=TRUE), 1),
    mean_age_dx    = round(mean(age_at_diagnosis, na.rm=TRUE), 1),
    mean_years_dx  = round(mean(years_to_diagnosis, na.rm=TRUE), 1),
    pct_hla        = round(100*mean(hla_carrier==1, na.rm=TRUE), 1),
    mean_prs       = round(mean(prs_score, na.rm=TRUE), 3)
), by=endophenotype][order(endophenotype)]
cat("\nCluster summary:\n"); print(clust_summary)
fwrite(clust_summary, file.path(OUT_DIR, "ms_endophenotype_clinical.csv"))

test_outcome <- function(var, label, type="chi") {
    sub <- umap_dt[!is.na(endophenotype) & !is.na(get(var))]
    if (type=="chi") {
        xt <- table(sub$endophenotype, sub[[var]])
        if (nrow(xt)>1 && ncol(xt)>1) {
            p <- tryCatch(chisq.test(xt)$p.value, warning=function(w) chisq.test(xt, simulate.p.value=TRUE)$p.value)
            cat(sprintf("  %-35s χ² p=%.3g\n", label, p))
        }
    } else {
        aov_f <- aov(get(var) ~ endophenotype, data=sub)
        p <- summary(aov_f)[[1]][["Pr(>F)"]][1]
        cat(sprintf("  %-35s ANOVA p=%.3g\n", label, p))
    }
}
test_outcome("hla_carrier",       "HLA-DRB1*15:01 carrier",    "chi")
test_outcome("prs_quartile",      "MS PRS quartile",            "chi")
test_outcome("prot_cluster",      "Proteomic cluster",          "chi")
test_outcome("comorbidity_cluster","Comorbidity cluster",       "chi")
test_outcome("onset_label",       "Pre/post-onset",             "chi")
test_outcome("prs_score",         "PRS (continuous)",           "anova")
test_outcome("age_at_diagnosis",  "Age at diagnosis",           "anova")
test_outcome("years_to_diagnosis","Years to diagnosis",         "anova")
test_outcome("sex",               "Sex",                        "chi")

cat("\n01_ms_multimodal_endophenotype.R complete.\n")
