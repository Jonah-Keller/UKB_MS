#!/usr/bin/env Rscript
# 03_ms_genetic_stratification.R  — Figure 5
#
# DESIGN: HLA × MS interaction proteomics + stratified analyses
#
#   Core analysis: full interaction model in the matched cohort (MS + HC)
#     protein ~ ms_status * hla_carrier + age_at_sampling + sex
#   Estimates:
#     ms_effect    — MS vs HC in HLA- background (reference)
#     hla_effect   — HLA+ vs HLA- in HC (main effect in controls)
#     interaction  — does the MS proteomic effect differ by HLA carrier status?
#   Derived contrasts via limma:
#     HLA- MS vs HLA- HC  (pure MS effect without HLA modulation)
#     HLA+ MS vs HLA+ HC  (MS effect in HLA+ background)
#
#   Supplementary analyses:
#     B. HLA+ vs HLA- within MS (within-disease stratification)
#     C. Prodromal trajectory × HLA (pre-onset cases, spline × carrier interaction)
#     D. PRS Q4 vs Q1 within MS (HLA-adjusted)
#
#   Panels:
#     A — Interaction volcano: proteins where MS effect is HLA-modified
#     B — LogFC scatter: HLA- MS-effect vs HLA+ MS-effect (per protein)
#     C — Trajectory: NEFL/MOG/GFAP predicted means ± HLA (pre-onset only)
#     D — PRS Q4 vs Q1 volcano (HLA-adjusted, within MS)
#     E — HLA vs PRS logFC concordance scatter

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(limma)
    library(splines)
})

args       <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))

QC_FILE  <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
HLA_FILE <- file.path(PROJ_DIR, "data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE <- file.path(PROJ_DIR, "data", "ukb", "genetics", "ms_prs_scores.csv")
OUT_DIR  <- file.path(PROJ_DIR, "results", "endophenotype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_THR  <- 0.05
N_MATCH  <- 3   # HC per MS case in matched cohort

# ── 1. Load cohort ────────────────────────────────────────────────────────────
cat("Loading data...\n")
qc  <- fread(QC_FILE,  showProgress=FALSE)
hla <- fread(HLA_FILE, showProgress=FALSE)
prs <- fread(PRS_FILE, showProgress=FALSE)

meta_cols <- c("eid","ms_status","age_at_sampling","age_at_diagnosis",
               "years_to_diagnosis","sex","olink_instance","qc_outlier",
               "UMAP1","UMAP2","mean_npx","sd_npx","lower_iqr","upper_iqr","iqr_npx")

ms_all <- qc[qc_outlier==FALSE &
             ms_status %in% c("pre_onset","post_onset") &
             !is.na(age_at_diagnosis)]
hc_all <- qc[qc_outlier==FALSE & ms_status=="control"]

prot_cols <- setdiff(names(qc), meta_cols)

# Attach HLA to everyone
ms_all <- merge(ms_all, hla[, .(eid, hla_carrier=drb1_1501_carrier)], by="eid", all.x=TRUE)
hc_all <- merge(hc_all, hla[, .(eid, hla_carrier=drb1_1501_carrier)], by="eid", all.x=TRUE)

ms_all <- merge(ms_all, prs[, .(eid, prs_score)], by="eid", all.x=TRUE)

cat(sprintf("  MS cases: %d (HLA+=%.0f, HLA-=%.0f)\n",
    nrow(ms_all), sum(ms_all$hla_carrier==1, na.rm=TRUE),
    sum(ms_all$hla_carrier==0, na.rm=TRUE)))
cat(sprintf("  HC available: %d (HLA data: %d)\n",
    nrow(hc_all), sum(!is.na(hc_all$hla_carrier))))

# ── 1b. Control-only PC1/PC2 (project all participants into control PC space) ─
# Control-derived PCs capture technical/batch variance without disease signal
# contamination. All participants (MS + HC) are projected into this space.
cat("Computing control-only PC1/PC2...\n")
ctrl_prot_mat  <- t(as.matrix(hc_all[, ..prot_cols]))
colnames(ctrl_prot_mat) <- hc_all$eid
ctrl_row_means <- rowMeans(ctrl_prot_mat, na.rm = TRUE)
for (i in seq_len(nrow(ctrl_prot_mat))) {
    nas <- is.na(ctrl_prot_mat[i, ])
    if (any(nas)) ctrl_prot_mat[i, nas] <- ctrl_row_means[i]
}
pc_fit_gs  <- prcomp(t(ctrl_prot_mat), center = TRUE, scale. = FALSE)
pc_var_gs  <- 100 * pc_fit_gs$sdev^2 / sum(pc_fit_gs$sdev^2)
cat(sprintf("  Control PC1=%.1f%%  PC2=%.1f%%\n", pc_var_gs[1], pc_var_gs[2]))

all_qc_gs     <- qc[qc_outlier == FALSE]
all_prot_mat  <- t(as.matrix(all_qc_gs[, ..prot_cols]))
colnames(all_prot_mat) <- all_qc_gs$eid
all_row_means <- rowMeans(all_prot_mat, na.rm = TRUE)
for (i in seq_len(nrow(all_prot_mat))) {
    nas <- is.na(all_prot_mat[i, ])
    if (any(nas)) all_prot_mat[i, nas] <- all_row_means[i]
}
all_pcs_gs <- predict(pc_fit_gs, newdata = t(all_prot_mat))
pc_dt_gs   <- data.table(eid = as.integer(colnames(all_prot_mat)),
                          PC1 = all_pcs_gs[, 1], PC2 = all_pcs_gs[, 2])
ms_all  <- merge(ms_all,  pc_dt_gs, by = "eid", all.x = TRUE)
hc_all  <- merge(hc_all,  pc_dt_gs, by = "eid", all.x = TRUE)
cat(sprintf("  PCs merged: %d MS, %d HC with valid PCs\n",
    sum(!is.na(ms_all$PC1)), sum(!is.na(hc_all$PC1))))

# ── 2. Build matched cohort for interaction model ─────────────────────────────
# Age+sex match: for each MS case, sample N_MATCH controls from HCs with
# |age_at_sampling - hc_age| <= 2 and same sex. Use with-replacement if needed.
cat(sprintf("Building %d:1 age+sex matched HC cohort...\n", N_MATCH))
set.seed(42)
matched_hc_ids <- unlist(lapply(seq_len(nrow(ms_all)), function(i) {
    row   <- ms_all[i]
    pool  <- hc_all[abs(age_at_sampling - row$age_at_sampling) <= 2 &
                    sex == row$sex, eid]
    if (length(pool) == 0)
        pool <- hc_all[abs(age_at_sampling - row$age_at_sampling) <= 5 &
                       sex == row$sex, eid]
    if (length(pool) == 0) return(integer(0))
    sample(pool, min(N_MATCH, length(pool)), replace=length(pool) < N_MATCH)
}))
matched_hc <- hc_all[eid %in% matched_hc_ids]
matched_hc <- matched_hc[!duplicated(eid)]  # keep unique HCs

full_cohort <- rbind(
    ms_all[, c(meta_cols[1:6], "hla_carrier", "PC1", "PC2", prot_cols), with=FALSE],
    matched_hc[, c(meta_cols[1:6], "hla_carrier", "PC1", "PC2", prot_cols), with=FALSE],
    fill=TRUE
)
full_cohort[, ms_case := as.integer(ms_status %in% c("pre_onset","post_onset"))]

cat(sprintf("  Matched cohort: %d MS + %d HC = %d total\n",
    nrow(ms_all), nrow(matched_hc), nrow(full_cohort)))
cat(sprintf("  HLA coverage: MS=%.0f%%, HC=%.0f%%\n",
    100*mean(!is.na(ms_all$hla_carrier)),
    100*mean(!is.na(matched_hc$hla_carrier))))

# ── Helper: run_limma ─────────────────────────────────────────────────────────
run_limma <- function(dat, design_mat, coef_name) {
    expr <- t(as.matrix(dat[, ..prot_cols]))
    colnames(expr) <- dat$eid
    row_means <- rowMeans(expr, na.rm=TRUE)
    for (i in seq_len(nrow(expr)))
        expr[i, is.na(expr[i,])] <- row_means[i]
    fit <- eBayes(lmFit(expr, design_mat), trend=TRUE)
    res <- as.data.table(topTable(fit, coef=coef_name, n=Inf, sort.by="P"))
    res[, protein := rownames(topTable(fit, coef=coef_name, n=Inf, sort.by="P"))]
    res[, fdr := p.adjust(P.Value, method="BH")]
    res
}

# ── 3. Part A: HLA × MS interaction model ────────────────────────────────────
cat("\n=== Part A: HLA × MS interaction model ===\n")

int_dat <- full_cohort[!is.na(hla_carrier) & !is.na(age_at_sampling) & !is.na(sex)]
# Explicit interaction term (avoids : in coefficient names, which breaks makeContrasts)
# Coefficients:
#   Intercept     = HLA- HC mean
#   ms_case       = MS effect in HLA- background
#   hla_carrier   = HLA+ vs HLA- in HC
#   ms_hla_int    = INTERACTION: extra MS effect conferred by HLA+ carrier status
#   age_at_sampling, sex = covariates
int_dat[, ms_hla_int  := ms_case * hla_carrier]
# ytd_design: years_to_diagnosis (NA→0 for HC); adjusts for disease phase/duration
int_dat[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
design_int <- model.matrix(~ ms_case + hla_carrier + ms_hla_int +
                               ytd_design + age_at_sampling + sex + PC1 + PC2,
                            data = int_dat)

# Run limma for each coefficient of interest
cat("  Running interaction limma (n=", nrow(int_dat), ")...\n", sep="")
expr_int <- t(as.matrix(int_dat[, ..prot_cols]))
colnames(expr_int) <- int_dat$eid
row_means <- rowMeans(expr_int, na.rm=TRUE)
for (i in seq_len(nrow(expr_int)))
    expr_int[i, is.na(expr_int[i,])] <- row_means[i]
fit_int <- eBayes(lmFit(expr_int, design_int), trend=TRUE)

# Interaction term
dep_interaction <- as.data.table(topTable(fit_int, coef="ms_hla_int", n=Inf, sort.by="P"))
dep_interaction[, protein := rownames(topTable(fit_int, coef="ms_hla_int", n=Inf, sort.by="P"))]
dep_interaction[, fdr     := p.adjust(P.Value, method="BH")]
n_int <- sum(dep_interaction$fdr < FDR_THR)
cat(sprintf("  Interaction DEPs (FDR<0.05): %d\n", n_int))

# MS effect in HLA- background
dep_ms_hlaneg <- as.data.table(topTable(fit_int, coef="ms_case", n=Inf, sort.by="P"))
dep_ms_hlaneg[, protein := rownames(topTable(fit_int, coef="ms_case", n=Inf, sort.by="P"))]
dep_ms_hlaneg[, fdr := p.adjust(P.Value, method="BH")]
cat(sprintf("  MS effect (HLA- background) DEPs: %d\n", sum(dep_ms_hlaneg$fdr < FDR_THR)))

# MS effect in HLA+ background = ms_case + ms_hla_int (contrast)
cont_mat <- makeContrasts("ms_case + ms_hla_int", levels=design_int)
fit_cont <- contrasts.fit(fit_int, cont_mat)
fit_cont <- eBayes(fit_cont, trend=TRUE)
dep_ms_hlapos <- as.data.table(topTable(fit_cont, coef=1, n=Inf, sort.by="P"))
dep_ms_hlapos[, protein := rownames(topTable(fit_cont, coef=1, n=Inf, sort.by="P"))]
dep_ms_hlapos[, fdr := p.adjust(P.Value, method="BH")]
cat(sprintf("  MS effect (HLA+ background) DEPs: %d\n", sum(dep_ms_hlapos$fdr < FDR_THR)))

# Merge for Panel B scatter
scatter_dt <- merge(
    dep_ms_hlaneg[, .(protein, logFC_hlaneg=logFC, p_hlaneg=P.Value, fdr_hlaneg=fdr)],
    dep_ms_hlapos[, .(protein, logFC_hlapos=logFC, p_hlapos=P.Value, fdr_hlapos=fdr)],
    by="protein"
)
scatter_dt <- merge(scatter_dt,
    dep_interaction[, .(protein, logFC_int=logFC, fdr_int=fdr)], by="protein")
scatter_dt[, sig_int := fdr_int < FDR_THR]
cat(sprintf("  Proteins with discordant MS effect (sig interaction): %d\n",
            sum(scatter_dt$sig_int)))

# Save
fwrite(dep_interaction, file.path(OUT_DIR, "ms_hla_interaction_deps.csv"))
fwrite(scatter_dt,      file.path(OUT_DIR, "ms_hla_stratified_logfc.csv"))

# ── 4. Part B: HLA+ vs HLA- within MS only ───────────────────────────────────
cat("\n=== Part B: HLA+ vs HLA- within MS ===\n")
ms_hla <- ms_all[!is.na(hla_carrier)]
ms_hla[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
design_hla <- model.matrix(~ hla_carrier + ytd_design + age_at_sampling + sex + PC1 + PC2,
                             data = as.data.frame(ms_hla[, .(
                                 hla_carrier     = factor(hla_carrier),
                                 ytd_design, age_at_sampling, sex, PC1, PC2)]))
dep_hla <- run_limma(ms_hla, design_hla, "hla_carrier1")
cat(sprintf("  HLA+ vs HLA- DEPs (FDR<0.05): %d\n", sum(dep_hla$fdr < FDR_THR)))
fwrite(dep_hla, file.path(OUT_DIR, "ms_hla_dep_results.csv"))

# ── 5. Part C: Prodromal trajectory × HLA ────────────────────────────────────
cat("\n=== Part C: Prodromal trajectory × HLA (pre-onset only) ===\n")
pre_onset <- ms_all[ms_status=="pre_onset" & !is.na(hla_carrier)]
cat(sprintf("  Pre-onset cases: %d (HLA+=%.0f, HLA-=%.0f)\n",
            nrow(pre_onset), sum(pre_onset$hla_carrier==1), sum(pre_onset$hla_carrier==0)))

target_prots <- intersect(c("nefl","mog","gfap"), prot_cols)
pre_onset[, yrs_before_onset := -years_to_diagnosis]

traj_results <- rbindlist(lapply(target_prots, function(prot) {
    sub <- pre_onset[!is.na(get(prot)) & !is.na(yrs_before_onset)]
    if (nrow(sub) < 30) return(NULL)
    # Interaction: does slope differ by HLA?
    fit_full <- lm(get(prot) ~ ns(yrs_before_onset, df=3) * factor(hla_carrier) +
                   age_at_sampling + sex, data=sub)
    fit_base <- lm(get(prot) ~ ns(yrs_before_onset, df=3) + factor(hla_carrier) +
                   age_at_sampling + sex, data=sub)
    int_p    <- anova(fit_base, fit_full)$"Pr(>F)"[2]

    # Predicted trajectories
    grid <- CJ(yrs_before_onset=seq(0, max(sub$yrs_before_onset, na.rm=TRUE), length.out=50),
               hla_carrier=c(0,1),
               age_at_sampling=median(sub$age_at_sampling, na.rm=TRUE),
               sex=0)
    pred_df <- data.frame(
        yrs_before_onset = grid$yrs_before_onset,
        hla_carrier      = grid$hla_carrier,
        age_at_sampling  = grid$age_at_sampling,
        sex              = grid$sex)
    pred_df[["factor(hla_carrier)"]] <- factor(pred_df$hla_carrier, levels=c(0,1))
    pred_vals <- tryCatch(
        predict(fit_full, newdata=pred_df),
        error=function(e) predict(fit_base, newdata=pred_df))
    grid[, npx_pred := pred_vals]
    grid[, protein := prot]
    grid[, interaction_p := int_p]
    grid
}))

cat(sprintf("  Trajectory proteins fitted: %s\n",
            paste(unique(traj_results$protein), collapse=", ")))
for (prot in unique(traj_results$protein)) {
    p <- unique(traj_results[protein==prot, interaction_p])
    cat(sprintf("    %s: HLA×time interaction p=%.3g\n", toupper(prot), p))
}

# ── 6. Part D: PRS Q4 vs Q1 within MS (HLA-adjusted) ────────────────────────
cat("\n=== Part D: PRS Q4 vs Q1 (HLA-adjusted) ===\n")
ms_prs <- ms_all[!is.na(prs_score) & !is.na(hla_carrier)]
ms_prs[, prs_q := cut(prs_score, quantile(prs_score, c(0,.25,.75,1)), include.lowest=TRUE,
                        labels=c("Q1","mid","Q4"))]
prs_sub <- ms_prs[prs_q %in% c("Q1","Q4")]
cat(sprintf("  PRS Q4: %d | Q1: %d\n", sum(prs_sub$prs_q=="Q4"), sum(prs_sub$prs_q=="Q1")))
design_prs <- model.matrix(~ prs_q + age_at_sampling + sex + hla_carrier + PC1 + PC2,
                             data = as.data.frame(prs_sub[, .(
                                 prs_q = factor(prs_q, levels = c("Q1","Q4")),
                                 age_at_sampling, sex, hla_carrier, PC1, PC2)]))
dep_prs <- run_limma(prs_sub, design_prs, "prs_qQ4")
cat(sprintf("  PRS Q4 vs Q1 DEPs (FDR<0.05): %d\n", sum(dep_prs$fdr < FDR_THR)))
fwrite(dep_prs, file.path(OUT_DIR, "ms_prs_dep_results.csv"))

# ── 7. Part E: PRS × MS proteome interaction (genome-wide) ───────────────────
# Does higher MS PRS amplify the MS vs HC protein signature?
# Model: protein ~ ms_case + prs_z + ms_case:prs_z + ytd + age + sex + PC1 + PC2
# Interaction coef ms_case:prs_z answers: per 1SD PRS, how much does the MS–HC
# gap change? Run in the full matched cohort (MS + matched HC) with PRS merged.
cat("\n=== Part E: PRS × MS proteome interaction (genome-wide) ===\n")

fc_prs <- merge(full_cohort, prs[, .(eid, prs_score)], by = "eid", all.x = TRUE)
fc_prs <- fc_prs[!is.na(prs_score) & !is.na(age_at_sampling) & !is.na(sex) & !is.na(PC1)]
fc_prs[, prs_z      := as.numeric(scale(prs_score))]
fc_prs[, ytd_fc     := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
cat(sprintf("  PRS×MS cohort: %d MS + %d HC\n",
    sum(fc_prs$ms_case == 1), sum(fc_prs$ms_case == 0)))

design_prs_int <- model.matrix(
    ~ ms_case + prs_z + ms_case:prs_z + ytd_fc + age_at_sampling + sex + PC1 + PC2,
    data = as.data.frame(fc_prs)
)
colnames(design_prs_int) <- make.names(colnames(design_prs_int))
prs_int_coef <- grep("ms_case.*prs_z|prs_z.*ms_case", colnames(design_prs_int), value = TRUE)
cat(sprintf("  Interaction coefficient: %s\n", prs_int_coef))

expr_prs_int <- t(as.matrix(fc_prs[, ..prot_cols]))
colnames(expr_prs_int) <- fc_prs$eid
rm_pi <- rowMeans(expr_prs_int, na.rm = TRUE)
for (i in seq_len(nrow(expr_prs_int)))
    expr_prs_int[i, is.na(expr_prs_int[i, ])] <- rm_pi[i]

fit_prs_int     <- eBayes(lmFit(expr_prs_int, design_prs_int), trend = TRUE)
dep_prs_int     <- as.data.table(topTable(fit_prs_int, coef = prs_int_coef, n = Inf, sort.by = "P"))
dep_prs_int[, protein    := rownames(topTable(fit_prs_int, coef = prs_int_coef, n = Inf, sort.by = "P"))]
dep_prs_int[, fdr        := p.adjust(P.Value, method = "BH")]
dep_prs_int[, sex_group  := "All"]
n_prs_int <- sum(dep_prs_int$fdr < FDR_THR)
cat(sprintf("  PRS×MS interaction FDR<0.05: %d | nominal p<0.05: %d\n",
    n_prs_int, sum(dep_prs_int$P.Value < 0.05)))
cat("  Top 10 by p:\n")
print(dep_prs_int[order(P.Value)][1:10, .(protein, logFC, P.Value, fdr)])
fwrite(dep_prs_int, file.path(OUT_DIR, "ms_prs_interaction_genome.csv"))

# ── 8. Part F: Sex-stratified HLA and PRS×MS analyses ────────────────────────
# Run HLA×MS interaction, within-MS HLA+ vs HLA−, and PRS×MS interaction
# separately in females (sex=0) and males (sex=1) to test for sex-specific effects.
cat("\n=== Part F: Sex-stratified genetic analyses ===\n")

run_sex_stratum <- function(sex_val, sex_label) {
    cat(sprintf("\n--- %s (sex=%d) ---\n", sex_label, sex_val))
    out <- list()

    # F1: HLA×MS interaction
    fc_s <- full_cohort[sex == sex_val & !is.na(hla_carrier) &
                        !is.na(age_at_sampling) & !is.na(PC1)]
    fc_s[, ms_hla_int := ms_case * hla_carrier]
    fc_s[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
    cat(sprintf("  HLA×MS cohort: %d MS + %d HC\n",
        sum(fc_s$ms_case == 1), sum(fc_s$ms_case == 0)))

    if (sum(fc_s$ms_case == 1) >= 20) {
        design_s <- model.matrix(
            ~ ms_case + hla_carrier + ms_hla_int + ytd_design + age_at_sampling + PC1 + PC2,
            data = as.data.frame(fc_s))
        expr_s  <- t(as.matrix(fc_s[, ..prot_cols]))
        colnames(expr_s) <- fc_s$eid
        rm_s <- rowMeans(expr_s, na.rm = TRUE)
        for (i in seq_len(nrow(expr_s))) expr_s[i, is.na(expr_s[i,])] <- rm_s[i]
        fit_s   <- eBayes(lmFit(expr_s, design_s), trend = TRUE)
        int_s   <- as.data.table(topTable(fit_s, coef = "ms_hla_int", n = Inf, sort.by = "P"))
        int_s[, protein   := rownames(topTable(fit_s, coef = "ms_hla_int", n = Inf, sort.by = "P"))]
        int_s[, fdr       := p.adjust(P.Value, method = "BH")]
        int_s[, sex_group := sex_label]
        int_s[, analysis  := "hla_ms_interaction"]
        cat(sprintf("  HLA×MS FDR<0.05: %d | nominal p<0.05: %d\n",
            sum(int_s$fdr < FDR_THR), sum(int_s$P.Value < 0.05)))
        out$hla_int <- int_s
    }

    # F2: Within-MS HLA+ vs HLA−
    ms_s <- ms_all[sex == sex_val & !is.na(hla_carrier) & !is.na(PC1)]
    if (nrow(ms_s) >= 20) {
        ms_s[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
        design_hs <- model.matrix(
            ~ hla_carrier + ytd_design + age_at_sampling + PC1 + PC2,
            data = as.data.frame(ms_s[, .(hla_carrier = factor(hla_carrier),
                                          ytd_design, age_at_sampling, PC1, PC2)]))
        hla_s     <- run_limma(ms_s, design_hs, "hla_carrier1")
        hla_s[, sex_group := sex_label]
        hla_s[, analysis  := "hla_within_ms"]
        cat(sprintf("  Within-MS HLA+ vs HLA− FDR<0.05: %d | nominal p<0.05: %d\n",
            sum(hla_s$fdr < FDR_THR), sum(hla_s$P.Value < 0.05)))
        out$hla_within <- hla_s
    }

    # F3: PRS×MS interaction (sex-stratified)
    fc_ps <- merge(full_cohort[sex == sex_val], prs[, .(eid, prs_score)], by = "eid", all.x = TRUE)
    fc_ps <- fc_ps[!is.na(prs_score) & !is.na(age_at_sampling) & !is.na(PC1)]
    cat(sprintf("  PRS×MS cohort: %d MS + %d HC\n",
        sum(fc_ps$ms_case == 1), sum(fc_ps$ms_case == 0)))

    if (sum(fc_ps$ms_case == 1) >= 20) {
        fc_ps[, prs_z      := as.numeric(scale(prs_score))]
        fc_ps[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
        design_ps <- model.matrix(
            ~ ms_case + prs_z + ms_case:prs_z + ytd_design + age_at_sampling + PC1 + PC2,
            data = as.data.frame(fc_ps))
        colnames(design_ps) <- make.names(colnames(design_ps))
        prs_ic_s <- grep("ms_case.*prs_z|prs_z.*ms_case", colnames(design_ps), value = TRUE)
        expr_ps  <- t(as.matrix(fc_ps[, ..prot_cols]))
        colnames(expr_ps) <- fc_ps$eid
        rm_ps <- rowMeans(expr_ps, na.rm = TRUE)
        for (i in seq_len(nrow(expr_ps))) expr_ps[i, is.na(expr_ps[i,])] <- rm_ps[i]
        fit_ps  <- eBayes(lmFit(expr_ps, design_ps), trend = TRUE)
        prs_s   <- as.data.table(topTable(fit_ps, coef = prs_ic_s, n = Inf, sort.by = "P"))
        prs_s[, protein   := rownames(topTable(fit_ps, coef = prs_ic_s, n = Inf, sort.by = "P"))]
        prs_s[, fdr       := p.adjust(P.Value, method = "BH")]
        prs_s[, sex_group := sex_label]
        prs_s[, analysis  := "prs_ms_interaction"]
        cat(sprintf("  PRS×MS FDR<0.05: %d | nominal p<0.05: %d\n",
            sum(prs_s$fdr < FDR_THR), sum(prs_s$P.Value < 0.05)))
        out$prs_int <- prs_s
    }

    out
}

female_res <- run_sex_stratum(0, "Female")
male_res   <- run_sex_stratum(1, "Male")

sex_all <- rbindlist(c(
    list(female_res$hla_int,  female_res$hla_within,  female_res$prs_int,
         male_res$hla_int,    male_res$hla_within,    male_res$prs_int)
), fill = TRUE)

if (nrow(sex_all) > 0) {
    fwrite(sex_all, file.path(OUT_DIR, "ms_sex_stratified_genetic.csv"))
    cat("\nSaved: ms_sex_stratified_genetic.csv\n")
    # Summary by analysis + sex
    cat("\nNominal hits (p<0.05) by analysis + sex:\n")
    print(sex_all[, .(nominal_p05 = sum(P.Value < 0.05, na.rm=TRUE),
                      fdr05       = sum(fdr < FDR_THR, na.rm=TRUE)),
                  by = .(analysis, sex_group)])
}

# ── 9. Assemble figure ────────────────────────────────────────────────────────
cat("\nAssembling figure...\n")

# Panel A: Interaction volcano
dep_interaction[, direction := fcase(
    fdr < FDR_THR & logFC > 0, "Stronger in HLA+",
    fdr < FDR_THR & logFC < 0, "Weaker in HLA+",
    default="NS")]
dep_interaction[, label := ifelse(fdr < 0.01, toupper(protein), NA_character_)]
pA <- ggplot(dep_interaction, aes(logFC, -log10(P.Value), colour=direction, label=label)) +
    geom_point(size=0.8, alpha=0.6) +
    geom_text_repel(size=2.2, max.overlaps=20, segment.size=0.2, na.rm=TRUE) +
    scale_colour_manual(values=c("Stronger in HLA+"="#E6A817",
                                  "Weaker in HLA+"  ="#56B4E9",
                                  "NS"="grey65"), guide="none") +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey50") +
    labs(x=expression(Interaction~log[2]~FC),
         y=expression(-log[10](p)),
         title=sprintf("A  HLA × MS interaction\n(FDR<0.05: %d proteins)", n_int)) +
    theme_ukb(base_size=9)

# Panel B: HLA- MS-effect vs HLA+ MS-effect scatter
scatter_dt[, label_pt := ifelse(sig_int, toupper(protein), NA_character_)]
pB <- ggplot(scatter_dt, aes(logFC_hlaneg, logFC_hlapos,
                               colour=sig_int, label=label_pt)) +
    geom_point(size=0.8, alpha=0.6) +
    geom_text_repel(size=2.2, max.overlaps=15, segment.size=0.2, na.rm=TRUE) +
    geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=0.4) +
    geom_hline(yintercept=0, linewidth=0.3, colour="grey70") +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
    scale_colour_manual(values=c("TRUE"="#E6A817","FALSE"="grey65"),
                        labels=c("TRUE"="Sig. interaction","FALSE"="NS"),
                        name=NULL) +
    labs(x=expression(MS~log[2]~FC~(HLA^"-"~background)),
         y=expression(MS~log[2]~FC~(HLA^"+"~background)),
         title="B  MS proteomic effect: HLA- vs HLA+ background") +
    theme_ukb(base_size=9)

# Panel C: Trajectory plots
if (!is.null(traj_results) && nrow(traj_results) > 0) {
    prot_labels <- c(nefl="NEFL", mog="MOG", gfap="GFAP")
    traj_results[, prot_label := prot_labels[protein]]
    traj_results[, hla_label  := ifelse(hla_carrier==1, "HLA+", "HLA-")]
    # Add raw data points
    raw_long <- melt(pre_onset[, c("eid","yrs_before_onset","hla_carrier",
                                    intersect(target_prots, names(pre_onset))), with=FALSE],
                     id.vars=c("eid","yrs_before_onset","hla_carrier"),
                     variable.name="protein", value.name="npx")
    raw_long <- raw_long[!is.na(npx)]
    raw_long[, hla_label  := ifelse(hla_carrier==1, "HLA+", "HLA-")]
    raw_long[, prot_label := prot_labels[as.character(protein)]]

    pC <- ggplot() +
        geom_point(data=raw_long,
                   aes(yrs_before_onset, npx, colour=hla_label),
                   size=0.6, alpha=0.35) +
        geom_line(data=traj_results,
                  aes(yrs_before_onset, npx_pred, colour=hla_label),
                  linewidth=1) +
        facet_wrap(~prot_label, scales="free_y", nrow=1) +
        scale_colour_manual(values=c("HLA+"="#E6A817","HLA-"="#2B4C7E"), name="HLA status") +
        geom_text(data=traj_results[, .(interaction_p=unique(interaction_p),
                                         prot_label=unique(prot_label)), by=protein],
                  aes(x=Inf, y=Inf,
                      label=sprintf("Int. p=%.2g", interaction_p)),
                  hjust=1.05, vjust=1.3, size=2.5, colour="#1A1A1A") +
        labs(x="Years before onset", y="NPX",
             title="C  Prodromal trajectories by HLA status (pre-onset cases)") +
        theme_ukb(base_size=9) +
        theme(strip.text=element_text(face="bold"))
} else {
    pC <- ggplot(data.frame(x=1,y=1,l="No trajectory data")) +
        geom_text(aes(x,y,label=l)) + theme_ukb() + labs(title="C")
}

# Panel D: PRS volcano
dep_prs[, direction := fcase(fdr < FDR_THR & logFC > 0, "Up in Q4",
                               fdr < FDR_THR & logFC < 0, "Down in Q4",
                               default="NS")]
dep_prs[, label := ifelse(fdr < 0.10, toupper(protein), NA_character_)]
pD <- ggplot(dep_prs, aes(logFC, -log10(P.Value), colour=direction, label=label)) +
    geom_point(size=0.8, alpha=0.6) +
    geom_text_repel(size=2.2, max.overlaps=15, segment.size=0.2, na.rm=TRUE) +
    scale_colour_manual(values=c("Up in Q4"="#CC0066","Down in Q4"="#56B4E9","NS"="grey65"),
                        guide="none") +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey50") +
    labs(x=expression(log[2]~FC~(PRS~Q4~vs~Q1)),
         y=expression(-log[10](p)),
         title=sprintf("D  PRS Q4 vs Q1 (HLA-adjusted)\n(FDR<0.05: %d)", sum(dep_prs$fdr < FDR_THR))) +
    theme_ukb(base_size=9)

# Panel E: HLA within-MS vs PRS concordance
conc_dt <- merge(dep_hla[,  .(protein, logFC_hla=logFC, fdr_hla=fdr)],
                 dep_prs[,  .(protein, logFC_prs=logFC, fdr_prs=fdr)],
                 by="protein")
conc_dt[, sig_either := (fdr_hla < FDR_THR) | (fdr_prs < FDR_THR)]
conc_dt[, label := ifelse(sig_either, toupper(protein), NA_character_)]
pE <- ggplot(conc_dt, aes(logFC_hla, logFC_prs, colour=sig_either, label=label)) +
    geom_point(size=0.8, alpha=0.6) +
    geom_text_repel(size=2.2, max.overlaps=12, segment.size=0.2, na.rm=TRUE) +
    geom_hline(yintercept=0, linewidth=0.3, colour="grey70") +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
    geom_smooth(method="lm", se=FALSE, colour="#E6A817", linewidth=0.5) +
    scale_colour_manual(values=c("TRUE"="#CC0066","FALSE"="grey65"),
                        labels=c("TRUE"="Sig. in either","FALSE"="NS"), name=NULL) +
    labs(x=expression(HLA+~vs~HLA-~log[2]~FC~(within~MS)),
         y=expression(PRS~Q4~vs~Q1~log[2]~FC),
         title="E  HLA vs PRS proteomic signature concordance") +
    theme_ukb(base_size=9)

fig5 <- (pA | pB) / pC / (pD | pE) +
    plot_layout(heights=c(1.2, 1.0, 1.2))

out_path <- file.path(OUT_DIR, "fig5_genetic_stratification.pdf")
ggsave(out_path, fig5, width=14, height=16, device=cairo_pdf)
cat(sprintf("  Saved: %s\n", basename(out_path)))

cat("\n03_ms_genetic_stratification.R complete.\n")
