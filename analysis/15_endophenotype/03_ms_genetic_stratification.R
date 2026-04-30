#!/usr/bin/env Rscript
# 03_genetic_stratification.R  — Figure 5
#
# DESIGN: HLA × disease interaction proteomics + stratified analyses
#
#   Core analysis: full interaction model in the matched cohort (cases + HC)
#     protein ~ status * hla_carrier + age_at_sampling + sex
#   Estimates:
#     case_effect  — case vs HC in HLA- background (reference)
#     hla_effect   — HLA+ vs HLA- in HC (main effect in controls)
#     interaction  — does the case proteomic effect differ by HLA carrier status?
#   Derived contrasts via limma:
#     HLA- case vs HLA- HC  (pure case effect without HLA modulation)
#     HLA+ case vs HLA+ HC  (case effect in HLA+ background)
#
#   Supplementary analyses:
#     B. HLA+ vs HLA- within cases (within-disease stratification)
#     C. Prodromal trajectory × HLA (pre-onset cases, spline × carrier interaction)
#     D. PRS Q4 vs Q1 within cases (HLA-adjusted)
#
#   Disease-, HLA-, PRS-, and trajectory-marker constants come from
#   configs/disease.yaml.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(limma)
    library(splines)
    library(here)
    library(glue)
})

source(here::here("analysis", "helpers", "disease_config.R"))
cfg <- load_disease_config()

STATUS_COL    <- cfg$cohort_status_col
PRE_ONSET     <- cfg$status_values$pre_onset
POST_ONSET    <- cfg$status_values$post_onset
CONTROL       <- cfg$status_values$control
COHORT        <- cfg$cohort_short
DISEASE_LBL   <- cfg$disease_short_caps
HLA_CARRIER   <- cfg$hla_carrier_col
PRS_COL       <- cfg$prs_combined_col
HLA_ALLELE    <- cfg$hla_allele
PRS_LABEL     <- cfg$prs_label
TRAJ_MARKERS  <- cfg$cns_injury_markers

source(here::here("analysis", "helpers", "ukb_theme.R"))

QC_FILE  <- here::here("data", "ukb", "olink", "processed",
                       glue::glue("{COHORT}_olink_qc.csv"))
HLA_FILE <- here::here("data", "ukb", "genetics", "hla_drb1_imputed.csv")
PRS_FILE <- here::here("data", "ukb", "genetics",
                       glue::glue("{COHORT}_prs_scores.csv"))
OUT_DIR  <- here::here("results", "endophenotype")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
FDR_THR  <- 0.05
N_MATCH  <- 3   # HC per case in matched cohort

# ── 1. Load cohort ────────────────────────────────────────────────────────────
# Genetic stratification needs both HLA and PRS data; skip cleanly if either
# is absent (matching the pattern in 07_hla_subgroup and 10_prs/01_ms_prs).
cat("Loading data...\n")
if (!file.exists(QC_FILE))  stop("QC file not found: ", QC_FILE)
if (!file.exists(HLA_FILE)) {
    cat(sprintf("Genetic-stratification stage skipped: %s not found.\n",
                basename(HLA_FILE)))
    cat("  See EXTRACT_TODO.md (HLA section) to enable.\n")
    quit(save = "no", status = 0)
}
if (!file.exists(PRS_FILE)) {
    cat(sprintf("Genetic-stratification stage skipped: %s not found.\n",
                basename(PRS_FILE)))
    cat("  Populate cfg$prs_pgs_ids and run 10_prs/00_merge_prs_scores.R first.\n")
    quit(save = "no", status = 0)
}
qc  <- fread(QC_FILE,  showProgress=FALSE)
hla <- fread(HLA_FILE, showProgress=FALSE)
prs <- fread(PRS_FILE, showProgress=FALSE)

meta_cols <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex", "olink_instance", "qc_outlier",
               "UMAP1", "UMAP2", "mean_npx", "sd_npx",
               "lower_iqr", "upper_iqr", "iqr_npx")

case_all <- qc[qc_outlier == FALSE &
               get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET) &
               !is.na(age_at_diagnosis)]
hc_all   <- qc[qc_outlier == FALSE & get(STATUS_COL) == CONTROL]

prot_cols <- setdiff(names(qc), meta_cols)

# Attach HLA to everyone
case_all <- merge(case_all,
                  hla[, .(eid, hla_carrier = get(HLA_CARRIER))],
                  by = "eid", all.x = TRUE)
hc_all   <- merge(hc_all,
                  hla[, .(eid, hla_carrier = get(HLA_CARRIER))],
                  by = "eid", all.x = TRUE)

case_all <- merge(case_all, prs[, .(eid, prs_score = get(PRS_COL))],
                  by = "eid", all.x = TRUE)

cat(sprintf("  %s cases: %d (HLA+=%.0f, HLA-=%.0f)\n",
    DISEASE_LBL, nrow(case_all),
    sum(case_all$hla_carrier == 1, na.rm = TRUE),
    sum(case_all$hla_carrier == 0, na.rm = TRUE)))
cat(sprintf("  HC available: %d (HLA data: %d)\n",
    nrow(hc_all), sum(!is.na(hc_all$hla_carrier))))

# ── 1b. Control-only PC1/PC2 (project all participants into control PC space) ─
# Control-derived PCs capture technical/batch variance without disease signal
# contamination. All participants (cases + HC) are projected into this space.
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
case_all <- merge(case_all, pc_dt_gs, by = "eid", all.x = TRUE)
hc_all   <- merge(hc_all,   pc_dt_gs, by = "eid", all.x = TRUE)
cat(sprintf("  PCs merged: %d cases, %d HC with valid PCs\n",
    sum(!is.na(case_all$PC1)), sum(!is.na(hc_all$PC1))))

# ── 2. Build matched cohort for interaction model ─────────────────────────────
# Age+sex match: for each case, sample N_MATCH controls from HCs with
# |age_at_sampling - hc_age| <= 2 and same sex. Use with-replacement if needed.
cat(sprintf("Building %d:1 age+sex matched HC cohort...\n", N_MATCH))
set.seed(42)
matched_hc_ids <- unlist(lapply(seq_len(nrow(case_all)), function(i) {
    row   <- case_all[i]
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

base_cols <- c("eid", STATUS_COL, "age_at_sampling", "age_at_diagnosis",
               "years_to_diagnosis", "sex")
full_cohort <- rbind(
    case_all[,   c(base_cols, "hla_carrier", "PC1", "PC2", prot_cols), with = FALSE],
    matched_hc[, c(base_cols, "hla_carrier", "PC1", "PC2", prot_cols), with = FALSE],
    fill = TRUE
)
full_cohort[, case_flag := as.integer(get(STATUS_COL) %in% c(PRE_ONSET, POST_ONSET))]

cat(sprintf("  Matched cohort: %d cases + %d HC = %d total\n",
    nrow(case_all), nrow(matched_hc), nrow(full_cohort)))
cat(sprintf("  HLA coverage: cases=%.0f%%, HC=%.0f%%\n",
    100 * mean(!is.na(case_all$hla_carrier)),
    100 * mean(!is.na(matched_hc$hla_carrier))))

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

# ── 3. Part A: HLA × case interaction model ──────────────────────────────────
cat(sprintf("\n=== Part A: HLA × %s interaction model ===\n", DISEASE_LBL))

int_dat <- full_cohort[!is.na(hla_carrier) & !is.na(age_at_sampling) & !is.na(sex)]
# Explicit interaction term (avoids : in coefficient names, which breaks makeContrasts)
# Coefficients:
#   Intercept       = HLA- HC mean
#   case_flag       = case effect in HLA- background
#   hla_carrier     = HLA+ vs HLA- in HC
#   case_hla_int    = INTERACTION: extra case effect conferred by HLA+ carrier
#   age_at_sampling, sex = covariates
int_dat[, case_hla_int := case_flag * hla_carrier]
# ytd_design: years_to_diagnosis (NA→0 for HC); adjusts for disease phase/duration
int_dat[, ytd_design   := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
design_int <- model.matrix(~ case_flag + hla_carrier + case_hla_int +
                               ytd_design + age_at_sampling + sex + PC1 + PC2,
                            data = int_dat)

cat("  Running interaction limma (n=", nrow(int_dat), ")...\n", sep="")
expr_int <- t(as.matrix(int_dat[, ..prot_cols]))
colnames(expr_int) <- int_dat$eid
row_means <- rowMeans(expr_int, na.rm=TRUE)
for (i in seq_len(nrow(expr_int)))
    expr_int[i, is.na(expr_int[i,])] <- row_means[i]
fit_int <- eBayes(lmFit(expr_int, design_int), trend=TRUE)

# Interaction term
dep_interaction <- as.data.table(topTable(fit_int, coef="case_hla_int", n=Inf, sort.by="P"))
dep_interaction[, protein := rownames(topTable(fit_int, coef="case_hla_int", n=Inf, sort.by="P"))]
dep_interaction[, fdr     := p.adjust(P.Value, method="BH")]
n_int <- sum(dep_interaction$fdr < FDR_THR)
cat(sprintf("  Interaction DEPs (FDR<0.05): %d\n", n_int))

# Case effect in HLA- background
dep_case_hlaneg <- as.data.table(topTable(fit_int, coef="case_flag", n=Inf, sort.by="P"))
dep_case_hlaneg[, protein := rownames(topTable(fit_int, coef="case_flag", n=Inf, sort.by="P"))]
dep_case_hlaneg[, fdr := p.adjust(P.Value, method="BH")]
cat(sprintf("  Case effect (HLA- background) DEPs: %d\n", sum(dep_case_hlaneg$fdr < FDR_THR)))

# Case effect in HLA+ background = case_flag + case_hla_int (contrast)
cont_mat <- makeContrasts("case_flag + case_hla_int", levels = design_int)
fit_cont <- contrasts.fit(fit_int, cont_mat)
fit_cont <- eBayes(fit_cont, trend=TRUE)
dep_case_hlapos <- as.data.table(topTable(fit_cont, coef=1, n=Inf, sort.by="P"))
dep_case_hlapos[, protein := rownames(topTable(fit_cont, coef=1, n=Inf, sort.by="P"))]
dep_case_hlapos[, fdr := p.adjust(P.Value, method="BH")]
cat(sprintf("  Case effect (HLA+ background) DEPs: %d\n", sum(dep_case_hlapos$fdr < FDR_THR)))

# Merge for Panel B scatter
scatter_dt <- merge(
    dep_case_hlaneg[, .(protein, logFC_hlaneg=logFC, p_hlaneg=P.Value, fdr_hlaneg=fdr)],
    dep_case_hlapos[, .(protein, logFC_hlapos=logFC, p_hlapos=P.Value, fdr_hlapos=fdr)],
    by="protein"
)
scatter_dt <- merge(scatter_dt,
    dep_interaction[, .(protein, logFC_int=logFC, fdr_int=fdr)], by="protein")
scatter_dt[, sig_int := fdr_int < FDR_THR]
cat(sprintf("  Proteins with discordant case effect (sig interaction): %d\n",
            sum(scatter_dt$sig_int)))

fwrite(dep_interaction, file.path(OUT_DIR, glue::glue("{COHORT}_hla_interaction_deps.csv")))
fwrite(scatter_dt,      file.path(OUT_DIR, glue::glue("{COHORT}_hla_stratified_logfc.csv")))

# ── 4. Part B: HLA+ vs HLA- within cases only ────────────────────────────────
cat(sprintf("\n=== Part B: HLA+ vs HLA- within %s ===\n", DISEASE_LBL))
case_hla <- case_all[!is.na(hla_carrier)]
case_hla[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
design_hla <- model.matrix(~ hla_carrier + ytd_design + age_at_sampling + sex + PC1 + PC2,
                             data = as.data.frame(case_hla[, .(
                                 hla_carrier     = factor(hla_carrier),
                                 ytd_design, age_at_sampling, sex, PC1, PC2)]))
dep_hla <- run_limma(case_hla, design_hla, "hla_carrier1")
cat(sprintf("  HLA+ vs HLA- DEPs (FDR<0.05): %d\n", sum(dep_hla$fdr < FDR_THR)))
fwrite(dep_hla, file.path(OUT_DIR, glue::glue("{COHORT}_hla_dep_results.csv")))

# ── 5. Part C: Prodromal trajectory × HLA ────────────────────────────────────
cat("\n=== Part C: Prodromal trajectory × HLA (pre-onset only) ===\n")
pre_onset <- case_all[get(STATUS_COL) == PRE_ONSET & !is.na(hla_carrier)]
cat(sprintf("  Pre-onset cases: %d (HLA+=%.0f, HLA-=%.0f)\n",
            nrow(pre_onset), sum(pre_onset$hla_carrier==1), sum(pre_onset$hla_carrier==0)))

target_prots <- intersect(TRAJ_MARKERS, prot_cols)
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

# ── 6. Part D: PRS Q4 vs Q1 within cases (HLA-adjusted) ──────────────────────
cat("\n=== Part D: PRS Q4 vs Q1 (HLA-adjusted) ===\n")
case_prs <- case_all[!is.na(prs_score) & !is.na(hla_carrier)]
case_prs[, prs_q := cut(prs_score, quantile(prs_score, c(0,.25,.75,1)), include.lowest=TRUE,
                        labels=c("Q1","mid","Q4"))]
prs_sub <- case_prs[prs_q %in% c("Q1","Q4")]
cat(sprintf("  PRS Q4: %d | Q1: %d\n", sum(prs_sub$prs_q=="Q4"), sum(prs_sub$prs_q=="Q1")))
design_prs <- model.matrix(~ prs_q + age_at_sampling + sex + hla_carrier + PC1 + PC2,
                             data = as.data.frame(prs_sub[, .(
                                 prs_q = factor(prs_q, levels = c("Q1","Q4")),
                                 age_at_sampling, sex, hla_carrier, PC1, PC2)]))
dep_prs <- run_limma(prs_sub, design_prs, "prs_qQ4")
cat(sprintf("  PRS Q4 vs Q1 DEPs (FDR<0.05): %d\n", sum(dep_prs$fdr < FDR_THR)))
fwrite(dep_prs, file.path(OUT_DIR, glue::glue("{COHORT}_prs_dep_results.csv")))

# ── 7. Part E: PRS × case proteome interaction (genome-wide) ─────────────────
# Does higher PRS amplify the case vs HC protein signature?
# Model: protein ~ case_flag + prs_z + case_flag:prs_z + ytd + age + sex + PC1 + PC2
# Interaction coef case_flag:prs_z answers: per 1SD PRS, how much does the
# case–HC gap change? Run in the full matched cohort with PRS merged.
cat(sprintf("\n=== Part E: PRS × %s proteome interaction (genome-wide) ===\n",
            DISEASE_LBL))

fc_prs <- merge(full_cohort, prs[, .(eid, prs_score = get(PRS_COL))],
                by = "eid", all.x = TRUE)
fc_prs <- fc_prs[!is.na(prs_score) & !is.na(age_at_sampling) & !is.na(sex) & !is.na(PC1)]
fc_prs[, prs_z      := as.numeric(scale(prs_score))]
fc_prs[, ytd_fc     := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
cat(sprintf("  PRS×case cohort: %d cases + %d HC\n",
    sum(fc_prs$case_flag == 1), sum(fc_prs$case_flag == 0)))

design_prs_int <- model.matrix(
    ~ case_flag + prs_z + case_flag:prs_z + ytd_fc + age_at_sampling + sex + PC1 + PC2,
    data = as.data.frame(fc_prs)
)
colnames(design_prs_int) <- make.names(colnames(design_prs_int))
prs_int_coef <- grep("case_flag.*prs_z|prs_z.*case_flag", colnames(design_prs_int), value = TRUE)
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
cat(sprintf("  PRS×case interaction FDR<0.05: %d | nominal p<0.05: %d\n",
    n_prs_int, sum(dep_prs_int$P.Value < 0.05)))
cat("  Top 10 by p:\n")
print(dep_prs_int[order(P.Value)][1:10, .(protein, logFC, P.Value, fdr)])
fwrite(dep_prs_int,
       file.path(OUT_DIR, glue::glue("{COHORT}_prs_interaction_genome.csv")))

# ── 8. Part F: Sex-stratified HLA and PRS×case analyses ──────────────────────
# Run HLA×case interaction, within-case HLA+ vs HLA−, and PRS×case interaction
# separately in females (sex=0) and males (sex=1) to test for sex-specific effects.
cat("\n=== Part F: Sex-stratified genetic analyses ===\n")

run_sex_stratum <- function(sex_val, sex_label) {
    cat(sprintf("\n--- %s (sex=%d) ---\n", sex_label, sex_val))
    out <- list()

    # F1: HLA×case interaction
    fc_s <- full_cohort[sex == sex_val & !is.na(hla_carrier) &
                        !is.na(age_at_sampling) & !is.na(PC1)]
    fc_s[, case_hla_int := case_flag * hla_carrier]
    fc_s[, ytd_design   := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
    cat(sprintf("  HLA×case cohort: %d cases + %d HC\n",
        sum(fc_s$case_flag == 1), sum(fc_s$case_flag == 0)))

    if (sum(fc_s$case_flag == 1) >= 20) {
        design_s <- model.matrix(
            ~ case_flag + hla_carrier + case_hla_int + ytd_design + age_at_sampling + PC1 + PC2,
            data = as.data.frame(fc_s))
        expr_s  <- t(as.matrix(fc_s[, ..prot_cols]))
        colnames(expr_s) <- fc_s$eid
        rm_s <- rowMeans(expr_s, na.rm = TRUE)
        for (i in seq_len(nrow(expr_s))) expr_s[i, is.na(expr_s[i,])] <- rm_s[i]
        fit_s   <- eBayes(lmFit(expr_s, design_s), trend = TRUE)
        int_s   <- as.data.table(topTable(fit_s, coef = "case_hla_int", n = Inf, sort.by = "P"))
        int_s[, protein   := rownames(topTable(fit_s, coef = "case_hla_int", n = Inf, sort.by = "P"))]
        int_s[, fdr       := p.adjust(P.Value, method = "BH")]
        int_s[, sex_group := sex_label]
        int_s[, analysis  := "hla_case_interaction"]
        cat(sprintf("  HLA×case FDR<0.05: %d | nominal p<0.05: %d\n",
            sum(int_s$fdr < FDR_THR), sum(int_s$P.Value < 0.05)))
        out$hla_int <- int_s
    }

    # F2: Within-case HLA+ vs HLA−
    case_s <- case_all[sex == sex_val & !is.na(hla_carrier) & !is.na(PC1)]
    if (nrow(case_s) >= 20) {
        case_s[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
        design_hs <- model.matrix(
            ~ hla_carrier + ytd_design + age_at_sampling + PC1 + PC2,
            data = as.data.frame(case_s[, .(hla_carrier = factor(hla_carrier),
                                            ytd_design, age_at_sampling, PC1, PC2)]))
        hla_s     <- run_limma(case_s, design_hs, "hla_carrier1")
        hla_s[, sex_group := sex_label]
        hla_s[, analysis  := "hla_within_case"]
        cat(sprintf("  Within-case HLA+ vs HLA− FDR<0.05: %d | nominal p<0.05: %d\n",
            sum(hla_s$fdr < FDR_THR), sum(hla_s$P.Value < 0.05)))
        out$hla_within <- hla_s
    }

    # F3: PRS×case interaction (sex-stratified)
    fc_ps <- merge(full_cohort[sex == sex_val],
                   prs[, .(eid, prs_score = get(PRS_COL))], by = "eid", all.x = TRUE)
    fc_ps <- fc_ps[!is.na(prs_score) & !is.na(age_at_sampling) & !is.na(PC1)]
    cat(sprintf("  PRS×case cohort: %d cases + %d HC\n",
        sum(fc_ps$case_flag == 1), sum(fc_ps$case_flag == 0)))

    if (sum(fc_ps$case_flag == 1) >= 20) {
        fc_ps[, prs_z      := as.numeric(scale(prs_score))]
        fc_ps[, ytd_design := fifelse(is.na(years_to_diagnosis), 0, years_to_diagnosis)]
        design_ps <- model.matrix(
            ~ case_flag + prs_z + case_flag:prs_z + ytd_design + age_at_sampling + PC1 + PC2,
            data = as.data.frame(fc_ps))
        colnames(design_ps) <- make.names(colnames(design_ps))
        prs_ic_s <- grep("case_flag.*prs_z|prs_z.*case_flag", colnames(design_ps), value = TRUE)
        expr_ps  <- t(as.matrix(fc_ps[, ..prot_cols]))
        colnames(expr_ps) <- fc_ps$eid
        rm_ps <- rowMeans(expr_ps, na.rm = TRUE)
        for (i in seq_len(nrow(expr_ps))) expr_ps[i, is.na(expr_ps[i,])] <- rm_ps[i]
        fit_ps  <- eBayes(lmFit(expr_ps, design_ps), trend = TRUE)
        prs_s   <- as.data.table(topTable(fit_ps, coef = prs_ic_s, n = Inf, sort.by = "P"))
        prs_s[, protein   := rownames(topTable(fit_ps, coef = prs_ic_s, n = Inf, sort.by = "P"))]
        prs_s[, fdr       := p.adjust(P.Value, method = "BH")]
        prs_s[, sex_group := sex_label]
        prs_s[, analysis  := "prs_case_interaction"]
        cat(sprintf("  PRS×case FDR<0.05: %d | nominal p<0.05: %d\n",
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
    fwrite(sex_all,
           file.path(OUT_DIR, glue::glue("{COHORT}_sex_stratified_genetic.csv")))
    cat(sprintf("\nSaved: %s_sex_stratified_genetic.csv\n", COHORT))
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
         title=sprintf("A  HLA × %s interaction\n(FDR<0.05: %d proteins)",
                       DISEASE_LBL, n_int)) +
    theme_ukb(base_size=9)

# Panel B: HLA- case-effect vs HLA+ case-effect scatter
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
    labs(x=bquote(.(DISEASE_LBL)~log[2]~FC~(HLA^"-"~background)),
         y=bquote(.(DISEASE_LBL)~log[2]~FC~(HLA^"+"~background)),
         title=sprintf("B  %s proteomic effect: HLA- vs HLA+ background", DISEASE_LBL)) +
    theme_ukb(base_size=9)

# Panel C: Trajectory plots
if (!is.null(traj_results) && nrow(traj_results) > 0) {
    prot_labels <- setNames(toupper(target_prots), target_prots)
    traj_results[, prot_label := prot_labels[protein]]
    traj_results[, hla_label  := ifelse(hla_carrier==1, "HLA+", "HLA-")]
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

# Panel E: HLA within-case vs PRS concordance
conc_dt <- merge(dep_hla[, .(protein, logFC_hla=logFC, fdr_hla=fdr)],
                 dep_prs[, .(protein, logFC_prs=logFC, fdr_prs=fdr)],
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
    labs(x=bquote(HLA+~vs~HLA-~log[2]~FC~(within~.(DISEASE_LBL))),
         y=expression(PRS~Q4~vs~Q1~log[2]~FC),
         title="E  HLA vs PRS proteomic signature concordance") +
    theme_ukb(base_size=9)

fig5 <- (pA | pB) / pC / (pD | pE) +
    plot_layout(heights=c(1.2, 1.0, 1.2))

out_path <- file.path(OUT_DIR, "fig5_genetic_stratification.pdf")
ggsave(out_path, fig5, width=14, height=16, device=cairo_pdf)
cat(sprintf("  Saved: %s\n", basename(out_path)))

cat("\n03_genetic_stratification.R complete.\n")
