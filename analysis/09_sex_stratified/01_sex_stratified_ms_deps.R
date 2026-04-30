#!/usr/bin/env Rscript
# 01_sex_stratified_ms_deps.R
# Sex-stratified differential proteomics: MS vs HC separately by sex
#
# Rationale: MS has a 3:1 female:male ratio; sex-dimorphic proteomic
# signatures may explain heterogeneous disease course (Abdelhak extension).
#
# Approach:
#   - Run limma on female-only and male-only PSM-matched cohorts
#   - Combined contrast only (powered by n)
#   - Compare DEP lists: shared, female-specific, male-specific
#   - Interaction test: identify proteins with sex-by-status interaction
#
# Output: results/sex_stratified/
#   ms_female_vs_hc.csv
#   ms_male_vs_hc.csv
#   ms_sex_dep_overlap.pdf
#   ms_sex_logfc_scatter.pdf
#   ms_sex_interaction.csv

suppressPackageStartupMessages({
    library(limma)
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(MatchIt)
    library(patchwork)
})

args      <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args, value = TRUE)
SCRIPT_DIR <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
PROJ_DIR   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."))
source(file.path(PROJ_DIR, "analysis", "helpers", "ukb_theme.R"))
source(file.path(PROJ_DIR, "analysis", "helpers", "limma_utils.R"))

QC_FILE     <- file.path(PROJ_DIR, "data", "ukb", "olink", "processed", "ms_olink_qc.csv")
CADASIL_DIR <- file.path(dirname(PROJ_DIR), "CADASIL_Proteome_ML_Keller_2024_Rebuttal", "data", "ukb")
OUT_DIR     <- file.path(PROJ_DIR, "results", "sex_stratified")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_THR    <- 0.05
PSM_RATIO  <- 10
PSM_CAL    <- 0.2
LOGFC_LABEL <- 0.25
set.seed(42)

load_covariates <- function(cadasil_dir) {
    misc <- file.path(cadasil_dir, "misc")
    bmi   <- fread(file.path(misc, "ukb_bmi_age.csv"), showProgress=FALSE)
    setnames(bmi, c("participant.eid","participant.p21001_i0"), c("eid","bmi"), skip_absent=TRUE)
    bmi <- bmi[, .(eid, bmi)]
    smoke <- fread(file.path(misc, "ukb_smoking_status.csv"), showProgress=FALSE)
    setnames(smoke, c("participant.eid","participant.p20116_i0"), c("eid","smoking_code"), skip_absent=TRUE)
    smoke[, ever_smoker := as.integer(smoking_code %in% c(1L,2L))]
    smoke <- smoke[smoking_code != -3L, .(eid, ever_smoker)]
    wmh <- fread(file.path(misc, "wmh_htn_dm.csv"), showProgress=FALSE)
    wmh[, diabetes := as.integer(p2443_i0=="Yes")]
    wmh <- wmh[p2443_i0 != "Prefer not to answer", .(eid, diabetes)]
    alc <- fread(file.path(misc, "ukb_alcohol_intake.csv"), showProgress=FALSE)
    setnames(alc, c("participant.eid","participant.p1558_i0"), c("eid","alcohol_freq"), skip_absent=TRUE)
    alc <- alc[alcohol_freq > 0L, .(eid, alcohol_freq)]
    Reduce(function(x,y) merge(x,y,by="eid",all.x=TRUE), list(bmi,smoke,wmh,alc))
}

# ── 1. Load & filter ──────────────────────────────────────────────────────────
cat("Loading QC'd MS Olink data...\n")
dt <- fread(QC_FILE, showProgress=FALSE)
covs <- load_covariates(CADASIL_DIR)
dt   <- merge(dt, covs, by="eid", all.x=TRUE)
dt[, sex_num := as.integer(factor(sex)) - 1L]

META_COLS <- c("eid","ms_status","age_at_sampling","age_at_diagnosis","years_to_diagnosis",
               "sex","olink_instance","qc_outlier","UMAP1","UMAP2","mean_npx",
               "bmi","ever_smoker","diabetes","alcohol_freq","sex_num",
               "PC1","PC2")
protein_cols <- setdiff(names(dt), META_COLS)

dt_filt <- dt[qc_outlier==FALSE & !is.na(age_at_sampling) & !is.na(sex)]
cat(sprintf("  After filter: %d participants\n", nrow(dt_filt)))

# Compute control-only PC1/PC2 once on the filtered cohort so every sex/stage
# sub-analysis uses the same PC space (matches 04_differential pipeline).
cat("Computing control-only PC1/PC2...\n")
pc_dt <- compute_control_pcs(dt_filt, protein_cols)
for (col in c("PC1", "PC2")) if (col %in% names(dt_filt)) set(dt_filt, j = col, value = NULL)
dt_filt <- merge(dt_filt, pc_dt, by = "eid", all.x = TRUE)
    # sex: 0=Female, 1=Male (UKB coding)
cat(sprintf("  Sex: Female(0)=%d, Male(1)=%d\n",
            sum(dt_filt$sex==0), sum(dt_filt$sex==1)))

# ── 2. PSM + limma per sex ────────────────────────────────────────────────────
run_sex_limma <- function(dt_sex, sex_label) {
    cat(sprintf("\n--- %s ---\n", sex_label))
    cat(sprintf("  N: %d cases, %d controls\n",
                sum(dt_sex$ms_status != "control"),
                sum(dt_sex$ms_status == "control")))

    PSM_VARS <- c("age_at_sampling","bmi","ever_smoker","diabetes","alcohol_freq")
    dt_sex[, is_case := as.integer(ms_status != "control")]
    dt_c <- dt_sex[complete.cases(as.data.frame(dt_sex[, c("is_case", ..PSM_VARS)]))]

    cat(sprintf("  PSM input: %d cases, %d controls\n",
                sum(dt_c$is_case), sum(!dt_c$is_case)))

    # No sex_num needed (already stratified by sex)
    psm_df <- as.data.frame(dt_c[, c("eid","is_case", ..PSM_VARS)])
    m <- matchit(is_case ~ age_at_sampling + bmi + ever_smoker + diabetes + alcohol_freq,
                 data=psm_df, method="nearest", ratio=PSM_RATIO,
                 caliper=PSM_CAL, std.caliper=TRUE)
    matched_eids <- as.integer(match.data(m)$eid)
    dt_m <- dt_c[eid %in% matched_eids]
    cat(sprintf("  After PSM: %d cases, %d controls\n",
                sum(dt_m$is_case), sum(!dt_m$is_case)))

    # Limma
    dt_m[, ms_status_f := factor(fifelse(ms_status=="control","control","ms"),
                                  levels=c("control","ms"))]
    meta <- as.data.frame(dt_m[, .(eid, ms_status_f, mean_npx, age_at_sampling, PC1, PC2)])
    prot_mat <- t(as.matrix(dt_m[, ..protein_cols]))
    colnames(prot_mat) <- dt_m$eid

    design <- model.matrix(~0 + ms_status_f + mean_npx + age_at_sampling + PC1 + PC2, meta)
    colnames(design) <- gsub("ms_status_f","",colnames(design))

    fit <- lmFit(prot_mat, design)
    cm  <- makeContrasts(ms - control, levels=design)
    fit2 <- eBayes(contrasts.fit(fit, cm))
    res  <- topTable_safe(fit2, coef=1)
    res[, sex := sex_label]
    n_dep <- sum(res$adj.P.Val < FDR_THR, na.rm=TRUE)
    cat(sprintf("  DEPs at FDR<0.05: %d\n", n_dep))
    res
}

dt_female <- dt_filt[sex == 0]   # 0 = Female in UKB
dt_male   <- dt_filt[sex == 1]   # 1 = Male in UKB

res_female <- run_sex_limma(dt_female, "Female")
res_male   <- run_sex_limma(dt_male,   "Male")

# Save FDR-significant results
fwrite(res_female[adj.P.Val < FDR_THR], file.path(OUT_DIR, "ms_female_vs_hc.csv"))
fwrite(res_male[adj.P.Val   < FDR_THR], file.path(OUT_DIR, "ms_male_vs_hc.csv"))
# Save full results (all proteins) for volcano plots
fwrite(res_female, file.path(OUT_DIR, "ms_female_vs_hc_all.csv"))
fwrite(res_male,   file.path(OUT_DIR, "ms_male_vs_hc_all.csv"))
cat(sprintf("\n  Saved: ms_female_vs_hc.csv (%d DEPs)\n", sum(res_female$adj.P.Val < FDR_THR)))
cat(sprintf("  Saved: ms_male_vs_hc.csv (%d DEPs)\n",   sum(res_male$adj.P.Val < FDR_THR)))

# ── 3. Overlap & scatter ───────────────────────────────────────────────────────
cat("\nComparing sex-specific DEP lists...\n")
f_sig <- res_female[adj.P.Val < FDR_THR, protein]
m_sig <- res_male[adj.P.Val   < FDR_THR, protein]
shared <- intersect(f_sig, m_sig)
f_only <- setdiff(f_sig, m_sig)
m_only <- setdiff(m_sig, f_sig)

cat(sprintf("  Female-only: %d | Male-only: %d | Shared: %d\n",
            length(f_only), length(m_only), length(shared)))

# logFC scatter for all proteins
scatter <- merge(
    res_female[, .(protein, female_logFC=logFC, female_sig=adj.P.Val<FDR_THR)],
    res_male[,   .(protein, male_logFC=logFC,   male_sig=adj.P.Val<FDR_THR)],
    by="protein"
)
scatter[, group := fcase(
    female_sig & male_sig  & (female_logFC>0)==(male_logFC>0), "Shared (concordant)",
    female_sig & male_sig  & (female_logFC>0)!=(male_logFC>0), "Shared (discordant)",
    female_sig & !male_sig, "Female only",
    !female_sig & male_sig, "Male only",
    default="NS"
)]
scatter[, label := fifelse(group!="NS", toupper(protein), NA_character_)]
scatter[, group_f := factor(group, levels=c("NS","Male only","Female only",
                                             "Shared (discordant)","Shared (concordant)"))]

SEX_COLS <- c("Shared (concordant)"="#CC0066","Shared (discordant)"="#E6A817",
               "Female only"="#56B4E9","Male only"="#2B4C7E","NS"="grey80")
SEX_SIZE  <- c("Shared (concordant)"=2.2,"Shared (discordant)"=2.2,
               "Female only"=1.6,"Male only"=1.6,"NS"=0.6)

p_scatter <- ggplot(scatter, aes(x=female_logFC, y=male_logFC,
                                   colour=group_f, size=group_f)) +
    geom_hline(yintercept=0, linewidth=0.3, colour="grey50", linetype="dashed") +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey50", linetype="dashed") +
    geom_point(alpha=0.6, shape=16) +
    geom_text_repel(aes(label=label), size=2.5, colour="black",
                    max.overlaps=30, box.padding=0.3, segment.size=0.25,
                    segment.colour="grey40", na.rm=TRUE) +
    scale_colour_manual(values=SEX_COLS, name=NULL,
                        guide=guide_legend(override.aes=list(size=2.2, alpha=1))) +
    scale_size_manual(values=SEX_SIZE, guide="none") +
    labs(x="Female logFC (MS vs HC)",
         y="Male logFC (MS vs HC)",
         title="Sex-stratified MS DEPs: logFC comparison") +
    theme_ukb()

ggsave(file.path(OUT_DIR, "ms_sex_logfc_scatter.pdf"),
       p_scatter, width=5.5, height=5, device=cairo_pdf)
cat("  Saved: ms_sex_logfc_scatter.pdf\n")

# DEP count bars
count_dt <- rbind(
    data.table(group="Female only", n=length(f_only)),
    data.table(group="Male only",   n=length(m_only)),
    data.table(group="Shared",      n=length(shared))
)
p_counts <- ggplot(count_dt, aes(x=group, y=n, fill=group)) +
    geom_col(width=0.6, colour="white") +
    geom_text(aes(label=n), vjust=-0.4, size=3.5) +
    scale_fill_manual(values=c("Female only"="#56B4E9","Male only"="#2B4C7E","Shared"="#CC0066"),
                      guide="none") +
    labs(x=NULL, y="Number of DEPs (FDR<0.05)",
         title="Sex-stratified MS DEP overlap") +
    theme_ukb()

ggsave(file.path(OUT_DIR, "ms_sex_dep_overlap.pdf"),
       p_counts, width=4, height=3.5, device=cairo_pdf)
cat("  Saved: ms_sex_dep_overlap.pdf\n")

# ── 4. Sex-by-status interaction (full cohort) ────────────────────────────────
cat("\nRunning sex-by-status interaction model...\n")
dt_int <- dt_filt[complete.cases(
    as.data.frame(dt_filt[, .(age_at_sampling, sex, PC1, PC2, mean_npx)])
)]
dt_int[, is_ms := as.integer(ms_status != "control")]
dt_int[, sex_f := factor(sex, levels=c(1,0), labels=c("Male","Female"))]

meta_int <- as.data.frame(dt_int[, .(is_ms, sex_f, mean_npx, age_at_sampling, PC1, PC2)])
prot_mat_int <- t(as.matrix(dt_int[, ..protein_cols]))
colnames(prot_mat_int) <- dt_int$eid

design_int <- model.matrix(~is_ms * sex_f + mean_npx + age_at_sampling + PC1 + PC2, meta_int)
fit_int  <- lmFit(prot_mat_int, design_int)
fit_int2 <- eBayes(fit_int)

# Interaction coefficient = is_ms:sex_fFemale
int_res <- topTable_safe(fit_int2, coef="is_ms:sex_fFemale")
n_int <- sum(int_res$adj.P.Val < FDR_THR, na.rm=TRUE)
cat(sprintf("  Sex-by-MS interaction DEPs at FDR<0.05: %d\n", n_int))
print(int_res[1:min(10,.N), .(protein, logFC=round(logFC,3), adj.P.Val=round(adj.P.Val,4))])

fwrite(int_res, file.path(OUT_DIR, "ms_sex_interaction.csv"))
cat("  Saved: ms_sex_interaction.csv\n")

# ── 5. Stage-stratified sex analysis (pre-onset & post-onset separately) ──────
# Answers: is sex dimorphism already present before diagnosis?
cat("\nRunning stage-stratified sex analyses...\n")

run_stage_sex <- function(dt_stage, stage_label) {
    cat(sprintf("\n=== %s ===\n", stage_label))
    results <- list()
    for (sx in c(0L, 1L)) {
        sx_lbl <- if (sx == 0L) "Female" else "Male"
        dt_sx <- dt_stage[sex == sx]
        n_case <- sum(dt_sx$ms_status != "control")
        if (n_case < 20) {
            cat(sprintf("  %s: only %d cases, skipping\n", sx_lbl, n_case))
            next
        }
        cat(sprintf("  %s: %d cases, %d controls\n", sx_lbl, n_case,
                    sum(dt_sx$ms_status == "control")))
        res <- run_sex_limma(dt_sx, sprintf("%s %s", stage_label, sx_lbl))
        results[[sx_lbl]] <- res
    }
    results
}

# Pre-onset: stage cases + all controls
dt_pre  <- dt_filt[ms_status %in% c("pre_onset",  "control")]
dt_post <- dt_filt[ms_status %in% c("post_onset", "control")]

res_pre  <- run_stage_sex(dt_pre,  "Pre-onset")
res_post <- run_stage_sex(dt_post, "Post-onset")

for (sx_lbl in names(res_pre)) {
    tag <- tolower(sx_lbl)
    fwrite(res_pre[[sx_lbl]],
           file.path(OUT_DIR, sprintf("ms_pre_%s_vs_hc_all.csv",  tag)))
    cat(sprintf("  Saved: ms_pre_%s_vs_hc_all.csv (%d DEPs)\n", tag,
                sum(res_pre[[sx_lbl]]$adj.P.Val < FDR_THR, na.rm=TRUE)))
}
for (sx_lbl in names(res_post)) {
    tag <- tolower(sx_lbl)
    fwrite(res_post[[sx_lbl]],
           file.path(OUT_DIR, sprintf("ms_post_%s_vs_hc_all.csv", tag)))
    cat(sprintf("  Saved: ms_post_%s_vs_hc_all.csv (%d DEPs)\n", tag,
                sum(res_post[[sx_lbl]]$adj.P.Val < FDR_THR, na.rm=TRUE)))
}

cat("\n01_sex_stratified_ms_deps.R complete.\n")
