# Disease replication agent — handoff document

## TL;DR — Quick-start for the next disease

```bash
# 1. Clone the template into a new disease folder
cd /Users/jonahkeller/ELab
git clone git@github.com:Jonah-Keller/UKB_MS.git UKB_<DISEASE>
cd UKB_<DISEASE>
ln -s /Users/jonahkeller/ELab/UKB_MS/data data            # symlink shared UKB extract
gh repo create Jonah-Keller/UKB_<DISEASE> --private --source=. --remote=origin --push
git remote set-url origin git@github.com:Jonah-Keller/UKB_<DISEASE>.git
git remote add upstream git@github.com:Jonah-Keller/UKB_MS.git

# 2. Send the user the 6 OpenEvidence/PGS-Catalog/UKB-Showcase queries (§1).
# 3. Edit configs/disease.yaml with their answers (§3).
# 4. Run preflight — it tells you exactly what's blocked, what runs, what skips:
Rscript scripts/preflight/preflight.R

# 5. Build the cohort:
python3 analysis/01_cohort/01_build_ms_cohort.py
# (or 02_build_<disease>_cohort.py if the disease needs a custom case definition)

# 6. Run the pipeline in the order in §5.  Watch for stages that skip
#    cleanly with a "see EXTRACT_TODO.md §N" pointer — those are not bugs.

# 7. Build the manuscript:
Rscript scripts/build/build_submission.R

# 8. Final report (§10).
```

If preflight fails, fix the underlying issue first.  Don't spawn worker
agents until preflight is clean — the failure modes you'd hit during
parallel execution are exactly the ones preflight catches up-front.

---


You are an autonomous agent assigned to **replicate the UKB plasma-proteomics
study on a new disease cohort** using the template + the patterns proven on
two cohorts so far:

- `github.com/Jonah-Keller/UKB_MS` — the original template (multiple sclerosis)
- `github.com/Jonah-Keller/UKB_ALS` — first replication (amyotrophic lateral sclerosis)

The user will tell you the new disease (e.g., stroke, type 2 diabetes,
Alzheimer's, Parkinson's, vascular dementia). Everything else is in this
document and in `configs/disease.yaml`.

This supersedes the earlier `AGENT_TASK_disease_replication.md` —
**use this file**.

---

## 0. What success looks like

A new repository `UKB_<DISEASE>` (e.g. `UKB_Stroke`) where:
- The pipeline ran end-to-end with **zero analysis-script edits**
- All disease-/HLA-/PRS-/Mendelian-/phenomic-specific values live in
  `configs/disease.yaml`
- A `EXTRACT_TODO.md` documents every UKB-RAP pull the user still has to run
- Figures are regenerated from the new cohort's data, with all protein labels
  data-driven (top-N FDR-aware via `analysis/helpers/top_n_selector.R`)
- A final report names: sample sizes, biological findings (top DEPs, top
  hazards, prodromal PheWAS, cluster solutions), and skipped stages with
  the data extracts they need

The user provides the disease + answers to a small set of OpenEvidence /
PGS Catalog / UKB Showcase queries. You do everything else.

---

## 1. Phase 1 — Pre-flight: ask the user 6 scoped questions

Use AskUserQuestion or write the queries inline; whatever the harness allows.
Do **not** guess these — every one of them maps directly into
`configs/disease.yaml`.

The user has paid access to **OpenEvidence** (medical Q&A LLM), **PGS
Catalog** (free), and the **UKB Showcase** (free). They will run the queries
on your behalf and paste the answers.

### Question template (copy into your message verbatim)

```
For the <DISEASE> replication, please run these queries:

1. Open Evidence — disease cohort:
   "What is the UK Biobank first-occurrence field ID for <DISEASE>?
    Confirm the primary ICD-10 code(s) and any subtype-specific code if a
    narrower definition is preferred (e.g. ischemic-only stroke)."

2. Open Evidence — control exclusions:
   "List ICD-10 codes that should disqualify a participant from the control
    pool for <DISEASE>. Include diseases that share clinical phenotype,
    treatment confounding, or biology with <DISEASE>."

3. Open Evidence — HLA association:
   "In the largest GWAS of <DISEASE> to date, did any HLA-region locus reach
    genome-wide significance? Name the lead SNP and tagged HLA allele if so.
    If no replicable HLA finding exists, say so explicitly. Otherwise list
    candidate alleles by replicated effect size (OR > 1.2, p < 5e-8 in n>1000
    cases)."

4. Open Evidence — Mendelian / large-effect variants:
   "Are there Mendelian or large-effect variants enriched in <DISEASE>
    cases that could be extracted from UKB whole-exome/whole-genome
    sequencing data? Report carrier rate, expected count among UKB
    <DISEASE> cases, and whether the variant is callable from the
    standard UKB-RAP pipeline (no additional bioinformatics required)."

5. PGS Catalog (https://www.pgscatalog.org/trait_search/?q=<DISEASE>):
   "List 1-2 best polygenic scores for <DISEASE>. For each: PGS ID, source
    GWAS, variant count, training-data ancestry, and validation AUC or
    pseudo-R² in independent EUR cohorts."

6. Open Evidence — phenomic / lifestyle risk factors:
   "Rank the top 3 environmental / lifestyle risk factors for <DISEASE> by
    replicated effect size in cohort studies (≥3 independent replications,
    n > 5000 cases). State which UKB fields capture them. Flag any that
    cannot be reliably captured (free-text only, US-only ICD subdivisions
    that international ICD-10 lacks, etc.)."
```

If **any answer is "no robust association"** for HLA, Mendelian, PRS, or a
specific phenomic factor, capture that explicitly — leave the corresponding
cfg field as the empty/sentinel value and the pipeline will skip that stage
with a clear message. **Do not invent values.**

---

## 2. Phase 2 — Setup

```bash
cd /Users/jonahkeller/ELab
git clone git@github.com:Jonah-Keller/UKB_MS.git UKB_<DISEASE>
cd UKB_<DISEASE>

# Symlink the shared UKB extract (data/ is gitignored)
ln -s /Users/jonahkeller/ELab/UKB_MS/data data

# Create the new GitHub repo and re-point origin
gh repo create Jonah-Keller/UKB_<DISEASE> --private --source=. --remote=origin --push
git remote set-url origin git@github.com:Jonah-Keller/UKB_<DISEASE>.git
git push -u origin main
```

Add the UKB_MS template as `upstream` so future template fixes flow in:

```bash
git remote add upstream git@github.com:Jonah-Keller/UKB_MS.git
git fetch upstream main
```

---

## 3. Phase 3 — Configure `configs/disease.yaml`

Edit only this file. The schema (with all fields used by the proven
template) is documented below. Fields marked **REQUIRED** must be set;
fields marked _optional_ can be left as sentinel/empty and the pipeline
skips the corresponding stage gracefully.

```yaml
# ── Identifying labels (REQUIRED) ──────────────────────────────────────────
cohort_short:        "stroke"                         # filename prefix
disease_long:        "ischemic stroke"                # caption substitution
disease_short_caps:  "Stroke"                         # display label
project_short:       "ukb_stroke_proteomics"
project_title:       "..."
comparison_cohort_short: "ms"                         # which cohort to cross-compare to

# ── Cohort definition (REQUIRED) ───────────────────────────────────────────
icd_codes:                ["I63"]                     # primary case codes
first_occurrence_field:   "p131360"                   # UKB date-of-first-occurrence field
cohort_status_col:        "stroke_status"             # column written by 01_cohort
status_values:
  pre_onset:  "pre_onset"
  post_onset: "post_onset"
  control:    "control"

control_exclusion_codes:                              # NESTED dict, group keys are docs
  cerebrovascular:  ["I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69"]
  related_overlap:  [...]                             # disease-specific from question 2

# ── HLA risk locus (optional — set hla_allele: "NONE" to skip stage 07) ────
hla_allele:           "DRB1*15:01"                    # or "NONE", or multi-subtype encoding
hla_search_patterns:  ["15:01", "1501", "DRB1*15"]
hla_carrier_col:      "drb1_1501_carrier"
hla_dosage_col:       "drb1_1501_dosage"
# ── Multi-subtype OR with explicit exclusion (optional, e.g. ALS DRB1*04) ──
hla_excluded_subtypes: []                             # populate to enable case 4 of carrier derivation
hla_dosage_threshold:  0.8                            # set dosage <X to 0 per UKB docs
hla_direction:         "risk"                         # or "protective" — affects plot legend wording
hla_field:             22182

# ── Mendelian carrier (optional — set null to skip stage 07b) ──────────────
# Single Mendelian variant or expansion stratification.  ALS uses C9orf72.
mendelian_carrier_col:        "c9orf72_expansion_carrier"
mendelian_threshold_repeats:  30                      # primary cutoff
mendelian_threshold_strict:   100                     # for sensitivity analysis
mendelian_source:             "ukb_wgs_expansionhunter_REPCN_long"
other_mendelian_carrier_col:  null                    # pool of rare variants if disease has no single dominant one

# ── Phenomic / exposure stratifiers (optional) ─────────────────────────────
# Each entry is either:
#   - a nested dict with field, levels (optional), and use_as: "stratifier" | "covariate" | "covariate_or_dichotomize_at_<N>"
#   - a documentation string (e.g. "NOT_RELIABLY_CAPTURED_IN_UKB") — skipped at runtime
phenomic_stratifiers:
  smoking_status:
    field:  20116
    levels: ["never", "former", "current"]
    use_as: "stratifier"
  bmi:
    field:  21001
    use_as: "covariate_or_dichotomize_at_25"
  head_trauma_history:
    source: "hes_icd10"
    codes_prefix: ["S00", "S01", "S02", "S06", "S09"]
    use_as: "stratifier"
  bulbar_onset: "NOT_RELIABLY_CAPTURED_IN_UKB"

# ── Sex handling (REQUIRED — controls whether stage 09 runs) ───────────────
# primary_use:
#   "stratifier" → stage 09_sex_stratified runs the standard split-by-sex limma
#   "covariate"  → stage 09 exits cleanly; sex is a fixed effect in all models;
#                  sex × <factor> interactions are tested by stages 07b and 17
# do_not_stratify_on_sex_alone: TRUE forces 09 to skip even when primary_use is "stratifier"
# interaction_terms_to_test: list of "sex_x_<column_name>"; the matching stage
#   runs the 3-way sex × case × <column> sub-model
sex_handling:
  primary_use: "covariate"                            # OR "stratifier"
  field:       31                                     # self-reported; or 22001 for genetic
  do_not_stratify_on_sex_alone: true
  interaction_terms_to_test:
    - "sex_x_smoking_status"
    - "sex_x_c9orf72_carrier"

# ── Polygenic risk score (optional — empty list = skip stage 10) ───────────
prs_pgs_ids:          []                              # populate from PGS Catalog query
prs_label:            "Stroke PRS"
prs_combined_col:     "prs_score"

# ── Programmatic top-N protein selection (template-wide, do not change) ────
top_n_volcano:    20
top_n_heatmap:    30
top_n_trajectory: 3
min_sig_for_fdr:  5

# ── CNS biology reference (optional, biology-driven) ───────────────────────
cns_injury_markers:   ["nefl", "mog", "gfap"]         # canonical CNS-injury panel
trajectory_proteins:  ["NEFL", "GFAP"]                # per-protein trajectory panels in fig 2
```

After editing, **validate it loads** before running anything else:

```bash
Rscript -e 'source(here::here("analysis","helpers","disease_config.R")); cfg <- load_disease_config(); str(cfg, max.level=1)'
python3 -c 'from analysis.helpers.disease_config import load_disease_config; print(load_disease_config())'
```

Both must succeed. Commit + push the YAML:

```bash
git add configs/disease.yaml
git commit -m "config: switch disease.yaml to <DISEASE> replication"
git push
```

---

## 4. Phase 4 — Run the cohort builder

If the disease has a special case definition (e.g. ALS uses Chia et al.
2025 G12.2-confirmed-at-death; Alzheimer's might require G30 + cognitive
test threshold; vascular dementia is debated), check whether a parallel
cohort builder already exists in `analysis/01_cohort/`. The repo already
ships:

- `01_build_ms_cohort.py` — simple ICD-first-occurrence (covers most diseases)
- `02_build_als_cohort.py` — ALS-specific (G12.2 + death confirmation)

If your disease fits the simple case definition, just run the templated:

```bash
python3 analysis/01_cohort/01_build_ms_cohort.py
```

Otherwise create a parallel `0X_build_<disease>_cohort.py` modeled on
`02_build_als_cohort.py`. Both write `data/ukb/cohort/<cohort_short>_cohort.csv`.

**Sanity-check sample sizes immediately after:**

```r
qc <- data.table::fread(sprintf("data/ukb/olink/processed/%s_olink_qc.csv", cfg$cohort_short))
print(table(qc[[cfg$cohort_status_col]], useNA = "ifany"))
```

**STOP if pre-onset N < 30.** Pre-onset is the highest-value contrast; tiny
samples make the replication uninformative. Tell the user, propose a fallback
(e.g. broader ICD definition, longer follow-up window), and wait.

---

## 5. Phase 5 — Run the pipeline

Order of execution. Each stage skips gracefully if its inputs are missing.
**Watch the output of every stage** for "skipped" messages — those are
real, but document the cause in your final report.

```bash
# Olink QC (probably already exists; check)
[ -f data/ukb/olink/processed/$(yq -r .cohort_short configs/disease.yaml)_olink_qc.csv ] || \
  Rscript analysis/02_olink_qc/01_olink_qc.R

# Differential expression
Rscript analysis/04_differential/01_limma_ms_vs_hc.R
Rscript analysis/04_differential/03_cns_trajectories.R

# Networks + cell-types
Rscript analysis/05_networks/01_ppi_network.R
Rscript analysis/05_networks/04_prepost_go.R
Rscript analysis/05_networks/05_hpa_celltype_heatmap.R
Rscript analysis/05_networks/06_celltype_overrepresentation.R

# ML classifiers (slow — 10-30 min). Run in background:
nohup Rscript analysis/06_ml/01_ms_ml_classifier.R > /tmp/ml.log 2>&1 &

# HLA stratification (skips if no allele data extracted)
Rscript analysis/07_hla_subgroup/01_hla_interaction_limma.R
Rscript analysis/07_hla_subgroup/01_hla_ms_subgroup.R

# Mendelian carrier stratification (NEW — skips if c9orf72_expansion.csv missing)
Rscript analysis/07b_mendelian_subgroup/01_c9orf72_subgroup.R

# MRI / cognition
Rscript analysis/08_mri/01_wmh_protein_correlation.R
Rscript analysis/08_mri/02_brain_mri_dti_ms.R    # if applicable
Rscript analysis/08_mri/03_cognitive_ms.R         # if applicable

# Sex-stratified (skips automatically if cfg$sex_handling$primary_use == "covariate")
Rscript analysis/09_sex_stratified/01_sex_stratified_ms_deps.R

# PRS
Rscript analysis/10_prs/00_merge_prs_scores.R    # only if prs_pgs_ids non-empty
Rscript analysis/10_prs/01_ms_prs.R

# Cross-disease comparison
Rscript analysis/11_cross_disease/01_ms_als_comparison.R

# Survival + comorbidity + endophenotype
Rscript analysis/13_survival/01_ms_protein_cox.R
Rscript analysis/14_comorbidity/02a_ms_prems_phewas.R
Rscript analysis/15_endophenotype/02_ms_prems_comorbidity_endophenotype.R
for s in analysis/15_endophenotype/0[3-9]_*.R analysis/15_endophenotype/1[0-3]_*.R; do
  Rscript "$s" 2>&1 | tail -5
done
Rscript analysis/16_cluster_gwas/01_build_phenotype.R     # builds pheno; REGENIE runs on RAP

# Phenomic stratifiers (NEW)
Rscript analysis/17_phenomic_subgroup/01_phenomic_stratifiers.R
```

---

## 6. Phase 6 — Manuscript figures + build

```bash
Rscript analysis/manuscript/figure1.R    # pre-onset signature
Rscript analysis/manuscript/figure2.R    # temporal window
Rscript analysis/manuscript/figure3.R    # HLA + PRS (skips panels gracefully if data missing)
Rscript analysis/manuscript/figure4.R    # sex differences (panels reflect cfg$sex_handling)
Rscript analysis/manuscript/figure5.R    # endophenotype clusters

Rscript scripts/build/build_submission.R
```

**Outputs land in `manuscript/submission/`:**
- `main_text/manuscript.docx`
- `main_materials.docx` (figures + tables in one)
- `supplementary/Supplementary_Material.docx`
- `supplementary/Supplementary_Tables.xlsx`

---

## 7. Phase 7 — Document missing data

Write `EXTRACT_TODO.md` (template in `UKB_ALS/EXTRACT_TODO.md`). For each
stage that skipped due to missing data, include:
- What field/extract is needed
- Exact RAP query (Cohort Browser path or dx-toolkit command)
- Where to write the output CSV
- Expected sample counts for sanity-checking the pull

---

## 8. Known bug patterns (already fixed; watch for analogues)

These were surfaced during MS→ALS replication. If you find a similar pattern
in a new context, fix it generically, not as a special case. Commit message:
`fix: <pattern> surfaced during <DISEASE> replication`.

| # | Pattern | Symptom | Fix |
|---|---|---|---|
| 1 | Hardcoded `ms_status` in helper functions | "Object 'ms_status' not found" data.table NSE error | Add `status_col` argument with default; thread through call sites |
| 2 | `run_go()` arity collision (3 positional args expecting `min_genes`) | `invalid format '%d'; use format %s for character objects` in sprintf | Use named arg `universe = bg_vector`; add `universe` parameter to helper |
| 3 | HLA stage column-not-found | data.table .SDcols error with the configured `hla_carrier_col` | Both 07 scripts now bail gracefully — copy that pattern when adding new optional-data stages |
| 4 | PRS empty list crashes script | "File 'als_prs_scores.csv' does not exist or is non-readable" | Guard at top: if `length(cfg$prs_pgs_ids) == 0` then skip with a message |
| 5 | Empty cor_dt with small cohorts | "object 'p_value' not found" downstream | Wrap downstream operations in `if (nrow(...) == 0) skip` |
| 6 | `fcase` glue-class mismatch | "Argument #6 has different class than argument #2" | Cast `glue("...")` outputs with `as.character()` before `fcase`. **Recurs** in 11, 13, and figure3 — always wrap glue results that flow into a fcase output column |
| 7 | Hardcoded sex-stratification | Stage 09 produced underpowered separate-sex limmas for ALS regardless of cfg | Read `cfg$sex_handling$primary_use`, exit cleanly when "covariate" |
| 8 | Hardcoded interaction term names | `"sex_x_c9orf72_carrier" %in% cfg_list` doesn't match `sex_x_<other_col>` for other diseases | Use regex `^sex_x_(.+)$` to extract the column name and match generically |
| 9 | **Hardcoded k=3 cluster count** | Stages 04, 08, 11, 13 silently dropped clusters beyond C0/C1/C2 when ALS produced k=5 | Discover cluster set via `sort(unique(dt$cluster))` at runtime; build factor levels, contrast specs, and panel-letter assignments from the discovered set |
| 10 | **Cluster contrast helper hardcoded for 3 clusters** | `run_cluster_psm_limma` had explicit `c("C0","C1","C2")` loop and 6-contrast `makeContrasts` call | Iterate over `setdiff(levels(cluster_f), "None")` and build contrasts via `do.call(makeContrasts, c(specs, list(levels=design)))` |
| 11 | **CLUST_COLS palette had 3 hardcoded entries** | Cluster colors only defined for C0/C1/C2 | New `cluster_palette(n)` helper anchors first 3 on canonical UKB palette, extends with Okabe-Ito and `scales::hue_pal()` for higher k |
| 12 | **GMM G_TARGET hardcoded to 3** | Stage 12 GMM components hardcoded to align with k=3 k-means | Read cluster count from cluster assignments CSV at runtime; expand `G_RANGE` to span the discovered k |
| 13 | **caret `resamples()` requires ≥2 models** | "at least two train objects are needed" when only glmnet trained | Guard with `if (length(model_list) >= 2L) resamples() else build comparison from single model's `$resample` data.frame manually` |
| 14 | **MS-specific biological labels hardcoded in panel titles** | Volcano panel titles like `"C0 vs None  ·  spine / connective tissue"` survive into ALS output | Drop curated biology subtitles; use generic `"Ck vs None"` titles. Biology-specific labels belong in manuscript captions, not figure scripts |

### Latent bugs that MS+ALS happened to dodge

These were surfaced by audit (not by an actual run yet), with fixes
landed in commits 5193... onwards.  They will likely surface on stroke
(very short pre-onset window), DM2 (decades-long prodrome, 15K+ cases),
or AD (different biomarker repertoire) but were silent on MS+ALS because
both happen to have multi-year prodrome + few-hundred case cohort + NEFL
as a top biomarker.  Watch for analogues if you find a regression.

| # | Pattern | Symptom | Fix |
|---|---|---|---|
| 15 | **Hardcoded YTD bins assume multi-year prodrome** | Rolling Cox bins fixed at 0–4y, 4–8y, 8–12y, 12+y. Acute stroke has zero pre-Dx; DM2 has 20+ years | `cfg$ytd_bins` (named list of low/high pairs) overrides the default; data-driven bin choice via `quantile()` if not specified |
| 16 | **Featured-protein subplot hardcoded to NEFL by name** | `13_survival` had `nefl_roll <- roll_dt[protein == "nefl"]`; if disease has no NEFL signal, panel silently empty | Take `top_prots[1]` instead — whichever protein ranked top in the overall Cox; output filename uses `tolower(feat_disp)` |
| 17 | **PSM ratio 1:10 ignores actual control supply** | When a cluster has 8 cases and 50 controls, demanding 1:10 quietly fails | Auto-scale: `eff_ratio = min(target, floor(n_ctrl / max(1, n_treat)))`, log when scaled down |
| 18 | **Empty rolling-Cox output crashes downstream GO** | `02_temporal_go_enrichment.R` ran on `ms_protein_cox_rolling.csv` regardless of cohort + crashed on empty input | Cohort-aware path `{COHORT}_protein_cox_rolling.csv`; explicit empty-file + missing-file guards both bail cleanly with pointers |
| 19 | **`figure2.R BIN_BREAKS` hardcoded -10..+10 yr** | Bins assume disease window matches MS/ALS  | `cfg$temporal_bins` (numeric vector) overrides; auto-generate `BIN_LABELS` from breaks |
| 20 | **MS-only / ALS-only Python differential scripts mixed in with templated path** | `06_ms_differential/01_ms_differential.py`, `07_als_differential/01_als_differential.py` hardcode protein highlight lists from Abdelhak / Chia papers | Marked DEPRECATED in module docstring; pointer to canonical `04_differential/01_limma_ms_vs_hc.R`; new replications must NOT run these |

---

## 9. Validation checklist before final report

```bash
# Step 1 — preflight (run BEFORE doing anything substantive on a new disease)
Rscript scripts/preflight/preflight.R

# Preflight prints a per-stage table showing run/skip/blocked status.
# Exit code 0 means safe to proceed; 1 means a required stage is blocked.

# Step 2 — run the pipeline
# (see §5 for the canonical order)

# Step 3 — final-report validation
```

```r
# Loader sanity
Rscript -e 'source(here::here("analysis","helpers","disease_config.R")); cfg <- load_disease_config(); str(cfg, max.level=1)'

# Cohort sample sizes
Rscript -e 'cfg <- yaml::read_yaml(here::here("configs","disease.yaml"))
qc <- data.table::fread(sprintf("data/ukb/olink/processed/%s_olink_qc.csv", cfg$cohort_short))
print(table(qc[[cfg$cohort_status_col]], qc$qc_outlier, useNA="ifany"))'

# Top-N selector behavior
Rscript -e 'source(here::here("analysis","helpers","top_n_selector.R"))
cfg <- yaml::read_yaml(here::here("configs","disease.yaml"))
diff <- data.table::fread(sprintf("results/differential/%s_pre_vs_hc.csv", cfg$cohort_short))
sel <- top_n_by_direction(diff, n_each = 10)
print(sel)
print(attr(sel, "ordered_by"))'

# All stages exited cleanly (no Execution halted in the logs)
grep -L "Execution halted" logs/*.log

# Output completeness
ls manuscript/submission/main_text/manuscript.docx
ls manuscript/submission/main_materials.docx
ls manuscript/submission/supplementary/Supplementary_Material.docx
ls manuscript/submission/supplementary/Supplementary_Tables.xlsx
```

---

## 10. Final report format

When you stop, push everything and post a single-message summary:

1. **Sample sizes**: total / pre-onset / post-onset / controls.
   Pre-onset stratification by years-to-Dx bin. Sex breakdown.
2. **Stage-by-stage status**: ✅ complete / ❌ skipped / 🟡 partial.
   For skipped: the EXTRACT_TODO section blocking it.
3. **Top biological findings**: top 5 pre-onset DEPs, top 5 prospective
   hazards, top 5 prodromal ICD enrichments, cross-disease shared/discordant
   proteins, cluster solution k + silhouette.
4. **Stratification results** (if data was available): HLA × disease, PRS ×
   disease, Mendelian carrier × disease, smoking × disease — top 5 interaction
   DEPs each, with effect sizes.
5. **Bug fixes** you pushed: any new patterns from §8.
6. **Open data extracts** still pending: link to `EXTRACT_TODO.md`.

---

## 11. Suggested replication queue (per user request)

1. **Stroke** (`I63` ischemic, `I61`+`I60` hemorrhagic, or combined `I60-I69`)
   — large UKB N (probably >2000 cases), strong cardiovascular prodrome
   expected, no dominant HLA, well-validated PRS available
2. **Type 2 diabetes** (`E11`) — very large UKB N (>15,000 cases), metabolic
   prodrome, no HLA, multiple validated PRS scores from DIAGRAM consortium
3. **Alzheimer's disease** (`G30`) — moderate UKB N, well-known APOE
   association (treat as Mendelian-style stratifier), pre-clinical biomarker
   literature exists, IGAP/EADB PRS available
4. **Parkinson's disease** (`G20`) — moderate UKB N, prodromal symptoms
   (hyposmia, RBD) reported, multiple GWAS-derived PRS, GBA/LRRK2 as
   Mendelian-style stratifiers
5. **Vascular dementia** (`F01`/`I67`) — small N, but interesting cross-
   comparison to AD

For each: this same document applies. Just substitute `<DISEASE>`.

---

## 12. Don'ts

- **Don't** modify analysis scripts to handle "this cohort is special"
  cases. If a script breaks under cfg variation, fix it generically and
  push the helper-level fix as a PR back to UKB_MS template.
- **Don't** commit the `data/` symlink. `.gitignore` covers `data/` (dir)
  but the bare `data` line (without slash) is required for the symlink itself.
- **Don't** push huge artifacts (>100MB). Already in `.gitignore`:
  `results/ml/models/`, `results/ml/models_sex/`, STRING cache.
- **Don't** invent PGS Catalog IDs, HLA alleles, or Mendelian variants. If
  the user can't supply them via OpenEvidence, leave the field empty and
  document the gap in the final report.
- **Don't** run UKB-RAP jobs locally. REGENIE GWAS belongs on RAP. Build
  the phenotype, document the RAP recipe, and stop.
- **Don't** modify the UKB_MS template (`upstream`). Helper-level fixes that
  affect all cohorts go as PRs to UKB_MS, raised from the replication repo.
- **Don't** leave hardcoded protein names, interaction names, or analysis
  decisions in scripts. Every value belongs in `disease.yaml`. If you find
  yourself writing `if (cfg$cohort_short == "stroke") ...`, **stop** — that's
  the wrong pattern. Generalise the logic and add a config field instead.

---

## Validated cohort capabilities

The pipeline has been validated end-to-end on:

- **MS** (k=3 endophenotype clusters; HLA-DRB1*15:01; 2 PGS Catalog scores)
- **ALS** (k=5 endophenotype clusters; HLA-DRB1*04 multi-subtype OR with
  *04:05 exclusion, currently skipped pending data extract; C9orf72 mendelian
  carrier scaffolded; smoking + head trauma phenomic scaffolds)

After MS→ALS:
- All 11 endophenotype-stage scripts (15/03–13) are k-agnostic
- Cluster contrast helper generates contrasts dynamically
- Sex stratification fully cfg-driven
- HLA carrier derivation supports single-allele AND multi-subtype-OR encodings
- Mendelian carrier stage parallel to HLA exists
- Phenomic stratifier stage iterates `cfg$phenomic_stratifiers`
- Top-N protein selection programmatic across all figure scripts

End of agent handoff document. Last revised after MS→ALS endophenotype
backport (UKB_MS commit 66f5518, UKB_ALS commit 9f82658).
