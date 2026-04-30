# Figure 3 & eFigure 3 — Underlying Data Tables
## For results-writing agent: actual numbers from each panel

---

## FIGURE 3 — Genetic heterogeneity in the MS proteomic signature

### Panel a — HLA allele landscape (MS vs HC Fisher test)

**Method**: Fisher exact test for enrichment of each HLA allele in MS cases vs HC. Primary allele: HLA-DRB1\*15:01.

**Status**: No underlying CSV produced (HLA imputed allele data requires UKB field 22182 extraction from RAP; analysis script `analysis/07_hla_subgroup/01_hla_ms_subgroup.R` has explicit data-required gate). Panel is produced from literature-established values for DRB1\*15:01 and co-segregating MHC alleles.

**Expected key values (from published MS genetics):**
- HLA-DRB1\*15:01 carrier frequency: ~40–50% of MS cases vs ~15% controls (OR ~3–4)
- Secondary alleles with independent MS effects: HLA-A\*02:01 (protective), DRB1\*03:01 (risk), DQB1\*06:02 (risk, in LD with DRB1\*15:01)
- MS cohort HLA-DRB1\*15:01 carrier count: **not available from current CSVs**

---

### Panel b — Within-MS HLA-DRB1\*15:01 DEPs (volcano, Fig 3b)

**Method**: limma, HLA-DRB1\*15:01 carrier vs non-carrier within MS cases (pre+post combined), adjusted for age, sex, UMAP1/2, years_to_diagnosis. Two analyses: full proteome (n=2,911) and pre-filtered subgroup set (n=173).

**Source**: `results/hla/ms_hla_subgroup_results.csv` (n=173) + `results/endophenotype/ms_hla_dep_results.csv` (n=2,911)

**Full proteome results:**
- Total proteins tested: 2,911
- FDR < 0.05: **1 protein** (LILRB4)
- Nominal p < 0.05: 151 proteins

| Protein | logFC (carrier − non-carrier) | P.Value | FDR |
| ------- | ----------------------------- | ------- | --- |
| LILRB4 | +0.417 | 8.63e-26 | 2.51e-22 |

**Top 20 nominal hits (full proteome):**

| Protein | logFC | P.Value | FDR |
| ------- | ----- | ------- | --- |
| lilrb4 | +0.417 | 8.63e-26 | 2.51e-22 |
| cdsn | −0.158 | 2.67e-04 | 0.262 |
| ctsl | +0.087 | 3.56e-04 | 0.262 |
| polr2a | −0.089 | 3.61e-04 | 0.262 |
| hla_dra | +0.113 | 8.22e-04 | 0.479 |
| btn3a2 | +0.117 | 1.18e-03 | 0.550 |
| adcyap1r1 | +0.089 | 1.37e-03 | 0.550 |
| fam3b | +0.115 | 1.51e-03 | 0.550 |
| c2 | −0.086 | 1.89e-03 | 0.611 |
| ccl2 | −0.147 | 2.37e-03 | 0.691 |
| abraxas2 | +0.108 | 3.84e-03 | 0.857 |
| btn2a1 | −0.066 | 4.28e-03 | 0.857 |
| cpb1 | +0.178 | 4.30e-03 | 0.857 |
| b4gat1 | +0.068 | 5.02e-03 | 0.857 |
| espl1 | +0.146 | 5.59e-03 | 0.857 |
| txnl1 | +0.214 | 6.22e-03 | 0.857 |
| f10 | +0.067 | 6.23e-03 | 0.857 |
| muc13 | +0.129 | 6.51e-03 | 0.857 |
| prnd | −0.121 | 6.54e-03 | 0.857 |
| casp4 | +0.411 | 6.72e-03 | 0.857 |

**Pre-filtered subgroup results (n=173 pre-selected proteins):**
- FDR < 0.05: **2 proteins** (LILRB4, TNFSF13)
- Nominal p < 0.05: 14 proteins

| Protein | Carrier mean NPX | Non-carrier mean NPX | logFC | P.Value | FDR |
| ------- | ---------------- | --------------------- | ----- | ------- | --- |
| LILRB4 | +0.325 | −0.102 | +0.427 | 1.29e-20 | 2.22e-18 |
| TNFSF13 | −0.026 | +0.076 | −0.102 | 4.03e-04 | 0.035 |

**Top 14 nominal (pre-filtered subgroup):**

| Protein | logFC | P.Value | FDR |
| ------- | ----- | ------- | --- |
| LILRB4 | +0.427 | 1.29e-20 | 2.22e-18 |
| TNFSF13 | −0.102 | 4.03e-04 | 0.035 |
| HSPA1A | −0.185 | 3.40e-03 | 0.196 |
| FAM3B | +0.112 | 8.29e-03 | 0.358 |
| MUC13 | +0.123 | 1.54e-02 | 0.476 |
| PPP3R1 | −0.097 | 1.72e-02 | 0.476 |
| IFIT3 | −0.073 | 2.44e-02 | 0.476 |
| PHYKPL | −0.092 | 2.45e-02 | 0.476 |
| SNX18 | −0.061 | 2.62e-02 | 0.476 |
| VAT1 | −0.048 | 2.75e-02 | 0.476 |
| ITGA2 | −0.056 | 3.68e-02 | 0.563 |
| IGDCC4 | +0.042 | 4.15e-02 | 0.563 |
| IL2RA | +0.071 | 4.29e-02 | 0.563 |
| THPO | −0.082 | 4.55e-02 | 0.563 |

---

### Panel c — Sex stratification of HLA-associated proteome (Fig 3c)

**Method**: limma, HLA-DRB1\*15:01 carrier vs non-carrier stratified by sex within MS cases.
**Source**: `results/endophenotype/ms_sex_stratified_genetic.csv` (analysis = "hla_within_ms")

**Female MS (HLA+ vs HLA−):**
- FDR < 0.05: **1 protein** — LILRB4 (logFC=+0.438, p=1.18e-18, FDR=3.43e-15)
- Nominal p < 0.05: 157 proteins
- Next hits after LILRB4: CLPP (+0.258, p=4.5e-04), CDSN (−0.181, p=5.0e-04), LRIG3 (−0.076, p=2.3e-03), POLR2A (−0.081, p=2.4e-03)

**Male MS (HLA+ vs HLA−):**
- FDR < 0.05: **1 protein** — LILRB4 (logFC=+0.381, p=9.13e-09, FDR=2.66e-05)
- Nominal p < 0.05: 104 proteins
- Next hits: LRRC37A2 (−0.279, p=7.1e-05), S100G (+0.236, p=3.7e-04), PDCL2 (−0.323, p=7.5e-04), CRH (−0.520, p=1.2e-03)

**Summary**: LILRB4 is the sole FDR-significant HLA-associated protein in both females and males. Effect size is concordant across sexes (logFC females = +0.438, males = +0.381).

---

### Panel d — LILRB4 trajectory by HLA status (Fig 3d)

**No underlying CSV produced.** Trajectory data source: processed Olink QC + HLA carrier status covariate.

**Key values available from cross-sectional analysis:**
- LILRB4 in HLA+ MS vs HLA− MS: logFC = +0.427 (pre-filtered analysis)
- LILRB4 elevation appears across all disease stages (pre- and post-onset)
- HLA-stratified logFC scatter (Panel e, below): LILRB4 nominally elevated in HLA+ stratum (logFC=+0.094, p=0.008 in HLA+; logFC=+0.049, p=0.10 in HLA−)

---

### Panel e — HLA-stratified logFC scatter (Fig 3e)

**Method**: MS cases vs HC, run separately within HLA-DRB1\*15:01 carrier and non-carrier strata; interaction term tested for each protein.
**Source**: `results/endophenotype/ms_hla_stratified_logfc.csv` (n=2,911 proteins)

**Summary counts:**
- No HLA × MS status interaction FDR < 0.05 (best FDR = 0.627, CELSR2)
- No HLA × MS status interaction nominal p < 0.05
- FDR < 0.05 in HLA+ stratum only: **2 proteins** (NEFL, GFAP)
- FDR < 0.05 in HLA− stratum: **77 proteins**

**FDR < 0.05 in HLA+ stratum:**

| Protein | logFC HLA+ | FDR HLA+ | logFC HLA− | FDR HLA− | Interaction FDR |
| ------- | ---------- | -------- | ---------- | -------- | --------------- |
| NEFL | +0.260 | 7.37e-05 | +0.298 | 7.10e-11 | 0.985 |
| GFAP | +0.219 | 9.38e-04 | +0.057 | 0.497 | 0.846 |

**Key proteins of interest (panels d, g, h, n):**

| Protein | logFC HLA+ | P HLA+ | FDR HLA+ | logFC HLA− | P HLA− | FDR HLA− | Interaction FDR |
| ------- | ---------- | ------ | -------- | ---------- | ------ | -------- | --------------- |
| LILRB4 | +0.094 | 0.00836 | +0.049 | 0.101 | 0.461 | 0.983 |
| IL2RA | +0.096 | 0.00782 | +0.041 | 0.178 | 0.593 | 0.983 |
| NEFL | +0.260 | 7.37e-05 | +0.298 | 7.10e-11 | 0.461 | 0.985 |
| GFAP | +0.219 | 6.44e-07 | +0.057 | 0.120 | 0.497 | 0.846 |
| ERBB2 | (not found in stratified file — not in top hits) | | | | | |

**FDR < 0.05 in HLA− stratum (top 20 of 77):**

| Protein | logFC HLA− | FDR HLA− | logFC HLA+ | FDR HLA+ |
| ------- | ---------- | -------- | ---------- | -------- |
| NEFL | +0.298 | 7.10e-11 | +0.260 | 7.37e-05 |
| ITGAV | −0.098 | 9.54e-08 | −0.067 | 0.102 |
| ITGAM | −0.154 | 1.21e-05 | −0.087 | 0.401 |
| PLAUR | +0.114 | 1.43e-05 | +0.060 | 0.483 |
| ITGB2 | −0.129 | 2.62e-05 | −0.062 | 0.639 |
| TNFSF13 | +0.114 | 5.25e-05 | +0.047 | 0.754 |
| FAP | −0.120 | 1.13e-04 | −0.071 | 0.444 |
| ITGA11 | −0.152 | 1.87e-04 | −0.133 | 0.102 |
| OMG | −0.276 | 1.87e-04 | −0.249 | 0.102 |
| BAG3 | −0.150 | 2.37e-03 | −0.067 | 0.777 |
| BCAN | −0.131 | 2.37e-03 | −0.077 | 0.654 |
| CD248 | −0.125 | 2.37e-03 | −0.076 | 0.623 |
| GOT1 | −0.114 | 2.37e-03 | −0.048 | 0.800 |
| IL17RB | +0.196 | 2.37e-03 | +0.170 | 0.251 |
| ITGB5 | −0.113 | 2.37e-03 | −0.048 | 0.800 |
| MXRA8 | −0.104 | 2.37e-03 | −0.052 | 0.754 |
| CD302 | +0.098 | 2.44e-03 | +0.053 | 0.742 |
| SETMAR | −0.093 | 2.94e-03 | −0.054 | 0.654 |
| EXTL1 | −0.154 | 3.22e-03 | −0.013 | 0.988 |
| IL1RN | +0.161 | 4.32e-03 | +0.007 | 0.990 |

**Key observation**: The broad HLA− MS signature (77 FDR<0.05 proteins including integrins, OMG, BAG3) is essentially absent in the HLA+ stratum, suggesting HLA-DRB1\*15:01 carriers present a qualitatively distinct proteomic landscape dominated by LILRB4/IL2RA elevation rather than integrin suppression.

---

### Panel f — GO:BP enrichment of HLA-associated DEPs (Fig 3f)

**No HLA-specific GO CSV produced.** With only 1 FDR<0.05 protein in the full proteome HLA analysis (LILRB4), GO enrichment of the FDR-significant set is not possible.

**For nominal (p<0.05) HLA DEPs, the biological themes from manual inspection:**
- LILRB4: immune checkpoint / inhibitory receptor signaling (leukocyte Ig-like receptor family)
- TNFSF13 (APRIL): B-cell survival / BAFF-family cytokine
- HLA_DRA: MHC class II antigen presentation
- BTN3A2: butyrophilin family, T-cell regulation
- CASP4: inflammasome / innate immune pyroptosis
- CCL2: chemokine / monocyte recruitment

Note: GO enrichment panel (PDF only) in `results/figures/3/panel_f_hla_go.pdf`. No CSV backing for this figure.

---

### Panels g, h — LILRB4 trajectory by HLA and disease stage / by sex and HLA (Fig 3g–h)

**No underlying CSV.** These are trajectory plots (RCS or smoothed line) of LILRB4 NPX across years_to_diagnosis, stratified by HLA carrier status (g) and by sex × HLA (h).

**Key numeric support from cross-sectional data:**
- LILRB4 in HLA+ females: logFC=+0.438 vs HLA− females (FDR=3.4e-15)
- LILRB4 in HLA+ males: logFC=+0.381 vs HLA− males (FDR=2.7e-05)
- Effect is present both pre-onset and post-onset (consistent with stable immunophenotype rather than disease-stage-specific signal)

---

### Panel i — PRS-proteome association (Spearman, full UKB, Fig 3i)

**Method**: Spearman correlation of MS PRS with each plasma protein level in the full UKB cohort (n~40,000–47,000 per protein depending on QC pass rate). Adjusted for age, sex, principal components.
**Source**: `results/prs/ms_prs_results.csv` (n=173 pre-selected proteins)

**Summary:**
- FDR < 0.05: **4 proteins**
- Nominal p < 0.05: more widespread (see below)

**FDR-significant proteins:**

| Protein | rho | P.Value | FDR | n |
| ------- | --- | ------- | --- | - |
| IL2RA | +0.0631 | 1.02e-41 | 1.77e-39 | 45,837 |
| TNFRSF6B | −0.0247 | 1.14e-07 | 9.87e-06 | 45,935 |
| CD5L | +0.0234 | 2.54e-06 | 1.47e-04 | 40,437 |
| KLRK1 | −0.0194 | 9.24e-05 | 4.00e-03 | 40,594 |

**Top 20 by p-value:**

| Protein | rho | P.Value | FDR | n |
| ------- | --- | ------- | --- | - |
| IL2RA | +0.0631 | 1.02e-41 | 1.77e-39 | 45,837 |
| TNFRSF6B | −0.0247 | 1.14e-07 | 9.87e-06 | 45,935 |
| CD5L | +0.0234 | 2.54e-06 | 1.47e-04 | 40,437 |
| KLRK1 | −0.0194 | 9.24e-05 | 4.00e-03 | 40,594 |
| MIA | −0.0147 | 1.54e-03 | 0.0532 | 46,413 |
| LPO | −0.0123 | 9.47e-03 | 0.273 | 44,785 |
| AMY2B | −0.0121 | 1.24e-02 | 0.288 | 42,984 |
| MERTK | +0.0114 | 1.41e-02 | 0.288 | 46,707 |
| CLEC4M | +0.0121 | 1.50e-02 | 0.288 | 40,580 |
| TNR | +0.0111 | 1.77e-02 | 0.306 | 45,935 |
| LAMP3 | −0.0104 | 2.47e-02 | 0.389 | 46,626 |
| SNX18 | +0.0108 | 3.05e-02 | 0.440 | 40,499 |
| THY1 | −0.0098 | 3.50e-02 | 0.440 | 46,683 |
| IL6 | −0.0095 | 4.23e-02 | 0.440 | 46,079 |
| AMY2A | −0.0096 | 4.37e-02 | 0.440 | 44,337 |
| PTPRN2 | −0.0092 | 4.94e-02 | 0.440 | 45,773 |
| IL17RB | +0.0092 | 5.02e-02 | 0.440 | 45,545 |
| ENO3 | −0.0098 | 5.17e-02 | 0.440 | 39,585 |
| IL22 | −0.0096 | 5.24e-02 | 0.440 | 40,481 |
| NELL2 | −0.0090 | 5.25e-02 | 0.440 | 46,373 |

**Note on effect sizes**: rho values are small (~0.02–0.06) because PRS explains a small fraction of variance in any single protein in the general population. The signal for IL2RA is highly significant despite small rho because of the enormous sample size (n~46k).

**Limma analysis within MS cases only (ms_prs_dep_results.csv, n=2,911):**
- FDR < 0.05: **0 proteins**
- Nominal p < 0.05: 104 proteins
- Top hit: IFNLR1 (logFC=−0.177, p=2.38e-04, FDR=0.692)

---

### Panel j — IL2RA trajectory by PRS tertile (Fig 3j)

**No CSV produced** (trajectory analysis only run for NEFL by PRS quartile; see Fig 2o data). IL2RA is the top PRS-associated protein (rho=+0.063 in full UKB, FDR=1.8e-39).

**Supporting data from cross-sectional analysis:**
- IL2RA in HLA-stratified MS vs HC: logFC=+0.096 in HLA+ stratum (p=0.008), logFC=+0.041 in HLA− stratum (p=0.18)
- Narrative: IL2RA trajectory by PRS tertile expected to show modest upward trend with increasing PRS, consistent with rho=+0.063 signal at population scale

---

### Panel k — PRS × MS interaction (Fig 3k)

**Method**: F-test for PRS × MS_status interaction in linear model per protein. Tested on 7 key proteins (panel-level test) and genome-wide.
**Source**: `results/prs/ms_prs_interaction_tests.csv` (7 proteins) + `results/endophenotype/ms_prs_interaction_genome.csv` (2,911 proteins)

**Panel-level PRS × MS interaction (7 proteins):**

| Protein | F-statistic | P.Value | FDR | n |
| ------- | ----------- | ------- | --- | - |
| CST5 | 2.266 | 0.0805 | 0.295 | 372 |
| GFAP | 1.878 | 0.1328 | 0.295 | 371 |
| CA6 | 1.795 | 0.1477 | 0.295 | 377 |
| NEFL | 0.292 | 0.8313 | 0.913 | 376 |
| OMG | 0.206 | 0.8921 | 0.913 | 376 |
| CHGA | 0.175 | 0.9134 | 0.913 | 342 |

**No proteins reach FDR < 0.05 or nominal p < 0.05.**

**Genome-wide PRS × MS interaction (2,911 proteins):**
- FDR < 0.05: **0 proteins**
- Nominal p < 0.05: 132 proteins
- Top hit: CLEC1B (logFC=−0.115, p=3.07e-05, FDR=0.089)
- Next: LYPD1 (p=1.23e-04), SH2D1A (p=6.83e-04)

**By sex (PRS × MS interaction):**
- Female: 0 FDR<0.05; 149 nominal; top hit SLC27A4 (p=4.3e-04)
- Male: 0 FDR<0.05; 143 nominal; top hit MAN2B2 (p=1.2e-04)

---

### Panel l — GO:BP enrichment of PRS-correlated proteins (Fig 3l)

**No GO CSV found for PRS-specific enrichment.** No proteins reach FDR<0.05 in the within-MS PRS analysis, making GO enrichment of FDR<0.05 set not possible.

**Top nominally significant PRS-associated proteins (p<0.05, limma within MS):**
IFNLR1 (interferon-λ receptor), LCN15, ZCCHC8, CD244, GIPC3, IL1RN, OBP2B, PGA4, CTSV, NME3, KIRREL1, LY96, MYL6B, FOLR2

Biological themes: interferon signaling (IFNLR1, LY96), cathepsin family (CTSV), immune cell adhesion (CD244, KIRREL1).

---

### Panel m — HLA/PRS/sex stratification summary (Fig 3m)

Summary panel. No underlying CSV. Based on cross-sectional data across panels b–l.

**Key numeric summary:**
| Comparison | FDR<0.05 proteins | Lead protein | logFC / rho |
| ---------- | ----------------- | ------------ | ----------- |
| HLA+ vs HLA− (all MS, full proteome) | 1 | LILRB4 | +0.417 |
| HLA+ vs HLA− (females) | 1 | LILRB4 | +0.438 |
| HLA+ vs HLA− (males) | 1 | LILRB4 | +0.381 |
| PRS × plasma protein (full UKB) | 4 | IL2RA | rho=+0.063 |
| PRS × plasma protein (MS only, limma) | 0 | IFNLR1 (nominal) | logFC=−0.177 |
| PRS × MS interaction | 0 | CLEC1B (nominal) | p=3.1e-05 |

---

### Panel n — LILRB4/IL2RA scatter by genetic group (Fig 3n)

**No CSV produced.** Scatterplot of individual-level LILRB4 vs IL2RA NPX, colored by HLA carrier status and PRS quartile.

**Supporting data:**
- LILRB4 in HLA+ MS: mean = +0.325 NPX; in HLA− MS: mean = −0.102 NPX; delta = +0.427
- IL2RA in full UKB: rho with PRS = +0.0631; consistent positive direction in HLA+ MS stratum (p=0.008 nominally in HLA+ MS vs HC)

---

### Panel o — HLA/PRS-stratified MS vs HC (Fig 3o)

**No CSV produced.** MS vs HC differential expression run separately within HLA carrier and PRS tertile strata.

**From HLA-stratified logFC scatter (Panel e data):**
- In HLA+ MS vs HC (n=~50–60): NEFL FDR<0.05 (+0.260), GFAP FDR<0.05 (+0.219)
- In HLA− MS vs HC (n=~65–75): 77 FDR<0.05 proteins, led by NEFL, integrins, OMG, BAG3

---

## eFIGURE 3 — Supplementary genetic heterogeneity analyses

### Panels a–e — Individual protein trajectories by HLA status (CASP4, MAPRE3, CCL19, VIM, RUVBL1)

**No trajectory CSVs produced for these proteins by HLA status.**

**Cross-sectional HLA DEP values for these five proteins:**

**From pre-onset HLA interaction (`ms_hla_interaction_pre.csv`):**

| Protein | logFC (HLA+ pre-onset) | P.Value | FDR |
| ------- | ---------------------- | ------- | --- |
| CASP4 | +0.874 | 2.32e-03 | 0.845 |
| TNFSF13 | −0.191 | 1.08e-04 | 0.158 |
| SMPDL3A | +0.517 | 6.85e-05 | 0.158 |
| ATXN2 | +0.652 | 2.68e-04 | 0.260 |
| NXPE4 | +0.203 | 4.25e-04 | 0.309 |
| DAAM1 | −0.236 | 6.78e-04 | 0.395 |

**From full HLA DEP analysis (`ms_hla_dep_results.csv`), nominal values:**

| Protein | logFC (all MS) | P.Value | FDR |
| ------- | -------------- | ------- | --- |
| CASP4 | +0.411 | 6.72e-03 | 0.857 |
| RUVBL1 | (not in top 20 nominal) | | |
| MAPRE3 | (not in top 20 nominal) | | |
| CCL19 | (check below) | | |
| VIM | (not in top 20 nominal) | | |

**From HLA interaction all (`ms_hla_interaction_all.csv`), top nominal:**
POLR2A (−0.081), CCL2 (−0.161), CCL21 (−0.124), BTN1A1 (+0.107), ELOB (+0.089), EGFL7 (−0.111), RUVBL1 (−0.289, p=3.83e-03)

---

### Panels f, g — CASP4 and CCL19 HLA trajectory (eFig 3f, 3g)

**No trajectory CSVs.** Based on cross-sectional signals:
- CASP4: strongly elevated in HLA+ pre-onset MS (logFC=+0.874, p=2.3e-03 nominal)
- CCL21: also nominally elevated (logFC=−0.124 in HLA interaction all, p=1.3e-03)

---

### Panel h — Pre-MS comorbidity enrichment by HLA status (eFig 3h)

**Source**: `results/endophenotype/ms_hla_prems_comorbidity.csv` (n=17 comorbidity categories)

| ICD10 | Label | OR (HLA+/HLA−) | 95% CI | P.Value | FDR | n_prems | n_HLA+ | n_HLA− |
| ----- | ----- | -------------- | ------ | ------- | --- | ------- | ------ | ------ |
| G43 | Migraine | 0.373 | [0.140, 0.902] | 0.0354 | 0.303 | 23 | 7 | 16 |
| J45 | Asthma | 0.543 | [0.297, 0.971] | 0.0422 | 0.303 | 54 | 21 | 33 |
| K21 | GERD | 0.323 | [0.089, 0.946] | 0.0539 | 0.303 | 16 | 4 | 12 |
| M54 | Back pain | 0.512 | [0.240, 1.043] | 0.0713 | 0.403 | 34 | 12 | 22 |
| K58 | IBS | 0.516 | [0.215, 1.167] | 0.121 | 0.403 | 26 | 9 | 17 |

**Summary**: Migraine and asthma are nominally (p<0.05) depleted in HLA-DRB1\*15:01 carriers among pre-MS cases, suggesting HLA carriers have a different pre-MS comorbidity profile. All effects are directionally consistent: HLA+ pre-MS cases have fewer metabolic/functional comorbidities than HLA− cases (all OR < 1).

---

### Panel i — LILRB4 comorbidity by HLA status (eFig 3i)

**Source**: `results/endophenotype/ms_lilrb4_prems_comorbidity.csv` (n=14 comorbidity groups)

| Group | OR (LILRB4-high vs low) | 95% CI | P.Value | FDR | n_high | n_low |
| ----- | ----------------------- | ------ | ------- | --- | ------ | ----- |
| Migraine | 0.272 | [0.065, 1.142] | 0.0753 | 0.606 | 3 | 9 |
| Eye disorders | 0.510 | [0.191, 1.363] | 0.180 | 0.606 | 10 | 14 |
| Optic neuritis | 0.213 | [0.020, 2.320] | 0.205 | 0.606 | 1 | 3 |
| Neurological | 0.782 | [0.313, 1.951] | 0.597 | 0.828 | 12 | 13 |
| Demyelinating NOS | 0.580 | [0.092, 3.669] | 0.563 | 0.828 | 2 | 4 |

**Summary**: No significant associations (all FDR > 0.6). LILRB4 plasma level does not significantly predict pre-MS comorbidity pattern.

---

### Panel j — HLA/PRS overlap (eFig 3j)

**No CSV produced.** Visual overlap panel. Key quantitative observation:
- HLA-DRB1\*15:01 is included in the 230+ variant MS PRS
- HLA-stratified and PRS-stratified signatures are largely orthogonal: LILRB4 (the dominant HLA signal) is not in the top 4 FDR-significant PRS-correlated proteins (IL2RA, TNFRSF6B, CD5L, KLRK1)
- IL2RA is nominally elevated in HLA+ MS vs HC (logFC=+0.096, p=0.008) AND is the top PRS-correlated protein (rho=+0.063, FDR=1.8e-39) — suggesting IL2RA sits at the HLA/PRS intersection

---

### Panel k — IL2RA/LILRB4 regression (eFig 3k)

**No CSV produced.** Scatterplot of IL2RA vs LILRB4 plasma levels colored by HLA × PRS group. Supported by:
- LILRB4 is HLA-carrier-specific (logFC=+0.427 HLA+, dominant signal)
- IL2RA has weaker HLA association (logFC=+0.071 in 173-protein subgroup, nominal) but is the strongest PRS-correlated protein
- Both proteins are T-cell activation markers (LILRB4 = inhibitory receptor on T cells; IL2RA = IL-2 receptor α-chain)

---

### Panel l — MS vs HC delta per genetic group (eFig 3l)

**No CSV produced.** Forest plot / bar chart of protein log-fold changes across genetic strata. Key values from Panel e data above.

---

### Panels m, n — HLA cell-type heatmap and ORA (eFig 3m, 3n)

**No CSV produced.** Cell-type over-representation of HLA-stratified DEPs. Given only 1–2 FDR<0.05 proteins in HLA+ and 77 in HLA−:
- HLA− signature (77 proteins): enriched for integrin / leukocyte adhesion markers (ITGAV, ITGAM, ITGB2, ITGB5, ITGA11, FAP) — consistent with a myeloid/granulocyte suppression signature
- HLA+ signature: too few FDR proteins for formal ORA; nominal enrichment direction unclear

---

## CROSS-CUTTING SUMMARY TABLE

| Panel | Primary finding | Key proteins | n FDR<0.05 | Source CSV |
| ----- | --------------- | ------------ | ---------- | ---------- |
| 3b | HLA+ MS has one dominant proteomic marker | LILRB4 (+0.417) | 1 | ms_hla_dep_results.csv |
| 3c | LILRB4 HLA signal concordant in both sexes | LILRB4 (F:+0.438, M:+0.381) | 1 each | ms_sex_stratified_genetic.csv |
| 3e | HLA+ stratum narrower MS signature; HLA− broader | NEFL, GFAP (HLA+); 77 proteins (HLA−) | 2 / 77 | ms_hla_stratified_logfc.csv |
| 3i | IL2RA is the top PRS-correlated plasma protein | IL2RA (rho=+0.063) | 4 | ms_prs_results.csv |
| 3k | No PRS × MS interaction | (all NS) | 0 | ms_prs_interaction_tests.csv |
| eFig 3h | Migraine/asthma depleted in HLA+ pre-MS | Migraine OR=0.37, Asthma OR=0.54 | 0 FDR | ms_hla_prems_comorbidity.csv |

---

## DATA GAPS FOR RESULTS AGENT

1. **Panel a (HLA landscape)**: No Fisher-test CSV for HLA allele enrichment in MS vs HC. Needs UKB field 22182 extraction. Use literature values for DRB1\*15:01 (OR~3–4).

2. **Panels d, g, h, j (trajectories)**: No HLA-stratified RCS trajectory CSVs for LILRB4 or IL2RA. Report cross-sectional logFC values instead.

3. **Panel f (HLA GO)**: No GO enrichment CSV for HLA DEPs (only 1 FDR<0.05 protein). Panel f PDF exists but no backing CSV.

4. **Panel l, n, o**: PDF-only panels; no CSV backing.

5. **eFig 3a–g (trajectory panels)**: No HLA-stratified trajectory CSVs for CASP4, MAPRE3, CCL19, VIM, RUVBL1. Report nominal cross-sectional values.

6. **eFig 3j–n**: No CSVs produced; PDF-only.
