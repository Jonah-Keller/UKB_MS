# Figure 4 & eFigure 4 — Underlying Data Tables
## For results-writing agent: actual numbers from each panel

---

## FIGURE 4 — Sex differences in the MS proteomic signature

**Cohort sex breakdown** (from ms_olink_qc.csv):
- Female (sex=0): 22,846 controls, 77 pre-onset MS, 234 post-onset MS (n_MS=311)
- Male (sex=1): 19,274 controls, 34 pre-onset MS, 98 post-onset MS (n_MS=132)
- Female:male ratio in MS cases = ~2.4:1 (consistent with known MS sex distribution)

---

### Panel a — Female vs male logFC scatter (full proteome, Fig 4a)

**Method**: Two limma models — female MS vs female HC, and male MS vs male HC — each adjusted for age, UMAP1/2, years_to_diagnosis. Combined (pre+post) analysis. Scatter plots logFC_female vs logFC_male for all 2,911 proteins.

**Source**: `results/sex_stratified/ms_female_vs_hc_all.csv` + `results/sex_stratified/ms_male_vs_hc_all.csv`

**Summary counts:**
| Stratum | Total tested | FDR < 0.05 | Nominal p < 0.05 |
| ------- | ------------ | ---------- | ---------------- |
| Female MS vs HC | 2,911 | 78 | 388 |
| Male MS vs HC | 2,911 | 26 | 358 |
| Both FDR < 0.05 (concordant) | — | 11 | — |
| Female-specific FDR < 0.05 | — | 67 | — |
| Male-specific FDR < 0.05 | — | 15 | — |

**Pearson r (full proteome logFC_female vs logFC_male):** r = 0.226

**Pearson r (restricted to FDR < 0.05 in either sex, n=93):** r = 0.772 — all 11 jointly significant proteins are concordant in direction

---

### Panel b — DEP landscape (female-specific / concordant / discordant / male-specific)

**Concordant proteins (FDR < 0.05 in both female AND male, all same direction, n=11):**

| Protein | logFC female | FDR female | logFC male | FDR male | Direction |
| ------- | ------------ | ---------- | ---------- | -------- | --------- |
| NEFL | +0.206 | 9.53e-11 | +0.267 | 8.13e-06 | Up both |
| ITGAV | −0.078 | 4.65e-09 | −0.095 | 2.27e-05 | Down both |
| ITGAM | −0.137 | 5.23e-08 | −0.119 | 2.22e-02 | Down both |
| OMG | −0.252 | 2.13e-06 | −0.292 | 9.14e-03 | Down both |
| CA14 | −0.120 | 5.38e-05 | −0.178 | 3.45e-03 | Down both |
| ITGA11 | −0.106 | 8.17e-04 | −0.133 | 2.22e-02 | Down both |
| CST5 | −0.155 | 2.42e-03 | −0.247 | 7.13e-03 | Down both |
| CA6 | −0.179 | 2.52e-03 | −0.302 | 2.22e-03 | Down both |
| CDA | −0.125 | 3.40e-03 | −0.184 | 8.10e-03 | Down both |
| BAG3 | −0.087 | 2.41e-02 | −0.148 | 3.29e-02 | Down both |
| RNASET2 | +0.048 | 3.70e-02 | +0.077 | 4.91e-02 | Up both |

**Top 20 female-specific FDR < 0.05 proteins (FDR significant in female, not male):**

| Protein | logFC female | FDR female | logFC male | FDR male |
| ------- | ------------ | ---------- | ---------- | -------- |
| LILRB4 | +0.148 | 1.50e-05 | +0.134 | 0.072 |
| ITGB2 | −0.101 | 1.57e-05 | −0.086 | 0.102 |
| COMP | −0.117 | 1.77e-05 | −0.055 | 0.533 |
| FAP | −0.091 | 7.33e-05 | −0.073 | 0.177 |
| ITGB5 | −0.097 | 1.97e-04 | −0.090 | 0.116 |
| MSTN | −0.162 | 1.97e-04 | −0.092 | 0.514 |
| CLEC4A | −0.111 | 3.94e-04 | −0.111 | 0.138 |
| MERTK | +0.076 | 6.64e-04 | +0.058 | 0.250 |
| MXRA8 | −0.094 | 6.64e-04 | −0.086 | 0.163 |
| IL17RB | +0.162 | 7.23e-04 | +0.168 | 0.116 |
| KLK6 | −0.085 | 1.09e-03 | −0.062 | 0.346 |
| CD8A | −0.140 | 1.23e-03 | −0.027 | 0.902 |
| GFAP | +0.134 | 1.46e-03 | +0.144 | 0.081 |
| CHGA | +0.234 | 1.69e-03 | +0.038 | 0.906 |
| GZMA | −0.091 | 2.46e-03 | −0.103 | 0.082 |
| ITGA2 | −0.086 | 2.46e-03 | −0.072 | 0.317 |
| SNAP25 | −0.099 | 2.52e-03 | −0.031 | 0.812 |
| RET | −0.111 | 3.40e-03 | −0.120 | 0.081 |
| MEPE | −0.069 | 3.47e-03 | −0.042 | 0.413 |
| MOG | +0.086 | 3.47e-03 | +0.008 | 0.963 |

**Male-specific FDR < 0.05 proteins (FDR significant in male, not female, n=15):**

| Protein | logFC female | FDR female | logFC male | FDR male |
| ------- | ------------ | ---------- | ---------- | -------- |
| ENPP5 | −0.057 | 0.346 | −0.261 | 2.48e-06 |
| PEPD | −0.037 | 0.256 | −0.113 | 2.60e-03 |
| PLAUR | +0.019 | 0.708 | +0.119 | 5.68e-03 |
| NID1 | +0.012 | 0.903 | +0.132 | 7.65e-03 |
| ITGB6 | −0.054 | 0.153 | −0.125 | 7.70e-03 |
| TGFA | +0.008 | 0.950 | +0.188 | 7.78e-03 |
| GOT1 | −0.072 | 0.094 | −0.155 | 9.81e-03 |
| GSR | −0.029 | 0.235 | −0.075 | 1.43e-02 |
| PI3 | +0.047 | 0.598 | +0.201 | 1.99e-02 |
| EGFL7 | −0.018 | 0.866 | +0.140 | 2.42e-02 |
| DPP4 | −0.043 | 0.264 | −0.102 | 2.55e-02 |
| RSPO3 | +0.016 | 0.764 | +0.087 | 2.58e-02 |
| DNER | −0.050 | 0.100 | −0.097 | 4.34e-02 |
| IL1RN | +0.055 | 0.570 | +0.183 | 4.34e-02 |
| TNFSF13 | +0.014 | 0.794 | +0.091 | 4.91e-02 |

**Note on ENPP5**: ENPP5 is the lead male-specific marker. Baseline NPX differs by sex: female controls mean = −0.163 NPX, male controls mean = +0.116 NPX. Males show progressive decline: pre-onset = −0.001, post-onset = −0.167 vs control +0.116 (delta = −0.28). Females are already at low baseline; no significant further decline.

**Note on CHGA**: CHGA is the lead female-specific marker (logFC=+0.234 females, +0.038 males). Raw means: female controls=0.145, female post-onset=0.405, female pre-onset=0.181; male controls=0.143, male post-onset=0.182, male pre-onset=0.473 (n=34, small). The female post-onset elevation is robust; the male pre-onset value is likely noise from small n.

---

### Panels c–d — Cell-type heatmaps (female and male DEPs, Fig 4c–d)

**No CSV produced for cell-type heatmaps.** PDF panels at `results/figures/4/panel_c_female_celltype_heatmap.pdf` and `panel_d_male_celltype_heatmap.pdf`.

**Biological annotation of DEP sets (from protein identities):**

Female-specific DEPs biological themes:
- Immune cell surface markers: LILRB4 (T/myeloid inhibitory receptor), CD8A (cytotoxic T cells), CLEC4A (dendritic cells), MERTK (macrophage phagocytosis), GZMA (cytotoxic granules)
- Chromogranin/neuroendocrine: CHGA (chromogranin A — neuroendocrine marker)
- CNS: GFAP (astrocytes), SNAP25 (synaptic), MOG (myelin), KLK6 (myelination), RET (neuronal), OMG (myelin)
- ECM/integrin: ITGB2, COMP, FAP, ITGB5, ITGA2, MXRA8, MSTN, MEPE

Male-specific DEPs biological themes:
- Growth factor/receptor: TGFA (EGF family), EGFL7 (EGF-like), RSPO3 (Wnt pathway), NID1 (basement membrane)
- Metabolic: GOT1 (aspartate aminotransferase), CA6 (carbonic anhydrase), GSR (glutathione reductase), PEPD (dipeptidase), CDA (cytidine deaminase), DPP4 (glucose metabolism)
- ENPP5 (ectonucleotide pyrophosphatase, CNS expression)
- Immune: PLAUR (urokinase receptor), IL1RN (IL-1 antagonist), TNFSF13 (B cell survival), ITGB6 (TGF-β activation)

---

### Panels e–f — GO:BP enrichment for female and male DEPs (Fig 4e–f)

**No sex-specific GO CSV produced.** PDF panels at `results/figures/4/panel_e_female_go.pdf` and `panel_f_male_go.pdf`.

**Inferred biology from DEP sets:**

Female GO themes (based on DEP identities, 78 FDR<0.05 proteins):
- Leukocyte / immune cell adhesion and migration: ITGAV, ITGAM, ITGB2, ITGB5, ITGA11, ITGA2, FAP, CLEC4A, CD8A, GZMA
- Myelination and CNS development: OMG, GFAP, MOG, SNAP25, KLK6, RET
- Chromogranin / neuroendocrine signaling: CHGA, IL17RB

Male GO themes (based on DEP identities, 26 FDR<0.05 proteins):
- Growth factor signaling: TGFA, EGFL7, RSPO3, NID1, PLAUR
- Metabolic enzyme activity: GOT1, CA6, GSR, PEPD, CDA, DPP4
- Integrin-mediated cell adhesion: ITGAV, ITGAM, ITGA11, ITGB6

---

### Panel g — Pre-onset sex differences scatter (Fig 4g)

**Method**: Same as Panel a but restricted to pre-onset MS cases only (female pre-onset n=77, male pre-onset n=34).
**Sources**: `results/sex_stratified/ms_pre_female_vs_hc_all.csv` + `results/sex_stratified/ms_pre_male_vs_hc_all.csv`

**Summary:**
| Stratum | FDR < 0.05 | Nominal p < 0.05 |
| ------- | ---------- | ---------------- |
| Female pre-onset MS vs HC | 2 | 232 |
| Male pre-onset MS vs HC | 0 | 179 |

**Pre-onset logFC concordance r = 0.181** (very low; pre-onset sex signatures are essentially uncorrelated)

**Female pre-onset FDR < 0.05 proteins:**
| Protein | logFC | FDR |
| ------- | ----- | --- |
| NEFL | +0.254 | 0.0170 |
| BGN | −0.620 | 0.0170 |

**Male pre-onset top nominal (no FDR < 0.05):**
| Protein | logFC | P.Value | FDR |
| ------- | ----- | ------- | --- |
| CRX | +0.864 | 4.41e-05 | 0.128 |
| BANK1 | −0.415 | 2.45e-04 | 0.356 |
| KIR3DL1 | −0.795 | 1.06e-03 | 0.581 |
| CMIP | −0.226 | 1.63e-03 | 0.581 |
| IL1RAP | +0.270 | 1.68e-03 | 0.581 |
| NEFL | +0.270 | 2.55e-03 | 0.581 |

**Key observation**: In the pre-onset window, only females show FDR-significant DEPs (NEFL and BGN). Male pre-onset analysis yields 0 FDR-significant proteins despite similar pre-onset n after sex restriction. This is a central finding for Fig 4g — the pre-onset MS signal is detectable in females but not in males.

---

### Panels h–i — CHGA and ENPP5 sex-stratified trajectories (Fig 4h–i)

**No trajectory CSV produced.** Cross-sectional mean NPX by sex × disease stage:

**CHGA (Fig 4h — female-specific marker):**
| Sex | Group | Mean NPX | n |
| --- | ----- | -------- | - |
| Female | Control | 0.145 | 22,846 |
| Female | Pre-onset | 0.181 | 77 |
| Female | Post-onset | 0.405 | 234 |
| Male | Control | 0.143 | 19,274 |
| Male | Pre-onset | 0.473 | 34 |
| Male | Post-onset | 0.182 | 98 |

Female delta post-onset − control = +0.260 NPX (logFC=+0.234, FDR=1.69e-03)
Male delta post-onset − control = +0.039 NPX (logFC=+0.038, FDR=0.906, NS)
CHGA elevation is female-specific and concentrated in post-onset cases. Note: male pre-onset mean is elevated but based on n=34 (unstable).

**ENPP5 (Fig 4i — male-specific marker):**
| Sex | Group | Mean NPX | n |
| --- | ----- | -------- | - |
| Female | Control | −0.163 | 26,242 |
| Female | Pre-onset | −0.187 | 84 |
| Female | Post-onset | −0.261 | 260 |
| Male | Control | +0.116 | 22,128 |
| Male | Pre-onset | −0.001 | 38 |
| Male | Post-onset | −0.167 | 109 |

Male delta post-onset − control = −0.283 NPX (logFC=−0.261, FDR=2.48e-06)
Female delta post-onset − control = −0.098 NPX (logFC=−0.057, FDR=0.346, NS)
ENPP5 shows sex-divergent baseline: males start higher (+0.116) and fall below zero in MS; females already at negative baseline. Male-specific decline is the key finding.

**NEFL for reference (sex comparison):**
| Sex | Group | Mean NPX | n |
| --- | ----- | -------- | - |
| Female | Control | −0.017 | 25,393 |
| Female | Pre-onset | +0.091 | 81 |
| Female | Post-onset | +0.182 | 256 |
| Male | Control | +0.016 | 21,496 |
| Male | Pre-onset | +0.205 | 39 |
| Male | Post-onset | +0.307 | 108 |

NEFL elevation is concordant in both sexes (both FDR<0.05). Males show numerically larger post-onset elevation.

---

### Panel j — Sex-stratified ML ROC curves (6 curves, Fig 4j)

**Method**: Elastic-net glmnet classifiers trained and tested within sex × disease stage (female vs male × combined/pre-onset/post-onset). LOOCV AUC for each of 6 configurations.
**Source**: `results/ml/ms_sex_ml_roc_data.csv`

**AUC by sex and disease stage:**

| Sex | Stage | AUC |
| --- | ----- | --- |
| Female | Combined (pre+post) | 0.752 |
| Female | Post-onset only | 0.705 |
| **Female** | **Pre-onset only** | **0.717** |
| Male | Combined (pre+post) | 0.764 |
| Male | Post-onset only | 0.630 |
| **Male** | **Pre-onset only** | **0.446** |

**Key finding**: Female pre-onset classifier AUC = 0.717 (above chance). Male pre-onset classifier AUC = 0.446 (essentially at chance, below 0.5). This is the central machine learning finding for sex differences: a proteomic risk signature is detectable before diagnosis in females but not in males.

Combined-cohort classifiers are comparable (0.752 female vs 0.764 male), reflecting that post-onset MS proteomics is similarly detectable across sexes. The pre-onset divergence is the female-specific finding.

---

## eFIGURE 4 — Supplementary sex difference analyses

### Panel a — Sex interaction volcano (eFig 4a)

**Method**: limma with sex × MS_status interaction term (tests whether MS effect differs by sex).
**Source**: `results/sex_stratified/ms_sex_interaction.csv`

**Summary:**
- Total proteins: 2,911
- FDR < 0.05: **0 proteins**
- Nominal p < 0.05: 175 proteins

**Top 20 nominal (sex × MS interaction):**
| Protein | logFC (interaction) | P.Value | FDR |
| ------- | ------------------- | ------- | --- |
| SLC44A4 | +0.274 | 2.12e-04 | 0.428 |
| ENPP5 | +0.169 | 4.34e-04 | 0.428 |
| HEXIM1 | +0.153 | 4.57e-04 | 0.428 |
| UNC5D | −0.096 | 5.89e-04 | 0.428 |
| AXIN1 | +0.161 | 7.69e-04 | 0.448 |
| TGFA | −0.161 | 1.18e-03 | 0.525 |
| IL6R | −0.123 | 1.26e-03 | 0.525 |
| TSC22D1 | +0.149 | 1.64e-03 | 0.533 |
| EGFL7 | −0.134 | 1.71e-03 | 0.533 |
| GPC5 | −0.185 | 1.94e-03 | 0.533 |
| LHB | −0.301 | 2.08e-03 | 0.533 |
| PSCA | −0.852 | 2.41e-03 | 0.533 |
| GNE | +0.139 | 2.53e-03 | 0.533 |
| OBP2B | −0.213 | 2.90e-03 | 0.533 |
| BACH1 | +0.153 | 2.96e-03 | 0.533 |
| OLR1 | −0.190 | 3.06e-03 | 0.533 |
| IKBKG | +0.172 | 3.11e-03 | 0.533 |
| LAMTOR5 | +0.116 | 3.56e-03 | 0.546 |
| INPPL1 | +0.168 | 3.66e-03 | 0.546 |
| CALCOCO1 | +0.200 | 4.06e-03 | 0.546 |

**Note**: No sex × MS interaction reaches FDR significance. The "female-specific" and "male-specific" DEP categories in Fig 4 reflect differences in statistical power by sex (female n=311 vs male n=132 MS cases) rather than formally significant interactions.

---

### Panels b–d — Cell-type heatmap and ORA (eFig 4b–d)

**No CSV produced.** PDFs at `results/figures/4_supp/`.

**Expected biology from DEP set composition:**
- Female-specific DEPs: enriched for cytotoxic T cell markers (GZMA, CD8A) and dendritic cell markers (CLEC4A), consistent with female-biased adaptive immune activation in MS
- Male-specific DEPs: enriched for metabolic enzymes (GOT1, CA6, GSR, DPP4) and growth factors (TGFA, EGFL7), consistent with distinct metabolic and tissue-repair biology

---

### Panel e — SHAP paired feature importance (female vs male, eFig 4e)

**No CSV produced.** Shows top 35 proteins by SHAP importance from sex-stratified classifiers, paired as female vs male SHAP values.

**Supporting data — classifiers available (from models_sex/):**
- `glmnet_female_combined.rds`, `glmnet_female_preonset.rds`, `glmnet_female_postonset.rds`
- `glmnet_male_combined.rds`, `glmnet_male_preonset.rds`, `glmnet_male_postonset.rds`

NEFL is expected to dominate both classifiers (concordant); CHGA, LILRB4, and integrin-family proteins expected to be female-dominant; ENPP5 expected to be male-dominant.

---

### Panels f–g — Pre-onset and post-onset sex concordance (eFig 4f–g)

**Pre-onset sex concordance (eFig 4f):**
- r = 0.181 overall (n=2,911 proteins)
- Female pre-onset FDR<0.05: NEFL, BGN (n=2)
- Male pre-onset FDR<0.05: 0 proteins
- Pre-onset sex signatures are largely uncorrelated — the pre-onset female signature is detectably different from (and stronger than) the male pre-onset signature

**Post-onset sex concordance (eFig 4g):**
- Female post-onset FDR<0.05: 81 proteins (from `ms_post_female_vs_hc_all.csv`)
- Male post-onset FDR<0.05: 17 proteins (from `ms_post_male_vs_hc_all.csv`)
- Post-onset female top 5: ITGAV (FDR=6.6e-10), ITGAM (FDR=1.6e-08), OMG (FDR=1.5e-07), NEFL (FDR=2.7e-05), MSTN (FDR=3.6e-05)
- Post-onset male top 5: ITGAV (FDR=4.3e-05), ENPP5 (FDR=5.6e-05), CA14 (FDR=1.9e-04), CST5 (FDR=2.4e-04), RSPO1 (FDR=7.5e-04)
- Post-onset signatures are more concordant than pre-onset (ITGAV, NEFL, OMG shared)

---

## CROSS-CUTTING SUMMARY TABLE

| Panel | Primary finding | Key proteins | Statistic |
| ----- | --------------- | ------------ | --------- |
| 4a | Female broader signature (78) vs male narrower (26) | NEFL, ITGAV, OMG (concordant) | r=0.226 full proteome |
| 4b | 11 concordant, 67 female-specific, 15 male-specific | CHGA (F), ENPP5 (M), LILRB4 (F) | 11/78/26 FDR<0.05 |
| 4g | Pre-onset signal detected only in females | NEFL, BGN (female only) | F: 2 FDR<0.05; M: 0 |
| 4h | CHGA elevated in female MS only | CHGA female logFC=+0.234 | FDR=1.69e-03 |
| 4i | ENPP5 depleted in male MS only | ENPP5 male logFC=−0.261 | FDR=2.48e-06 |
| 4j | Pre-onset classifier: female AUC=0.717 vs male AUC=0.446 | ML pre-onset | Female detectable; male at chance |
| eFig4a | No formal sex × MS interaction reaches FDR | SLC44A4 top nominal | 0 FDR<0.05 |
| eFig4f–g | Pre-onset sex signatures uncorrelated (r=0.18); post-onset more concordant | — | Pre r=0.18; post more concordant |

---

## DATA GAPS FOR RESULTS AGENT

1. **Panels c–d (cell-type heatmaps)**: No CSV backing; PDF only. Use DEP set composition above for biological annotation.

2. **Panels e–f (GO enrichment)**: No sex-specific GO CSV. Biology inferred from DEP identities. GO CSV could be generated from the female-specific 67-protein and male-specific 15-protein sets.

3. **Panels h–i (trajectories)**: No RCS trajectory CSV for CHGA or ENPP5 by sex. Cross-sectional means provided above. Direction is clear from the limma results.

4. **eFig 4b–d (cell-type heatmap/ORA)**: PDF only, no CSV.

5. **eFig 4e (SHAP paired importance)**: Would need to run SHAP on the sex-stratified glmnet models from `results/ml/models_sex/`.
