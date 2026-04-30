# Figure 2 & eFigure 2 — Underlying Data Tables
## For results-writing agent: actual numbers from each panel

---

## FIGURE 2 — The temporal proteomic window in pre-diagnostic MS

### Panel a — Sampling window distribution

| Group | n | Blood draw window |
|---|---|---|
| Pre-onset MS | 125 | Up to ~8 years before diagnosis |
| Post-onset MS | 113 | Up to ~12 years after diagnosis |

Spline models span −8 to +12 years relative to diagnosis. Pre-onset cases cluster within 4 years of diagnosis (n=8 in 0–2 yr bin, n=25 in 2–5 yr bin, n=56 in 5–10 yr bin, n=36 >10 yr). No underlying CSV; visual histogram only.

---

### Panels b–d — CNS biomarker RCS trajectories (NEFL, GFAP, MOG)

**Method**: Restricted cubic spline (df=3) NPX vs years_to_diagnosis. Divergence = when prediction lower 95% CI (pred_lci) first exceeds HC mean (≈0 in NPX units). Source: `results/differential/cns_trajectories.csv`

**NEFL (Panel b) — Spline predictions at key timepoints:**

| years_to_dx | pred_npx | pred_lci | pred_uci | hc_mean | hc_lci | hc_uci |
|---|---|---|---|---|---|---|
| −8.00 | +0.176 | −0.112 | +0.463 | −0.006 | −0.485 | +0.473 |
| −6.00 | +0.180 | −0.007 | +0.368 | −0.006 | −0.485 | +0.473 |
| **−5.75** | **+0.181** | **+0.002** | **+0.360** | −0.006 | −0.485 | +0.473 |
| −5.00 | +0.182 | +0.021 | +0.344 | −0.006 | −0.485 | +0.473 |
| −4.00 | +0.184 | +0.026 | +0.342 | −0.006 | −0.485 | +0.473 |
| 0.00 | +0.185 | +0.011 | +0.359 | −0.006 | −0.485 | +0.473 |
| +5.00 | +0.283 | +0.153 | +0.413 | −0.006 | −0.485 | +0.473 |

**Divergence**: pred_lci first crosses above HC mean at **−5.75 years** before diagnosis (pred_lci = +0.002). NEFL rises continuously through the post-onset window. Replicates Abdelhak 2026 (~−6 yr).

**GFAP (Panel c):**
From `ms_pre_temporal.csv`: logFC_ms = +0.037, FDR = 0.951 (not significant pre-onset). No pre-onset divergence. Post-onset FDR<0.05: logFC = +0.175, FDR = 0.010. First divergence: **+5.0 years after diagnosis**. Astrocytic marker of established disease only.

**MOG (Panel d):**
logFC_ms = +0.072, FDR = 0.850; logFC_ytd = +0.0009, *P*_ytd = 0.908. **No divergence** from HC across the −8 to +12 year window.

---

### eFigure 2 / Panel context — ERBB2 trajectory

**ERBB2** — Spline predictions at key timepoints (from `cns_trajectories.csv`):

| years_to_dx | pred_npx | pred_lci | pred_uci | hc_mean |
|---|---|---|---|---|
| −8.00 | −0.174 | −0.289 | −0.060 | +0.001 |
| −7.00 | −0.177 | −0.269 | −0.084 | +0.001 |
| −5.00 | −0.180 | −0.250 | −0.109 | +0.001 |
| −3.00 | −0.175 | −0.248 | −0.103 | +0.001 |
| 0.00 | −0.155 | −0.236 | −0.074 | +0.001 |
| +3.00 | −0.095 | −0.189 | −0.001 | +0.001 |
| +5.00 | +0.006 | −0.093 | +0.106 | +0.001 |

**Key finding**: ERBB2 is suppressed throughout the entire pre-onset window (~−0.17 NPX below HC). Suppression is present at the −8 yr window boundary — the true onset likely predates the measurement window. Post-onset: ERBB2 logFC = +0.037 (not FDR-significant), near-complete reversal from pre-onset logFC = −0.107 (FDR = 0.010). Stage-specific biomarker.

---

### Panels e–h — Immune/inflammatory protein trajectories (BGN, KIR3DL1, OSM, OMG)

From `ms_pre_temporal.csv`:

| Protein | Panel | logFC_ms | FDR_ms | logFC_ytd (per yr) | *P*_ytd | gap_5yr |
|---|---|---|---|---|---|---|
| BGN | e | −0.242 | 0.850 | −0.022 | 0.414 | −0.351 |
| KIR3DL1 | f | −0.412 | 0.732 | +0.003 | 0.932 | −0.399 |
| OSM | g | +0.312 | 0.571 | −0.007 | 0.678 | +0.275 |
| OMG | h | −0.391 | 0.455 | +0.022 | 0.201 | −0.283 |

None show significant time × MS interaction (*P*_ytd > 0.10): their elevation/suppression relative to HC is relatively constant across the pre-diagnostic window. gap_5yr = estimated MS−HC difference at 5 years before diagnosis.

---

### Panel i — Pre-onset vs post-onset logFC concordance scatter

Source: `ms_pre_vs_hc.csv` and `ms_post_vs_hc.csv` (n=2,911 proteins each).

| Subset | n proteins | Pearson *r* | Same direction |
|---|---|---|---|
| All 2,911 proteins | 2,911 | 0.005 | 49.8% |
| Nominal *P*<0.05 in either stage | 670 | 0.052 | 54.3% |
| FDR<0.05 in either stage | 40 | **0.694** | 80.0% |
| Nominal *P*<0.05 in both stages | 50 | 0.647 | 78.0% |

**Key discordant proteins** (opposite direction, nominal *P*<0.05 in both stages):

| Protein | logFC pre-onset | logFC post-onset |
|---|---|---|
| HAGH | +0.147 | −0.171 |
| AK1 | +0.117 | −0.126 |
| KLK14 | +0.107 | −0.124 |
| KYNU | −0.129 | +0.087 |
| CD22 | −0.086 | +0.083 |
| ERBB2 | −0.107 | +0.037 |

**Interpretation**: Among the biologically meaningful (FDR-significant) proteins, the two stages are moderately concordant (r=0.694, 80% same direction). At the full-proteome level, pre- and post-onset signatures are essentially orthogonal (r=0.005), confirming that the pre-diagnostic proteome is not a diluted version of established disease. ERBB2 is the most biologically interpretable discordant marker.

---

### Panel k — ML ROC curves

Source: `ms_ml_results.csv`, `ms_preonset_ml_results.csv`

| Classifier | Model | AUC |
|---|---|---|
| Full cohort (all MS vs HC) | glmnet | **0.769** |
| Full cohort | GBM | 0.677 |
| Full cohort | Random forest | 0.639 |
| Pre-onset only (pre-MS vs HC) | glmnet | **0.569** |
| Pre-onset only | nnet | 0.658 |
| Pre-onset only | GBM | 0.559 |
| Pre-onset only | Random forest | 0.501 |

---

### Panels l — Pre-diagnostic risk score by time bin

Source: `results/ml/ms_preonset_risk_score_timebins.csv` (newly extracted)

HC reference: mean predicted probability = **0.447**, SD = 0.186 (n=47,801 controls).

| Time bin | n cases | Mean MS prob. | SE | SD |
|---|---|---|---|---|
| >10 yr before Dx | 36 | 0.550 | 0.029 | 0.172 |
| 5–10 yr before Dx | 56 | 0.585 | 0.025 | 0.191 |
| 2–5 yr before Dx | 25 | 0.566 | 0.039 | 0.197 |
| 0–2 yr before Dx | 8 | 0.596 | 0.095 | 0.269 |

**Key finding**: Pre-onset MS cases show elevated mean predicted probability above the HC mean (0.447) in all time bins, including those >10 years before diagnosis (mean 0.550). The risk score does not show a strong monotonic increase toward diagnosis — the 5–10 yr bin (0.585) is numerically higher than the 2–5 yr bin (0.566), though the 0–2 yr bin has only n=8 cases with a wide SE. The signal is detectable at population scale even in the most distal pre-diagnostic window.

---

### Panel m — SHAP / feature importance (pre-onset classifier)

Source: `ms_preonset_ml_feature_importance.csv`, `ms_preonset_boruta_features.csv`

**Boruta-confirmed features (10 proteins + 2 clinical):**

| Rank | Feature | Mean importance |
|---|---|---|
| 1 | NEFL | 88.3 |
| 2 | BAG3 | 69.4 |
| 3 | CD248 | 55.9 |
| 4 | LILRB4 | 51.2 |
| 5 | ARL13B | 49.5 |
| 6 | KLK6 | 39.9 |
| 7 | age_at_sampling | 39.3 |
| 8 | GALNT7 | 36.2 |
| 9 | OMG | 33.0 |
| 10 | CD302 | 25.1 |
| 11 | BCAN | 22.7 |
| 12 | sex | 22.2 |

NEFL dominates the pre-onset classifier. BAG3 (anti-apoptotic chaperone), CD248 (endosialin, vascular/fibroblast marker), and LILRB4 (inhibitory immune receptor) are the next most important features — none of which are the most prominent proteins in the simple volcano ranking.

---

### Panel n — Cox rolling hazard ratio

Source: `ms_protein_cox_results.csv`, `ms_protein_cox_rolling.csv`

**Overall Cox (full pre-onset window; 37 proteins FDR<0.05; n_events=119–122):**

| Protein | HR | 95% CI | *P*-value | FDR | Direction |
|---|---|---|---|---|---|
| NEFL | 1.567 | [1.349, 1.821] | 4.49×10⁻⁹ | 1.31×10⁻⁵ | Risk ↑ |
| ERBB2 | 0.715 | [0.637, 0.804] | 1.75×10⁻⁸ | 2.55×10⁻⁵ | Protective ↓ |
| SPOCK1 | 0.676 | [0.562, 0.813] | 3.32×10⁻⁵ | 0.021 | Protective ↓ |
| SETMAR | 0.739 | [0.641, 0.853] | 3.62×10⁻⁵ | 0.021 | Protective ↓ |
| MMP9 | 1.469 | [1.222, 1.766] | 4.19×10⁻⁵ | 0.021 | Risk ↑ |
| BGN | 0.651 | [0.528, 0.802] | 5.40×10⁻⁵ | 0.021 | Protective ↓ |
| FGL1 | 1.473 | [1.220, 1.779] | 5.69×10⁻⁵ | 0.021 | Risk ↑ |
| CLEC5A | 1.437 | [1.204, 1.714] | 5.79×10⁻⁵ | 0.021 | Risk ↑ |
| OSM | 1.429 | [1.198, 1.704] | 7.24×10⁻⁵ | 0.023 | Risk ↑ |
| FASLG | 0.708 | [0.597, 0.840] | 7.71×10⁻⁵ | 0.023 | Protective ↓ |
| GZMA | 0.706 | [0.591, 0.843] | 1.16×10⁻⁴ | 0.026 | Protective ↓ |
| RNASE6 | 1.250 | [1.116, 1.400] | 1.18×10⁻⁴ | 0.026 | Risk ↑ |
| PGLYRP1 | 1.342 | [1.152, 1.564] | 1.56×10⁻⁴ | 0.032 | Risk ↑ |
| ITGAV | 0.718 | [0.605, 0.853] | 1.66×10⁻⁴ | 0.032 | Protective ↓ |

**Rolling Cox HR by time bin (n cases: 0–4 yr=24, 4–8 yr=46, 8–12 yr=38, 12+ yr=17):**

| Protein | 0–4 yr HR | 4–8 yr HR | 8–12 yr HR | 12+ yr HR |
|---|---|---|---|---|
| NEFL | 1.70 (*P*=5.8×10⁻⁴) | 1.71 (*P*=1.1×10⁻⁶) | 1.41 (*P*=0.029) | 1.20 (*P*=0.53) |
| ERBB2 | 0.66 (*P*=1.1×10⁻⁶) | 0.83 (*P*=0.21) | 0.69 (*P*=5.1×10⁻⁵) | 0.76 (*P*=0.18) |
| OSM | 1.54 (*P*=0.032) | 1.41 (*P*=0.021) | 1.29 (*P*=0.13) | 1.73 (*P*=0.026) |
| FGL1 | 2.14 (*P*=4.3×10⁻⁵) | 1.14 (*P*=0.41) | 1.47 (*P*=0.034) | 1.54 (*P*=0.088) |
| CLEC5A | 1.51 (*P*=0.044) | 1.44 (*P*=0.012) | 1.36 (*P*=0.058) | 1.51 (*P*=0.11) |
| BGN | 0.81 (*P*=0.37) | 0.60 (*P*=0.004) | 0.77 (*P*=0.17) | 0.39 (*P*=0.002) |
| MMP9 | 1.47 (*P*=0.066) | 1.48 (*P*=0.011) | 1.39 (*P*=0.054) | 1.66 (*P*=0.051) |
| SPOCK1 | 0.51 (*P*=0.001) | 0.71 (*P*=0.025) | 0.69 (*P*=0.031) | 0.87 (*P*=0.59) |

**Key finding**: NEFL HR is significant and consistent across 0–4, 4–8, and 8–12 yr bins — the strongest and most sustained prospective hazard signal. ERBB2 is protective in the 0–4 yr and 8–12 yr bins (non-monotonic, possibly reflecting different biological phases). FGL1 is strongest at 0–4 yr. BGN protective effect is strongest at 12+ yr. OSM risk HR is consistent across all bins.

---

### Panel o — PRS-stratified NEFL trajectory

Source: `results/prs/ms_prs_nefl_trajectories.csv`, `results/prs/ms_prs_nefl_keytimes.csv` (newly extracted)

**Sample sizes per PRS quartile:** Q1 (low)=94, Q2=94, Q3=94, Q4 (high)=94 (n=376 total MS cases with PRS+NEFL).

**PRS × time interaction test (NEFL):** F=0.292, *P*=0.831 (not significant; FDR=0.913). PRS does not significantly modify the temporal trajectory of NEFL.

**Mean NEFL (NPX) at key timepoints by PRS quartile** (±1.5 yr window around each timepoint):

| Years to Dx | Q1 low (n) | Q2 (n) | Q3 (n) | Q4 high (n) |
|---|---|---|---|---|
| −8 yr | 0.170 (11) | 0.316 (14) | −0.019 (10) | 0.124 (11) |
| −5 yr | 0.115 (14) | 0.307 (7) | 0.141 (9) | −0.100 (8) |
| −2 yr | 0.562 (5) | 0.046 (15) | 0.659 (9) | 0.037 (4) |
| 0 yr | 0.289 (9) | 0.263 (9) | 0.281 (3) | −0.263 (5) |
| +2 yr | 0.564 (6) | −0.198 (4) | 0.242 (5) | 0.131 (9) |
| +5 yr | 0.409 (8) | 0.018 (10) | 0.061 (8) | 0.391 (10) |
| +10 yr | 0.083 (8) | −0.313 (5) | −0.055 (8) | −0.126 (5) |

**Key finding**: No consistent PRS-quartile ordering of NEFL levels across time bins. The PRS×time interaction is not statistically significant (*P*=0.831). Small per-cell n (3–15) limits power to detect an interaction. PRS genetic risk does not appear to modify the NEFL temporal trajectory in this cohort.

---

### Temporal protein × time interaction (Panel j heatmap context)

Top 17 proteins by significance of *years_to_diagnosis* interaction (*P*_ytd), from `ms_pre_temporal.csv`:

| Protein | logFC_ms | FDR_ms | logFC_ytd (per yr) | *P*_ytd | gap_5yr |
|---|---|---|---|---|---|
| PLAU | −0.218 | 0.183 | +0.029 | 1.11×10⁻⁴ | −0.074 |
| ARSB | −0.246 | 0.390 | +0.037 | 1.15×10⁻⁴ | −0.060 |
| APOH | +0.257 | 0.305 | −0.035 | 1.72×10⁻⁴ | +0.081 |
| GADD45B | −0.259 | 0.559 | +0.045 | 4.37×10⁻⁴ | −0.035 |
| IFNLR1 | −0.212 | 0.455 | +0.032 | 7.35×10⁻⁴ | −0.053 |
| PPP2R5A | +0.245 | 0.522 | −0.037 | 1.45×10⁻³ | +0.061 |
| B3GAT3 | +0.115 | 0.305 | −0.013 | 1.60×10⁻³ | +0.049 |
| RNASE6 | −0.083 | 0.765 | +0.020 | 1.63×10⁻³ | +0.018 |
| FGF23 | +0.443 | 0.183 | −0.044 | 2.21×10⁻³ | +0.223 |
| CD3E | +0.146 | 0.458 | −0.020 | 2.66×10⁻³ | +0.048 |
| MAGED1 | +0.220 | 0.559 | −0.033 | 2.74×10⁻³ | +0.058 |
| PARD3 | −0.195 | 0.735 | +0.042 | 3.55×10⁻³ | +0.015 |
| IL21R | +0.280 | 0.455 | −0.034 | 3.63×10⁻³ | +0.110 |
| IGLON5 | +0.487 | 0.455 | −0.060 | 4.18×10⁻³ | +0.185 |
| SLAMF7 | −0.274 | 0.570 | +0.040 | 6.08×10⁻³ | −0.075 |
| ECM1 | −0.182 | 0.455 | +0.022 | 6.10×10⁻³ | −0.071 |
| NEFL | +0.387 | 0.013 | −0.023 | 2.87×10⁻² | +0.274 |

NEFL has the only FDR<0.05 main MS effect (FDR=0.013) among these top interaction proteins, with a significant negative ytd interaction (*P*=0.029): levels rise as diagnosis approaches (gap narrows from 5 yr to 0 yr). All other proteins show temporal dynamics without reaching FDR significance for the baseline MS effect.

---

### Temporal GO enrichment (Panel j context)

From `ms_temporal_go_enrichment.csv`, **nominal *P*<0.05** terms by time bin:

| Time bin | GO term | *P*-value | query_size | intersect |
|---|---|---|---|---|
| 12+ yr | Regulation of EGFR signaling pathway | 0.024 | 6 | 2 |
| 12+ yr | Regulation of ERBB signaling pathway | 0.024 | 6 | 2 |
| 12+ yr | Positive reg. of ERBB signaling | 0.024 | 6 | 2 |
| 12+ yr | Release of cytochrome c from mitochondria | 0.024 | 6 | 2 |
| 12+ yr | Regulation of release of cytochrome c | 0.024 | 6 | 2 |
| 12+ yr | Positive reg. release of cytochrome c | 0.024 | 6 | 2 |
| 12+ yr | Regulation of DNA binding | 0.024 | 6 | 2 |
| 12+ yr | Positive regulation of DNA binding | 0.024 | 6 | 2 |
| 12+ yr | Apoptotic mitochondrial changes | 0.024 | 6 | 2 |

**Note**: Only the **12+ years** pre-diagnostic bin shows any enriched GO:BP terms at nominal *P*<0.05 (all driven by ERBB2-associated signaling and apoptotic pathways, 2/6 proteins in the gene set). Bins 0–4 yr, 4–8 yr, and 8–12 yr show no nominally significant GO enrichment — the temporal signal is too diffuse at the pathway level except in the most distal window, where ERBB2 suppression is the dominant biology.

---

---

## eFIGURE 2 — ERBB2 trajectory and additional CNS biomarker context

*(Reconciled per figures_detail.md: eFig2 covers ERBB2 trajectory + GFAP and MOG as supplementary context. figures.yaml updated accordingly.)*

### Panel a — ERBB2 trajectory (full disease course)

See ERBB2 spline table above (Panels b–d section). Key numbers:

| Metric | Value |
|---|---|
| Pre-onset logFC vs HC | −0.107 (FDR = 0.010) |
| Post-onset logFC vs HC | +0.037 (not FDR sig) |
| Direction reversal | Pre: suppressed → Post: neutral |
| Window boundary | −8 yr: pred_npx = −0.174, still clearly suppressed |
| True onset | Likely before −8 yr measurement boundary |

### Panel b — GFAP trajectory (extended post-onset window)

| Metric | Value |
|---|---|
| Pre-onset logFC vs HC | +0.037 (FDR = 0.951) |
| Post-onset logFC vs HC | +0.175 (FDR = 0.010) |
| First divergence | +5.0 yr after diagnosis |
| Pre-onset divergence | None detected |

### Panel c — MOG trajectory

| Metric | Value |
|---|---|
| logFC_ms | +0.072 (FDR = 0.850) |
| logFC_ytd (interaction) | +0.001 (*P* = 0.908) |
| Any divergence | None across −8 to +12 yr |

---

## Key numbers summary

| Quantity | Value | Source |
|---|---|---|
| NEFL first divergence | −5.75 yr before diagnosis | cns_trajectories.csv |
| GFAP first divergence | +5.0 yr after diagnosis | figures_detail |
| MOG divergence | None | cns_trajectories.csv |
| ERBB2 pre logFC | −0.107 (FDR 0.010) | ms_pre_vs_hc.csv |
| ERBB2 post logFC | +0.037 (NS) | ms_post_vs_hc.csv |
| Pre vs post r (FDR-sig proteins, n=40) | 0.694 | computed |
| Pre vs post r (nom sig both stages, n=50) | 0.647 | computed |
| Pre vs post same direction (FDR-sig) | 80% | computed |
| Pre-onset glmnet AUC | 0.569 | ms_preonset_ml_results.csv |
| Full-cohort glmnet AUC | 0.769 | ms_ml_results.csv |
| Boruta proteins | 10 (NEFL, BAG3, CD248, LILRB4, ARL13B, KLK6, GALNT7, OMG, CD302, BCAN) | ms_preonset_boruta_features.csv |
| Risk score HC mean prob | 0.447 (SD 0.186) | extracted |
| Risk score >10 yr bin | 0.550 (n=36) | ms_preonset_risk_score_timebins.csv |
| Risk score 5–10 yr bin | 0.585 (n=56) | ms_preonset_risk_score_timebins.csv |
| Risk score 2–5 yr bin | 0.566 (n=25) | ms_preonset_risk_score_timebins.csv |
| Risk score 0–2 yr bin | 0.596 (n=8) | ms_preonset_risk_score_timebins.csv |
| NEFL Cox HR (overall) | 1.567 [1.349–1.821], FDR=1.3×10⁻⁵ | ms_protein_cox_results.csv |
| ERBB2 Cox HR (overall) | 0.715 [0.637–0.804], FDR=2.6×10⁻⁵ | ms_protein_cox_results.csv |
| Cox FDR<0.05 proteins | 37 | ms_protein_cox_results.csv |
| PRS × NEFL interaction *P* | 0.831 (NS) | ms_prs_interaction_tests.csv |
| Temporal GO enrichment | Nominal P<0.05 in 12+ yr bin only (ERBB2/apoptosis) | ms_temporal_go_enrichment.csv |
