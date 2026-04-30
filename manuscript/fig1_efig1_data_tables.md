# Figure 1 & eFigure 1 — Underlying Data Tables
## For results-writing agent: actual numbers from each panel

---

## FIGURE 1 — The pre-onset MS plasma proteomic signature

### Panel b — UMAP cohort overview
| Group | n |
|---|---|
| Pre-onset MS | 125 |
| Post-onset MS | 113 |
| Healthy controls (HC) | 47,801 |
| Olink proteins tested | 2,911 |

UMAP1/2 coordinates carried forward as limma covariates (age, sex also included).

---

### Panel c — Pre-onset volcano (pre-onset MS n=125 vs HC n=47,801)

**Summary**: 2,911 proteins tested. **FDR<0.05: 2 proteins** (NEFL↑, ERBB2↓). Nominal *P*<0.05: 342 proteins.

**Top 25 proteins by *P*-value:**

| Protein | logFC | *P*-value | FDR (adj.*P*) | Direction |
|---|---|---|---|---|
| NEFL | +0.2274 | 6.09×10⁻⁸ | **0.000177** | Up (FDR<0.05) |
| ERBB2 | −0.1074 | 7.06×10⁻⁶ | **0.0103** | Down (FDR<0.05) |
| PTPRC | −0.0747 | 1.11×10⁻⁴ | 0.0893 | Down (nominal) |
| SCLY | −0.1535 | 1.63×10⁻⁴ | 0.0893 | Down (nominal) |
| SETMAR | −0.0915 | 2.17×10⁻⁴ | 0.0893 | Down (nominal) |
| PI3 | +0.1922 | 2.52×10⁻⁴ | 0.0893 | Up (nominal) |
| BGN | −0.3956 | 2.52×10⁻⁴ | 0.0893 | Down (nominal) |
| OSM | +0.2597 | 3.07×10⁻⁴ | 0.0893 | Up (nominal) |
| KCTD5 | −0.1072 | 3.36×10⁻⁴ | 0.0893 | Down (nominal) |
| ARL13B | +0.1733 | 3.40×10⁻⁴ | 0.0893 | Up (nominal) |
| TGFA | +0.1614 | 3.62×10⁻⁴ | 0.0893 | Up (nominal) |
| MMP9 | +0.2243 | 3.86×10⁻⁴ | 0.0893 | Up (nominal) |
| SPOCK1 | −0.0937 | 3.99×10⁻⁴ | 0.0893 | Down (nominal) |
| OMG | −0.2384 | 4.92×10⁻⁴ | 0.0974 | Down (nominal) |
| CLEC5A | +0.0987 | 5.02×10⁻⁴ | 0.0974 | Up (nominal) |
| GBP1 | −0.1713 | 5.37×10⁻⁴ | 0.0977 | Down (nominal) |
| PKD1 | −0.0797 | 7.68×10⁻⁴ | 0.116 | Down (nominal) |
| CDA | −0.1344 | 7.70×10⁻⁴ | 0.116 | Down (nominal) |
| GZMA | −0.1164 | 7.86×10⁻⁴ | 0.116 | Down (nominal) |
| IL1RAP | +0.1541 | 8.59×10⁻⁴ | 0.116 | Up (nominal) |
| PGLYRP1 | +0.1394 | 8.79×10⁻⁴ | 0.116 | Up (nominal) |
| RNF41 | −0.1383 | 9.01×10⁻⁴ | 0.116 | Down (nominal) |
| HGFAC | −0.0793 | 9.42×10⁻⁴ | 0.116 | Down (nominal) |
| ADGRE5 | −0.0975 | 9.68×10⁻⁴ | 0.116 | Down (nominal) |
| FGL1 | +0.1928 | 9.97×10⁻⁴ | 0.116 | Up (nominal) |

---

### Panel d — Ranked barplot (top 20 by nominal *P*, pre-onset)

Same as panel c rows 1–20 above. Lead signals: NEFL (logFC +0.23, FDR 0.000177) and ERBB2 (logFC −0.11, FDR 0.0103).

---

### Panel e — Pre-onset GO:BP enrichment

**Input**: nominal *P*<0.05 pre-onset DEPs (342 proteins). Top enriched GO:BP terms:

| GO term | Direction | GeneRatio | Fold Enrichment | FDR | Proteins |
|---|---|---|---|---|---|
| Positive regulation of cell division | Up | 5/25 | 38.10 | 0.000143 | OSM/TGFA/PPBP/PDGFA/PDGFB |
| Regulation of cell division | Up | 5/25 | 19.85 | 0.00181 | OSM/TGFA/PPBP/PDGFA/PDGFB |
| Positive reg. of protein phosphorylation | Up | 5/25 | 16.33 | 0.00312 | OSM/MMP9/PDGFA/PDGFB/PLAUR |
| Positive regulation of phosphorylation | Up | 5/25 | 14.85 | 0.00370 | OSM/MMP9/PDGFA/PDGFB/PLAUR |
| Negative reg. of blood coagulation | Up | 3/25 | 43.52 | 0.00496 | PDGFA/PDGFB/PLAUR |
| Negative regulation of hemostasis | Up | 3/25 | 42.70 | 0.00496 | PDGFA/PDGFB/PLAUR |
| Negative regulation of coagulation | Up | 3/25 | 40.41 | 0.00496 | PDGFA/PDGFB/PLAUR |
| Positive reg. of phosphorus metabolic process | Up | 5/25 | 11.79 | 0.00496 | OSM/MMP9/PDGFA/PDGFB/PLAUR |
| Antimicrobial humoral response | Up | 4/25 | 16.40 | 0.00639 | PI3/PPBP/PGLYRP1/WFDC2 |
| Defense response to bacterium | Up | 5/25 | 10.45 | 0.00639 | PI3/PPBP/PGLYRP1/CLEC4D/WFDC2 |
| Regulation of blood coagulation | Up | 3/25 | 32.33 | 0.00639 | PDGFA/PDGFB/PLAUR |
| Regulation of hemostasis | Up | 3/25 | 31.43 | 0.00639 | PDGFA/PDGFB/PLAUR |
| Positive reg. of protein modification | Up | 5/25 | 10.06 | 0.00639 | OSM/MMP9/PDGFA/PDGFB/PLAUR |
| Reg. of protein autophosphorylation | Up | 2/25 | 116.06 | 0.00639 | PDGFA/PDGFB |
| PDGF receptor-beta signaling pathway | Up | 2/25 | 100.59 | 0.00752 | PDGFA/PDGFB |
| Negative regulation of wound healing | Up | 3/25 | 29.39 | 0.00639 | PDGFA/PDGFB/PLAUR |
| Regulation of phosphorylation | Up | 5/25 | 8.75 | 0.00904 | OSM/MMP9/PDGFA/PDGFB/PLAUR |

**Note**: All enriched terms are in **up-regulated** pre-onset DEPs. Pathway theme: growth factor/cytokine signaling (OSM, PDGFA/B), innate immune/antimicrobial (PI3, PGLYRP1), coagulation regulation (PDGFA/B, PLAUR). No FDR<0.05 GO terms in down-regulated pre-onset DEPs.

---

### Panel h — Pre-onset cell-type over-representation (Walchli/HPA)

**FDR<0.05 results only (all in down-regulated DEPs):**

| Cell type | Direction | Odds Ratio | FDR |
|---|---|---|---|
| Granulocytes | Down | 0.197 | 0.00275 |
| Macrophages | Down | 0.287 | 0.00733 |
| B-cells | Down | 0.270 | 0.00764 |
| Dendritic cells | Down | 0.308 | 0.0257 |

**Interpretation**: OR < 1 = depletion of cell-type markers among pre-onset down-regulated DEPs. Nominal (not FDR<0.05): monocytes (OR=0.437, FDR=0.105), microglial cells (OR=0.488, FDR=0.171). No significant enrichment (OR > 1) in up-regulated pre-onset DEPs.

---

### Panel i — Pre vs post-onset GO:BP pathway comparison heatmap

**Pre-onset enriched pathways (FDR<0.05)**: growth factor signaling (OSM, PDGF), innate immune/antimicrobial (PI3, PGLYRP1), coagulation regulation — all in **up-regulated** DEPs only.

**Post-onset enriched pathways (FDR<0.05, top terms)** for comparison:

| GO term | Direction | GeneRatio | Fold Enrichment | FDR | Proteins |
|---|---|---|---|---|---|
| Cell adhesion mediated by integrin | Down | 9/88 | 25.38 | 1.39×10⁻⁷ | ITGAV/ITGB2/ITGB5/ITGA11/ITGA2/RET/HRG/ITGB6/DPP4 |
| Cell-matrix adhesion | Down | 11/88 | 9.58 | 1.80×10⁻⁵ | ITGAV/ITGAM/ITGB2/ITGB5/ITGA11/ITGA2/TNXB/HRG/ITGB6/THY1/MYOC |
| Integrin-mediated signaling | Down | 8/88 | 14.65 | 4.47×10⁻⁵ | ITGAV/ITGAM/ITGB2/ITGB5/ITGA11/ITGA2/ITGB6/THY1 |
| Collagen metabolic process | Down | 7/88 | 12.93 | 2.99×10⁻⁴ | ITGA2/FAP/KLK6/TNXB/PEPD/NPPC/GOT1 |
| Reg. of activated T cell proliferation | Up | 4/36 | 45.56 | 1.34×10⁻³ | LILRB4/LGALS9/IGFBP2/IL2RA |
| Neuron projection regeneration | Down | 5/88 | 17.86 | 1.61×10⁻³ | OMG/RGMA/THY1/EPHA4/APOA4 |
| Negative reg. of TNF production | Up | 4/36 | 25.87 | 4.72×10⁻³ | LILRB4/LGALS9/ACP5/CX3CL1 |
| Positive regulation of cell development | Up | 7/36 | 8.13 | 4.72×10⁻³ | NEFL/GFAP/LILRB4/LGALS9/IL15/IL2RA/CX3CL1 |

**Key contrast pre vs post**: Upregulated DEPs shift from growth-factor/proliferation (pre: OSM, PDGF) to immune regulatory/T-cell (post: LILRB4, IL2RA, LGALS9). Downregulated DEPs shift from absent (pre: no FDR GO terms) to integrin/ECM adhesion (post: ITGAV, ITGAM, ITGB2 family).

---

---

## eFIGURE 1 — Full differential analysis results

### Panels A–E — Combined MS (pre+post, n=238) vs HC (n=47,801)

**Summary**: 2,911 proteins tested. **FDR<0.05: 184 proteins**. Nominal *P*<0.05: 502 proteins.

**Top 20 proteins by *P*-value:**

| Protein | logFC | *P*-value | FDR | Direction |
|---|---|---|---|---|
| ITGAV | −0.0859 | 4.84×10⁻²² | 1.41×10⁻¹⁸ | Down |
| NEFL | +0.2165 | 2.38×10⁻²¹ | 3.47×10⁻¹⁸ | Up |
| ITGAM | −0.1396 | 5.94×10⁻¹⁷ | 5.76×10⁻¹⁴ | Down |
| OMG | −0.2780 | 9.46×10⁻¹⁶ | 6.88×10⁻¹³ | Down |
| CA14 | −0.1408 | 5.51×10⁻¹⁴ | 3.21×10⁻¹¹ | Down |
| ITGB2 | −0.1067 | 7.29×10⁻¹³ | 3.54×10⁻¹⁰ | Down |
| CA6 | −0.2283 | 2.36×10⁻¹¹ | 9.61×10⁻⁹ | Down |
| CDA | −0.1375 | 2.67×10⁻¹¹ | 9.61×10⁻⁹ | Down |
| LILRB4 | +0.1383 | 2.97×10⁻¹¹ | 9.61×10⁻⁹ | Up |
| ITGA11 | −0.1174 | 8.50×10⁻¹¹ | 2.47×10⁻⁸ | Down |
| CST5 | −0.1808 | 4.19×10⁻¹⁰ | 1.11×10⁻⁷ | Down |
| MXRA8 | −0.0846 | 7.01×10⁻¹⁰ | 1.70×10⁻⁷ | Down |
| FAP | −0.0899 | 1.09×10⁻⁹ | 2.45×10⁻⁷ | Down |
| GFAP | +0.1350 | 2.87×10⁻⁹ | 5.97×10⁻⁷ | Up |
| ITGB5 | −0.0927 | 5.61×10⁻⁹ | 1.09×10⁻⁶ | Down |
| COMP | −0.0977 | 1.06×10⁻⁸ | 1.94×10⁻⁶ | Down |
| RET | −0.1184 | 1.88×10⁻⁸ | 3.21×10⁻⁶ | Down |
| IL17RB | +0.1580 | 2.06×10⁻⁸ | 3.33×10⁻⁶ | Up |
| MERTK | +0.0729 | 2.41×10⁻⁸ | 3.69×10⁻⁶ | Up |
| CLEC4A | −0.1070 | 2.65×10⁻⁸ | 3.86×10⁻⁶ | Down |

**Combined GO:BP top terms (panel C):**

| GO term | Direction | GeneRatio | Fold Enrich. | FDR | Key proteins |
|---|---|---|---|---|---|
| Cell adhesion mediated by integrin | Down | 11/108 | 25.28 | 1.18×10⁻⁹ | ITGAV/ITGB2/ITGA11/ITGB5/RET/ITGA2/ITGB6/ADA/HRG/DPP4/ITGB1 |
| Cell-matrix adhesion | Down | 14/108 | 9.94 | 1.46×10⁻⁷ | ITGAV/ITGAM/ITGB2/ITGA11/ITGB5/ITGA2/TNXB/ITGB6/THY1/MYOC/HRG/CD34/ITGB1/L1CAM |
| Symbiont entry into host cell | Down | 11/108 | 14.12 | 2.52×10⁻⁷ | ITGAV/ITGB5/ITGA2/ITGB6/CLEC4M/HSPA1A/DPP4/ITGB1/EGFR/CR1/CD209 |
| Cell-substrate adhesion | Down | 15/108 | 7.10 | 1.30×10⁻⁶ | ITGAV/ITGAM/ITGB2/ITGA11/ITGB5/ITGA2/TNXB/ITGB6/THY1/MYOC/HRG/CD34/ITGB1/L1CAM/SPOCK1 |
| Collagen metabolic process | Down | 8/108 | 12.04 | 8.21×10⁻⁵ | FAP/PEPD/ITGA2/KLK6/GOT1/TNXB/NPPC/ITGB1 |
| Positive reg. of leukocyte activation | Up | 10/57 | 8.68 | 1.66×10⁻⁴ | LILRB4/IL15/IGFBP2/LGALS9/CD83/IL2RA/TNFSF13/TREM2/CTSC/IL6 |
| Regulation of microglial cell activation | Up | 4/57 | 69.66 | 1.66×10⁻⁴ | TREM2/CTSC/CX3CL1/IL6 |
| Positive reg. of cell activation | Up | 10/57 | 8.25 | 1.66×10⁻⁴ | LILRB4/IL15/IGFBP2/LGALS9/CD83/IL2RA/TNFSF13/TREM2/CTSC/IL6 |
| Myeloid leukocyte activation | Up | 8/57 | 10.38 | 2.47×10⁻⁴ | LILRB4/CHGA/IL15/LGALS9/TREM2/CTSC/CX3CL1/IL6 |
| Interleukin-1β production | Up | 6/57 | 17.11 | 2.47×10⁻⁴ | LILRB4/LGALS9/ACP5/TREM2/CX3CL1/IL6 |

**Combined cell-type ORA FDR<0.05 results (panel E):**

| Cell type | Direction | Odds Ratio | FDR |
|---|---|---|---|
| Dendritic cells | Down | 0.252 | 1.39×10⁻⁷ |
| B-cells | Down | 0.290 | 2.31×10⁻⁷ |
| Plasma cells | Down | 0.282 | 1.39×10⁻⁷ |
| Granulocytes | Down | 0.294 | 1.39×10⁻⁷ |
| Monocytes | Down | 0.316 | 6.69×10⁻⁷ |
| Macrophages | Down | 0.364 | 2.40×10⁻⁶ |
| NK-cells | Down | 0.358 | 1.29×10⁻⁵ |
| T-cells | Down | 0.407 | 4.28×10⁻⁵ |
| Endothelial | Down | 0.383 | 6.87×10⁻⁵ |
| Smooth muscle | Down | 0.412 | 2.19×10⁻³ |

---

### Panels F–J — Post-onset MS (n=113) vs HC (n=47,801)

**Summary**: 2,911 proteins tested. **FDR<0.05: 39 proteins**. Nominal *P*<0.05: 378 proteins.

**Top 20 proteins by *P*-value:**

| Protein | logFC | *P*-value | FDR | Direction |
|---|---|---|---|---|
| NEFL | +0.3182 | 5.09×10⁻¹³ | 1.48×10⁻⁹ | Up |
| ITGAV | −0.0952 | 1.03×10⁻⁸ | 1.50×10⁻⁵ | Down |
| ITGAM | −0.1589 | 2.46×10⁻⁷ | 2.38×10⁻⁴ | Down |
| OMG | −0.3103 | 1.79×10⁻⁶ | 1.30×10⁻³ | Down |
| TNFSF13 | +0.1083 | 7.33×10⁻⁶ | 4.27×10⁻³ | Up |
| LGALS9 | +0.1137 | 1.92×10⁻⁵ | 9.31×10⁻³ | Up |
| MXRA8 | −0.1064 | 2.78×10⁻⁵ | 1.02×10⁻² | Down |
| CA6 | −0.2622 | 3.27×10⁻⁵ | 1.02×10⁻² | Down |
| GFAP | +0.1748 | 3.36×10⁻⁵ | 1.02×10⁻² | Up |
| COL15A1 | −0.0871 | 3.51×10⁻⁵ | 1.02×10⁻² | Down |
| PRKAB1 | +0.1532 | 4.31×10⁻⁵ | 1.14×10⁻² | Up |
| PRR4 | −0.2391 | 4.87×10⁻⁵ | 1.18×10⁻² | Down |
| IL6 | +0.2872 | 7.12×10⁻⁵ | 1.50×10⁻² | Up |
| LTO1 | +0.1614 | 7.23×10⁻⁵ | 1.50×10⁻² | Up |
| CA14 | −0.1392 | 8.01×10⁻⁵ | 1.55×10⁻² | Down |
| HPGDS | −0.1424 | 9.10×10⁻⁵ | 1.66×10⁻² | Down |
| TREM2 | +0.1896 | 1.02×10⁻⁴ | 1.67×10⁻² | Up |
| CLEC7A | +0.2188 | 1.03×10⁻⁴ | 1.67×10⁻² | Up |
| CRKL | +0.1302 | 1.10×10⁻⁴ | 1.69×10⁻² | Up |
| DBNL | +0.1007 | 1.18×10⁻⁴ | 1.71×10⁻² | Up |

**Post-onset GO:BP top terms (panel H):**

| GO term | Direction | GeneRatio | Fold Enrich. | FDR | Key proteins |
|---|---|---|---|---|---|
| Cell adhesion mediated by integrin | Down | 9/88 | 25.38 | 1.39×10⁻⁷ | ITGAV/ITGB2/ITGB5/ITGA11/ITGA2/RET/HRG/ITGB6/DPP4 |
| Cell-matrix adhesion | Down | 11/88 | 9.58 | 1.80×10⁻⁵ | ITGAV/ITGAM/ITGB2/ITGB5/ITGA11/ITGA2/TNXB/HRG/ITGB6/THY1/MYOC |
| Integrin-mediated signaling | Down | 8/88 | 14.65 | 4.47×10⁻⁵ | ITGAV/ITGAM/ITGB2/ITGB5/ITGA11/ITGA2/ITGB6/THY1 |
| Collagen metabolic process | Down | 7/88 | 12.93 | 2.99×10⁻⁴ | ITGA2/FAP/KLK6/TNXB/PEPD/NPPC/GOT1 |
| Reg. of activated T cell proliferation | Up | 4/36 | 45.56 | 1.34×10⁻³ | LILRB4/LGALS9/IGFBP2/IL2RA |
| Activated T cell proliferation | Up | 4/36 | 42.77 | 1.34×10⁻³ | LILRB4/LGALS9/IGFBP2/IL2RA |
| Neuron projection regeneration | Down | 5/88 | 17.86 | 1.61×10⁻³ | OMG/RGMA/THY1/EPHA4/APOA4 |
| Negative reg. of TNF production | Up | 4/36 | 25.87 | 4.72×10⁻³ | LILRB4/LGALS9/ACP5/CX3CL1 |
| Positive regulation of cell development | Up | 7/36 | 8.13 | 4.72×10⁻³ | NEFL/GFAP/LILRB4/LGALS9/IL15/IL2RA/CX3CL1 |

**Post-onset cell-type ORA FDR<0.05 (panel J):**

| Cell type | Direction | Odds Ratio | FDR |
|---|---|---|---|
| Plasma cells | Down | 0.157 | 3.08×10⁻¹⁰ |
| B-cells | Down | 0.180 | 2.72×10⁻⁹ |
| Monocytes | Down | 0.193 | 2.72×10⁻⁹ |
| Granulocytes | Down | 0.215 | 9.17×10⁻⁹ |
| Macrophages | Down | 0.241 | 1.00×10⁻⁸ |
| NK-cells | Down | 0.240 | 1.78×10⁻⁷ |
| Dendritic cells | Down | 0.206 | 2.48×10⁻⁸ |
| T-cells | Down | 0.296 | 1.01×10⁻⁶ |
| Endothelial | Down | 0.289 | 8.97×10⁻⁶ |
| Smooth muscle | Down | 0.273 | 1.42×10⁻⁴ |
| Oligodendrocytes | Down | 0.522 | 1.40×10⁻² |

---

### Panels K–N — Pre-onset vs post-onset (within MS, direct contrast)

**Summary**: 2,912 proteins tested. **FDR<0.05: 3 proteins**. Nominal *P*<0.05: 232 proteins.
Positive logFC = higher in pre-onset relative to post-onset (i.e., pre-elevated or post-depleted).

**Top 20 proteins by *P*-value:**

| Protein | logFC | *P*-value | FDR | Interpretation |
|---|---|---|---|---|
| FGF23 | +0.2767 | 2.45×10⁻⁶ | **0.00713** | Higher pre-onset (FDR<0.05) |
| CLC | +0.3329 | 1.28×10⁻⁵ | **0.01870** | Higher pre-onset (FDR<0.05) |
| CLEC5A | +0.1361 | 3.61×10⁻⁵ | **0.0350** | Higher pre-onset (FDR<0.05) |
| FST | +0.2232 | 2.53×10⁻⁴ | 0.145 | Higher pre-onset (nominal) |
| S100A11 | +0.1777 | 3.10×10⁻⁴ | 0.145 | Higher pre-onset (nominal) |
| PPY | +0.4506 | 3.42×10⁻⁴ | 0.145 | Higher pre-onset (nominal) |
| ACP5 | −0.1354 | 4.79×10⁻⁴ | 0.145 | Higher post-onset (nominal) |
| OSM | +0.3398 | 6.03×10⁻⁴ | 0.145 | Higher pre-onset (nominal) |
| MMP8 | +0.3202 | 7.22×10⁻⁴ | 0.145 | Higher pre-onset (nominal) |
| RETN | +0.1666 | 7.20×10⁻⁴ | 0.145 | Higher pre-onset (nominal) |
| GFAP | −0.1893 | 9.43×10⁻⁴ | 0.145 | Higher post-onset (nominal) |
| CLEC7A | −0.1710 | 1.28×10⁻³ | 0.167 | Higher post-onset (nominal) |
| TREM2 | −0.1609 | 1.75×10⁻³ | 0.203 | Higher post-onset (nominal) |
| LTO1 | −0.1397 | 2.30×10⁻³ | 0.237 | Higher post-onset (nominal) |

**Key stage-specific signals**: GFAP (logFC = −0.189, *P* = 9.43×10⁻⁴) is nominally higher in post-onset, consistent with its post-diagnosis trajectory. NEFL does not appear in the pre vs post contrast top hits (similar elevation in both stages). FGF23, CLC, CLEC5A are significantly higher pre-onset than post-onset.
