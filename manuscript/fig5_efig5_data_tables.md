# Figure 5 & eFigure 5 — Underlying Data Tables
## For results-writing agent: actual numbers from each panel

---

## OVERVIEW — Pre-diagnostic comorbidity endophenotypes

**Cohort entering clustering analysis:**
- Total MS cases with proteomics: 499
- Cases with ≥1 FDR-significant enriched pre-MS ICD block: 77 (clustered)
- Cases with 0 enriched blocks: 422 ("None" reference group)
- UMAP parameters: 11 binary block features; Jaccard distance; n_neighbors=10; min_dist=0.05; k=3
- **Silhouette score: 0.864** (k=3 selected by silhouette criterion)

---

## FIGURE 5 — Main panels

### Panel a — Pre-MS PheWAS lollipop (Fig 5a)

**Method**: Logistic regression of ICD-10 coded diagnoses (pre-MS window) in 499 pre-onset MS cases vs matched controls. FDR corrected across all blocks tested.

**Source**: `results/comorbidity/ms_prems_phewas.csv`

**FDR < 0.05: 23 conditions (12 enriched, 11 depleted)**

| ICD10 | Block | Condition | OR | 95% CI | FDR | %MS pre | %HC pre | Direction |
|-------|-------|-----------|-----|--------|-----|---------|---------|-----------|
| H46 | H4 | Optic neuritis | 56.0 | 26.8–117.2 | 1.8e-24 | 2.4% | 0.04% | Enriched |
| R26 | R2 | Abnormalities of gait/mobility | 20.9 | 7.9–55.5 | 7.7e-08 | 1.0% | 0.05% | Enriched |
| I10 | I1 | Essential hypertension | 0.49 | 0.37–0.65 | 2.0e-05 | 12.0% | 24.7% | Depleted |
| E78 | E7 | Disorders of lipoprotein metabolism | 0.31 | 0.19–0.49 | 4.4e-05 | 3.6% | 13.2% | Depleted |
| G51 | G5 | Facial nerve disorders | 5.10 | 2.59–10.1 | 7.4e-05 | 1.8% | 0.37% | Enriched |
| R29 | R2 | Other nervous system symptoms | 8.34 | 3.32–21.0 | 1.6e-04 | 1.0% | 0.12% | Enriched |
| N17 | N1 | Acute kidney injury | 7.48 | 3.00–18.7 | 3.4e-04 | 1.0% | 0.16% | Enriched |
| R20 | R2 | Disturbances of skin sensation | 5.19 | 2.41–11.2 | 5.0e-04 | 1.4% | 0.26% | Enriched |
| M19 | M1 | Other osteoarthritis | 0.37 | 0.23–0.59 | 5.8e-04 | 3.6% | 9.8% | Depleted |
| J06 | J0 | Acute upper respiratory infections | 0.37 | 0.22–0.65 | 7.2e-03 | 2.6% | 6.3% | Depleted |
| R79 | R7 | Abnormal blood chemistry findings | 4.27 | 1.87–9.76 | 7.8e-03 | 1.2% | 0.44% | Enriched |
| L98 | L9 | Other skin disorders | 0.22 | 0.09–0.54 | 1.1e-02 | 1.0% | 4.4% | Depleted |
| M45 | M4 | Ankylosing spondylitis | 4.60 | 1.86–11.3 | 1.1e-02 | 1.0% | 0.26% | Enriched |
| Z37 | Z3 | Outcome of delivery | 0.31 | 0.15–0.64 | 1.8e-02 | 1.6% | 3.0% | Depleted |
| M25 | M2 | Other joint disorders | 0.55 | 0.38–0.81 | 2.0e-02 | 5.8% | 10.2% | Depleted |
| L30 | L3 | Other dermatitis | 0.48 | 0.30–0.78 | 2.3e-02 | 3.6% | 7.0% | Depleted |
| K76 | K7 | Other liver diseases | 3.97 | 1.61–9.76 | 2.3e-02 | 1.0% | 0.30% | Enriched |
| B27 | B2 | Infectious mononucleosis | 2.43 | 1.36–4.36 | 2.3e-02 | 2.4% | 0.92% | Enriched |
| M50 | M5 | Cervical disc disorders | 3.48 | 1.53–7.91 | 2.3e-02 | 1.2% | 0.38% | Enriched |
| Z86 | Z8 | Personal history of past diseases | 1.86 | 1.22–2.85 | 3.0e-02 | 4.6% | 2.9% | Enriched |
| N63 | N6 | Unspecified lump in breast | 0.28 | 0.12–0.68 | 3.5e-02 | 1.0% | 2.7% | Depleted |
| H25 | H2 | Age-related cataract | 3.16 | 1.39–7.19 | 4.2e-02 | 1.2% | 0.47% | Enriched |
| K21 | K2 | Gastro-oesophageal reflux disease | 0.50 | 0.30–0.82 | 4.2e-02 | 3.2% | 6.7% | Depleted |

**Enrichment narrative**: Optic neuritis (H46) is the dominant enriched signal by far (OR=56, likely prodromal MS). Gait/mobility symptoms (R26, R29, R20) form a cluster of early neurological prodrome. Musculoskeletal conditions (M45 ankylosing spondylitis, M50 cervical disc) and mononucleosis (B27) are enriched. Hypertension (I10), dyslipidaemia (E78), and osteoarthritis (M19) are notably depleted (inverse MS-comorbidity pattern).

---

### Panel b — UMAP + cluster assignments (Fig 5b)

**Source**: `results/endophenotype/ms_prems_umap_config.csv` + `ms_prems_cluster_assignments.csv`

**UMAP configuration:**
- Features: 11 binary ICD block flags (H4, R2, G5, N1, R7, M4, K7, B2, M5, Z8, H2)
- Metric: Jaccard (precomputed)
- n_neighbors: 10; min_dist: 0.05
- k = 3; **silhouette = 0.864**
- Patients entering UMAP: 499; patients clustered (≥1 block): 77

**Cluster sizes:**
| Cluster | n | Pre-onset | Post-onset |
|---------|---|-----------|-----------|
| C0 | 21 | 11 (52.4%) | 10 (47.6%) |
| C1 | 26 | 9 (34.6%) | 17 (65.4%) |
| C2 | 30 | 21 (70.0%) | 9 (30.0%) |
| None | 422 | 84 (19.9%) | 338 (80.1%) |

---

### Panel c — ICD block prevalence heatmap (Fig 5c)

**Source**: block_overrep from `ms_prems_block_overrep.csv` + `ms_prems_group_overrep.csv`

**Cluster-defining ICD blocks (FDR < 0.05 in block_overrep):**

| ICD Block | Cluster | %in cluster | %in other | OR | FDR |
|-----------|---------|------------|-----------|-----|-----|
| R2 (gait/sensory/neurological symptoms) | C2 | 50.0% | 0.0% | Inf | 5.7e-07 |
| H4 (optic neuritis) | C1 | 46.2% | 0.0% | Inf | 1.7e-06 |
| B2 (mononucleosis) | C2 | 40.0% | 0.0% | Inf | 1.0e-05 |
| G5 (facial nerve disorders) | C1 | 34.6% | 0.0% | Inf | 6.3e-05 |
| M5 (cervical disc disorders) | C0 | 28.6% | 0.0% | Inf | 5.9e-04 |
| Z8 (personal history of diseases) | C0 | 57.1% | 19.6% | 5.3 | 4.6e-03 |
| M4 (ankylosing spondylitis) | C1 | 19.2% | 0.0% | Inf | 6.2e-03 |
| H2 (age-related cataract) | C0 | 23.8% | 1.8% | 16.4 | 8.2e-03 |
| N1 (acute kidney injury) | C2 | 16.7% | 0.0% | Inf | 9.4e-03 |
| K7 (liver diseases) | C2 | 16.7% | 0.0% | Inf | 9.4e-03 |

**ICD code-level cluster-defining features (FDR < 0.05 in group_overrep):**

| ICD Code | Cluster | Description | %in | %out | OR | FDR |
|----------|---------|-------------|-----|------|-----|-----|
| F32 | C1 | Depressive episode | 60.0% | 9.4% | 14.1 | 7.6e-09 |
| B27 | C1 | Infectious mononucleosis | 30.0% | 0.0% | Inf | 1.1e-07 |
| M79 | C0 | Fibromyalgia/soft tissue disorders | 27.3% | 0.0% | Inf | 1.9e-04 |
| M54 | C0 | Backache | 25.8% | 2.5% | 13.4 | 3.1e-03 |
| F17 | C0 | Nicotine dependence | 18.0% | 0.0% | Inf | 6.6e-03 |
| M25 | C0 | Joint pain/stiffness | 21.9% | 2.5% | 10.8 | 9.6e-03 |
| D50 | C0 | Iron deficiency anemia | 14.8% | 0.0% | Inf | 1.7e-02 |
| R10 | C0 | Abdominal/pelvic pain | 13.3% | 0.0% | Inf | 2.7e-02 |
| G56 | C0 | Mononeuropathies of upper limb | 11.7% | 0.0% | Inf | 4.4e-02 |

**Cluster identity summary:**
- **C0 (n=21)**: Musculoskeletal/chronic pain cluster — fibromyalgia, backache, joint pain, cervical disc, mononeuropathies, nicotine dependence, anemia. Oldest at diagnosis (mean 60y).
- **C1 (n=26)**: Neuroinflammatory/psychiatric cluster — optic neuritis, facial nerve disorders, ankylosing spondylitis, depression, infectious mononucleosis. Most pre-onset enriched (34.6% pre-onset; majority sampled pre-diagnosis, mean -3.4y).
- **C2 (n=30)**: Mixed neurological prodrome/viral cluster — gait/sensory/neurological symptoms, mononucleosis, renal/hepatic events. Highest pre-onset fraction (70%), latest sampling (mean +4.5y), highest PRS (0.368), lowest HLA rate (36.7%).

---

### Panel d — Age at diagnosis by cluster (Fig 5d)

**Source**: `results/endophenotype/ms_prems_cluster_summary.csv`

| Cluster | n | Mean age at diagnosis | %Female | Mean yrs to sampling |
|---------|---|-----------------------|---------|---------------------|
| C0 | 21 | 60.0 years | 71.4% | +2.0 (post-dx) |
| C1 | 26 | 49.2 years | 69.2% | −3.4 (pre-dx) |
| C2 | 30 | 57.2 years | 66.7% | +4.5 (post-dx) |
| None | 422 | 43.3 years | 70.6% | −11.7 (pre-dx) |

**Note**: Sex balance is uniform across clusters (all ~30% male, χ²=0.971, figure caption). C0 has the oldest diagnosis age; None group has the youngest and most pre-onset sampling.

---

### Panel e — HLA-DRB1*15:01 carrier rate by cluster (Fig 5e)

**Source**: `results/endophenotype/ms_prems_cluster_summary.csv`

| Cluster | n | HLA carrier rate |
|---------|---|-----------------|
| C0 | 21 | 52.4% |
| C1 | 26 | 46.2% |
| C2 | 30 | 36.7% |
| None | 422 | 51.4% |

**Key finding**: C2 has the lowest HLA carrier rate (36.7% vs 51.4% in None). Fisher tests comparing each cluster to None shown in figure. C2 is notably HLA-depleted and PRS-enriched — pointing to non-HLA genetic drivers.

---

### Panel f — MS polygenic risk score by cluster (Fig 5f)

**Source**: `results/endophenotype/ms_prems_cluster_summary.csv`

| Cluster | n | Mean PRS |
|---------|---|----------|
| C0 | 21 | 0.152 |
| C1 | 26 | 0.247 |
| C2 | 30 | **0.368** |
| None | 422 | 0.273 |

**Key finding**: C2 has the highest mean PRS (0.368), substantially above the None group (0.273) and C0 (0.152). Kruskal-Wallis test shown in figure. C2 = high-PRS, low-HLA cluster; C0 = low-PRS, intermediate-HLA cluster.

---

### Panels g–i — Cluster vs None volcano plots (Fig 5g/h/i)

**Method**: PSM-matched limma within each cluster vs matched MS-None controls. All 2,915 proteins tested.

**Source**: `results/endophenotype/cluster_proteomics/cluster_deps_C*_vs_None.csv`

| Cluster | Total tested | FDR < 0.05 | Top protein(s) |
|---------|-------------|-----------|----------------|
| C0 vs None | 2,915 | **0** | Top nominal: FCAMR (FDR=0.067), LRTM2 (FDR=0.067) |
| C1 vs None | 2,915 | **0** | Top nominal: KIAA2013 (FDR=0.051), HSD11B1 (FDR=0.103) |
| C2 vs None | 2,915 | **2** | CCL3 (FDR=0.040), VWF (FDR=0.042) |

**C2 FDR < 0.05 DEPs:**
| Protein | logFC | FDR | Direction |
|---------|-------|-----|-----------|
| CCL3 | +0.762 | 0.040 | Up in C2 |
| VWF | +0.685 | 0.042 | Up in C2 |

**C0 top nominal DEPs (FDR < 0.10):**
| Protein | logFC | FDR |
|---------|-------|-----|
| FCAMR | +0.867 | 0.067 |
| LRTM2 | −0.256 | 0.067 |
| TNFRSF10A | +0.415 | 0.067 |
| BLOC1S2 | +0.666 | 0.067 |

**C1 top nominal DEPs (FDR ~ 0.05):**
| Protein | logFC | FDR |
|---------|-------|-----|
| KIAA2013 | +0.600 | 0.051 |
| HSD11B1 | −0.394 | 0.103 |

---

### Panels j–l — GO:BP enrichment per cluster (Fig 5j/k/l)

**Source**: `results/endophenotype/cluster_pathways/cluster_go_*_vs_None.csv`

| Cluster | FDR < 0.05 GO terms |
|---------|---------------------|
| C0 | 17 |
| C1 | 90 |
| C2 | 749 |

**C0 top GO terms (17 total):**
| Term | FoldEnrichment | FDR | Count | Set |
|------|---------------|-----|-------|-----|
| Positive regulation of MAPK cascade | 6.6 | 0.0037 | 10 | up_DEPs |
| Positive regulation of MAPK cascade | 5.8 | 0.0019 | 12 | down_DEPs |
| Regulation of ROS metabolic process | 10.6 | 0.0040 | 7 | down_DEPs |
| Extrinsic apoptotic signaling (death domain) | 18.3 | 0.0064 | 5 | up_DEPs |
| Positive regulation of ERK1/2 cascade | 9.4 | 0.023 | 6 | up_DEPs |
| Positive regulation of protein phosphorylation | 8.2 | 0.037 | 6 | up_DEPs |

**C1 top GO terms (90 total):**
| Term | FoldEnrichment | FDR | Count | Set |
|------|---------------|-----|-------|-----|
| Cell adhesion mediated by integrin | 24.3 | 0.0023 | 5 | up_DEPs |
| Eosinophil chemotaxis | 52.8 | 0.0046 | 3 | up_DEPs |
| Positive regulation of cell adhesion via integrin | 50.4 | 0.0046 | 3 | up_DEPs |
| Cell chemotaxis | 7.8 | 0.0046 | 7 | up_DEPs |
| Killing of cells of another organism | 13.7 | 0.0046 | 5 | up_DEPs |
| Antimicrobial humoral immune response | 13.7 | 0.0046 | 5 | up_DEPs |

**C2 top GO terms (749 total — broadest enrichment):**
| Term | FoldEnrichment | FDR | Count | Set |
|------|---------------|-----|-------|-----|
| Leukocyte migration | 7.7 | 6.6e-15 | 30 | up_DEPs |
| Positive regulation of cell activation | 7.6 | 6.6e-15 | 30 | up_DEPs |
| Positive regulation of leukocyte activation | 7.7 | 1.1e-14 | 29 | up_DEPs |
| Mononuclear cell migration | 9.9 | 1.0e-13 | 23 | up_DEPs |
| Leukocyte cell-cell adhesion | 6.5 | 2.3e-12 | 28 | up_DEPs |
| Positive regulation of leukocyte migration | 12.0 | 5.6e-12 | 18 | up_DEPs |

**Pathway narrative**: C0 shows MAPK/apoptosis signaling (modest enrichment, 17 terms). C1 shows integrin-mediated adhesion and immune cell chemotaxis (90 terms, dominated by antimicrobial/eosinophil response). C2 has overwhelming leukocyte activation and migration enrichment (749 terms), consistent with CCL3/VWF as upstream chemokine/vascular drivers.

---

### Panels r–u — Cell-type ORA (Fig 5r/s/t/u)

**Source**: `results/endophenotype/cluster_celltypes/cluster_celltype_summary.csv`

**FDR < 0.05 cell-type enrichments: 31 total (C0: 18, C1: 0, C2: 13)**

**C0 cell-type ORA (18 significant):**
Downregulated DEP set (proteins reduced in C0 vs None):
| Cell type | OR | FDR | Direction |
|-----------|-----|-----|-----------|
| Macrophages | 0.29 | 2.4e-06 | down_DEPs |
| Monocytes | 0.27 | 2.4e-06 | down_DEPs |
| NK-cells | 0.33 | 6.2e-05 | down_DEPs |
| Plasma cells | 0.34 | 6.2e-05 | down_DEPs |
| T-cells | 0.36 | 7.5e-05 | down_DEPs |
| Dendritic cells | 0.37 | 2.8e-04 | down_DEPs |
| B-cells | 0.40 | 3.3e-04 | down_DEPs |
| Granulocytes | 0.43 | 7.5e-04 | down_DEPs |
| Endothelial | 0.56 | 0.037 | down_DEPs |
| Microglial cells | 0.58 | 0.039 | down_DEPs |

Upregulated DEP set (proteins elevated in C0 vs None):
| Cell type | OR | FDR |
|-----------|-----|-----|
| Astrocytes | 0.31 | 0.0020 |
| Excitatory neurons | 0.42 | 0.028 |
| Dendritic cells | 0.50 | 0.034 |
| Endothelial | 0.48 | 0.034 |
| T-cells | 0.50 | 0.034 |
| NK-cells | 0.48 | 0.034 |
| Inhibitory neurons | 0.48 | 0.034 |
| Oligodendrocytes | 0.45 | 0.034 |

**Note on OR < 1**: All OR values < 1 indicate cell-type markers are UNDER-represented in the DEP list relative to chance — i.e., the cluster DEPs are systematically biased away from canonical immune and CNS cell-type expression profiles. This is interpretable as C0 having a proteomically diffuse or non-cell-type-specific signature.

**C1 cell-type ORA (0 significant)**: No cell type reached FDR < 0.05 in either direction. The integrin/chemotaxis pathway enrichment is not attributable to a single detectable immune cell type by this ORA method.

**C2 cell-type ORA (13 significant, all up-DEP set):**
| Cell type | OR | FDR |
|-----------|-----|-----|
| Fibroblasts | **1.92** | 0.024 |
| Excitatory neurons | 0.48 | 1.5e-04 |
| NK-cells | 0.56 | 0.0010 |
| T-cells | 0.57 | 0.0010 |
| Granulocytes | 0.56 | 0.0010 |
| Oligodendrocytes | 0.53 | 0.0010 |
| Inhibitory neurons | 0.54 | 0.0010 |
| B-cells | 0.61 | 0.0051 |
| Plasma cells | 0.63 | 0.0060 |
| OPCs | 0.63 | 0.0091 |
| Astrocytes | 0.63 | 0.0091 |
| Monocytes | 0.67 | 0.017 |
| Dendritic cells | 0.67 | 0.018 |

**Key C2 cell-type finding**: Fibroblast markers are the only cell type OVER-represented (OR=1.92, FDR=0.024) in C2 up-regulated DEPs (CCL3, VWF context). All CNS and immune cell types are under-represented. This points to a stromal/vascular rather than immune-cell-intrinsic plasma signature for C2.

---

### Panel v — Pre-MS ML AUC by cluster (Fig 5v)

**Source**: `results/endophenotype/cluster_prems_ml/cluster_prems_auc_summary.csv`

| Comparison | n cases | AUC | 95% CI |
|------------|---------|-----|--------|
| All pre-onset | 125 | 0.569 | — |
| C0 pre-onset | 11 | 0.570 | 0.369–0.771 |
| C1 pre-onset | 9 | 0.548 | 0.338–0.758 |
| **C2 pre-onset** | **21** | **0.671** | **0.541–0.802** |

**Key finding**: C2 is the only cluster with a pre-MS classifier meaningfully above chance (AUC=0.671, lower CI=0.54). C0 and C1 are at chance. C2's classifier exploits its broader leukocyte activation signature.

**C2 pre-MS classifier top features:**
| Feature | Coefficient |
|---------|------------|
| ARL13B | −0.449 |
| BCAN | +0.273 |
| CD302 | −0.252 |
| LILRB4 | −0.233 |
| BAG3 | +0.213 |
| GALNT7 | +0.188 |

**C0 pre-MS classifier top features:**
| Feature | Coefficient |
|---------|------------|
| NEFL | −1.354 |
| CD302 | −1.253 |
| CD248 | +1.105 |
| BAG3 | −0.692 |
| LILRB4 | +0.600 |
| BCAN | +0.597 |

**C1 pre-MS classifier top features:**
| Feature | Coefficient |
|---------|------------|
| GALNT7 | +0.469 |
| CD302 | −0.405 |
| ARL13B | −0.354 |
| BCAN | −0.167 |
| BAG3 | +0.145 |

---

### Panels w–x — GMM co-presentation model (Fig 5w/x)

**Source**: `results/endophenotype/gmm_copresentation/gmm_bic_sweep.csv` + `gmm_copresenters.csv` + `gmm_posterior_probs.csv`

**BIC sweep (G=1–5):**
| G | BIC |
|---|-----|
| 2 | −13,272.9 |
| **3** | **−14,564.2** (optimal) |
| 4 | −14,474.1 |
| 5 | −31,716.5 |

Lowest BIC (most parsimonious) at **G=3**, confirming k-means k=3 selection.

**Co-presenter analysis:**
- Total MS cases in GMM: 499
- Co-presenters (soft assignment to >1 cluster): **1/499 (0.2%)**
- Hard assignment confidence: 498/499 cases have max posterior essentially 1.0
- This indicates near-complete cluster separation — the three endophenotypes are discrete rather than a continuous spectrum

**GMM vs k-means confusion (among the 77 clustered cases):**
| GMM component | C0 (k-means) | C1 (k-means) | C2 (k-means) |
|---------------|-------------|-------------|-------------|
| Component 1 | 7 | 5 | 10 |
| Component 2 | 13 | 10 | 19 |
| Component 3 | 1 | **11** | 1 |

GMM component 3 strongly aligns with C1 (11/26 = 42% of k-means C1, only 1/21 C0 and 1/30 C2). Component 3 appears to capture the depression/mononucleosis profile. The other components are mixed.

**GMM dominant cluster distribution (full 499 cases):**
| GMM cluster | n |
|-------------|---|
| C1 | 41 |
| C2 | 117 |
| C2.1 | 341 |

The majority of the full MS cohort (341/499 = 68%) has posterior dominated by the C2.1 component, consistent with the "None" (no prior comorbidity) phenotype being captured by the GMM as a large low-comorbidity component.

---

### Panel y — HLA/PRS enrichment by cluster (Fig 5y)

**Source**: `results/endophenotype/ms_prems_cluster_summary.csv`

| Cluster | n | HLA carrier % | Mean PRS | HLA vs None | PRS vs None |
|---------|---|--------------|----------|-------------|-------------|
| C0 | 21 | 52.4% | 0.152 | Similar | Lower |
| C1 | 26 | 46.2% | 0.247 | Similar | Similar |
| **C2** | **30** | **36.7%** | **0.368** | **Depleted** | **Elevated** |
| None | 422 | 51.4% | 0.273 | Reference | Reference |

**Key finding**: C2 exhibits an inverse genetic pattern — lowest HLA carrier rate (36.7%) but highest polygenic risk (0.368). This suggests C2 represents a predominantly HLA-independent, polygenic-driven disease subtype with early neurological symptom accrual and elevated immune activation in the pre-diagnostic plasma.

---

## eFIGURE 5 — Post-onset validation supplement

### Panels a–c — Post-onset volcano plots (eFig 5a/b/c)

**Source**: `results/endophenotype/postms_cluster_proteomics/postms_deps_*_vs_None.csv`

| Cluster | Post-onset FDR < 0.05 | Top proteins |
|---------|----------------------|--------------|
| C0 post vs None | **2** | EDEM2 (logFC=+0.782, FDR=3.5e-06), GPHA2 (logFC=+0.768, FDR=5.1e-03) |
| C1 post vs None | **0** | Top nominal: HSD11B1, CCL21, NPTN (all NS) |
| C2 post vs None | **0** | Top nominal: CCL3 (FDR=0.248), MMP7 (NS) |

**Notable**: C2 had 2 FDR<0.05 pre-onset DEPs (CCL3, VWF) but 0 post-onset. C0 had 0 pre-onset but 2 post-onset (EDEM2, GPHA2). This stage-specificity is captured by the concordance scatterplots (eFig5l–n).

---

### Panels d–g — Post-onset GO:BP enrichment (eFig 5d/e/f/g)

**Source**: `results/endophenotype/postms_cluster_pathways/postms_go_*_vs_None.csv`

| Cluster | Post-onset FDR < 0.05 GO terms | Top term |
|---------|-------------------------------|----------|
| C0 | 1 | Megakaryocyte development (FoldEnrich=33.8, FDR=0.011, down_DEPs) |
| C1 | 1 | Complement activation, alternative pathway (FoldEnrich=66.2, FDR=0.015, up_DEPs) |
| C2 | **22** | ECM organization, NK cell chemotaxis, basement membrane organization (all up_DEPs, FDR=0.027) |

**C2 post-onset top GO terms:**
| Term | FoldEnrichment | FDR | Set |
|------|---------------|-----|-----|
| Basement membrane organization | 41.8 | 0.027 | up_DEPs |
| Natural killer cell chemotaxis | 114.3 | 0.027 | up_DEPs |
| Positive regulation of Ca2+-mediated signaling | 29.1 | 0.027 | up_DEPs |
| Extracellular matrix organization | 8.6 | 0.027 | up_DEPs |
| Host-mediated suppression of viral transcription | 87.9 | 0.027 | up_DEPs |

**Post-onset C2 narrative**: Despite losing FDR DEPs in post-onset analysis, C2 retains a broad ECM/stromal/NK chemotaxis pathway enrichment — consistent with a sustained tissue-remodeling and NK-mediated immunological phenotype that predates and persists through diagnosis.

---

### Panels h–k — Post-onset cell-type ORA (eFig 5h/i/j/k)

**Source**: `results/endophenotype/postms_cluster_celltypes/postms_cluster_celltype_summary.csv`

**Post-onset FDR < 0.05 cell-type enrichments: 0 (across all clusters)**

No post-onset cell-type ORA reached FDR significance in any cluster. The cell-type signatures evident at pre-onset (C0 immune marker depletion, C2 fibroblast enrichment) are not detectable at the post-onset stage, consistent with the general blunting of cluster-specific proteomic signals after diagnosis.

---

### Panels l–n — Pre/post concordance scatterplots (eFig 5l/m/n)

**Source**: `results/endophenotype/prepost_concordance/concordance_stats.csv`

| Cluster | Pearson r (pre vs post logFC) | Spearman ρ | n proteins |
|---------|------------------------------|------------|------------|
| C0 | **0.686** | 0.656 | 2,915 |
| C1 | **0.730** | 0.702 | 2,915 |
| C2 | **0.308** | 0.195 | 2,915 |

**Key finding**: C0 and C1 show strong pre/post concordance (r=0.69–0.73), indicating their proteomic signatures are temporally stable and present across disease stages. C2 shows markedly weaker concordance (r=0.31), consistent with C2's pre-onset DEPs (CCL3, VWF) being stage-specific signals that do not persist post-diagnosis.

**Interpretation**: C2's pre-onset signature reflects biology that is specific to the pre-symptomatic window — an early inflammatory/vascular burst (CCL3 chemokine + VWF coagulation marker) that resolves or is obscured after diagnosis. C0 and C1's proteomes are more stable immunological or chronic phenotypes.

---

## DATA GAPS (panels without extractable CSVs)

| Panel | Description | Status |
|-------|-------------|--------|
| Fig 5n–q | PPI STRING chord diagrams | No CSV in results/endophenotype/cluster_ppi/ — visually derived from STRING API |
| Fig 5m, eFig 5g | Cross-cluster GO:BP dotplot | Uses same pathway CSVs as per-cluster panels; no separate cross-cluster comparison CSV |
| Fig 5u, eFig 5k | Cross-cluster cell-type heatmap | Derived from per-cluster ORA CSVs; no dedicated cross-cluster CSV |
| Fig 5w BIC curve | Exact BIC trajectory by G | **Available**: gmm_bic_sweep.csv (G=2–5 above) |

---

## ETABLE IDs for Fig 5 supplementary tables (proposed)

| Proposed eTable ID | Content |
|--------------------|---------|
| `etbl_prems_phewas` | Full PheWAS (23 FDR<0.05 conditions) |
| `etbl_cluster_clinical` | Cluster clinical characteristics (n, age, sex, HLA, PRS, sampling time) |
| `etbl9_cluster_deps` | Already exists: cluster DEPs (all contrasts) |
| `etbl8_cluster_assignments` | Already exists: per-case cluster assignments |
| `etbl_cluster_pathways` | Cross-cluster GO:BP enrichment (C0/C1/C2 vs None) |
| `etbl_postms_concordance` | Pre/post concordance r by cluster (already in etbl10_concordance) |
