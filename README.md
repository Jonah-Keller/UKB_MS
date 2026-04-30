# UK Biobank Disease-Proteomics Replication Template

End-to-end pipeline for a UK Biobank plasma-proteomics study of a single
disease cohort: cohort engineering, Olink QC, differential expression,
network/cell-type enrichment, ML classifiers, HLA × disease interaction,
PRS interaction, MRI association, sex-stratified analyses,
pre-diagnostic comorbidity endophenotyping, cluster GWAS, and the
manuscript build (5 main figures, 9 extended-data figures, 33
supplementary tables).

The repository ships configured for **multiple sclerosis (MS)** as the
default cohort. To replicate the same study on a different disease, edit
`configs/disease.yaml` — no analysis script needs to be modified.

## Repository layout

```
configs/
├── disease.yaml          ← swap config: ICD codes, HLA allele, PRS IDs
└── figures.yaml          ← figure/table registry (captions use {{placeholders}})

scripts/
├── build/                ← manuscript build orchestrators
│   ├── build_submission.R       ← top-level: knit + render + assemble docx
│   ├── build_tables_docx.R      ← stage CSVs + render main tables
│   ├── build_supplementary.R    ← supplementary docx + xlsx
│   ├── compose_figure.R         ← composite multi-panel figure rendering
│   └── build_unified_data.py    ← consolidate RAP-extracted data
└── extraction/           ← UK Biobank RAP extraction utilities

analysis/
├── 00_rap_extraction/    ← extraction from UKB RAP
├── 01_cohort/            ← case/control engineering (reads disease.yaml)
├── 02_olink_qc/          ← Olink NPX QC (cohort-agnostic)
├── 03_covariates/        ← demographics, PCs
├── 04_differential/      ← limma DEPs (pre/post/combined)
├── 05_networks/          ← PPI, GO, cell-type enrichment
├── 06_ml/                ← elastic net, GBM, RF classifiers
├── 07_hla_subgroup/      ← HLA × disease interaction (reads hla_allele)
├── 08_mri/               ← protein × WMH/DTI/cognition
├── 09_sex_stratified/    ← per-sex limma + interaction
├── 10_prs/               ← PRS scoring + interaction (reads prs_pgs_ids)
├── 11_cross_disease/     ← cross-cohort comparison
├── 12_clustering/        ← unsupervised patient clustering
├── 13_survival/          ← Cox PH on protein vs time-to-Dx
├── 14_comorbidity/       ← pre-diagnostic PheWAS
├── 15_endophenotype/     ← comorbidity-driven cluster proteomics
├── 16_cluster_gwas/      ← REGENIE GWAS on cluster phenotype
├── manuscript/           ← figureN.R + manuscript.Rmd
└── helpers/
    ├── disease_config.{py,R}        ← loads configs/disease.yaml (cached)
    ├── ukb_disease_outcomes.py      ← case/control + timing classifiers
    └── ukb_ms_outcomes.py           ← back-compat shim → ukb_disease_outcomes

code/
├── config/00_config.R    ← project paths + disease.yaml metadata
├── registry/             ← figures.yaml loader + Jinja-style render
└── tables/               ← flextable styling, CSV staging, supp xlsx packer

data/                     ← UKB-restricted (gitignored). Required subdirs:
results/                  ← analysis outputs (most committed; participant-level CSVs gitignored)
manuscript/               ← reference.docx + submission/
```

## Replicating on a new disease cohort

1. **Clone** this repo as a template.
2. **Edit `configs/disease.yaml`** — at minimum:
   - `cohort_short` (e.g., `"als"`) — used as filename prefix
   - `disease_long` / `disease_short_caps` — caption substitution
   - `icd_codes` (e.g., `["G12"]`) — primary case codes
   - `first_occurrence_field` — UKB field for date-of-first-occurrence
   - `control_exclusion_codes` — codes that disqualify from controls
   - `hla_allele` + `hla_search_patterns` + `hla_carrier_col` — HLA risk locus
   - `prs_pgs_ids` — PGS Catalog IDs for the disease
3. **Drop UKB extracts into `data/ukb/`** (see existing `analysis/00_rap_extraction/` scripts).
4. **Run the pipeline** (each stage in order):
   ```bash
   python -m analysis.01_cohort.01_build_cohort
   Rscript analysis/02_olink_qc/01_olink_qc.R
   # … (see analysis/ stage numbers)
   ```
5. **Build the manuscript**:
   ```bash
   Rscript scripts/build/build_submission.R
   ```

The figure/table registry (`configs/figures.yaml`) substitutes `{{disease_long}}`,
`{{disease_short_caps}}`, `{{hla_allele}}`, `{{prs_label}}`, and
`{{cohort_short}}` placeholders at build time using values from
`configs/disease.yaml`.

## Citation

Pre-diagnostic comorbidity endophenotypes reveal proteomic and genetic
heterogeneity in multiple sclerosis. Keller JN, Kantor JA, Elahi FM. (2026).

## License

See `LICENSE`.
