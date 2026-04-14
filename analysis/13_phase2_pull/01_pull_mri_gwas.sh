#!/bin/bash
# Phase 2 UKB Data Extraction Stub
# To be executed on the UKB Research Analysis Platform (RAP)

echo "=== UKB Phase 2 Multi-Modal Extraction Stub ==="

# 1. GWAS SNPs for MS (HLA-DRB1*15:01)
# Field 22182: HLA imputation values
echo "1. Run dx extract_dataset to get Field 22182 (HLA) for MS/ALS participants"
# dx extract_dataset project-xxxx:record-xxxx -fields="participant.eid,participant.p22182" -o ukb_hla.csv --eids ms_als_eids.txt

# 2. Genetics PCs (to fix the empty genetics in covariates table)
# Field 22009: Genetic principal components
echo "2. Run dx extract_dataset for Field 22009 (Genetic PCs 1-40)"
# dx extract_dataset project-xxxx:record-xxxx -fields="participant.eid,participant.p22009_a*" -o ukb_genetic_pcs.csv --eids ms_als_eids.txt

# 3. Brain MRI Data (Phase 2 extension)
echo "3. Run dx extract_dataset for Brain MRI imaging derived phenotypes (IDPs)"
# Example: T1 structural IDPs (Category 110, 111, 112)
# dx extract_dataset project-xxxx:record-xxxx -fields="participant.eid,participant.p25000_i2,..." -o ukb_brain_mri_idps.csv --eids ms_als_eids.txt

echo "Note: Download outputs via dx download and place them in:"
echo " - data/ukb/genetics/"
echo " - data/ukb/brain_mri/"
