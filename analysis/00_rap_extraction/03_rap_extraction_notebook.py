"""RAP JupyterLab Extraction Notebook — MS/ALS Project

Upload this script to RAP and run in a JupyterLab session.
Extracts ALL fields needed for the complete MS + ALS proteomics replication.

Categories:
  1. Genetic PCs (field 22009) — population stratification
  2. HLA imputation (field 22182) — MS immune cluster
  3. HES ICD-10 diagnoses (41270, 41202, 41204) — G12.2 ALS refinement
  4. Death cause ICD-10 (40001, 40002) — already in demographics but we want clean extract
  5. Genetic QC (22006, 22019, 22021, 22027) — QC filtering
  6. Assessment centre + dates (54, 53) — batch effects
  7. Ethnicity (21000) — covariate
  8. Townsend deprivation (189) — socioeconomic covariate
  9. MS PRS / genetic risk fields if available

Output: CSVs saved to /opt/notebooks/data/ms_als_extraction/ → then upload to project

To run:
  1. Start a JupyterLab session on RAP (Spark cluster, ≥8GB RAM)
  2. Upload this file as a notebook or run cells in sequence
  3. After completion, download outputs with:
     dx download "project-GxX43xjJpjp2G7XK6ffz0qb1:/data/ms_als_extraction/*"
"""

# ============================================================
# Cell 1: Imports and Spark Setup
# ============================================================
import pyspark
import dxpy
import dxdata
import numpy as np
import pandas as pd

import os
import re
import multiprocessing
import psutil
import math
from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession

# Get available system resources
num_cores = multiprocessing.cpu_count()
available_memory_bytes = psutil.virtual_memory().available
available_memory_gb = math.floor(available_memory_bytes / (1024**3))

# Reserve some memory for the OS
spark_memory = max(available_memory_gb - 2, 4)

spark = (
    SparkSession.builder
    .appName("MS_ALS_Extraction")
    .config("spark.driver.memory", f"{spark_memory}g")
    .config("spark.executor.memory", f"{spark_memory}g")
    .config("spark.driver.maxResultSize", "4g")
    .getOrCreate()
)

sc = spark.sparkContext
engine = dxdata.connect()
print(f"Spark initialized: {num_cores} cores, {spark_memory}GB RAM")

# ============================================================
# Cell 2: Connect to dispensed dataset
# ============================================================
dispensed_database = dxpy.find_one_data_object(
    classname='database',
    name='app*',
    folder='/',
    name_mode='glob',
    describe=True
)
dispensed_database_name = dispensed_database['describe']['name']

dispensed_dataset = dxpy.find_one_data_object(
    typename='Dataset',
    name='app*.dataset',
    folder='/',
    name_mode='glob'
)
dispensed_dataset_id = dispensed_dataset['id']

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
participant = dataset['participant']
print(f"Dataset: {dispensed_database_name}")
print(f"Dataset ID: {dispensed_dataset_id}")

# ============================================================
# Cell 3: Helper functions
# ============================================================
def fields_for_id(field_id):
    """Get all field objects for a given UKB field ID."""
    from distutils.version import LooseVersion
    field_id = str(field_id)
    fields = participant.find_fields(
        name_regex=r'^p{}(?!\d)(?:_i\d+)?(?:_a\d+)?$'.format(field_id)
    )
    return sorted(fields, key=lambda f: LooseVersion(f.name))

def field_names_for_id(field_id):
    """Get field column names for a given UKB field ID."""
    return [f.name for f in fields_for_id(field_id)]

def extract_fields(field_ids, output_name, max_arrays=None):
    """Extract fields from participant table and save as CSV.
    
    Parameters
    ----------
    field_ids : list of int
        UKB field IDs to extract
    output_name : str
        Output filename prefix
    max_arrays : int, optional
        Maximum number of array indices to include per field.
        If None, include all arrays.
    """
    all_field_names = ['eid']
    missing_fields = []
    
    for fid in field_ids:
        try:
            names = field_names_for_id(fid)
            if max_arrays is not None:
                # Filter to only keep limited array indices
                filtered = []
                for n in names:
                    # Extract array index if present
                    array_match = re.search(r'_a(\d+)', n)
                    if array_match:
                        aidx = int(array_match.group(1))
                        if aidx < max_arrays:
                            filtered.append(n)
                    else:
                        filtered.append(n)
                names = filtered
            all_field_names.extend(names)
            print(f"  Field {fid}: {len(names)} columns")
        except Exception as e:
            missing_fields.append(fid)
            print(f"  Field {fid}: NOT FOUND ({e})")
    
    all_field_names = list(dict.fromkeys(all_field_names))  # deduplicate
    print(f"Total columns: {len(all_field_names)}")
    
    if missing_fields:
        print(f"Missing fields: {missing_fields}")
    
    # Retrieve data
    print("Retrieving data from Spark...")
    df = participant.retrieve_fields(names=all_field_names, engine=engine)
    pdf = df.toPandas()
    print(f"Retrieved: {pdf.shape[0]} rows × {pdf.shape[1]} columns")
    
    # Save
    out_dir = '/opt/notebooks/data/ms_als_extraction/'
    os.makedirs(out_dir, exist_ok=True)
    out_path = f'{out_dir}{output_name}.csv'
    pdf.to_csv(out_path, index=False)
    print(f"Saved: {out_path} ({os.path.getsize(out_path) / 1e6:.1f} MB)")
    
    return pdf

# ============================================================
# Cell 4: Extract Genetic PCs (field 22009)
# Purpose: Population stratification covariates
# ============================================================
print("=" * 60)
print("EXTRACTION 1: Genetic Principal Components")
print("=" * 60)
genetic_pc_fields = [22009]  # 40 PCs as array indices a1-a40
pc_df = extract_fields(genetic_pc_fields, 'genetic_pcs')

# ============================================================
# Cell 5: Extract HLA Imputation (field 22182)
# Purpose: HLA-DRB1*15:01 dosage for MS immune cluster
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 2: HLA Imputation Values")
print("=" * 60)
hla_fields = [22182]  # HLA imputation values (362 alleles as arrays)
# Limit to first 100 arrays to avoid memory issues — DRB1*15:01 should be in first few
hla_df = extract_fields(hla_fields, 'hla_imputation', max_arrays=100)

# ============================================================
# Cell 6: Extract HES ICD-10 Diagnoses (41270, 41202, 41204)
# Purpose: Full 4-char ICD-10 codes (G12.2 for ALS)
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 3: HES Diagnoses (ICD-10)")
print("=" * 60)
# 41270 = all ICD-10 diagnoses combined (huge array field)
# 41202 = main ICD-10
# 41204 = secondary ICD-10
hes_fields = [41270, 41202, 41204]
# These are array fields with potentially hundreds of entries per person
# Limit to first 100 array indices
hes_df = extract_fields(hes_fields, 'hes_icd10_diagnoses', max_arrays=100)

# ============================================================
# Cell 7: Extract Death Causes (40001, 40002)
# Purpose: Cause of death ICD-10 codes for ALS confirmation
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 4: Death Cause ICD-10")
print("=" * 60)
death_fields = [
    40001,  # Primary cause of death (ICD-10)
    40002,  # Contributory causes of death (ICD-10)
    40000,  # Date of death
    40007,  # Age at death
]
death_df = extract_fields(death_fields, 'death_causes')

# ============================================================
# Cell 8: Extract Genetic QC Fields
# Purpose: Sample/kinship filtering
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 5: Genetic QC")
print("=" * 60)
genetic_qc_fields = [
    22006,  # Genetic ethnic grouping
    22019,  # Sex chromosome aneuploidy
    22021,  # Genetic kinship to other participants
    22027,  # Outliers for heterozygosity or missing rate
    22001,  # Genetic sex
    22000,  # Genotype measurement batch
    22020,  # Used in genetic principal components
]
qc_df = extract_fields(genetic_qc_fields, 'genetic_qc')

# ============================================================
# Cell 9: Extract Covariates & Demographics
# Purpose: BMI, ethnicity, Townsend, smoking, alcohol, centre
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 6: Covariates & Demographics")
print("=" * 60)
covariate_fields = [
    31,     # Sex
    34,     # Year of birth
    52,     # Month of birth
    21003,  # Age at recruitment
    21001,  # BMI
    21000,  # Ethnic background
    189,    # Townsend deprivation index
    20116,  # Smoking status
    20117,  # Alcohol intake frequency
    54,     # Assessment centre
    53,     # Date of attending assessment centre
    22189,  # Townsend supplementary
    20022,  # Birth year (estimated)
]
cov_df = extract_fields(covariate_fields, 'covariates_demographics')

# ============================================================
# Cell 10: Extract MS/ALS-specific exclusion ICD codes
# Purpose: Neuropathy/myopathy exclusion for Chia control def
# G60-G64 neuropathy, G70-G73 myopathy (Chia control exclusion)
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 7: Neurological Exclusion Codes")
print("=" * 60)
neuro_excl_fields = [
    131082,  # G60 Hereditary/idiopathic neuropathy
    131084,  # G61 Inflammatory polyneuropathy
    131086,  # G62 Other polyneuropathies
    131088,  # G63 Polyneuropathy in diseases classified elsewhere
    131090,  # G64 Other disorders of peripheral nervous system
    131092,  # G70 Myasthenia gravis
    131094,  # G71 Primary disorders of muscles
    131096,  # G72 Other myopathies
    131098,  # G73 Myoneural junction/muscle disorders
    131016,  # G12 Motor neuron disease (our ALS field)
    131042,  # G35 Multiple sclerosis (our MS field)
]
excl_df = extract_fields(neuro_excl_fields, 'neurological_exclusion_codes')

# ============================================================
# Cell 11: Extract Blood Sample Dates
# Purpose: Timing calculations for pre/post-onset classification
# ============================================================
print("\n" + "=" * 60)
print("EXTRACTION 8: Blood Sample & Assessment Dates")
print("=" * 60)
date_fields = [
    3166,   # Date/time of blood sample collection
    53,     # Date of attending assessment centre (already in covariates)
]
date_df = extract_fields(date_fields, 'blood_sample_dates')

# ============================================================
# Cell 12: Upload outputs to RAP project
# ============================================================
print("\n" + "=" * 60)
print("UPLOADING TO PROJECT")
print("=" * 60)

# Create output folder on RAP
import subprocess

out_dir = '/opt/notebooks/data/ms_als_extraction/'
rap_dest = '/data/ms_als_extraction/'

# Upload each CSV
for f in os.listdir(out_dir):
    if f.endswith('.csv'):
        local_path = os.path.join(out_dir, f)
        print(f"Uploading {f} ({os.path.getsize(local_path) / 1e6:.1f} MB)...")
        subprocess.run(['dx', 'upload', local_path, '--dest', f'{rap_dest}{f}', '--brief'], check=True)

print("\nDone! All files uploaded to:", rap_dest)
print("\nTo download locally, run:")
print(f'  dx download "project-GxX43xjJpjp2G7XK6ffz0qb1:{rap_dest}*" -o data/ukb/rap_extraction/')
