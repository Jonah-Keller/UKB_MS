"""RAP Extraction — ICD-10 E/F First-Occurrence Chapters

Upload and run in a RAP JupyterLab session (Spark cluster, ≥8 GB RAM).

Extracts UKB first-occurrence date fields for ICD-10 chapters:
  E — Endocrine, nutritional and metabolic diseases (70 codes)
  F — Mental and behavioural disorders (77 codes)

These chapters were absent from the original ukb_first_occurrence_FINAL
export which only covered G–Q (p131xxx field IDs).

Output saved to:
  /opt/notebooks/data/ms_als_extraction/icd_ef_first_occurrence.csv

After running:
  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_ef_first_occurrence.csv"
Then run locally:
  python analysis/00_rap_extraction/06_parse_ef_chapters.py
"""

import os
import re
import math
import multiprocessing
import psutil

import dxpy
import dxdata
import pyspark
from pyspark.sql import SparkSession

num_cores = multiprocessing.cpu_count()
available_memory_gb = math.floor(psutil.virtual_memory().available / (1024 ** 3))
spark_memory = max(available_memory_gb - 2, 4)

spark = (
    SparkSession.builder
    .appName("UKB_EF_ICD_Extraction")
    .config("spark.driver.memory", f"{spark_memory}g")
    .config("spark.executor.memory", f"{spark_memory}g")
    .config("spark.driver.maxResultSize", "4g")
    .getOrCreate()
)
engine = dxdata.connect()
print(f"Spark ready: {num_cores} cores, {spark_memory} GB RAM")

# ── Connect to dispensed dataset ──────────────────────────────────────────────
dispensed_dataset = dxpy.find_one_data_object(
    typename="Dataset",
    name="app*.dataset",
    folder="/",
    name_mode="glob",
)
dataset = dxdata.load_dataset(id=dispensed_dataset["id"])
participant = dataset["participant"]

# ── E/F field IDs → ICD-10 code mapping ──────────────────────────────────────
# field_id (int) : (icd10_3char, description)
EF_FIELD_MAP = {
    130692: ("E01", "iodine-deficiency-related thyroid disorders"),
    130694: ("E02", "subclinical iodine-deficiency hypothyroidism"),
    130696: ("E03", "other hypothyroidism"),
    130698: ("E04", "other non-toxic goitre"),
    130700: ("E05", "thyrotoxicosis [hyperthyroidism]"),
    130702: ("E06", "thyroiditis"),
    130704: ("E07", "other disorders of thyroid"),
    130706: ("E10", "insulin-dependent diabetes mellitus"),
    130708: ("E11", "non-insulin-dependent diabetes mellitus"),
    130710: ("E12", "malnutrition-related diabetes mellitus"),
    130712: ("E13", "other specified diabetes mellitus"),
    130714: ("E14", "unspecified diabetes mellitus"),
    130716: ("E15", "nondiabetic hypoglycaemic coma"),
    130718: ("E16", "other disorders of pancreatic internal secretion"),
    130720: ("E20", "hypoparathyroidism"),
    130722: ("E21", "hyperparathyroidism and other disorders of parathyroid gland"),
    130724: ("E22", "hyperfunction of pituitary gland"),
    130726: ("E23", "hypofunction and other disorders of pituitary gland"),
    130728: ("E24", "cushings syndrome"),
    130730: ("E25", "adrenogenital disorders"),
    130732: ("E26", "hyperaldosteronism"),
    130734: ("E27", "other disorders of adrenal gland"),
    130736: ("E28", "ovarian dysfunction"),
    130738: ("E29", "testicular dysfunction"),
    130740: ("E30", "disorders of puberty not elsewhere classified"),
    130742: ("E31", "polyglandular dysfunction"),
    130744: ("E32", "diseases of thymus"),
    130746: ("E34", "other endocrine disorders"),
    130748: ("E35", "disorders of endocrine glands in diseases classified elsewhere"),
    130750: ("E40", "kwashiorkor"),
    130752: ("E41", "nutritional marasmus"),
    130756: ("E43", "unspecified severe protein-energy malnutrition"),
    130758: ("E44", "protein-energy malnutrition of moderate and mild degree"),
    130760: ("E45", "retarded development following protein-energy malnutrition"),
    130762: ("E46", "unspecified protein-energy malnutrition"),
    130764: ("E50", "vitamin a deficiency"),
    130766: ("E51", "thiamine deficiency"),
    130768: ("E52", "niacin deficiency [pellagra]"),
    130770: ("E53", "deficiency of other b group vitamins"),
    130772: ("E54", "ascorbic acid deficiency"),
    130774: ("E55", "vitamin d deficiency"),
    130776: ("E56", "other vitamin deficiencies"),
    130778: ("E58", "dietary calcium deficiency"),
    130780: ("E59", "dietary selenium deficiency"),
    130782: ("E60", "dietary zinc deficiency"),
    130784: ("E61", "deficiency of other nutrient elements"),
    130786: ("E63", "other nutritional deficiencies"),
    130788: ("E64", "sequelae of malnutrition and other nutritional deficiencies"),
    130790: ("E65", "localised adiposity"),
    130792: ("E66", "obesity"),
    130794: ("E67", "other hyperalimentation"),
    130796: ("E68", "sequelae of hyperalimentation"),
    130798: ("E70", "disorders of aromatic amino-acid metabolism"),
    130800: ("E71", "disorders of branched-chain amino-acid metabolism"),
    130802: ("E72", "other disorders of amino-acid metabolism"),
    130804: ("E73", "lactose intolerance"),
    130806: ("E74", "other disorders of carbohydrate metabolism"),
    130808: ("E75", "disorders of sphingolipid metabolism and other lipid storage"),
    130810: ("E76", "disorders of glycosaminoglycan metabolism"),
    130812: ("E77", "disorders of glycoprotein metabolism"),
    130814: ("E78", "disorders of lipoprotein metabolism and other lipidaemias"),
    130816: ("E79", "disorders of purine and pyrimidine metabolism"),
    130818: ("E80", "disorders of porphyrin and bilirubin metabolism"),
    130820: ("E83", "disorders of mineral metabolism"),
    130822: ("E84", "cystic fibrosis"),
    130824: ("E85", "amyloidosis"),
    130826: ("E86", "volume depletion"),
    130828: ("E87", "other disorders of fluid electrolyte and acid-base balance"),
    130830: ("E88", "other metabolic disorders"),
    130832: ("E89", "postprocedural endocrine and metabolic disorders"),
    130836: ("F00", "dementia in alzheimers disease"),
    130838: ("F01", "vascular dementia"),
    130840: ("F02", "dementia in other diseases classified elsewhere"),
    130842: ("F03", "unspecified dementia"),
    130844: ("F04", "organic amnesic syndrome not induced by alcohol"),
    130846: ("F05", "delirium not induced by alcohol and other psychoactive substances"),
    130848: ("F06", "other mental disorders due to brain damage and dysfunction"),
    130850: ("F07", "personality and behavioural disorders due to brain disease"),
    130852: ("F09", "unspecified organic or symptomatic mental disorder"),
    130854: ("F10", "mental and behavioural disorders due to use of alcohol"),
    130856: ("F11", "mental and behavioural disorders due to use of opioids"),
    130858: ("F12", "mental and behavioural disorders due to use of cannabinoids"),
    130860: ("F13", "mental and behavioural disorders due to use of sedatives or hypnotics"),
    130862: ("F14", "mental and behavioural disorders due to use of cocaine"),
    130864: ("F15", "mental and behavioural disorders due to use of other stimulants"),
    130866: ("F16", "mental and behavioural disorders due to use of hallucinogens"),
    130868: ("F17", "mental and behavioural disorders due to use of tobacco"),
    130870: ("F18", "mental and behavioural disorders due to use of volatile solvents"),
    130872: ("F19", "mental and behavioural disorders due to multiple drug use"),
    130874: ("F20", "schizophrenia"),
    130876: ("F21", "schizotypal disorder"),
    130878: ("F22", "persistent delusional disorders"),
    130880: ("F23", "acute and transient psychotic disorders"),
    130882: ("F24", "induced delusional disorder"),
    130884: ("F25", "schizoaffective disorders"),
    130886: ("F28", "other nonorganic psychotic disorders"),
    130888: ("F29", "unspecified nonorganic psychosis"),
    130890: ("F30", "manic episode"),
    130892: ("F31", "bipolar affective disorder"),
    130894: ("F32", "depressive episode"),
    130896: ("F33", "recurrent depressive disorder"),
    130898: ("F34", "persistent mood affective disorders"),
    130900: ("F38", "other mood affective disorders"),
    130902: ("F39", "unspecified mood affective disorder"),
    130904: ("F40", "phobic anxiety disorders"),
    130906: ("F41", "other anxiety disorders"),
    130908: ("F42", "obsessive-compulsive disorder"),
    130910: ("F43", "reaction to severe stress and adjustment disorders"),
    130912: ("F44", "dissociative conversion disorders"),
    130914: ("F45", "somatoform disorders"),
    130916: ("F48", "other neurotic disorders"),
    130918: ("F50", "eating disorders"),
    130920: ("F51", "nonorganic sleep disorders"),
    130922: ("F52", "sexual dysfunction not caused by organic disorder"),
    130924: ("F53", "mental and behavioural disorders associated with the puerperium"),
    130926: ("F54", "psychological and behavioural factors associated with disorders"),
    130928: ("F55", "abuse of non-dependence-producing substances"),
    130930: ("F59", "unspecified behavioural syndromes associated with physiological"),
    130932: ("F60", "specific personality disorders"),
    130934: ("F61", "mixed and other personality disorders"),
    130936: ("F62", "enduring personality changes not attributable to brain damage"),
    130938: ("F63", "habit and impulse disorders"),
    130940: ("F64", "gender identity disorders"),
    130942: ("F65", "disorders of sexual preference"),
    130944: ("F66", "psychological and behavioural disorders associated with sexuality"),
    130946: ("F68", "other disorders of adult personality and behaviour"),
    130948: ("F69", "unspecified disorder of adult personality and behaviour"),
    130950: ("F70", "mild mental retardation"),
    130952: ("F71", "moderate mental retardation"),
    130954: ("F72", "severe mental retardation"),
    130958: ("F78", "other mental retardation"),
    130960: ("F79", "unspecified mental retardation"),
    130962: ("F80", "specific developmental disorders of speech and language"),
    130964: ("F81", "specific developmental disorders of scholastic skills"),
    130966: ("F82", "specific developmental disorder of motor function"),
    130968: ("F83", "mixed specific developmental disorders"),
    130970: ("F84", "pervasive developmental disorders"),
    130972: ("F88", "other disorders of psychological development"),
    130974: ("F89", "unspecified disorder of psychological development"),
    130976: ("F90", "hyperkinetic disorders"),
    130978: ("F91", "conduct disorders"),
    130980: ("F92", "mixed disorders of conduct and emotions"),
    130982: ("F93", "emotional disorders with onset specific to childhood"),
    130984: ("F94", "disorders of social functioning with onset specific to childhood"),
    130986: ("F95", "tic disorders"),
    130988: ("F98", "other behavioural and emotional disorders with onset usually in childhood"),
    130990: ("F99", "mental disorder not otherwise specified"),
}

# ── Build column name list ────────────────────────────────────────────────────
print("Resolving field names...")
col_names = ["eid"]
field_to_col = {}  # field_id → column name in dataset

for fid in EF_FIELD_MAP:
    col = f"p{fid}"
    try:
        fields = participant.find_fields(name_regex=rf"^{col}$")
        if fields:
            col_names.append(col)
            field_to_col[fid] = col
            print(f"  {col} -> {EF_FIELD_MAP[fid][0]} OK")
        else:
            print(f"  {col} -> NOT FOUND")
    except Exception as e:
        print(f"  {col} -> ERROR: {e}")

print(f"\nFound {len(field_to_col)} / {len(EF_FIELD_MAP)} E/F fields in dataset")

# ── Extract ───────────────────────────────────────────────────────────────────
print("Extracting from participant table...")
df = participant.retrieve_fields(names=col_names, engine=engine)
pdf = df.toPandas()
print(f"Extracted: {pdf.shape[0]} participants x {pdf.shape[1]} columns")

# Add year_birth for age computation
try:
    yob_df = participant.retrieve_fields(names=["eid", "p34"], engine=engine).toPandas()
    pdf = pdf.merge(yob_df.rename(columns={"p34": "year_birth"}), on="eid", how="left")
    print("Added year_birth (field 34)")
except Exception as e:
    print(f"year_birth not available: {e}")

# ── Save ──────────────────────────────────────────────────────────────────────
out_dir = "/opt/notebooks/data/ms_als_extraction/"
os.makedirs(out_dir, exist_ok=True)
out_path = f"{out_dir}icd_ef_first_occurrence.csv"
pdf.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
print(f"  Size: {os.path.getsize(out_path) / 1e6:.1f} MB")
print(f"  Rows: {len(pdf)}")
print(f"  Cols: {list(pdf.columns[:10])} ...")

# Save field mapping for local parser
import json
mapping_path = f"{out_dir}ef_field_mapping.json"
with open(mapping_path, "w") as fh:
    json.dump({str(k): v for k, v in EF_FIELD_MAP.items()}, fh, indent=2)
print(f"  Mapping: {mapping_path}")

print("\nDone! Download with:")
print(f'  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_ef_first_occurrence.csv"')
print(f'  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/ef_field_mapping.json"')
