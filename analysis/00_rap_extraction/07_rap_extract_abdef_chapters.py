"""RAP Extraction — Missing ICD-10 Chapters A, B, D, E, F

Uses dx extract_dataset (no Java/Spark needed — works in dxjupyterlab cmd mode).

Chapters:
  A  Infectious and parasitic diseases    (81 fields: p130000-p130170)
  B  Viral / parasitic diseases           (79 fields: p130174-p130344)
  D  Neoplasms first-occurrence           (34 fields: p130622-p130688)
  E  Endocrine, nutritional, metabolic    (70 fields: p130692-p130832)
  F  Mental and behavioural disorders     (77 fields: p130836-p130990)

Output: project:/data/ms_als_extraction/icd_abdef_first_occurrence.csv
        project:/data/ms_als_extraction/abdef_field_mapping.json

After downloading:
  python analysis/00_rap_extraction/08_parse_abdef_chapters.py
"""

import os, json, subprocess
import dxpy

RAP_PROJECT = "project-J2P0fqjJ1Gg4gByJ6BJ51pYj"
RAP_FOLDER  = "/data/ms_als_extraction"

# Find dispensed dataset record ID
dispensed_dataset = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)
dataset_id = dispensed_dataset["id"]
print(f"Dataset: {dataset_id}")

# Field mapping: UKB field ID -> (ICD-10 3-char code, description)
FIELD_MAP = {
    130000: ("A00", "cholera"),
    130002: ("A01", "typhoid and paratyphoid fevers"),
    130004: ("A02", "other salmonella infections"),
    130006: ("A03", "shigellosis"),
    130008: ("A04", "other bacterial intestinal infections"),
    130010: ("A05", "other bacterial foodborne intoxications"),
    130012: ("A06", "amoebiasis"),
    130014: ("A07", "other protozoal intestinal diseases"),
    130016: ("A08", "viral and other specified intestinal infections"),
    130018: ("A09", "diarrhoea and gastro-enteritis of presumed infectious origin"),
    130020: ("A15", "respiratory tuberculosis, bacteriologically and histologically confirmed"),
    130022: ("A16", "respiratory tuberculosis, not confirmed bacteriologically or histologically"),
    130024: ("A17", "tuberculosis of nervous system"),
    130026: ("A18", "tuberculosis of other organs"),
    130028: ("A19", "miliary tuberculosis"),
    130030: ("A20", "plague"),
    130034: ("A22", "anthrax"),
    130036: ("A23", "brucellosis"),
    130038: ("A24", "glanders and melioidosis"),
    130040: ("A25", "rat-bite fevers"),
    130042: ("A26", "erysipeloid"),
    130044: ("A27", "leptospirosis"),
    130046: ("A28", "other zoonotic bacterial diseases, not elsewhere classified"),
    130048: ("A30", "leprosy [hansen's disease]"),
    130050: ("A31", "infection due to other mycobacteria"),
    130052: ("A32", "listeriosis"),
    130054: ("A33", "tetanus neonatorum"),
    130058: ("A35", "other tetanus"),
    130060: ("A36", "diphtheria"),
    130062: ("A37", "whooping cough"),
    130064: ("A38", "scarlet fever"),
    130066: ("A39", "meningococcal infection"),
    130068: ("A40", "streptococcal septicaemia"),
    130070: ("A41", "other septicaemia"),
    130072: ("A42", "actinomycosis"),
    130074: ("A43", "nocardiosis"),
    130076: ("A44", "bartonellosis"),
    130078: ("A46", "erysipelas"),
    130080: ("A48", "other bacterial diseases, not elsewhere classified"),
    130082: ("A49", "bacterial infection of unspecified site"),
    130084: ("A50", "congenital syphilis"),
    130086: ("A51", "early syphilis"),
    130088: ("A52", "late syphilis"),
    130090: ("A53", "other and unspecified syphilis"),
    130092: ("A54", "gonococcal infection"),
    130094: ("A55", "chlamydial lymphogranuloma (venereum"),
    130096: ("A56", "other sexually transmitted chlamydial diseases"),
    130100: ("A58", "granuloma inguinale"),
    130102: ("A59", "trichomoniasis"),
    130104: ("A60", "anogenital herpesviral [herpes simplex] infections"),
    130106: ("A63", "other predominantly sexually transmitted diseases, not elsewhere classified"),
    130108: ("A64", "unspecified sexually transmitted disease"),
    130112: ("A66", "yaws"),
    130114: ("A67", "pinta [carate]"),
    130116: ("A68", "relapsing fevers"),
    130118: ("A69", "other spirochaetal infections"),
    130120: ("A70", "chlamydia psittaci infection"),
    130122: ("A71", "trachoma"),
    130124: ("A74", "other diseases caused by chlamydiae"),
    130126: ("A75", "typhus fever"),
    130128: ("A77", "spotted fever [tick-borne rickettsioses]"),
    130130: ("A78", "q fever"),
    130132: ("A79", "other rickettsioses"),
    130134: ("A80", "acute poliomyelitis"),
    130136: ("A81", "atypical virus infections of central nervous system"),
    130138: ("A82", "rabies"),
    130140: ("A83", "mosquito-borne viral encephalitis"),
    130142: ("A84", "tick-borne viral encephalitis"),
    130144: ("A85", "other viral encephalitis, not elsewhere classified"),
    130146: ("A86", "unspecified viral encephalitis"),
    130148: ("A87", "viral meningitis"),
    130150: ("A88", "other viral infections of central nervous system, not elsewhere classified"),
    130152: ("A89", "unspecified viral infection of central nervous system"),
    130154: ("A90", "dengue fever [classical dengue]"),
    130156: ("A91", "dengue haemorrhagic fever"),
    130158: ("A92", "other mosquito-borne viral fevers"),
    130160: ("A93", "other arthropod-borne viral fevers, not elsewhere classified"),
    130162: ("A94", "unspecified arthropod-borne viral fever"),
    130164: ("A95", "yellow fever"),
    130168: ("A97", "dengue"),
    130170: ("A98", "other viral haemorrhagic fevers, not elsewhere classified"),
    130174: ("B00", "herpesviral [herpes simplex] infections"),
    130176: ("B01", "varicella [chickenpox]"),
    130178: ("B02", "zoster [herpes zoster]"),
    130180: ("B03", "smallpox"),
    130184: ("B05", "measles"),
    130186: ("B06", "rubella [german measles]"),
    130188: ("B07", "viral warts"),
    130190: ("B08", "other viral infections characterised by skin and mucous membrane lesions, not el"),
    130192: ("B09", "unspecified viral infection characterised by skin and mucous membrane lesions"),
    130194: ("B15", "acute hepatitis a"),
    130196: ("B16", "acute hepatitis b"),
    130198: ("B17", "other acute viral hepatitis"),
    130200: ("B18", "chronic viral hepatitis"),
    130202: ("B19", "unspecified viral hepatitis"),
    130204: ("B20", "human immunodeficiency virus [hiv] disease resulting in infectious and parasitic"),
    130206: ("B21", "human immunodeficiency virus [hiv] disease resulting in malignant neoplasms"),
    130208: ("B22", "human immunodeficiency virus [hiv] disease resulting in other specified diseases"),
    130210: ("B23", "human immunodeficiency virus [hiv] disease resulting in other conditions"),
    130212: ("B24", "unspecified human immunodeficiency virus [hiv] disease"),
    130214: ("B25", "cytomegaloviral disease"),
    130216: ("B26", "mumps"),
    130218: ("B27", "infectious mononucleosis"),
    130220: ("B30", "viral conjunctivitis"),
    130222: ("B33", "other viral diseases, not elsewhere classified"),
    130224: ("B34", "viral infection of unspecified site"),
    130226: ("B35", "dermatophytosis"),
    130228: ("B36", "other superficial mycoses"),
    130230: ("B37", "candidiasis"),
    130232: ("B38", "coccidioidomycosis"),
    130234: ("B39", "histoplasmosis"),
    130236: ("B40", "blastomycosis"),
    130240: ("B42", "sporotrichosis"),
    130242: ("B43", "chromomycosis and phaeomycotic abscess"),
    130244: ("B44", "aspergillosis"),
    130246: ("B45", "cryptococcosis"),
    130248: ("B46", "zygomycosis"),
    130250: ("B47", "mycetoma"),
    130252: ("B48", "other mycoses, not elsewhere classified"),
    130254: ("B49", "unspecified mycosis"),
    130256: ("B50", "plasmodium falciparum malaria"),
    130258: ("B51", "plasmodium vivax malaria"),
    130260: ("B52", "plasmodium malariae malaria"),
    130262: ("B53", "other parasitologically confirmed malaria"),
    130264: ("B54", "unspecified malaria"),
    130266: ("B55", "leishmaniasis"),
    130270: ("B57", "chagas' disease"),
    130272: ("B58", "toxoplasmosis"),
    130274: ("B59", "pneumocystosis"),
    130276: ("B60", "other protozoal diseases, not elsewhere classified"),
    130280: ("B65", "schistosomiasis [bilharziasis]"),
    130282: ("B66", "other fluke infections"),
    130284: ("B67", "echinococcosis"),
    130286: ("B68", "taeniasis"),
    130288: ("B69", "cysticercosis"),
    130292: ("B71", "other cestode infections"),
    130296: ("B73", "onchocerciasis"),
    130298: ("B74", "filariasis"),
    130300: ("B75", "trichinellosis"),
    130302: ("B76", "hookworm diseases"),
    130304: ("B77", "ascariasis"),
    130306: ("B78", "strongyloidiasis"),
    130308: ("B79", "trichuriasis"),
    130310: ("B80", "enterobiasis"),
    130312: ("B81", "other intestinal helminthiases, not elsewhere classified"),
    130314: ("B82", "unspecified intestinal parasitism"),
    130316: ("B83", "other helminthiases"),
    130318: ("B85", "pediculosis and phthiriasis"),
    130320: ("B86", "scabies"),
    130322: ("B87", "myiasis"),
    130324: ("B88", "other infestations"),
    130326: ("B89", "unspecified parasitic disease"),
    130328: ("B90", "sequelae of tuberculosis"),
    130330: ("B91", "sequelae of poliomyelitis"),
    130334: ("B94", "sequelae of other and unspecified infectious and parasitic diseases"),
    130336: ("B95", "streptococcus and staphylococcus as the cause of diseases classified to other ch"),
    130338: ("B96", "other bacterial agents as the cause of diseases classified to other chapters"),
    130340: ("B97", "viral agents as the cause of diseases classified to other chapters"),
    130342: ("B98", "other specified infectious agents as the cause of diseases classified to other c"),
    130344: ("B99", "other and unspecified infectious diseases"),
    130622: ("D50", "iron deficiency anaemia"),
    130624: ("D51", "vitamin b12 deficiency anaemia"),
    130626: ("D52", "folate deficiency anaemia"),
    130628: ("D53", "other nutritional anaemias"),
    130630: ("D55", "anaemia due to enzyme disorders"),
    130632: ("D56", "thalassaemia"),
    130634: ("D57", "sickle-cell disorders"),
    130636: ("D58", "other hereditary haemolytic anaemias"),
    130638: ("D59", "acquired haemolytic anaemia"),
    130640: ("D60", "acquired pure red cell aplasia [erythroblastopenia]"),
    130642: ("D61", "other aplastic anaemias"),
    130644: ("D62", "acute posthaemorrhagic anaemia"),
    130646: ("D63", "anaemia in chronic diseases classified elsewhere"),
    130648: ("D64", "other anaemias"),
    130650: ("D65", "disseminated intravascular coagulation [defibrination syndrome]"),
    130652: ("D66", "hereditary factor viii deficiency"),
    130654: ("D67", "hereditary factor ix deficiency"),
    130656: ("D68", "other coagulation defects"),
    130658: ("D69", "purpura and other haemorrhagic conditions"),
    130660: ("D70", "agranulocytosis"),
    130662: ("D71", "functional disorders of polymorphonuclear neutrophils"),
    130664: ("D72", "other disorders of white blood cells"),
    130666: ("D73", "diseases of spleen"),
    130668: ("D74", "methaemoglobinaemia"),
    130670: ("D75", "other diseases of blood and blood-forming organs"),
    130672: ("D76", "certain diseases involving lymphoreticular tissue and reticulohistiocytic system"),
    130674: ("D77", "other disorders of blood and blood-forming organs in diseases classified elsewhe"),
    130676: ("D80", "immunodeficiency with predominantly antibody defects"),
    130678: ("D81", "combined immunodeficiencies"),
    130680: ("D82", "immunodeficiency associated with other major defects"),
    130682: ("D83", "common variable immunodeficiency"),
    130684: ("D84", "other immunodeficiencies"),
    130686: ("D86", "sarcoidosis"),
    130688: ("D89", "other disorders involving the immune mechanism, not elsewhere classified"),
    130692: ("E01", "iodine-deficiency-related thyroid disorders and allied conditions"),
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
    130728: ("E24", "cushing's syndrome"),
    130730: ("E25", "adrenogenital disorders"),
    130732: ("E26", "hyperaldosteronism"),
    130734: ("E27", "other disorders of adrenal gland"),
    130736: ("E28", "ovarian dysfunction"),
    130738: ("E29", "testicular dysfunction"),
    130740: ("E30", "disorders of puberty, not elsewhere classified"),
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
    130800: ("E71", "disorders of branched-chain amino-acid metabolism and fatty-acid metabolism"),
    130802: ("E72", "other disorders of amino-acid metabolism"),
    130804: ("E73", "lactose intolerance"),
    130806: ("E74", "other disorders of carbohydrate metabolism"),
    130808: ("E75", "disorders of sphingolipid metabolism and other lipid storage disorders"),
    130810: ("E76", "disorders of glycosaminoglycan metabolism"),
    130812: ("E77", "disorders of glycoprotein metabolism"),
    130814: ("E78", "disorders of lipoprotein metabolism and other lipidaemias"),
    130816: ("E79", "disorders of purine and pyrimidine metabolism"),
    130818: ("E80", "disorders of porphyrin and bilirubin metabolism"),
    130820: ("E83", "disorders of mineral metabolism"),
    130822: ("E84", "cystic fibrosis"),
    130824: ("E85", "amyloidosis"),
    130826: ("E86", "volume depletion"),
    130828: ("E87", "other disorders of fluid, electrolyte and acid-base balance"),
    130830: ("E88", "other metabolic disorders"),
    130832: ("E89", "postprocedural endocrine and metabolic disorders, not elsewhere classified"),
    130836: ("F00", "dementia in alzheimer's disease"),
    130838: ("F01", "vascular dementia"),
    130840: ("F02", "dementia in other diseases classified elsewhere"),
    130842: ("F03", "unspecified dementia"),
    130844: ("F04", "organic amnesic syndrome, not induced by alcohol and other psychoactive substanc"),
    130846: ("F05", "delirium, not induced by alcohol and other psychoactive substances"),
    130848: ("F06", "other mental disorders due to brain damage and dysfunction and to physical disea"),
    130850: ("F07", "personality and behavioural disorders due to brain disease, damage and dysfuncti"),
    130852: ("F09", "unspecified organic or symptomatic mental disorder"),
    130854: ("F10", "mental and behavioural disorders due to use of alcohol"),
    130856: ("F11", "mental and behavioural disorders due to use of opioids"),
    130858: ("F12", "mental and behavioural disorders due to use of cannabinoids"),
    130860: ("F13", "mental and behavioural disorders due to use of sedatives or hypnotics"),
    130862: ("F14", "mental and behavioural disorders due to use of cocaine"),
    130864: ("F15", "mental and behavioural disorders due to use of other stimulants, including caffe"),
    130866: ("F16", "mental and behavioural disorders due to use of hallucinogens"),
    130868: ("F17", "mental and behavioural disorders due to use of tobacco"),
    130870: ("F18", "mental and behavioural disorders due to use of volatile solvents"),
    130872: ("F19", "mental and behavioural disorders due to multiple drug use and use of other psych"),
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
    130898: ("F34", "persistent mood [affective] disorders"),
    130900: ("F38", "other mood [affective] disorders"),
    130902: ("F39", "unspecified mood [affective] disorder"),
    130904: ("F40", "phobic anxiety disorders"),
    130906: ("F41", "other anxiety disorders"),
    130908: ("F42", "obsessive-compulsive disorder"),
    130910: ("F43", "reaction to severe stress, and adjustment disorders"),
    130912: ("F44", "dissociative [conversion] disorders"),
    130914: ("F45", "somatoform disorders"),
    130916: ("F48", "other neurotic disorders"),
    130918: ("F50", "eating disorders"),
    130920: ("F51", "nonorganic sleep disorders"),
    130922: ("F52", "sexual dysfunction, not caused by organic disorder or disease"),
    130924: ("F53", "mental and behavioural disorders associated with the puerperium, not elsewhere c"),
    130926: ("F54", "psychological and behavioural factors associated with disorders or diseases clas"),
    130928: ("F55", "abuse of non-dependence-producing substances"),
    130930: ("F59", "unspecified behavioural syndromes associated with physiological disturbances and"),
    130932: ("F60", "specific personality disorders"),
    130934: ("F61", "mixed and other personality disorders"),
    130936: ("F62", "enduring personality changes, not attributable to brain damage and disease"),
    130938: ("F63", "habit and impulse disorders"),
    130940: ("F64", "gender identity disorders"),
    130942: ("F65", "disorders of sexual preference"),
    130944: ("F66", "psychological and behavioural disorders associated with sexual development and o"),
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
    130984: ("F94", "disorders of social functioning with onset specific to childhood and adolescence"),
    130986: ("F95", "tic disorders"),
    130988: ("F98", "other behavioural and emotional disorders with onset usually occurring in childh"),
    130990: ("F99", "mental disorder, not otherwise specified"),
}

import math, pandas as pd
from functools import reduce

# dx extract_dataset has a 120s query timeout — must batch field requests
BATCH_SIZE = 20   # ICD date fields per call (safe margin below timeout)
icd_fields = [f"participant.p{fid}" for fid in sorted(FIELD_MAP.keys())]
n_batches = math.ceil(len(icd_fields) / BATCH_SIZE)
print(f"Extracting {len(icd_fields)} ICD fields in {n_batches} batches of {BATCH_SIZE}")

out_dir = "/tmp/ms_als_extraction/"
os.makedirs(out_dir, exist_ok=True)

def extract_fields(dataset_id, fields, out_path):
    cmd = ["dx", "extract_dataset", dataset_id,
           "--fields", ",".join(fields),
           "-o", out_path, "--delimiter", ","]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=200)
    if r.returncode != 0:
        raise RuntimeError(f"dx extract_dataset failed (rc={r.returncode}): {r.stderr[:500]}")
    df = pd.read_csv(out_path, low_memory=False)
    df.columns = [c.replace("participant.", "") for c in df.columns]
    os.remove(out_path)
    return df

# Step 1: get eid + year_birth
print("Fetching eid + year_birth (p34)...")
df_base = extract_fields(dataset_id, ["participant.eid", "participant.p34"],
                          out_dir + "batch_base.csv")
print(f"  {len(df_base)} participants")

# Step 2: fetch ICD fields in batches, merging onto df_base
df_merged = df_base
for batch_num, start in enumerate(range(0, len(icd_fields), BATCH_SIZE)):
    batch = icd_fields[start:start + BATCH_SIZE]
    batch_path = out_dir + f"batch_{batch_num}.csv"
    print(f"Batch {batch_num+1}/{n_batches}: {len(batch)} fields ({batch[0]}..{batch[-1]})...")
    df_batch = extract_fields(dataset_id, ["participant.eid"] + batch, batch_path)
    df_merged = df_merged.merge(df_batch, on="eid", how="left")
    print(f"  merged shape: {df_merged.shape}")

print(f"\nFinal shape: {df_merged.shape}")

# Rename p34 → year_birth for compatibility with local parser
df_merged = df_merged.rename(columns={"p34": "year_birth"})

# Save CSV
out_csv  = out_dir + "icd_abdef_first_occurrence.csv"
out_json = out_dir + "abdef_field_mapping.json"
df_merged.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}  ({os.path.getsize(out_csv)/1e6:.1f} MB)")

# Save field mapping JSON for local parser
with open(out_json, "w") as fh:
    json.dump({str(k): list(v) for k, v in FIELD_MAP.items()}, fh, indent=2)
print(f"Mapping saved: {out_json}")

# Upload both to project
print(f"\nUploading to {RAP_PROJECT}:{RAP_FOLDER} ...")
dxpy.upload_local_file(out_csv,  project=RAP_PROJECT, folder=RAP_FOLDER, wait_on_close=True)
dxpy.upload_local_file(out_json, project=RAP_PROJECT, folder=RAP_FOLDER, wait_on_close=True)
print("Upload complete!")

print("\n=== DONE ===\nDownload with:")
print('  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_abdef_first_occurrence.csv"')
print('  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/abdef_field_mapping.json"')
