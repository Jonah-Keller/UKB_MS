"""RAP Extraction — HES R/S/T/Z ICD-10 chapters via fields 41270/41280

These chapters (R=Symptoms, S=Injuries, T=Poisoning/Adverse, U=COVID, V=Transport,
Y=External causes, Z=Health contact factors) are NOT in first-occurrence fields
(p130xxx). They must be pulled from Hospital Episode Statistics:
  - p41270_a* : ICD-10 diagnosis codes (up to 213 per person, DataCoding=19)
  - p41280_a* : Corresponding first in-patient diagnosis dates (259 arrays)

Clinically relevant for MS research:
  R: R20 (sensory), R25 (motor), R26 (gait), R32-R39 (bladder/bowel),
     R41 (cognitive), R44 (visual), R51 (headache), R53 (fatigue),
     R73 (glycaemia), R90 (CNS imaging abnormal)
  T: T36-T50 (drug adverse effects, incl. corticosteroids T45)
  Z: Z51 (treatment), Z82/Z83/Z84/Z86/Z87 (family/personal history),
     Z96 (implants - DBS for MS tremor)
  S/U/V/Y included but not filtered - written only for codes with ≥1 case

Output:
  /data/ms_als_extraction/icd_hes_rstz_first_occurrence.csv
  /data/ms_als_extraction/hes_rstz_summary.json

After downloading:
  python analysis/00_rap_extraction/10_parse_hes_rstz.py

Upload & run:
  dx upload analysis/00_rap_extraction/09_rap_extract_hes_rstz.py \\
      --destination "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/code/09_rap_extract_hes_rstz.py"
  FILE_ID=$(dx find data --name 09_rap_extract_hes_rstz.py --project project-J2P0fqjJ1Gg4gByJ6BJ51pYj --brief)
  dx run dxjupyterlab \\
      --instance-type mem1_ssd1_v2_x16 \\
      --project project-J2P0fqjJ1Gg4gByJ6BJ51pYj \\
      --destination /data/ms_als_extraction \\
      -iin="${FILE_ID}" \\
      -icmd="python /home/dnanexus/in/09_rap_extract_hes_rstz.py" \\
      --name "hes_rstz_extraction" --brief --yes
"""

import os, json, math, subprocess
import numpy as np
import pandas as pd
import dxpy

RAP_PROJECT = "project-J2P0fqjJ1Gg4gByJ6BJ51pYj"
RAP_FOLDER  = "/data/ms_als_extraction"
TARGET_CHAPTERS = set("RSTUVYZ")

# ── Find dataset ──────────────────────────────────────────────────────────────
dispensed_dataset = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)
dataset_id = dispensed_dataset["id"]
print(f"Dataset: {dataset_id}")

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

# ── Pass 1: Extract p41270 (JSON list of ALL codes) + year_birth ──────────────
# p41270 is a single JSON-array field: ["G35","I10","R53",...]
# p41280_a{i} is the date for the code at sorted position i
print("Pass 1: Extracting p41270 (all diagnosis codes as JSON) + p34...")
df_codes = extract_fields(dataset_id,
                           ["participant.eid", "participant.p41270", "participant.p34"],
                           out_dir + "codes.csv")
print(f"  {len(df_codes)} participants")

# Build lookups from p41270
yob_dict: dict = {}
eid_pos_to_icd: dict[tuple, str] = {}  # (eid, position) -> icd3 for target codes
max_pos = 0

for row in df_codes.itertuples(index=False):
    eid, raw_codes, yob = row
    yob_dict[eid] = yob
    if not raw_codes or str(raw_codes) in ("", "[]", "nan", "None"):
        continue
    try:
        codes = json.loads(str(raw_codes))
    except Exception:
        continue
    for pos, code in enumerate(codes):
        if not code or code[0] not in TARGET_CHAPTERS:
            continue
        eid_pos_to_icd[(eid, pos)] = code[:3]
        if pos > max_pos:
            max_pos = pos

del df_codes
print(f"  {len(eid_pos_to_icd)} (eid, position) pairs with R/S/T/U/V/Y/Z codes")
print(f"  Max code position needed: {max_pos} (will pull p41280_a0..a{max_pos})")

# ── Pass 2: Extract p41280 date columns in batches, match by (eid, position) ──
N_DATE_SLOTS = min(max_pos + 1, 259)   # p41280 has 259 arrays (a0..a258)
BATCH_SIZE   = 15                       # date columns per dx extract_dataset call
n_batches    = math.ceil(N_DATE_SLOTS / BATCH_SIZE)

earliest: dict[tuple, str] = {}  # (eid, icd3) -> earliest date string
total_seen = 0

print(f"\nPass 2: Extracting {N_DATE_SLOTS} date columns in {n_batches} batches...")
for batch_num in range(n_batches):
    slot_start = batch_num * BATCH_SIZE
    slot_end   = min(slot_start + BATCH_SIZE, N_DATE_SLOTS)
    slots      = range(slot_start, slot_end)

    date_fields = [f"participant.p41280_a{i}" for i in slots]
    print(f"  Batch {batch_num+1}/{n_batches}: p41280_a{slot_start}..a{slot_end-1}...")
    df_dates = extract_fields(dataset_id, ["participant.eid"] + date_fields,
                               out_dir + f"dates_batch_{batch_num}.csv")

    for pos in slots:
        dt_col = f"p41280_a{pos}"
        if dt_col not in df_dates.columns:
            continue
        sub = df_dates[["eid", dt_col]].dropna(subset=[dt_col])
        for row in sub.itertuples(index=False):
            eid, date_val = row
            icd3 = eid_pos_to_icd.get((eid, pos))
            if icd3 is None:
                continue
            date_str = str(date_val)[:10]
            key = (eid, icd3)
            total_seen += 1
            if key not in earliest or date_str < earliest[key]:
                earliest[key] = date_str

    del df_dates
    print(f"    Cumulative (eid, icd3) pairs: {len(earliest)}")

print(f"\nTotal date matches: {total_seen}")
print(f"Unique (eid, icd3) pairs: {len(earliest)}")

# ── Build output dataframe ────────────────────────────────────────────────────
records = []
for (eid, icd3), date_str in earliest.items():
    yob = yob_dict.get(eid)
    try:
        year = int(date_str[:4])
        age = (year - int(yob)) if yob and not np.isnan(float(yob)) else np.nan
        if age < 0:
            age = np.nan
    except Exception:
        age = np.nan

    records.append({
        "eid": eid,
        "icd3": icd3,
        "chapter": icd3[0],
        "first_occurrence_date": date_str,
        "age_at_diagnosis": age,
    })

out_df = pd.DataFrame(records)
print(f"Output shape: {out_df.shape}")

chapter_counts = out_df.groupby("chapter")["eid"].count().to_dict()
code_counts = out_df.groupby("icd3")["eid"].count().sort_values(ascending=False)
print("\nPer-chapter case counts:")
for ch in sorted(chapter_counts):
    print(f"  Chapter {ch}: {chapter_counts[ch]} (eid, code) pairs")
print("\nTop 30 codes:")
print(code_counts.head(30).to_string())

# ── Save and upload ───────────────────────────────────────────────────────────
out_csv  = out_dir + "icd_hes_rstz_first_occurrence.csv"
json_path = out_dir + "hes_rstz_summary.json"

out_df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}  ({os.path.getsize(out_csv)/1e6:.1f} MB)")

summary = {
    "n_participants_any": int(out_df["eid"].nunique()),
    "n_unique_icd3_codes": int(out_df["icd3"].nunique()),
    "n_total_records": len(out_df),
    "per_chapter": {ch: int(n) for ch, n in chapter_counts.items()},
    "top_50_codes": {k: int(v) for k, v in code_counts.head(50).items()},
}
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nUploading to {RAP_PROJECT}:{RAP_FOLDER} ...")
dxpy.upload_local_file(out_csv,   project=RAP_PROJECT, folder=RAP_FOLDER, wait_on_close=True)
dxpy.upload_local_file(json_path, project=RAP_PROJECT, folder=RAP_FOLDER, wait_on_close=True)
print("Upload complete!")

print("\nDone! Download with:")
print(f'  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/icd_hes_rstz_first_occurrence.csv"')
print(f'  dx download "project-J2P0fqjJ1Gg4gByJ6BJ51pYj:/data/ms_als_extraction/hes_rstz_summary.json"')
