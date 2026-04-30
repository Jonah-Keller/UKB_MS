import pandas as pd
from pathlib import Path

# Paths
RAP_DIR = Path("/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction")
GEN_DIR = Path("/Users/jonahkeller/ELab/UKB_MS/data/ukb/genetics")
COV_DIR = Path("/Users/jonahkeller/ELab/UKB_MS/data/ukb/covariates")
GEN_DIR.mkdir(parents=True, exist_ok=True)

print("Loading local RAP extracts...")
# fast_extract contains PCs, demographics
df_fast = pd.read_csv(RAP_DIR / 'covariates_pcs.csv')
rename_fast = {"participant.eid": "eid", "participant.p31": "sex", "participant.p34": "year_of_birth"}
for i in range(1, 11):
    rename_fast[f"participant.p22009_a{i}"] = f"PC{i}"
df_fast = df_fast.rename(columns=rename_fast)

# Build genetic PCs table
pc_cols = ['eid'] + [f'PC{i}' for i in range(1, 11)]
df_pcs = df_fast[pc_cols].dropna(subset=['PC1'])
df_pcs.to_csv(GEN_DIR / 'genetic_pcs_top10.csv', index=False)
print(f"Saved {len(df_pcs)} genetic PCs")

# death_causes contains 40001
df_death = pd.read_csv(RAP_DIR / 'death_causes.csv', low_memory=False)
df_death = df_death.rename(columns={"participant.eid": "eid"})
g12_death_eids = set()
for col in df_death.columns:
    if col != 'eid':
        vals = df_death[col].astype(str)
        hits = df_death[vals.str.contains('G12.2|G122|G12', regex=True, na=False)]['eid']
        g12_death_eids.update(hits.values)

df_g12_death = pd.DataFrame({'eid': sorted(list(g12_death_eids)), 'g12_death': True})
df_g12_death.to_csv(RAP_DIR / 'g12_mnd_death_confirmed.csv', index=False)
print(f"Identified {len(g12_death_eids)} MND/ALS death confirmations")

print("All unified local datasets securely on disk!")
