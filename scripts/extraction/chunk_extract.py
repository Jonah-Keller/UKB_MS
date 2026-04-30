import subprocess
import os

def run_extract(fields, output_name):
    fields_file = f"/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/{output_name}_fields.txt"
    with open(fields_file, "w") as f:
        f.write('\n'.join(fields))
    
    cmd = [
        '/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx', 
        'extract_dataset', 'project-GxX43xjJpjp2G7XK6ffz0qb1:record-GxX62v8J88vBkPBk3pXXv0J8', 
        '--fields-file', fields_file, 
        '-o', f'/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/{output_name}.csv'
    ]
    print(f"Extracting {output_name} ({len(fields)} fields)...")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0:
        print(f"✅ Success: {output_name}.csv")
    else:
        print(f"❌ Failed {output_name}: {res.stderr}")

valid = set()
with open('/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/all_fields.csv') as f:
    for line in f:
        valid.add(line.split('\t')[0])

def get_fields_for_prefix(prefix, max_arrays=None):
    res = []
    for v in valid:
        if v.startswith(prefix):
            if '_a' in v and max_arrays is not None:
                idx = int(v.split('_a')[1])
                if idx <= max_arrays:
                    res.append(v)
            else:
                res.append(v)
    return res

# Chunk 1: Covariates and Genetic QC
qc_cov = []
for p in ['participant.p22006', 'participant.p22019', 'participant.p22021', 'participant.p22027', 'participant.p22001', 'participant.p22000', 'participant.p22020', 'participant.p21003_i0', 'participant.p31', 'participant.p34']:
    qc_cov.extend(get_fields_for_prefix(p))
run_extract(qc_cov, "covariates_qc")

# Chunk 2: Death Causes
death = get_fields_for_prefix('participant.p40001_i0', max_arrays=10) + get_fields_for_prefix('participant.p40002_i0', max_arrays=10)
run_extract(death, "death_causes")

# Chunk 3: Genetic PCs (Top 10)
pcs = get_fields_for_prefix('participant.p22009', max_arrays=10)
run_extract(pcs, "genetic_pcs")

# Chunk 4: HLA
hla = get_fields_for_prefix('participant.p22182')
run_extract(hla, "hla_imputation")
