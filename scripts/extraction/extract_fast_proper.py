import subprocess
fields = [
    'participant.eid',
    'participant.p22009_a1', 'participant.p22009_a2', 'participant.p22009_a3',
    'participant.p22009_a4', 'participant.p22009_a5', 'participant.p22009_a6',
    'participant.p22009_a7', 'participant.p22009_a8', 'participant.p22009_a9', 'participant.p22009_a10',
    'participant.p31', 'participant.p34', 'participant.p131016', 'participant.p131042', 'participant.p40007_i0'
]
with open('/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/fast_proper.txt', 'w') as f:
    f.write('\n'.join(fields))

print("Extracting proper covariates & PCs with EID...")
cmd = [
    '/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx', 
    'extract_dataset', 'project-GxX43xjJpjp2G7XK6ffz0qb1:record-GxX62v8J88vBkPBk3pXXv0J8', 
    '--fields-file', '/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/fast_proper.txt', 
    '-o', '/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/covariates_pcs.csv'
]
subprocess.run(cmd, check=True)
print("Extracted covariates_pcs.csv")

# Death causes
death_fields = ['participant.eid']
with open('/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/all_fields.csv') as f:
    for line in f:
        col = line.strip().split('\t')[0]
        if col.startswith('participant.p40001_i0') or col.startswith('participant.p40002_i0'):
            death_fields.append(col)

with open('/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/death_proper.txt', 'w') as f:
    f.write('\n'.join(death_fields[:15])) # Limit arrays to bypass timeout

print("Extracting proper death causes with EID...")
cmd2 = [
    '/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx', 
    'extract_dataset', 'project-GxX43xjJpjp2G7XK6ffz0qb1:record-GxX62v8J88vBkPBk3pXXv0J8', 
    '--fields-file', '/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/death_proper.txt', 
    '-o', '/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/death_causes.csv'
]
subprocess.run(cmd2, check=True)
print("Extracted death_causes.csv")
