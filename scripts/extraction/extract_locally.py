import subprocess
import os

fields = set()
with open('/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/all_fields.csv') as f:
    for line in f:
        col = line.strip().split('\t')[0]
        # Keep genetic PCs up to PC10
        if col.startswith('participant.p22009_a') and int(col.split('_a')[1]) <= 10:
            fields.add(col)
        # HES codes (all of 41270, 41202, 41204)
        if col.startswith('participant.p41270') or col.startswith('participant.p41202') or col.startswith('participant.p41204'):
            fields.add(col)
        # Death causes
        if col.startswith('participant.p40001_') or col.startswith('participant.p40002_') or col.startswith('participant.p40000_') or col.startswith('participant.p40007_'):
            fields.add(col)
        # Basic fields
        if col in ['participant.p31', 'participant.p34', 'participant.p131016', 'participant.p131042', 'participant.p21003_i0', 'participant.p189_i0', 'participant.p54_i0']:
            fields.add(col)
        # Genetic QC
        if col in ['participant.p22006', 'participant.p22019', 'participant.p22021', 'participant.p22027', 'participant.p22001', 'participant.p22000', 'participant.p22020']:
            fields.add(col)
        # HLA
        if col == "participant.p22182":
            fields.add(col)

fields_file = '/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/my_fields.txt'
with open(fields_file, 'w') as f:
    f.write('\n'.join(list(fields)))

print(f"Total fields to request: {len(fields)}")

cmd = [
    '/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx', 
    'extract_dataset', 'project-GxX43xjJpjp2G7XK6ffz0qb1:record-GxX62v8J88vBkPBk3pXXv0J8', 
    '--fields-file', fields_file, 
    '-o', '/Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/local_massive_extract.csv'
]

print('Running local extraction...')
subprocess.run(cmd, check=True)
print('Done!')
