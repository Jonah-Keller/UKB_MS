#!/bin/bash
set -e

echo "Building fields list..."
cat << 'FIELDS' > /Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/my_fields.txt
participant.p22009_a1
participant.p22009_a2
participant.p22009_a3
participant.p22009_a4
participant.p22009_a5
participant.p22009_a6
participant.p22009_a7
participant.p22009_a8
participant.p22009_a9
participant.p22009_a10
participant.p41270
participant.p41202
participant.p41204
participant.p40001_i0
participant.p40002_i0
participant.p40000_i0
participant.p40007_i0
participant.p22182
participant.p22006
participant.p22019
participant.p22021
participant.p22027
participant.p22001
participant.p22000
participant.p22020
FIELDS

echo "Running full dx extract_dataset locally (zero cloud cost)..."
/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx extract_dataset \
  project-GxX43xjJpjp2G7XK6ffz0qb1:record-GxX62v8J88vBkPBk3pXXv0J8 \
  --fields-file /Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/my_fields.txt \
  -o /Users/jonahkeller/ELab/UKB_MS/data/ukb/rap_extraction/local_massive_extract.csv

echo "Done!"
