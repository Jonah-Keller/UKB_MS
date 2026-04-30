#!/bin/bash
/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx run table-exporter \
  -idataset_or_cohort_or_dashboard=project-GxX43xjJpjp2G7XK6ffz0qb1:record-GxX62v8J88vBkPBk3pXXv0J8 \
  -ientity=participant \
  -ifield_titles="participant.p22009_a1,participant.p41270" \
  -ioutput_name="test_export" \
  -icoding_option="RAW" \
  -iheader_style="UKB-Format" \
  --yes
