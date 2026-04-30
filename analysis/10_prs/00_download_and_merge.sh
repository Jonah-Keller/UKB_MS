#!/bin/bash
# 00_download_and_merge.sh
# Run once RAP jobs complete — downloads .profile files and merges into PRS scores
#
# Check job status first:
#   dx find jobs --created-after=-6h | grep "Swiss Army"
#
# When all 44 jobs show (done), run this script from UKB_MS project root

set -e
PROJ_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PROFILE_DIR="${PROJ_ROOT}/data/ukb/genetics/prs_profiles"
mkdir -p "$PROFILE_DIR"

echo "Downloading .profile files from RAP /data/prs/results/..."
dx download "/data/prs/results/*.profile" -o "$PROFILE_DIR/" --no-progress

echo "Downloaded $(ls ${PROFILE_DIR}/*.profile 2>/dev/null | wc -l) profile files"

# Verify expected counts
for pgid in PGS000809 PGS004699; do
    n=$(ls ${PROFILE_DIR}/prs_${pgid}_c*.profile 2>/dev/null | wc -l)
    echo "  ${pgid}: ${n}/22 chromosomes"
done

echo ""
echo "Running merge script..."
/usr/local/bin/Rscript "${PROJ_ROOT}/analysis/10_prs/00_merge_prs_scores.R"
