#!/usr/bin/env bash
# ============================================================================
# bootstrap_new_disease.sh — set up a new UKB_<DISEASE> replication repo
# ============================================================================
# Clones the UKB_MS template into /Users/jonahkeller/ELab/UKB_<DISEASE>,
# creates a private GitHub repo at github.com/Jonah-Keller/UKB_<DISEASE>,
# symlinks the shared UKB data extract, scaffolds a disease.yaml stub with
# placeholders for every required field, and prints the OpenEvidence /
# PGS-Catalog / UKB-Showcase queries the user needs to run for the agent
# to fill in the config.
#
# Usage:
#   bash scripts/preflight/bootstrap_new_disease.sh <disease_short>
#
# Example:
#   bash scripts/preflight/bootstrap_new_disease.sh stroke
#
# After this exits, the agent (or you) can:
#   1. Send the user the printed query block.
#   2. Edit configs/disease.yaml with the user's answers.
#   3. cd into the new repo and run scripts/preflight/preflight.R.
#   4. Run the pipeline per AGENT_HANDOFF_disease_replication.md §5.
# ============================================================================

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <disease_short>"
    echo "Example: $0 stroke"
    exit 1
fi

DISEASE_SHORT="$1"
DISEASE_UPPER="$(echo "$DISEASE_SHORT" | tr '[:lower:]' '[:upper:]')"
TEMPLATE_REPO="Jonah-Keller/UKB_MS"
NEW_REPO="Jonah-Keller/UKB_${DISEASE_UPPER}"
ELAB_DIR="/Users/jonahkeller/ELab"
NEW_DIR="${ELAB_DIR}/UKB_${DISEASE_UPPER}"
DATA_SOURCE="${ELAB_DIR}/UKB_MS/data"

echo "══ Bootstrapping UKB_${DISEASE_UPPER} ══"
echo "Template:    ${TEMPLATE_REPO}"
echo "New repo:    ${NEW_REPO}"
echo "Local path:  ${NEW_DIR}"
echo

# ── Guards ──────────────────────────────────────────────────────────────────
if [[ -d "${NEW_DIR}" ]]; then
    echo "ERROR: ${NEW_DIR} already exists. Pick a different short name or"
    echo "       remove the existing directory first."
    exit 1
fi

if ! command -v gh >/dev/null 2>&1; then
    echo "ERROR: gh (GitHub CLI) not found in PATH."
    exit 1
fi

if [[ ! -d "${DATA_SOURCE}" ]]; then
    echo "ERROR: ${DATA_SOURCE} does not exist; cannot symlink shared data."
    exit 1
fi

# ── 1. GitHub repo: try template-clone first, fall back to manual ───────────
echo "[1/4] Creating GitHub repo + cloning..."
if gh repo create "${NEW_REPO}" --template "${TEMPLATE_REPO}" --private \
                                --clone --description "UKB ${DISEASE_UPPER} replication of UKB_MS plasma proteomics template" \
                                >/dev/null 2>&1; then
    echo "      template-clone succeeded (UKB_MS is a GitHub template repo)"
    mv "UKB_${DISEASE_UPPER}" "${NEW_DIR}"
else
    echo "      template-clone failed (UKB_MS not flagged as a template); falling back to git clone"
    git clone "git@github.com:${TEMPLATE_REPO}.git" "${NEW_DIR}"
    cd "${NEW_DIR}"
    gh repo create "${NEW_REPO}" --private --source=. --remote=origin --push \
       --description "UKB ${DISEASE_UPPER} replication of UKB_MS plasma proteomics template"
    git remote set-url origin "git@github.com:${NEW_REPO}.git"
    cd "${ELAB_DIR}"
fi

cd "${NEW_DIR}"

# Add UKB_MS as upstream so future template fixes can be pulled in
git remote add upstream "git@github.com:${TEMPLATE_REPO}.git" 2>/dev/null || true
git fetch upstream main >/dev/null 2>&1 || true

# ── 2. Symlink the shared UKB data extract ──────────────────────────────────
echo "[2/4] Symlinking shared UKB data extract..."
ln -s "${DATA_SOURCE}" data
echo "      data -> ${DATA_SOURCE}"

# ── 3. Scaffold a disease.yaml stub with placeholders ───────────────────────
echo "[3/4] Writing configs/disease.yaml stub with TBD placeholders..."
cat > configs/disease.yaml <<YAML
# ============================================================================
# disease.yaml — UKB_${DISEASE_UPPER} replication of the UKB_MS proteomics template
# ============================================================================
# This is a STUB.  Fill in every TBD field before running the pipeline.
# See AGENT_HANDOFF_disease_replication.md §1 for the OpenEvidence / PGS-Catalog
# / UKB-Showcase queries that produce the values for these fields.
# ============================================================================

# ── Identifying labels ─────────────────────────────────────────────────────
cohort_short:        "${DISEASE_SHORT}"
disease_long:        "TBD — full disease name"
disease_short_caps:  "${DISEASE_UPPER}"
project_short:       "ukb_${DISEASE_SHORT}_proteomics"
project_title:       "TBD — replication paper title"
comparison_cohort_short: "ms"

# ── Cohort definition (REQUIRED) ───────────────────────────────────────────
icd_codes:                ["TBD"]                       # primary case ICD-10 codes
first_occurrence_field:   "TBD_p1XXXXX"                 # UKB date-of-first-occurrence field
cohort_status_col:        "${DISEASE_SHORT}_status"
status_values:
  pre_onset:  "pre_onset"
  post_onset: "post_onset"
  control:    "control"

control_exclusion_codes:
  TBD_group:   ["TBD"]

# ── HLA risk locus (set hla_allele: "NONE" if no replicable HLA association) ─
hla_allele:           "TBD"                             # OR "NONE" to skip stage 07
hla_search_patterns:  ["TBD"]
hla_carrier_col:      "TBD_carrier"
hla_dosage_col:       "TBD_dosage"
hla_excluded_subtypes: []
hla_dosage_threshold:  0.8
hla_direction:         "risk"                           # or "protective"
hla_field:             22182

# ── Mendelian carrier (optional — set null to skip stage 07b) ──────────────
mendelian_carrier_col:        null                      # e.g. "apoe_e4_carrier", "c9orf72_expansion_carrier"
mendelian_threshold_repeats:  null
mendelian_threshold_strict:   null
mendelian_source:             null
other_mendelian_carrier_col:  null

# ── Phenomic / exposure stratifiers (optional; empty list to skip stage 17) ─
phenomic_stratifiers: {}

# ── Sex handling (REQUIRED) ────────────────────────────────────────────────
sex_handling:
  primary_use: "stratifier"                             # OR "covariate" if N too small
  field:       31
  do_not_stratify_on_sex_alone: false
  interaction_terms_to_test: []

# ── Polygenic risk score (empty list = skip stage 10) ──────────────────────
prs_pgs_ids:          []
prs_label:            "${DISEASE_UPPER} PRS"
prs_combined_col:     "prs_score"

# ── Programmatic top-N protein selection (template-wide; do not change) ────
top_n_volcano:    20
top_n_heatmap:    30
top_n_trajectory: 3
min_sig_for_fdr:  5

# ── CNS biology reference (only applies to CNS diseases — drop otherwise) ──
cns_injury_markers:   ["nefl", "mog", "gfap"]
trajectory_proteins:  ["NEFL", "GFAP"]
YAML
echo "      configs/disease.yaml staged with TBD placeholders"

# ── 4. Print the agent-ready prompt block ───────────────────────────────────
echo "[4/4] Printing agent prompt block..."
echo

cat <<PROMPT

══ AGENT-READY PROMPT (paste into a fresh agent) ══

You are replicating the UKB plasma-proteomics study on ${DISEASE_UPPER}.
The repo is cloned and configured at ${NEW_DIR} (origin: github.com/${NEW_REPO}).

Read AGENT_HANDOFF_disease_replication.md cover to cover before doing
anything substantive. The TL;DR at the top has the quick-start commands.
Pattern table in §8 documents 20 known bug patterns we've already fixed
once — recognise analogues immediately.

Send the user these 6 queries and wait for answers before editing
configs/disease.yaml:

----------------------------------------------------------------------
For the ${DISEASE_UPPER} replication, please run these queries:

1. Open Evidence: "What is the UK Biobank first-occurrence field ID for
   ${DISEASE_UPPER}? Confirm the primary ICD-10 code(s) and any subtype-
   specific code (e.g. ischemic-only stroke)."

2. Open Evidence: "List ICD-10 codes that should disqualify a participant
   from the control pool for ${DISEASE_UPPER}. Include diseases that share
   clinical phenotype, treatment confounding, or biology with ${DISEASE_UPPER}."

3. Open Evidence: "In the largest GWAS of ${DISEASE_UPPER} to date, did any
   HLA-region locus reach genome-wide significance? Name the lead SNP and
   tagged HLA allele if so. If no replicable HLA finding exists, say so
   explicitly. Otherwise list candidate alleles by replicated effect size."

4. Open Evidence: "Are there Mendelian or large-effect variants enriched
   in ${DISEASE_UPPER} cases that could be extracted from UKB whole-exome
   or whole-genome sequencing? Report carrier rate, expected count among
   UKB ${DISEASE_UPPER} cases, and whether the variant is callable from the
   standard UKB-RAP pipeline."

5. PGS Catalog (https://www.pgscatalog.org/trait_search/?q=${DISEASE_SHORT}):
   "List 1–2 best polygenic scores. For each: PGS ID, source GWAS, variant
   count, training-data ancestry, validation AUC or pseudo-R² in
   independent EUR cohorts."

6. Open Evidence: "Rank the top 3 environmental / lifestyle risk factors
   for ${DISEASE_UPPER} by replicated effect size in cohort studies. State
   which UKB fields capture them. Flag any that cannot be reliably captured."
----------------------------------------------------------------------

After you have answers and have edited disease.yaml:

  cd ${NEW_DIR}
  Rscript scripts/preflight/preflight.R   # MUST pass before pipeline runs

Then proceed per AGENT_HANDOFF §5.

══ END PROMPT ══

PROMPT

echo "Bootstrap complete."
echo "Next: cd ${NEW_DIR} and either drop the printed prompt into a fresh"
echo "agent, or fill in disease.yaml and run preflight yourself."
