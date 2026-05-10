#!/usr/bin/env bash

# =========================================================
# Step 09. Select PASS variants from joint-genotyped VCF
# Project: Integrated MRONJ risk prediction model
# Input:
#   Joint-genotyped cohort VCF
# Output:
#   Filtered VCF containing variants that passed previous filters
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 09] Selecting PASS variants using GATK SelectVariants"
echo "Input VCF : ${RAW_JOINT_VCF}"
echo "Output VCF: ${SELECTED_PASS_VCF}"

if [ ! -f "${RAW_JOINT_VCF}" ]; then
    echo "ERROR: Input VCF not found: ${RAW_JOINT_VCF}" >&2
    echo "Please run 08_joint_genotyping_genotypegvcfs.sh first." >&2
    exit 1
fi

gatk SelectVariants \
    -R "${REFERENCE}" \
    -V "${RAW_JOINT_VCF}" \
    -O "${SELECTED_PASS_VCF}" \
    --exclude-filtered true \
    > "${LOG_DIR}/select_pass_variants.log" 2>&1

gatk IndexFeatureFile \
    -I "${SELECTED_PASS_VCF}" \
    > "${LOG_DIR}/selected_pass_variants.index.log" 2>&1

echo "[Step 09] PASS variant selection completed successfully."
