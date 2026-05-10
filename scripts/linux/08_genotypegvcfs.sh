#!/usr/bin/env bash

# =========================================================
# Step 08. Joint genotyping using GATK GenotypeGVCFs
# Project: Integrated MRONJ risk prediction model
# Input:
#   GenomicsDB workspace
# Output:
#   Joint-genotyped cohort VCF
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 08] Joint genotyping using GATK GenotypeGVCFs"
echo "Reference genome : ${REFERENCE}"
echo "GenomicsDB input : ${GENOMICSDB_DIR}"
echo "Output VCF       : ${RAW_JOINT_VCF}"

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference genome not found: ${REFERENCE}" >&2
    exit 1
fi

if [ ! -d "${GENOMICSDB_DIR}" ]; then
    echo "ERROR: GenomicsDB workspace not found: ${GENOMICSDB_DIR}" >&2
    echo "Please run 07_genomicsdb_import.sh first." >&2
    exit 1
fi

gatk --java-options "-Xmx${GATK_JAVA_MEM}" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V "gendb://${GENOMICSDB_DIR}" \
    -O "${RAW_JOINT_VCF}" \
    > "${LOG_DIR}/genotypegvcfs.log" 2>&1

gatk IndexFeatureFile \
    -I "${RAW_JOINT_VCF}" \
    > "${LOG_DIR}/final_variants.index.log" 2>&1

echo "[Step 08] Joint genotyping completed successfully."
