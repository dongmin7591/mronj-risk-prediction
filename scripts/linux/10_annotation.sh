#!/usr/bin/env bash

# =========================================================
# Step 10. Functional annotation using GATK Funcotator
# Project: Integrated MRONJ risk prediction model
# Input:
#   Selected PASS variants VCF
# Output:
#   MAF file annotated by Funcotator
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 10] Functional annotation using GATK Funcotator"
echo "Input VCF            : ${SELECTED_PASS_VCF}"
echo "Funcotator datasource: ${FUNCOTATOR_DATASOURCE}"
echo "Output directory     : ${MAF_DIR}"

if [ ! -f "${SELECTED_PASS_VCF}" ]; then
    echo "ERROR: Input VCF not found: ${SELECTED_PASS_VCF}" >&2
    echo "Please run 09_select_pass_variants.sh first." >&2
    exit 1
fi

if [ ! -d "${FUNCOTATOR_DATASOURCE}" ]; then
    echo "ERROR: Funcotator data source directory not found: ${FUNCOTATOR_DATASOURCE}" >&2
    exit 1
fi

output_maf="${MAF_DIR}/selected_pass_variants.funcotator.maf.gz"

gatk --java-options "-Xmx${GATK_JAVA_MEM}" Funcotator \
    -R "${REFERENCE}" \
    -V "${SELECTED_PASS_VCF}" \
    -O "${output_maf}" \
    --output-file-format MAF \
    --data-sources-path "${FUNCOTATOR_DATASOURCE}" \
    --ref-version "${FUNCOTATOR_REF_VERSION}" \
    --annotation-default tumor_barcode:MRONJ_WGS_COHORT \
    > "${LOG_DIR}/funcotator_annotation_maf.log" 2>&1

echo "[Step 10] Funcotator annotation completed successfully."
echo "Output MAF: ${output_maf}"
