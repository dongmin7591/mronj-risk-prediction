#!/usr/bin/env bash

# =========================================================
# Step 04. Base quality score recalibration using GATK BQSR
# Project: Integrated MRONJ risk prediction model
# Input:
#   Duplicate-marked BAM files
# Output:
#   Final recalibrated BAM files
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 04] Base quality score recalibration using GATK"
echo "Input directory : ${DEDUP_BAM_DIR}"
echo "Output directory: ${FINAL_BAM_DIR}"

required_files=(
    "${REFERENCE}"
    "${DBSNP}"
    "${HAPMAP}"
    "${MILLS}"
    "${OMNI}"
)

for file in "${required_files[@]}"; do
    if [ ! -f "${file}" ]; then
        echo "ERROR: Required reference file not found: ${file}" >&2
        exit 1
    fi
done

shopt -s nullglob

dedup_bam_files=("${DEDUP_BAM_DIR}"/*.dedup.bam)

if [ ${#dedup_bam_files[@]} -eq 0 ]; then
    echo "ERROR: No duplicate-marked BAM files found in ${DEDUP_BAM_DIR}" >&2
    exit 1
fi

for input_bam in "${dedup_bam_files[@]}"; do
    sample_name="$(basename "${input_bam}" .dedup.bam)"

    echo "Processing sample: ${sample_name}"

    recal_table="${BQSR_DIR}/${sample_name}.recal_data.table"
    final_bam="${FINAL_BAM_DIR}/${sample_name}.final.bam"

    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${input_bam}" \
        --known-sites "${DBSNP}" \
        --known-sites "${HAPMAP}" \
        --known-sites "${MILLS}" \
        --known-sites "${OMNI}" \
        -O "${recal_table}" \
        > "${LOG_DIR}/${sample_name}.base_recalibrator.log" 2>&1

    gatk ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${input_bam}" \
        --bqsr-recal-file "${recal_table}" \
        -O "${final_bam}" \
        > "${LOG_DIR}/${sample_name}.apply_bqsr.log" 2>&1

    samtools index -@ "${THREADS}" "${final_bam}"

    echo "Completed BQSR: ${sample_name}"
done

echo "[Step 04] Base quality score recalibration completed."
