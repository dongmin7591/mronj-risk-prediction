#!/usr/bin/env bash

# =========================================================
# Step 03. Duplicate marking using GATK MarkDuplicatesSpark
# Project: Integrated MRONJ risk prediction model
# Input:
#   Sorted BAM files
# Output:
#   Duplicate-marked BAM files
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 03] Duplicate marking using GATK MarkDuplicatesSpark"
echo "Input directory : ${BAM_DIR}"
echo "Output directory: ${DEDUP_BAM_DIR}"

shopt -s nullglob

bam_files=("${BAM_DIR}"/*.sorted.bam)

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "ERROR: No sorted BAM files found in ${BAM_DIR}" >&2
    exit 1
fi

for input_bam in "${bam_files[@]}"; do
    sample_name="$(basename "${input_bam}" .sorted.bam)"

    echo "Processing sample: ${sample_name}"

    output_bam="${DEDUP_BAM_DIR}/${sample_name}.dedup.bam"

    gatk MarkDuplicatesSpark \
        -I "${input_bam}" \
        -O "${output_bam}" \
        --duplicate-tagging-policy OpticalOnly \
        --conf "spark.executor.instances=10" \
        --conf "spark.executor.cores=10" \
        --conf "spark.executor.memory=8G" \
        > "${LOG_DIR}/${sample_name}.markduplicates.log" 2>&1

    samtools index -@ "${THREADS}" "${output_bam}"

    echo "Completed duplicate marking: ${sample_name}"
done

echo "[Step 03] Duplicate marking completed."
