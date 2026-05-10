#!/usr/bin/env bash

# =========================================================
# Step 05. GVCF calling using GATK HaplotypeCaller
# Project: Integrated MRONJ risk prediction model
# Input:
#   Final recalibrated BAM files
# Output:
#   Per-sample GVCF files
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 05] GVCF calling using GATK HaplotypeCaller"
echo "Input directory : ${FINAL_BAM_DIR}"
echo "Output directory: ${GVCF_DIR}"

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference genome not found: ${REFERENCE}" >&2
    exit 1
fi

shopt -s nullglob

final_bam_files=("${FINAL_BAM_DIR}"/*.final.bam)

if [ ${#final_bam_files[@]} -eq 0 ]; then
    echo "ERROR: No final BAM files found in ${FINAL_BAM_DIR}" >&2
    exit 1
fi

for input_bam in "${final_bam_files[@]}"; do
    sample_name="$(basename "${input_bam}" .final.bam)"

    echo "Processing sample: ${sample_name}"

    output_gvcf="${GVCF_DIR}/${sample_name}.g.vcf.gz"
    bamout="${GVCF_DIR}/${sample_name}.bamout.bam"

    gatk --java-options "-Xmx${GATK_JAVA_MEM}" HaplotypeCaller \
        -R "${REFERENCE}" \
        -I "${input_bam}" \
        -ERC GVCF \
        -O "${output_gvcf}" \
        -bamout "${bamout}" \
        > "${LOG_DIR}/${sample_name}.haplotypecaller.log" 2>&1

    echo "Completed GVCF calling: ${sample_name}"
done

echo "[Step 05] GVCF calling completed."
