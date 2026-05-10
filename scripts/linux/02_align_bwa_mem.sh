#!/usr/bin/env bash

# =========================================================
# Step 02. Alignment to reference genome using BWA-MEM
# Project: Integrated MRONJ risk prediction model
# Input:
#   Trimmed paired-end FASTQ files
# Output:
#   Sorted BAM files and BAM index files
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 02] Alignment using BWA-MEM"
echo "Input directory : ${TRIMMED_FASTQ_DIR}"
echo "Output directory: ${BAM_DIR}"
echo "Reference       : ${REFERENCE}"

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference genome not found: ${REFERENCE}" >&2
    exit 1
fi

shopt -s nullglob

r1_files=("${TRIMMED_FASTQ_DIR}"/*_R1.trimmed.fastq.gz)

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "ERROR: No trimmed R1 FASTQ files found in ${TRIMMED_FASTQ_DIR}" >&2
    exit 1
fi

for r1_file in "${r1_files[@]}"; do
    sample_name="$(basename "${r1_file}" _R1.trimmed.fastq.gz)"
    r2_file="${TRIMMED_FASTQ_DIR}/${sample_name}_R2.trimmed.fastq.gz"

    if [ ! -f "${r2_file}" ]; then
        echo "WARNING: R2 file not found for sample ${sample_name}. Skipping." >&2
        continue
    fi

    echo "Processing sample: ${sample_name}"

    output_bam="${BAM_DIR}/${sample_name}.sorted.bam"

    bwa mem \
        -M \
        -t "${THREADS}" \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${sample_name}\tPL:ILLUMINA" \
        "${REFERENCE}" \
        "${r1_file}" \
        "${r2_file}" \
        2> "${LOG_DIR}/${sample_name}.bwa_mem.log" \
    | samtools view -@ "${THREADS}" -b - \
    | samtools sort -@ "${THREADS}" -O BAM -o "${output_bam}" -

    samtools index -@ "${THREADS}" "${output_bam}"

    echo "Completed alignment: ${sample_name}"
done

echo "[Step 02] Alignment completed."
