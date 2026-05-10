#!/usr/bin/env bash

# =========================================================
# Step 01. Adapter trimming using cutadapt
# Project: Integrated MRONJ risk prediction model
# Input:
#   Paired-end FASTQ files: *_R1.fastq.gz and *_R2.fastq.gz
# Output:
#   Trimmed FASTQ files: *_R1.trimmed.fastq.gz and *_R2.trimmed.fastq.gz
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 01] Adapter trimming using cutadapt"
echo "Input directory : ${RAW_FASTQ_DIR}"
echo "Output directory: ${TRIMMED_FASTQ_DIR}"

shopt -s nullglob

r1_files=("${RAW_FASTQ_DIR}"/*_R1.fastq.gz)

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "ERROR: No R1 FASTQ files found in ${RAW_FASTQ_DIR}" >&2
    exit 1
fi

for r1_file in "${r1_files[@]}"; do
    sample_name="$(basename "${r1_file}" _R1.fastq.gz)"
    r2_file="${RAW_FASTQ_DIR}/${sample_name}_R2.fastq.gz"

    if [ ! -f "${r2_file}" ]; then
        echo "WARNING: R2 file not found for sample ${sample_name}. Skipping." >&2
        continue
    fi

    echo "Processing sample: ${sample_name}"

    cutadapt \
        -a "${ADAPTER_R1}" \
        -A "${ADAPTER_R2}" \
        -j "${THREADS}" \
        -o "${TRIMMED_FASTQ_DIR}/${sample_name}_R1.trimmed.fastq.gz" \
        -p "${TRIMMED_FASTQ_DIR}/${sample_name}_R2.trimmed.fastq.gz" \
        "${r1_file}" \
        "${r2_file}" \
        > "${LOG_DIR}/${sample_name}.cutadapt.log" 2>&1

    echo "Completed trimming: ${sample_name}"
done

echo "[Step 01] Adapter trimming completed."
