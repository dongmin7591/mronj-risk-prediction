#!/usr/bin/env bash

# =========================================================
# Step 07. Import per-sample GVCFs into GenomicsDB
# Project: Integrated MRONJ risk prediction model
# Input:
#   GVCF sample map generated in Step 06
# Output:
#   GenomicsDB workspace
# =========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

echo "[Step 07] GenomicsDBImport"
echo "Sample map        : ${SAMPLE_MAP}"
echo "Interval list     : ${INTERVAL_LIST}"
echo "GenomicsDB output : ${GENOMICSDB_DIR}"
echo "Temporary directory: ${TMP_DIR}"

if [ ! -f "${SAMPLE_MAP}" ]; then
    echo "ERROR: Sample map not found: ${SAMPLE_MAP}" >&2
    echo "Please run 06_create_gvcf_sample_map.sh first." >&2
    exit 1
fi

if [ ! -f "${INTERVAL_LIST}" ]; then
    echo "ERROR: Interval list not found: ${INTERVAL_LIST}" >&2
    exit 1
fi

if [ -d "${GENOMICSDB_DIR}" ]; then
    echo "WARNING: Existing GenomicsDB workspace found: ${GENOMICSDB_DIR}" >&2
    echo "Removing existing workspace before import." >&2
    rm -rf "${GENOMICSDB_DIR}"
fi

mkdir -p "${TMP_DIR}"

gatk GenomicsDBImport \
    --genomicsdb-workspace-path "${GENOMICSDB_DIR}" \
    --sample-name-map "${SAMPLE_MAP}" \
    -L "${INTERVAL_LIST}" \
    --reader-threads "${GENOMICSDB_READER_THREADS}" \
    --batch-size "${GENOMICSDB_BATCH_SIZE}" \
    --consolidate \
    --tmp-dir "${TMP_DIR}" \
    --max-num-intervals-to-import-in-parallel "${GENOMICSDB_PARALLEL_INTERVALS}" \
    > "${LOG_DIR}/genomicsdb_import.log" 2>&1

echo "[Step 07] GenomicsDBImport completed successfully."
