#!/usr/bin/env bash

# =========================================================
# Configuration file for WGS preprocessing pipeline
# Project: Integrated MRONJ risk prediction model
# Description:
#   This file defines paths and computational parameters used
#   in the WGS preprocessing pipeline.
# =========================================================

set -euo pipefail

# -----------------------------
# Project directories
# -----------------------------
PROJECT_DIR="/path/to/project"

RAW_FASTQ_DIR="${PROJECT_DIR}/data/raw_fastq"
TRIMMED_FASTQ_DIR="${PROJECT_DIR}/data/trimmed_fastq"
BAM_DIR="${PROJECT_DIR}/results/bam"
DEDUP_BAM_DIR="${PROJECT_DIR}/results/bam/markduplicates"
BQSR_DIR="${PROJECT_DIR}/results/bam/bqsr"
FINAL_BAM_DIR="${PROJECT_DIR}/results/bam/final"
GVCF_DIR="${PROJECT_DIR}/results/gvcf"
LOG_DIR="${PROJECT_DIR}/logs"

# -----------------------------
# Reference files
# -----------------------------
REF_DIR="/path/to/reference/gatk_bundle"
REFERENCE="${REF_DIR}/Homo_sapiens_assembly38.fasta"

DBSNP="${REF_DIR}/dbsnp_146.hg38.vcf.gz"
HAPMAP="${REF_DIR}/hapmap_3.3_grch38_pop_stratified_af.vcf.gz"
MILLS="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
OMNI="${REF_DIR}/1000G_omni2.5.hg38.vcf.gz"

# -----------------------------
# Computational resources
# -----------------------------
THREADS=20
GATK_JAVA_MEM="4g"

# -----------------------------
# Adapter sequences for cutadapt
# -----------------------------
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# -----------------------------
# Create output directories
# -----------------------------
mkdir -p \
    "${TRIMMED_FASTQ_DIR}" \
    "${BAM_DIR}" \
    "${DEDUP_BAM_DIR}" \
    "${BQSR_DIR}" \
    "${FINAL_BAM_DIR}" \
    "${GVCF_DIR}" \
    "${LOG_DIR}"
