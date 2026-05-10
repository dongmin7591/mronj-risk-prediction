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


# -----------------------------
# Joint genotyping directories
# -----------------------------
GENOMICSDB_DIR="${PROJECT_DIR}/results/gvcf/genomicsdb"
JOINT_VCF_DIR="${PROJECT_DIR}/results/vcf"
SELECTED_VCF_DIR="${PROJECT_DIR}/results/vcf/selected_variants"
MAF_DIR="${PROJECT_DIR}/results/maf"
TMP_DIR="${PROJECT_DIR}/tmp"

# -----------------------------
# Joint genotyping files
# -----------------------------
SAMPLE_MAP="${GVCF_DIR}/gvcf_samples.list"
INTERVAL_LIST="${REF_DIR}/wgs_calling_regions.hg38.interval_list"

RAW_JOINT_VCF="${JOINT_VCF_DIR}/final_variants.vcf.gz"
SELECTED_PASS_VCF="${SELECTED_VCF_DIR}/selected_pass_variants.vcf.gz"

# -----------------------------
# Funcotator
# -----------------------------
FUNCOTATOR_DATASOURCE="${REF_DIR}/FuncotatorDataSource/germline_v2/test/funcotator_dataSources.v1.8.hg38.20230908g"
FUNCOTATOR_REF_VERSION="hg38"

# -----------------------------
# Computational resources for joint genotyping
# -----------------------------
GENOMICSDB_READER_THREADS=64
GENOMICSDB_BATCH_SIZE=50
GENOMICSDB_PARALLEL_INTERVALS=10

mkdir -p \
    "${GENOMICSDB_DIR}" \
    "${JOINT_VCF_DIR}" \
    "${SELECTED_VCF_DIR}" \
    "${MAF_DIR}" \
    "${TMP_DIR}"
