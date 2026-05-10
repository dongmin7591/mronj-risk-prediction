# mronj-risk-prediction

Code for the development and validation of an integrated medication-related osteonecrosis of the jaw (MRONJ) risk prediction model based on whole-genome sequencing (WGS) data and clinical variables.

## Overview

This repository contains Linux shell scripts and R scripts used for developing an integrated MRONJ risk prediction model.

The WGS preprocessing pipeline includes adapter trimming, read alignment, duplicate marking, base quality score recalibration, per-sample GVCF generation, joint genotyping, variant selection, and functional annotation.

The downstream R workflow includes clinical data preprocessing, WGS-derived variant matrix generation, genomic risk score estimation, integrated model training, internal validation, external validation, and publication-ready figure generation.

## Repository structure

```text
mronj-risk-prediction/
├── scripts/
│   ├── linux/
│   │   ├── 00_config.sh
│   │   ├── 01_trim_adapters.sh
│   │   ├── 02_align_bwa_mem.sh
│   │   ├── 03_mark_duplicates_gatk.sh
│   │   ├── 04_base_recalibration_bqsr.sh
│   │   ├── 05_call_gvcf_haplotypecaller.sh
│   │   ├── 06_create_gvcf.sh
│   │   ├── 07_genomicsdb.sh
│   │   ├── 08_genotypegvcfs.sh
│   │   ├── 09_select_variants.sh
│   │   └── 10_annotation.sh
│   └── R/
│       ├── 00_config.R
│       ├── 01_prepare_modeling_dataset.R
│       ├── 02_train_integrated_model.R
│       ├── 03_internal_validation.R
│       ├── 04_external_validation.R
│       └── 05_generate_publication_figures.R
├── data/
│   ├── README.md
│   ├── candidate_variants_template.csv
│   ├── clinical_variables_template.csv
│   └── rs_gene_mapping_template.csv
├── results/
├── figures/
├── README.md
└── .gitignore
