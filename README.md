# mronj-risk-prediction

Code for the development and validation of an integrated MRONJ risk prediction model based on whole-genome sequencing (WGS) data and clinical variables.

## Overview

This repository contains Linux shell scripts and R scripts used for developing an integrated risk prediction model for medication-related osteonecrosis of the jaw (MRONJ).

The WGS preprocessing pipeline includes adapter trimming, read alignment, duplicate marking, base quality score recalibration, and per-sample GVCF generation.

## Repository structure

```text
mronj-risk-prediction/
├── scripts/
│   ├── linux/
│   │   ├── 00_config.sh
│   │   ├── 01_trim_adapters_cutadapt.sh
│   │   ├── 02_align_bwa_mem.sh
│   │   ├── 03_mark_duplicates_gatk.sh
│   │   ├── 04_base_recalibration_bqsr.sh
│   │   └── 05_call_gvcf_haplotypecaller.sh
│   └── R/
├── data/
├── results/
├── figures/
└── README.md
