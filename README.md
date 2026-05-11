# mronj-risk-prediction

Code for the development and validation of an integrated medication-related osteonecrosis of the jaw (MRONJ) risk prediction model based on whole-genome sequencing (WGS) data and clinical variables.

## 🧬 Overview

This repository contains Linux shell scripts and R scripts used for developing an integrated MRONJ risk prediction model.

The WGS preprocessing pipeline includes adapter trimming, read alignment, duplicate marking, base quality score recalibration, per-sample GVCF generation, joint genotyping, variant selection, and functional annotation.

The downstream R workflow includes clinical data preprocessing, WGS-derived variant matrix generation, genomic risk score estimation, integrated model training, internal validation, external validation, and visualization.

## 📁 Repository structure

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
├── README.md
└── .gitignore
```

## ⚙️ WGS preprocessing and variant calling

Before running the WGS pipeline, modify the paths and computational parameters in:

```bash
scripts/linux/00_config.sh
```

Run the WGS preprocessing and variant calling scripts in the following order:

```bash
bash scripts/linux/01_trim_adapters.sh
bash scripts/linux/02_align_bwa_mem.sh
bash scripts/linux/03_mark_duplicates_gatk.sh
bash scripts/linux/04_base_recalibration_bqsr.sh
bash scripts/linux/05_call_gvcf_haplotypecaller.sh
bash scripts/linux/06_create_gvcf.sh
bash scripts/linux/07_genomicsdb.sh
bash scripts/linux/08_genotypegvcfs.sh
bash scripts/linux/09_select_variants.sh
bash scripts/linux/10_annotation.sh
```

## 🧪 WGS pipeline description

### Step 01. Adapter trimming

Paired-end FASTQ files are processed using `cutadapt` to remove sequencing adapter sequences.

### Step 02. Read alignment

Trimmed reads are aligned to the GRCh38 reference genome using `BWA-MEM`. The resulting alignments are converted to BAM format, sorted, and indexed using `SAMtools`.

### Step 03. Duplicate marking

Duplicate reads are marked using `GATK MarkDuplicatesSpark`.

### Step 04. Base quality score recalibration

Base quality score recalibration is performed using `GATK BaseRecalibrator` and `GATK ApplyBQSR` with known variant resources.

### Step 05. Per-sample GVCF generation

Per-sample GVCF files are generated using `GATK HaplotypeCaller` in GVCF mode.

### Step 06. GVCF sample map generation

A sample map file is generated for `GATK GenomicsDBImport`.

### Step 07. GenomicsDB import

Per-sample GVCF files are imported into a GenomicsDB workspace using `GATK GenomicsDBImport`.

### Step 08. Joint genotyping

Joint genotyping is performed using `GATK GenotypeGVCFs`.

### Step 09. Variant selection

Variants passing predefined filtering criteria are selected using `GATK SelectVariants`.

### Step 10. Functional annotation

Selected variants are functionally annotated using `GATK Funcotator` and exported in MAF format.

## 📊 Integrated MRONJ risk prediction model

Before running the R workflow, modify the input paths and analysis parameters in:

```bash
scripts/R/00_config.R
```

Run the R scripts in the following order:

```bash
Rscript scripts/R/01_prepare_modeling_dataset.R
Rscript scripts/R/02_train_integrated_model.R
Rscript scripts/R/03_internal_validation.R
Rscript scripts/R/04_external_validation.R
Rscript scripts/R/05_generate_publication_figures.R
```

## 🔬 R workflow description

### Step 01. Modeling dataset preparation

This step loads WGS-derived MAF data and clinical data, generates a binary candidate variant matrix, harmonizes clinical variables, defines MRONJ case/control labels, and splits the dataset into training and internal test sets.

Candidate variants should be predefined in:

```text
data/candidate_variants_template.csv
```

Clinical variables should be predefined in:

```text
data/clinical_variables_template.csv
```

### Step 02. Integrated model training

This step performs preprocessing using training data only, derives a WGS-based genomic risk score from predefined candidate variants, standardizes the score using training-set statistics, and trains prediction models.

The following model types are evaluated:

```text
Clinical only model
PRS only model
Integrated model: clinical variables + WGS-derived genomic risk score
```

Additional benchmark algorithms may include:

```text
Logistic regression
Elastic Net
Random Forest
XGBoost
```

### Step 03. Internal validation

This step evaluates model performance in the held-out internal test set.

The main evaluation metrics include:

```text
AUC
95% confidence interval
Sensitivity
Specificity
Accuracy
```

### Step 04. External validation

This step applies the trained model to an independent external validation cohort.

The external validation workflow uses the same candidate variants, clinical variables, model coefficients, imputation parameters, and PRS scaling parameters derived from the training set.

No post-hoc sample exclusion based on model-predicted probabilities should be performed in the main validation analysis.

### Step 05. Publication-ready figure generation

This step generates publication-ready visualizations, including:

```text
Internal test set ROC curve
External validation ROC curve
Risk score boxplot
Risk score density plot
Variant frequency summary
```

## 📦 Software requirements

The WGS preprocessing pipeline requires:

```text
cutadapt
BWA
SAMtools
GATK
```

The R analysis workflow requires:

```text
R
maftools
dplyr
tidyr
reshape2
caret
glmnet
pROC
ggplot2
ggpubr
randomForest
missForest
```

## 🗂️ Data availability

Raw sequencing data, FASTQ files, BAM files, GVCF files, VCF files, MAF files containing individual-level variants, and clinical datasets are not included in this repository.

This repository provides scripts and template files required to reproduce the WGS preprocessing, variant calling, variant annotation, and downstream MRONJ risk prediction workflow.

## ✅ Notes on reproducibility

All user-specific file paths and computational parameters should be configured in:

```text
scripts/linux/00_config.sh
scripts/R/00_config.R
```

Candidate variants and clinical variables should be predefined before model training.

Feature selection, imputation, PRS scaling, and model optimization should be performed using training data only to avoid information leakage.

The internal test set and external validation cohort should be used only for model evaluation.

## 📝 Citation

If you use this code, please cite the associated manuscript.
