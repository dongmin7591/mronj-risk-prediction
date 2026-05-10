# ==============================================================================
# 00_config.R
# Project: MRONJ Risk Prediction
# Description:
#   Configuration file for the integrated MRONJ risk prediction model.
# ==============================================================================

set.seed(2025)

# ------------------------------------------------------------------------------
# Directory configuration
# ------------------------------------------------------------------------------
project_dir <- "."

data_dir    <- file.path(project_dir, "data")
result_dir  <- file.path(project_dir, "results")
figure_dir  <- file.path(project_dir, "figures")
model_dir   <- file.path(result_dir, "models")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
wgs_maf_rds_file       <- file.path(data_dir, "WGS_Germline.maf.rds")
wgs_clinical_file      <- file.path(data_dir, "WGS_clinical_input.csv")

external_maf_rds_file  <- file.path(data_dir, "Sanger.maf.rds")
external_clinical_file <- file.path(data_dir, "Sanger_clinical_input.csv")

candidate_variant_file <- file.path(data_dir, "candidate_variants_template.csv")
clinical_variable_file <- file.path(data_dir, "clinical_variables_template.csv")
rs_gene_mapping_file   <- file.path(data_dir, "rs_gene_mapping_template.csv")

# ------------------------------------------------------------------------------
# Outcome and modeling parameters
# ------------------------------------------------------------------------------
case_label_raw <- "true 1"

train_fraction <- 0.70
random_seed <- 2025

# Remove variants showing extreme class imbalance in the training set.
near_zero_variance_cutoff <- 0.90

# Fixed prediction threshold for sensitivity/specificity reporting.
# This should be pre-specified or derived only from training data.
fixed_cutoff <- 0.50

# Elastic Net parameters for WGS-derived genomic risk score.
prs_alpha <- 0.10
prs_nfolds <- 5

# Cross-validation setting for model training.
cv_folds <- 5
cv_repeats <- 100
