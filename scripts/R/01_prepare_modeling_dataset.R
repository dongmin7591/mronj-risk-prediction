# ==============================================================================
# 01_prepare_modeling_dataset.R
# Project: MRONJ Risk Prediction
# Description:
#   Prepare WGS-derived variant matrix and clinical variables for integrated
#   MRONJ risk prediction modeling.
# ==============================================================================

source("scripts/R/00_config.R")

suppressPackageStartupMessages({
  library(maftools)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(caret)
})

message("Step 01: Preparing modeling dataset")

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

load_candidate_variants <- function(candidate_variant_file) {
  if (!file.exists(candidate_variant_file)) {
    stop("Candidate variant file not found: ", candidate_variant_file)
  }
  
  candidate_df <- read.csv(candidate_variant_file, stringsAsFactors = FALSE)
  
  if (!"dbSNP_RS" %in% colnames(candidate_df)) {
    stop("candidate_variants_template.csv must contain a column named 'dbSNP_RS'.")
  }
  
  candidate_rs <- unique(as.character(candidate_df$dbSNP_RS))
  candidate_rs <- gsub("^rs", "", candidate_rs)
  
  return(candidate_rs)
}

load_clinical_variables <- function(clinical_variable_file) {
  if (!file.exists(clinical_variable_file)) {
    stop("Clinical variable file not found: ", clinical_variable_file)
  }
  
  clinical_var_df <- read.csv(clinical_variable_file, stringsAsFactors = FALSE)
  
  if (!"Variable" %in% colnames(clinical_var_df)) {
    stop("clinical_variables_template.csv must contain a column named 'Variable'.")
  }
  
  clinical_vars <- unique(as.character(clinical_var_df$Variable))
  
  return(clinical_vars)
}

make_variant_matrix <- function(maf_object, candidate_rs) {
  maf_data <- rbind(maf_object@data, maf_object@maf.silent)
  
  if (!"dbSNP_RS" %in% colnames(maf_data)) {
    stop("The MAF object must contain a 'dbSNP_RS' column.")
  }
  
  maf_data$dbSNP_RS <- gsub("^rs", "", as.character(maf_data$dbSNP_RS))
  
  subset_mut <- maf_data %>%
    filter(dbSNP_RS %in% candidate_rs) %>%
    select(Tumor_Sample_Barcode, dbSNP_RS) %>%
    distinct() %>%
    mutate(Variant = 1)
  
  sample_ids <- unique(maf_object@clinical.data$Tumor_Sample_Barcode)
  
  if (nrow(subset_mut) == 0) {
    warning("No candidate variants were found in the MAF object.")
    
    mut_matrix <- data.frame(Tumor_Sample_Barcode = sample_ids)
    
    for (rs in candidate_rs) {
      mut_matrix[[rs]] <- 0
    }
    
    return(mut_matrix)
  }
  
  mut_matrix <- as.data.frame(
    reshape2::dcast(
      subset_mut,
      Tumor_Sample_Barcode ~ dbSNP_RS,
      value.var = "Variant",
      fill = 0
    )
  )
  
  missing_rs <- setdiff(candidate_rs, colnames(mut_matrix))
  
  for (rs in missing_rs) {
    mut_matrix[[rs]] <- 0
  }
  
  variant_cols <- setdiff(colnames(mut_matrix), "Tumor_Sample_Barcode")
  
  mut_matrix[, variant_cols] <- lapply(mut_matrix[, variant_cols, drop = FALSE], function(x) {
    ifelse(as.numeric(x) > 0, 1, 0)
  })
  
  # Ensure all samples are represented.
  missing_samples <- setdiff(sample_ids, mut_matrix$Tumor_Sample_Barcode)
  
  if (length(missing_samples) > 0) {
    missing_df <- data.frame(Tumor_Sample_Barcode = missing_samples)
    
    for (rs in candidate_rs) {
      missing_df[[rs]] <- 0
    }
    
    mut_matrix <- bind_rows(mut_matrix, missing_df)
  }
  
  return(mut_matrix)
}

clean_clinical_data <- function(clinical_df, clinical_vars) {
  clinical_df <- as.data.frame(clinical_df)
  
  if (!"Tumor_Sample_Barcode" %in% colnames(clinical_df)) {
    stop("Clinical data must contain 'Tumor_Sample_Barcode'.")
  }
  
  if (!"Group" %in% colnames(clinical_df)) {
    stop("Clinical data must contain 'Group'.")
  }
  
  clinical_df <- clinical_df %>%
    mutate(
      ONJ_Fac = factor(
        ifelse(Group == case_label_raw, "Case", "Control"),
        levels = c("Case", "Control")
      )
    )
  
  if ("Sex" %in% colnames(clinical_df)) {
    clinical_df <- clinical_df %>%
      mutate(
        Sex_Code = case_when(
          Sex == "M" ~ 1,
          Sex == "F" ~ 2,
          TRUE ~ as.numeric(NA)
        )
      )
  }
  
  existing_clin_vars <- intersect(clinical_vars, colnames(clinical_df))
  
  if (length(existing_clin_vars) == 0) {
    stop("None of the requested clinical variables were found in the clinical data.")
  }
  
  clinical_subset <- clinical_df %>%
    select(Tumor_Sample_Barcode, ONJ_Fac, all_of(existing_clin_vars))
  
  numeric_cols <- setdiff(colnames(clinical_subset), c("Tumor_Sample_Barcode", "ONJ_Fac"))
  
  clinical_subset <- clinical_subset %>%
    mutate(across(all_of(numeric_cols), ~ as.numeric(as.character(.)))) %>%
    mutate(across(all_of(numeric_cols), ~ ifelse(. %in% c(9999, 8888), NA, .)))
  
  return(clinical_subset)
}

# ------------------------------------------------------------------------------
# Load input files
# ------------------------------------------------------------------------------

candidate_rs <- load_candidate_variants(candidate_variant_file)
clinical_vars <- load_clinical_variables(clinical_variable_file)

if (!file.exists(wgs_maf_rds_file)) {
  stop("WGS MAF RDS file not found: ", wgs_maf_rds_file)
}

if (!file.exists(wgs_clinical_file)) {
  stop("WGS clinical input file not found: ", wgs_clinical_file)
}

wgs_maf <- readRDS(wgs_maf_rds_file)
wgs_clinical <- read.csv(wgs_clinical_file, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Generate modeling dataset
# ------------------------------------------------------------------------------

variant_matrix <- make_variant_matrix(wgs_maf, candidate_rs)
clinical_matrix <- clean_clinical_data(wgs_clinical, clinical_vars)

model_data <- clinical_matrix %>%
  left_join(variant_matrix, by = "Tumor_Sample_Barcode")

variant_cols <- intersect(candidate_rs, colnames(model_data))

if (length(variant_cols) > 0) {
  model_data[, variant_cols] <- lapply(model_data[, variant_cols, drop = FALSE], function(x) {
    x[is.na(x)] <- 0
    as.numeric(x)
  })
}

model_data <- as.data.frame(model_data)
colnames(model_data) <- make.names(colnames(model_data))

candidate_rs_safe <- make.names(candidate_rs)
clinical_vars_safe <- make.names(intersect(clinical_vars, colnames(clinical_matrix)))

candidate_rs_safe <- intersect(candidate_rs_safe, colnames(model_data))
clinical_vars_safe <- intersect(clinical_vars_safe, colnames(model_data))

# ------------------------------------------------------------------------------
# Train-test split
# ------------------------------------------------------------------------------

set.seed(random_seed)

train_idx <- createDataPartition(
  y = model_data$ONJ_Fac,
  p = train_fraction,
  list = FALSE
)

train_data_raw <- model_data[train_idx, ]
test_data_raw  <- model_data[-train_idx, ]

# ------------------------------------------------------------------------------
# Save objects
# ------------------------------------------------------------------------------

saveRDS(
  list(
    model_data = model_data,
    train_data_raw = train_data_raw,
    test_data_raw = test_data_raw,
    candidate_rs = candidate_rs_safe,
    clinical_vars = clinical_vars_safe
  ),
  file = file.path(result_dir, "01_prepared_modeling_dataset.rds")
)

write.csv(
  model_data,
  file.path(result_dir, "prepared_modeling_dataset.csv"),
  row.names = FALSE
)

message("Step 01 completed: prepared modeling dataset saved.")
