# ==============================================================================
# 04_external_validation.R
# Project: MRONJ Risk Prediction
# Description:
#   External validation of the integrated MRONJ risk prediction model.
#   No post-hoc sample exclusion based on model-predicted probabilities is applied.
# ==============================================================================

source("scripts/R/00_config.R")

suppressPackageStartupMessages({
  library(maftools)
  library(dplyr)
  library(reshape2)
  library(caret)
  library(pROC)
})

message("Step 04: External validation")

model_file <- file.path(model_dir, "integrated_model_training_output.rds")

if (!file.exists(model_file)) {
  stop("Model file not found. Please run 02_train_integrated_model.R first.")
}

obj <- readRDS(model_file)

if (!file.exists(external_maf_rds_file)) {
  stop("External MAF RDS file not found: ", external_maf_rds_file)
}

if (!file.exists(external_clinical_file)) {
  stop("External clinical file not found: ", external_clinical_file)
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

apply_simple_imputer <- function(df, imputer) {
  for (feature in names(imputer)) {
    if (!feature %in% colnames(df)) {
      df[[feature]] <- NA
    }
    
    fill_value <- imputer[[feature]]$value
    df[[feature]][is.na(df[[feature]])] <- fill_value
    df[[feature]] <- as.numeric(df[[feature]])
  }
  
  return(df)
}

make_external_variant_matrix <- function(maf_object, target_variants, rs_gene_mapping_file = NULL) {
  maf_data <- rbind(maf_object@data, maf_object@maf.silent)
  
  if (!"dbSNP_RS" %in% colnames(maf_data)) {
    maf_data$dbSNP_RS <- NA
  }
  
  if (!is.null(rs_gene_mapping_file) && file.exists(rs_gene_mapping_file)) {
    mapping_df <- read.csv(rs_gene_mapping_file, stringsAsFactors = FALSE)
    
    if (all(c("dbSNP_RS", "Hugo_Symbol") %in% colnames(mapping_df))) {
      mapping_df$dbSNP_RS <- gsub("^rs", "", as.character(mapping_df$dbSNP_RS))
      
      maf_data <- maf_data %>%
        select(-any_of("dbSNP_RS")) %>%
        left_join(mapping_df, by = "Hugo_Symbol", relationship = "many-to-one")
    }
  }
  
  maf_data$dbSNP_RS <- gsub("^rs", "", as.character(maf_data$dbSNP_RS))
  target_variants_raw <- gsub("^X", "", target_variants)
  
  subset_mut <- maf_data %>%
    filter(dbSNP_RS %in% target_variants_raw) %>%
    select(Tumor_Sample_Barcode, dbSNP_RS) %>%
    distinct() %>%
    mutate(Variant = 1)
  
  sample_ids <- unique(maf_object@clinical.data$Tumor_Sample_Barcode)
  
  if (nrow(subset_mut) == 0) {
    mut_matrix <- data.frame(Tumor_Sample_Barcode = sample_ids)
  } else {
    mut_matrix <- as.data.frame(
      reshape2::dcast(
        subset_mut,
        Tumor_Sample_Barcode ~ dbSNP_RS,
        value.var = "Variant",
        fill = 0
      )
    )
  }
  
  for (variant in target_variants_raw) {
    if (!variant %in% colnames(mut_matrix)) {
      mut_matrix[[variant]] <- 0
    }
  }
  
  variant_cols <- setdiff(colnames(mut_matrix), "Tumor_Sample_Barcode")
  
  mut_matrix[, variant_cols] <- lapply(mut_matrix[, variant_cols, drop = FALSE], function(x) {
    ifelse(as.numeric(x) > 0, 1, 0)
  })
  
  missing_samples <- setdiff(sample_ids, mut_matrix$Tumor_Sample_Barcode)
  
  if (length(missing_samples) > 0) {
    missing_df <- data.frame(Tumor_Sample_Barcode = missing_samples)
    
    for (variant in target_variants_raw) {
      missing_df[[variant]] <- 0
    }
    
    mut_matrix <- bind_rows(mut_matrix, missing_df)
  }
  
  colnames(mut_matrix) <- make.names(colnames(mut_matrix))
  
  return(mut_matrix)
}

calculate_metrics <- function(observed, predicted_prob, model_name, cutoff = 0.50) {
  roc_obj <- roc(
    observed,
    predicted_prob,
    levels = c("Control", "Case"),
    direction = "<",
    quiet = TRUE
  )
  
  ci_obj <- ci(roc_obj)
  
  pred_class <- factor(
    ifelse(predicted_prob >= cutoff, "Case", "Control"),
    levels = c("Case", "Control")
  )
  
  cm <- confusionMatrix(pred_class, observed, positive = "Case")
  
  metrics <- data.frame(
    Model = model_name,
    AUC = as.numeric(auc(roc_obj)),
    AUC_CI_Lower = as.numeric(ci_obj[1]),
    AUC_CI_Upper = as.numeric(ci_obj[3]),
    Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    Specificity = as.numeric(cm$byClass["Specificity"]),
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Cutoff = cutoff,
    stringsAsFactors = FALSE
  )
  
  return(metrics)
}

# ------------------------------------------------------------------------------
# Load external data
# ------------------------------------------------------------------------------

external_maf <- readRDS(external_maf_rds_file)
external_clinical <- read.csv(external_clinical_file, stringsAsFactors = FALSE)

external_clinical <- external_clinical %>%
  mutate(
    ONJ_Fac = factor(
      ifelse(Group == case_label_raw, "Case", "Control"),
      levels = c("Case", "Control")
    )
  )

if ("Sex" %in% colnames(external_clinical)) {
  external_clinical <- external_clinical %>%
    mutate(
      Sex_Code = case_when(
        Sex == "M" ~ 1,
        Sex == "F" ~ 2,
        TRUE ~ as.numeric(NA)
      )
    )
}

external_variant_matrix <- make_external_variant_matrix(
  maf_object = external_maf,
  target_variants = obj$final_genomic_features,
  rs_gene_mapping_file = rs_gene_mapping_file
)

external_clinical_subset <- external_clinical %>%
  select(Tumor_Sample_Barcode, ONJ_Fac, any_of(obj$clinical_vars))

external_data <- external_clinical_subset %>%
  left_join(external_variant_matrix, by = "Tumor_Sample_Barcode")

external_data <- as.data.frame(external_data)
colnames(external_data) <- make.names(colnames(external_data))

# ------------------------------------------------------------------------------
# Apply training-set preprocessing parameters
# ------------------------------------------------------------------------------

external_data <- apply_simple_imputer(external_data, obj$imputer)

for (feature in obj$final_genomic_features) {
  if (!feature %in% colnames(external_data)) {
    external_data[[feature]] <- 0
  }
}

external_geno <- as.matrix(external_data[, obj$final_genomic_features, drop = FALSE])

external_data$PRS <- as.numeric(external_geno %*% obj$final_weights)
external_data$PRS_Std <- (external_data$PRS - obj$prs_mean) / obj$prs_sd

external_data_final <- external_data[, c(
  "Tumor_Sample_Barcode",
  "ONJ_Fac",
  obj$clinical_vars,
  "PRS_Std"
), drop = FALSE]

# ------------------------------------------------------------------------------
# External prediction without post-hoc filtering
# ------------------------------------------------------------------------------

prob_clinical <- predict(
  obj$fit_enet_clinical,
  newdata = external_data_final,
  type = "prob"
)[, "Case"]

prob_prs <- predict(
  obj$fit_prs_only,
  newdata = external_data_final,
  type = "prob"
)[, "Case"]

prob_integrated <- predict(
  obj$fit_enet_integrated,
  newdata = external_data_final,
  type = "prob"
)[, "Case"]

external_metrics <- bind_rows(
  calculate_metrics(external_data_final$ONJ_Fac, prob_clinical, "Clinical only", fixed_cutoff),
  calculate_metrics(external_data_final$ONJ_Fac, prob_prs, "PRS only", fixed_cutoff),
  calculate_metrics(external_data_final$ONJ_Fac, prob_integrated, "Integrated", fixed_cutoff)
)

external_predictions <- data.frame(
  Tumor_Sample_Barcode = external_data_final$Tumor_Sample_Barcode,
  Outcome = external_data_final$ONJ_Fac,
  Prob_Clinical = prob_clinical,
  Prob_PRS = prob_prs,
  Prob_Integrated = prob_integrated,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

write.csv(
  external_metrics,
  file.path(result_dir, "external_validation_metrics.csv"),
  row.names = FALSE
)

write.csv(
  external_predictions,
  file.path(result_dir, "external_validation_predictions.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    external_metrics = external_metrics,
    external_predictions = external_predictions,
    external_data_final = external_data_final
  ),
  file.path(result_dir, "04_external_validation_output.rds")
)

print(external_metrics)

message("Step 04 completed: external validation results saved without post-hoc probability-based sample filtering.")
