# ==============================================================================
# 02_train_integrated_model.R
# Project: MRONJ Risk Prediction
# Description:
#   Train an integrated MRONJ risk prediction model using clinical variables and
#   WGS-derived candidate variants.
# ==============================================================================

source("scripts/R/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(glmnet)
  library(pROC)
  library(randomForest)
  library(xgboost)
})

message("Step 02: Training integrated MRONJ risk prediction model")

prepared_file <- file.path(result_dir, "01_prepared_modeling_dataset.rds")

if (!file.exists(prepared_file)) {
  stop("Prepared dataset not found. Please run 01_prepare_modeling_dataset.R first.")
}

prepared <- readRDS(prepared_file)

train_data_raw <- prepared$train_data_raw
test_data_raw  <- prepared$test_data_raw
candidate_rs   <- prepared$candidate_rs
clinical_vars  <- prepared$clinical_vars

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

get_mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

fit_simple_imputer <- function(train_df, feature_cols) {
  imputer <- list()
  
  for (feature in feature_cols) {
    x <- train_df[[feature]]
    unique_values <- unique(x[!is.na(x)])
    is_binary <- all(unique_values %in% c(0, 1))
    
    if (is_binary) {
      imputer[[feature]] <- list(
        type = "mode",
        value = as.numeric(get_mode_value(x))
      )
    } else {
      imputer[[feature]] <- list(
        type = "median",
        value = median(as.numeric(x), na.rm = TRUE)
      )
    }
  }
  
  return(imputer)
}

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

filter_near_zero_variants <- function(train_df, variant_cols, cutoff = 0.90) {
  kept <- c()
  
  stats <- data.frame(
    Variant = character(),
    Major_Class_Frequency = numeric(),
    Status = character(),
    stringsAsFactors = FALSE
  )
  
  for (variant in variant_cols) {
    x <- train_df[[variant]]
    x <- x[!is.na(x)]
    
    if (length(x) == 0) next
    
    freq <- table(x)
    major_freq <- max(freq) / sum(freq)
    
    if (major_freq >= cutoff) {
      status <- "Removed_near_zero_variance"
    } else {
      status <- "Kept"
      kept <- c(kept, variant)
    }
    
    stats <- rbind(
      stats,
      data.frame(
        Variant = variant,
        Major_Class_Frequency = major_freq,
        Status = status,
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(list(kept = kept, stats = stats))
}

# ------------------------------------------------------------------------------
# Feature preparation
# ------------------------------------------------------------------------------

candidate_rs <- intersect(candidate_rs, colnames(train_data_raw))
clinical_vars <- intersect(clinical_vars, colnames(train_data_raw))

if (length(candidate_rs) == 0) {
  stop("No candidate variant columns were found in the training dataset.")
}

if (length(clinical_vars) == 0) {
  stop("No clinical variable columns were found in the training dataset.")
}

feature_cols <- c(clinical_vars, candidate_rs)

imputer <- fit_simple_imputer(train_data_raw, feature_cols)

train_data <- apply_simple_imputer(train_data_raw, imputer)
test_data  <- apply_simple_imputer(test_data_raw, imputer)

# ------------------------------------------------------------------------------
# Variant filtering using training data only
# ------------------------------------------------------------------------------

variant_filter <- filter_near_zero_variants(
  train_df = train_data,
  variant_cols = candidate_rs,
  cutoff = near_zero_variance_cutoff
)

selected_variants <- variant_filter$kept
variant_filter_stats <- variant_filter$stats

if (length(selected_variants) == 0) {
  stop("No candidate variants remained after near-zero variance filtering.")
}

write.csv(
  variant_filter_stats,
  file.path(result_dir, "variant_filtering_summary.csv"),
  row.names = FALSE
)

message("Selected variants after near-zero variance filtering:")
print(selected_variants)

# ------------------------------------------------------------------------------
# WGS-derived genomic risk score using Elastic Net
# ------------------------------------------------------------------------------

x_prs <- as.matrix(train_data[, selected_variants, drop = FALSE])
y_prs <- ifelse(train_data$ONJ_Fac == "Case", 1, 0)

set.seed(random_seed)

cv_prs <- cv.glmnet(
  x = x_prs,
  y = y_prs,
  family = "binomial",
  alpha = prs_alpha,
  type.measure = "auc",
  nfolds = prs_nfolds
)

coef_prs <- coef(cv_prs, s = cv_prs$lambda.min)

coef_df <- data.frame(
  Feature = rownames(coef_prs),
  Coefficient = as.numeric(coef_prs[, 1]),
  stringsAsFactors = FALSE
) %>%
  filter(Feature != "(Intercept)") %>%
  filter(Coefficient != 0)

if (nrow(coef_df) == 0) {
  warning("Elastic Net selected no non-zero variants. Using all selected variants with equal weights.")
  
  coef_df <- data.frame(
    Feature = selected_variants,
    Coefficient = rep(1 / length(selected_variants), length(selected_variants)),
    stringsAsFactors = FALSE
  )
}

final_genomic_features <- coef_df$Feature
final_weights <- setNames(coef_df$Coefficient, coef_df$Feature)

train_geno <- as.matrix(train_data[, final_genomic_features, drop = FALSE])
test_geno  <- as.matrix(test_data[, final_genomic_features, drop = FALSE])

train_data$PRS <- as.numeric(train_geno %*% final_weights)
test_data$PRS  <- as.numeric(test_geno %*% final_weights)

prs_mean <- mean(train_data$PRS, na.rm = TRUE)
prs_sd   <- sd(train_data$PRS, na.rm = TRUE)

if (is.na(prs_sd) || prs_sd == 0) {
  stop("PRS standard deviation is zero or NA in the training set.")
}

train_data$PRS_Std <- (train_data$PRS - prs_mean) / prs_sd
test_data$PRS_Std  <- (test_data$PRS - prs_mean) / prs_sd

coef_df <- coef_df %>%
  mutate(
    Odds_Ratio = exp(Coefficient),
    Effect = ifelse(Coefficient >= 0, "Risk-increasing", "Risk-decreasing")
  )

write.csv(
  coef_df,
  file.path(result_dir, "genomic_score_variant_weights.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Model training
# ------------------------------------------------------------------------------

vars_clin_only  <- clinical_vars
vars_prs_only   <- "PRS_Std"
vars_integrated <- c(clinical_vars, "PRS_Std")

train_data_final <- train_data[, c("ONJ_Fac", vars_clin_only, "PRS_Std"), drop = FALSE]
test_data_final  <- test_data[, c("ONJ_Fac", vars_clin_only, "PRS_Std"), drop = FALSE]

train_control <- trainControl(
  method = "repeatedcv",
  number = cv_folds,
  repeats = cv_repeats,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  sampling = "up",
  allowParallel = TRUE,
  savePredictions = "final"
)

# Clinical only model
set.seed(random_seed)

fit_enet_clinical <- train(
  x = train_data_final[, vars_clin_only, drop = FALSE],
  y = train_data_final$ONJ_Fac,
  method = "glmnet",
  metric = "ROC",
  trControl = train_control,
  tuneLength = 5
)

# PRS only model
set.seed(random_seed)

fit_prs_only <- train(
  x = train_data_final[, vars_prs_only, drop = FALSE],
  y = train_data_final$ONJ_Fac,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = train_control
)

# Integrated Elastic Net model
set.seed(random_seed)

fit_enet_integrated <- train(
  x = train_data_final[, vars_integrated, drop = FALSE],
  y = train_data_final$ONJ_Fac,
  method = "glmnet",
  metric = "ROC",
  trControl = train_control,
  tuneLength = 5
)

# Optional benchmark: GLM
set.seed(random_seed)

fit_glm <- train(
  x = train_data_final[, vars_integrated, drop = FALSE],
  y = train_data_final$ONJ_Fac,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = train_control
)

# Optional benchmark: Random Forest
set.seed(random_seed)

fit_rf <- train(
  x = train_data_final[, vars_integrated, drop = FALSE],
  y = train_data_final$ONJ_Fac,
  method = "rf",
  metric = "ROC",
  trControl = train_control,
  tuneLength = 3
)

# Optional benchmark: XGBoost
xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(2, 3),
  eta = c(0.05, 0.10),
  gamma = 0,
  colsample_bytree = 0.8,
  subsample = 0.8,
  min_child_weight = 3
)

set.seed(random_seed)

fit_xgb <- train(
  x = train_data_final[, vars_integrated, drop = FALSE],
  y = train_data_final$ONJ_Fac,
  method = "xgbTree",
  metric = "ROC",
  trControl = train_control,
  tuneGrid = xgb_grid,
  verbosity = 0
)

# ------------------------------------------------------------------------------
# Save training output
# ------------------------------------------------------------------------------

training_output <- list(
  train_data_final = train_data_final,
  test_data_final = test_data_final,
  train_data_full = train_data,
  test_data_full = test_data,
  imputer = imputer,
  selected_variants = selected_variants,
  final_genomic_features = final_genomic_features,
  final_weights = final_weights,
  prs_mean = prs_mean,
  prs_sd = prs_sd,
  clinical_vars = clinical_vars,
  vars_clin_only = vars_clin_only,
  vars_prs_only = vars_prs_only,
  vars_integrated = vars_integrated,
  fit_enet_clinical = fit_enet_clinical,
  fit_prs_only = fit_prs_only,
  fit_enet_integrated = fit_enet_integrated,
  fit_glm = fit_glm,
  fit_rf = fit_rf,
  fit_xgb = fit_xgb
)

saveRDS(
  training_output,
  file.path(model_dir, "integrated_model_training_output.rds")
)

message("Step 02 completed: integrated model training output saved.")
