# ==============================================================================
# 03_internal_validation.R
# Project: MRONJ Risk Prediction
# Description:
#   Internal test-set validation of MRONJ risk prediction models.
# ==============================================================================

source("scripts/R/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(pROC)
})

message("Step 03: Internal validation")

model_file <- file.path(model_dir, "integrated_model_training_output.rds")

if (!file.exists(model_file)) {
  stop("Model file not found. Please run 02_train_integrated_model.R first.")
}

obj <- readRDS(model_file)

test_data_final <- obj$test_data_final

# ------------------------------------------------------------------------------
# Helper function
# ------------------------------------------------------------------------------

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
# Prediction
# ------------------------------------------------------------------------------

prob_clinical <- predict(
  obj$fit_enet_clinical,
  newdata = test_data_final,
  type = "prob"
)[, "Case"]

prob_prs <- predict(
  obj$fit_prs_only,
  newdata = test_data_final,
  type = "prob"
)[, "Case"]

prob_integrated <- predict(
  obj$fit_enet_integrated,
  newdata = test_data_final,
  type = "prob"
)[, "Case"]

prob_glm <- predict(
  obj$fit_glm,
  newdata = test_data_final,
  type = "prob"
)[, "Case"]

prob_rf <- predict(
  obj$fit_rf,
  newdata = test_data_final,
  type = "prob"
)[, "Case"]

prob_xgb <- predict(
  obj$fit_xgb,
  newdata = test_data_final,
  type = "prob"
)[, "Case"]

# ------------------------------------------------------------------------------
# Metrics
# ------------------------------------------------------------------------------

internal_metrics <- bind_rows(
  calculate_metrics(test_data_final$ONJ_Fac, prob_clinical, "Clinical only", fixed_cutoff),
  calculate_metrics(test_data_final$ONJ_Fac, prob_prs, "PRS only", fixed_cutoff),
  calculate_metrics(test_data_final$ONJ_Fac, prob_integrated, "Integrated", fixed_cutoff),
  calculate_metrics(test_data_final$ONJ_Fac, prob_glm, "GLM integrated", fixed_cutoff),
  calculate_metrics(test_data_final$ONJ_Fac, prob_rf, "Random forest integrated", fixed_cutoff),
  calculate_metrics(test_data_final$ONJ_Fac, prob_xgb, "XGBoost integrated", fixed_cutoff)
)

internal_predictions <- data.frame(
  Sample_ID = rownames(test_data_final),
  Outcome = test_data_final$ONJ_Fac,
  Prob_Clinical = prob_clinical,
  Prob_PRS = prob_prs,
  Prob_Integrated = prob_integrated,
  Prob_GLM = prob_glm,
  Prob_RF = prob_rf,
  Prob_XGB = prob_xgb,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

write.csv(
  internal_metrics,
  file.path(result_dir, "internal_validation_metrics.csv"),
  row.names = FALSE
)

write.csv(
  internal_predictions,
  file.path(result_dir, "internal_validation_predictions.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    internal_metrics = internal_metrics,
    internal_predictions = internal_predictions
  ),
  file.path(result_dir, "03_internal_validation_output.rds")
)

print(internal_metrics)

message("Step 03 completed: internal validation results saved.")
