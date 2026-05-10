# ==============================================================================
# 05_generate_publication_figures.R
# Project: MRONJ Risk Prediction
# Description:
#   Generate publication-ready figures for internal and external validation.
# ==============================================================================

source("scripts/R/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(scales)
})

message("Step 05: Generating publication figures")

internal_file <- file.path(result_dir, "03_internal_validation_output.rds")
external_file <- file.path(result_dir, "04_external_validation_output.rds")

if (!file.exists(internal_file)) {
  stop("Internal validation output not found. Please run 03_internal_validation.R first.")
}

if (!file.exists(external_file)) {
  stop("External validation output not found. Please run 04_external_validation.R first.")
}

internal_obj <- readRDS(internal_file)
external_obj <- readRDS(external_file)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

plot_roc_comparison <- function(pred_df, outcome_col, prob_cols, model_labels, plot_title, output_file) {
  roc_list <- list()
  
  for (i in seq_along(prob_cols)) {
    prob_col <- prob_cols[i]
    label <- model_labels[i]
    
    roc_obj <- roc(
      pred_df[[outcome_col]],
      pred_df[[prob_col]],
      levels = c("Control", "Case"),
      direction = "<",
      quiet = TRUE
    )
    
    ci_obj <- ci(roc_obj)
    
    roc_label <- paste0(
      label,
      " (AUC = ", sprintf("%.3f", auc(roc_obj)),
      " [", sprintf("%.2f", ci_obj[1]),
      "-", sprintf("%.2f", ci_obj[3]), "])"
    )
    
    roc_list[[roc_label]] <- smooth(roc_obj, method = "density")
  }
  
  p <- ggroc(roc_list, legacy.axes = TRUE, linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, color = "grey50", linetype = "dashed") +
    theme_bw(base_size = 14) +
    labs(
      title = plot_title,
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      legend.position = c(0.65, 0.20),
      legend.background = element_rect(
        fill = alpha("white", 0.9),
        color = "black",
        linewidth = 0.4
      ),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black")
    )
  
  ggsave(output_file, p, width = 7, height = 6, dpi = 300)
  
  return(p)
}

plot_risk_score_boxplot <- function(pred_df, prob_col, plot_title, output_file) {
  plot_df <- data.frame(
    Group = factor(pred_df$Outcome, levels = c("Control", "Case")),
    Risk_Score = pred_df[[prob_col]]
  )
  
  p <- ggplot(plot_df, aes(x = Group, y = Risk_Score, fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = Group), width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = c("Control" = "#4DBBD5FF", "Case" = "#E64B35FF")) +
    scale_color_manual(values = c("Control" = "#4DBBD5FF", "Case" = "#E64B35FF")) +
    theme_classic(base_size = 14) +
    labs(
      title = plot_title,
      x = "",
      y = "MRONJ risk score"
    ) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Control", "Case")),
      label = "p.format",
      tip.length = 0.02
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
    )
  
  ggsave(output_file, p, width = 4, height = 5, dpi = 300)
  
  return(p)
}

plot_risk_score_density <- function(pred_df, prob_col, plot_title, output_file) {
  plot_df <- data.frame(
    Group = factor(pred_df$Outcome, levels = c("Control", "Case")),
    Risk_Score = pred_df[[prob_col]]
  )
  
  p <- ggplot(plot_df, aes(x = Risk_Score, fill = Group)) +
    geom_density(alpha = 0.5, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Control" = "#4DBBD5FF", "Case" = "#E64B35FF")) +
    theme_classic(base_size = 14) +
    labs(
      title = plot_title,
      x = "MRONJ risk score",
      y = "Density",
      fill = "Actual outcome"
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
    )
  
  ggsave(output_file, p, width = 6, height = 4.5, dpi = 300)
  
  return(p)
}

# ------------------------------------------------------------------------------
# Internal validation figures
# ------------------------------------------------------------------------------

internal_predictions <- internal_obj$internal_predictions

p_internal_roc <- plot_roc_comparison(
  pred_df = internal_predictions,
  outcome_col = "Outcome",
  prob_cols = c("Prob_Clinical", "Prob_PRS", "Prob_Integrated"),
  model_labels = c("Clinical only", "PRS only", "Integrated"),
  plot_title = "Internal test set",
  output_file = file.path(figure_dir, "internal_test_set_roc.png")
)

p_internal_box <- plot_risk_score_boxplot(
  pred_df = internal_predictions,
  prob_col = "Prob_Integrated",
  plot_title = "Internal test set",
  output_file = file.path(figure_dir, "internal_integrated_risk_score_boxplot.png")
)

p_internal_density <- plot_risk_score_density(
  pred_df = internal_predictions,
  prob_col = "Prob_Integrated",
  plot_title = "Internal test set",
  output_file = file.path(figure_dir, "internal_integrated_risk_score_density.png")
)

# ------------------------------------------------------------------------------
# External validation figures
# ------------------------------------------------------------------------------

external_predictions <- external_obj$external_predictions

p_external_roc <- plot_roc_comparison(
  pred_df = external_predictions,
  outcome_col = "Outcome",
  prob_cols = c("Prob_Clinical", "Prob_PRS", "Prob_Integrated"),
  model_labels = c("Clinical only", "PRS only", "Integrated"),
  plot_title = "External validation cohort",
  output_file = file.path(figure_dir, "external_validation_roc.png")
)

p_external_box <- plot_risk_score_boxplot(
  pred_df = external_predictions,
  prob_col = "Prob_Integrated",
  plot_title = "External validation cohort",
  output_file = file.path(figure_dir, "external_integrated_risk_score_boxplot.png")
)

p_external_density <- plot_risk_score_density(
  pred_df = external_predictions,
  prob_col = "Prob_Integrated",
  plot_title = "External validation cohort",
  output_file = file.path(figure_dir, "external_integrated_risk_score_density.png")
)

message("Step 05 completed: publication figures saved.")
