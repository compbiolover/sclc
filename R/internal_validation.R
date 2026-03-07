# Name: internal_validation.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Comprehensive internal validation suite for a gene signature
# when external cohorts are unavailable. Produces multiple complementary
# analyses that each test a different aspect of the model's validity.

suppressMessages({
  library(glmnet)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(tidyverse)
})


#' Run a comprehensive internal validation suite for a gene signature
#'
#' Produces: Kaplan-Meier curves with log-rank test, calibration analysis,
#' permutation testing, coefficient stability via bootstrap, and
#' time-dependent AUC (if timeROC is available).
#'
#' @param cox_df Data frame with patient.id, vital.status, time, and gene
#'   expression columns.
#' @param signature_genes Character vector of gene names in the signature.
#' @param signature_coefs Named numeric vector of Cox coefficients for each
#'   gene. If NULL, coefficients are estimated from the data via ridge
#'   regression (alpha = 0).
#' @param n_permutations Number of permutations for the null distribution
#'   test (default 1000).
#' @param n_boot Number of bootstrap resamples for coefficient stability
#'   (default 1000).
#' @param seed Master seed for reproducibility.
#' @param progress_free Logical; if TRUE, use time2 column.
#' @param output_dir Directory to save figures and results.
#' @param verbose Logical; print progress.
#'
#' @return A list containing results from each validation analysis.
internal_validation <- function(cox_df,
                                signature_genes,
                                signature_coefs = NULL,
                                n_permutations = 1000,
                                n_boot = 1000,
                                seed = 42,
                                progress_free = FALSE,
                                output_dir = "Outputs/internal_validation",
                                verbose = TRUE) {
  set.seed(seed)

  # Input validation
  if (is.null(cox_df) || !is.data.frame(cox_df)) {
    stop("cox_df must be a data.frame.")
  }

  signature_genes <- tolower(signature_genes)
  available_genes <- intersect(signature_genes, colnames(cox_df))
  if (length(available_genes) == 0) {
    stop("No signature genes found in cox_df columns.")
  }
  if (length(available_genes) < length(signature_genes) && verbose) {
    missing <- setdiff(signature_genes, available_genes)
    message("Missing genes: ", paste(missing, collapse = ", "))
  }
  signature_genes <- available_genes

  time_col <- if (progress_free) "time2" else "time"

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  results <- list()

  # ===================================================================
  # 1. FIT THE MODEL (if coefficients not provided)
  # ===================================================================
  if (verbose) message("\n[1/6] Fitting Cox model...")

  formula_str <- paste0("~", paste0(make.names(signature_genes), collapse = "+"))
  df_clean <- cox_df
  colnames(df_clean) <- make.names(colnames(df_clean))

  model_x <- model.matrix(as.formula(formula_str), df_clean)
  model_y <- Surv(
    time = df_clean[[make.names(time_col)]],
    event = df_clean$vital.status
  )

  if (is.null(signature_coefs)) {
    # Fit ridge to get coefficients (no feature selection, just shrinkage)
    cv_ridge <- cv.glmnet(model_x, model_y, family = "cox", alpha = 0,
                          type.measure = "C", nfolds = 10)
    signature_coefs <- as.vector(coef(cv_ridge, s = "lambda.1se"))
    names(signature_coefs) <- colnames(model_x)
  }

  # Compute risk scores
  risk_scores <- as.vector(model_x %*% signature_coefs[colnames(model_x)])
  df_clean$risk_score <- risk_scores

  results$coefficients <- signature_coefs

  # ===================================================================
  # 2. KAPLAN-MEIER ANALYSIS (median split + tertile split)
  # ===================================================================
  if (verbose) message("[2/6] Kaplan-Meier survival analysis...")

  # --- Median split ---
  median_risk <- median(risk_scores)
  df_clean$risk_group <- ifelse(risk_scores > median_risk, "High", "Low")
  df_clean$risk_group <- factor(df_clean$risk_group, levels = c("Low", "High"))

  km_fit <- survfit(
    Surv(df_clean[[make.names(time_col)]], df_clean$vital.status) ~ risk_group,
    data = df_clean
  )
  logrank <- survdiff(
    Surv(df_clean[[make.names(time_col)]], df_clean$vital.status) ~ risk_group,
    data = df_clean
  )
  lr_pval <- 1 - pchisq(logrank$chisq, df = 1)

  km_plot <- ggsurvplot(
    fit = km_fit,
    data = df_clean,
    conf.int = TRUE,
    pval = TRUE,
    pval.size = 5,
    risk.table = TRUE,
    palette = c("blue", "red"),
    censor.shape = "|",
    surv.median.line = "hv",
    size = 1.5,
    censor.size = 4.5,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    legend.title = "Risk Group",
    title = "Kaplan-Meier: Median Risk Score Split",
    legend.labs = c("Low risk", "High risk"),
    ggtheme = theme(
      axis.line = element_line(color = "black", linewidth = 1.25),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.background = element_blank()
    )
  )

  ggsave(km_plot$plot,
         filename = "km_median_split.svg",
         path = output_dir, width = 8, height = 7)

  # --- Tertile split ---
  tertile_cuts <- quantile(risk_scores, probs = c(1/3, 2/3))
  df_clean$risk_tertile <- cut(
    risk_scores,
    breaks = c(-Inf, tertile_cuts[1], tertile_cuts[2], Inf),
    labels = c("Low", "Medium", "High")
  )

  km_fit_tert <- survfit(
    Surv(df_clean[[make.names(time_col)]], df_clean$vital.status) ~ risk_tertile,
    data = df_clean
  )
  logrank_tert <- survdiff(
    Surv(df_clean[[make.names(time_col)]], df_clean$vital.status) ~ risk_tertile,
    data = df_clean
  )
  lr_pval_tert <- 1 - pchisq(logrank_tert$chisq, df = 2)

  km_plot_tert <- ggsurvplot(
    fit = km_fit_tert,
    data = df_clean,
    conf.int = TRUE,
    pval = TRUE,
    pval.size = 5,
    risk.table = TRUE,
    palette = c("blue", "darkgreen", "red"),
    censor.shape = "|",
    size = 1.5,
    censor.size = 4.5,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    legend.title = "Risk Group",
    title = "Kaplan-Meier: Tertile Risk Score Split",
    legend.labs = c("Low risk", "Medium risk", "High risk"),
    ggtheme = theme(
      axis.line = element_line(color = "black", linewidth = 1.25),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.background = element_blank()
    )
  )

  ggsave(km_plot_tert$plot,
         filename = "km_tertile_split.svg",
         path = output_dir, width = 8, height = 7)

  results$km <- list(
    logrank_p_median = lr_pval,
    logrank_p_tertile = lr_pval_tert,
    km_fit_median = km_fit,
    km_fit_tertile = km_fit_tert
  )

  if (verbose) {
    message(sprintf("  Median split log-rank p = %.4e", lr_pval))
    message(sprintf("  Tertile split log-rank p = %.4e", lr_pval_tert))
  }

  # ===================================================================
  # 3. CALIBRATION: predicted vs observed survival at fixed timepoints
  # ===================================================================
  if (verbose) message("[3/6] Calibration analysis...")

  # Fit a Cox model on the risk score to get predicted survival
  cox_continuous <- coxph(
    Surv(df_clean[[make.names(time_col)]], df_clean$vital.status) ~ risk_score,
    data = df_clean
  )

  # C-index from continuous risk score
  cindex_continuous <- summary(cox_continuous)$concordance[1]
  cindex_se <- summary(cox_continuous)$concordance[2]

  # Create risk quantile groups for calibration
  n_groups <- 4
  df_clean$risk_quantile <- cut(
    risk_scores,
    breaks = quantile(risk_scores, probs = seq(0, 1, length.out = n_groups + 1)),
    include.lowest = TRUE,
    labels = paste0("Q", 1:n_groups)
  )

  # For each quantile group, compute observed vs predicted survival
  # at the median follow-up time
  median_time <- median(df_clean[[make.names(time_col)]])
  calibration_data <- data.frame(
    group = character(),
    n = integer(),
    observed_surv = numeric(),
    predicted_surv = numeric(),
    stringsAsFactors = FALSE
  )

  for (q in levels(df_clean$risk_quantile)) {
    q_data <- df_clean[df_clean$risk_quantile == q, ]

    # Observed: KM estimate at median time
    km_q <- survfit(
      Surv(q_data[[make.names(time_col)]], q_data$vital.status) ~ 1
    )
    # Find the KM estimate at or near the median time
    km_times <- km_q$time
    km_surv <- km_q$surv
    closest_idx <- which.min(abs(km_times - median_time))
    obs_surv <- if (length(closest_idx) > 0) km_surv[closest_idx] else NA

    # Predicted: from Cox model baseline survival
    pred_surv <- tryCatch({
      base_surv <- survfit(cox_continuous, newdata = q_data)
      # Average predicted survival at median time across group members
      surv_at_time <- summary(base_surv, times = median_time)$surv
      mean(surv_at_time, na.rm = TRUE)
    }, error = function(e) NA)

    calibration_data <- rbind(calibration_data, data.frame(
      group = q,
      n = nrow(q_data),
      observed_surv = round(obs_surv, 3),
      predicted_surv = round(pred_surv, 3),
      stringsAsFactors = FALSE
    ))
  }

  # Calibration plot
  if (all(!is.na(calibration_data$observed_surv)) &&
      all(!is.na(calibration_data$predicted_surv))) {
    cal_plot <- ggplot(calibration_data,
                       aes(x = predicted_surv, y = observed_surv)) +
      geom_point(size = 4) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_text(aes(label = group), vjust = -1, size = 4) +
      xlim(0, 1) + ylim(0, 1) +
      labs(
        x = "Predicted Survival Probability",
        y = "Observed Survival Probability (KM)",
        title = sprintf("Calibration at Median Follow-up (%.0f days)", median_time)
      ) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12)
      )

    ggsave(cal_plot,
           filename = "calibration_plot.svg",
           path = output_dir, width = 7, height = 7)
  }

  results$calibration <- list(
    c_index = cindex_continuous,
    c_index_se = cindex_se,
    calibration_data = calibration_data,
    cox_model = cox_continuous
  )

  if (verbose) {
    message(sprintf("  C-index (continuous): %.4f (SE: %.4f)",
                    cindex_continuous, cindex_se))
  }

  # ===================================================================
  # 4. PERMUTATION TEST: is the signature better than random?
  # ===================================================================
  if (verbose) message(sprintf("[4/6] Permutation test (%d permutations)...",
                               n_permutations))

  # Observed C-index
  observed_cindex <- cindex_continuous

  # Permute survival outcomes and recompute C-index
  perm_cindices <- numeric(n_permutations)
  for (p in 1:n_permutations) {
    perm_idx <- sample(nrow(df_clean))
    perm_y <- Surv(
      time = df_clean[[make.names(time_col)]][perm_idx],
      event = df_clean$vital.status[perm_idx]
    )
    perm_cox <- tryCatch({
      coxph(perm_y ~ risk_scores)
    }, error = function(e) NULL, warning = function(w) {
      suppressWarnings(coxph(perm_y ~ risk_scores))
    })

    if (!is.null(perm_cox)) {
      perm_cindices[p] <- summary(perm_cox)$concordance[1]
    } else {
      perm_cindices[p] <- 0.5
    }
  }

  perm_pval <- mean(perm_cindices >= observed_cindex)

  # Permutation distribution plot
  perm_df <- data.frame(cindex = perm_cindices)
  perm_plot <- ggplot(perm_df, aes(x = cindex)) +
    geom_histogram(bins = 50, fill = "gray70", color = "gray40") +
    geom_vline(xintercept = observed_cindex, color = "red",
               linewidth = 1.2, linetype = "solid") +
    annotate("text",
             x = observed_cindex, y = Inf,
             label = sprintf("Observed = %.4f\np = %.4f",
                             observed_cindex, perm_pval),
             vjust = 2, hjust = -0.1, color = "red", size = 4.5,
             fontface = "bold") +
    labs(
      x = "C-index (permuted data)",
      y = "Count",
      title = sprintf("Permutation Test (n = %d)", n_permutations)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12)
    )

  ggsave(perm_plot,
         filename = "permutation_test.svg",
         path = output_dir, width = 7, height = 5)

  results$permutation <- list(
    observed_cindex = observed_cindex,
    perm_cindices = perm_cindices,
    perm_pval = perm_pval
  )

  if (verbose) {
    message(sprintf("  Observed C-index: %.4f", observed_cindex))
    message(sprintf("  Permutation p-value: %.4f", perm_pval))
  }

  # ===================================================================
  # 5. BOOTSTRAP COEFFICIENT STABILITY
  # ===================================================================
  if (verbose) message(sprintf("[5/6] Bootstrap coefficient stability (%d resamples)...",
                               n_boot))

  boot_coefs <- matrix(NA, nrow = n_boot, ncol = length(signature_genes))
  colnames(boot_coefs) <- make.names(signature_genes)
  boot_cindices <- numeric(n_boot)

  for (b in 1:n_boot) {
    boot_idx <- sample(nrow(df_clean), replace = TRUE)
    boot_data <- df_clean[boot_idx, ]

    boot_x <- tryCatch(
      model.matrix(as.formula(formula_str), boot_data),
      error = function(e) NULL
    )
    if (is.null(boot_x)) next

    boot_y <- Surv(
      time = boot_data[[make.names(time_col)]],
      event = boot_data$vital.status
    )

    boot_cv <- tryCatch(
      cv.glmnet(boot_x, boot_y, family = "cox", alpha = 0,
                type.measure = "C", nfolds = 10),
      error = function(e) NULL
    )
    if (is.null(boot_cv)) next

    boot_c <- coef(boot_cv, s = "lambda.1se")
    for (g in colnames(boot_coefs)) {
      if (g %in% rownames(boot_c)) {
        boot_coefs[b, g] <- boot_c[g, 1]
      }
    }

    # C-index on the bootstrap sample
    boot_pred <- predict(boot_cv, newx = boot_x, s = "lambda.1se", type = "link")
    boot_conc <- tryCatch({
      survival::concordance(boot_y ~ as.vector(boot_pred))$concordance
    }, error = function(e) NA)
    boot_cindices[b] <- boot_conc
  }

  # Coefficient stability plot
  boot_coefs_long <- as.data.frame(boot_coefs) %>%
    pivot_longer(everything(), names_to = "gene", values_to = "coefficient") %>%
    filter(!is.na(coefficient))

  coef_plot <- ggplot(boot_coefs_long, aes(x = gene, y = coefficient)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Gene",
      y = "Bootstrap Coefficient",
      title = sprintf("Coefficient Stability (%d Bootstrap Resamples)", n_boot)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(coef_plot,
         filename = "coefficient_stability.svg",
         path = output_dir, width = 7, height = 5)

  # Summary of coefficient sign consistency
  sign_consistency <- apply(boot_coefs, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    max(mean(x > 0), mean(x < 0))
  })

  results$bootstrap <- list(
    boot_coefs = boot_coefs,
    boot_cindices = boot_cindices[!is.na(boot_cindices)],
    sign_consistency = sign_consistency,
    coef_summary = apply(boot_coefs, 2, function(x) {
      x <- x[!is.na(x)]
      c(mean = mean(x), sd = sd(x),
        q025 = quantile(x, 0.025), q975 = quantile(x, 0.975))
    })
  )

  if (verbose) {
    message("  Coefficient sign consistency:")
    for (g in names(sign_consistency)) {
      message(sprintf("    %s: %.1f%%", g, sign_consistency[g] * 100))
    }
    message(sprintf("  Bootstrap C-index: %.4f (SD: %.4f)",
                    mean(boot_cindices, na.rm = TRUE),
                    sd(boot_cindices, na.rm = TRUE)))
  }

  # ===================================================================
  # 6. TIME-DEPENDENT AUC (if timeROC is available)
  # ===================================================================
  has_timeROC <- requireNamespace("timeROC", quietly = TRUE)

  if (has_timeROC) {
    if (verbose) message("[6/6] Time-dependent AUC...")

    max_time <- max(df_clean[[make.names(time_col)]], na.rm = TRUE)

    # Set reasonable timepoints (1, 2, 3, 5 year marks)
    if (max_time > 1000) {
      # Time in days
      roc_times <- c(365, 730, 1095, 1825)
      roc_labels <- c("1-year", "2-year", "3-year", "5-year")
    } else {
      # Time in months
      roc_times <- c(12, 24, 36, 60)
      roc_labels <- c("1-year", "2-year", "3-year", "5-year")
    }
    # Only use timepoints within follow-up
    valid <- roc_times < max_time * 0.9
    roc_times <- roc_times[valid]
    roc_labels <- roc_labels[valid]

    if (length(roc_times) > 0) {
      troc <- tryCatch({
        timeROC::timeROC(
          T = df_clean[[make.names(time_col)]],
          delta = df_clean$vital.status,
          marker = risk_scores,
          cause = 1,
          times = roc_times,
          iid = TRUE
        )
      }, error = function(e) {
        if (verbose) message("  timeROC failed: ", e$message)
        NULL
      })

      if (!is.null(troc)) {
        # Extract AUCs (first element is t=0, skip it)
        auc_values <- troc$AUC[-1]
        names(auc_values) <- roc_labels

        # Time-dependent AUC plot
        auc_df <- data.frame(
          time_label = factor(roc_labels, levels = roc_labels),
          time = roc_times,
          auc = auc_values
        )

        auc_plot <- ggplot(auc_df, aes(x = time_label, y = auc)) +
          geom_col(fill = "steelblue", width = 0.6) +
          geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
          geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.5, size = 4.5) +
          ylim(0, 1) +
          labs(
            x = "Timepoint",
            y = "AUC",
            title = "Time-Dependent AUC"
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.text = element_text(size = 12)
          )

        ggsave(auc_plot,
               filename = "time_dependent_auc.svg",
               path = output_dir, width = 6, height = 5)

        results$time_auc <- list(
          auc_values = auc_values,
          timeROC_obj = troc
        )

        if (verbose) {
          for (i in seq_along(roc_labels)) {
            message(sprintf("  %s AUC: %.4f", roc_labels[i], auc_values[i]))
          }
        }
      }
    }
  } else {
    if (verbose) message("[6/6] Skipping time-dependent AUC (install timeROC)")
  }

  # ===================================================================
  # SUMMARY
  # ===================================================================
  if (verbose) {
    message("\n===== Internal Validation Summary =====")
    message(sprintf("Signature genes: %s",
                    paste(signature_genes, collapse = ", ")))
    message(sprintf("C-index:         %.4f (SE: %.4f)",
                    cindex_continuous, cindex_se))
    message(sprintf("KM log-rank p:   %.4e (median split)", lr_pval))
    message(sprintf("KM log-rank p:   %.4e (tertile split)", lr_pval_tert))
    message(sprintf("Permutation p:   %.4f", perm_pval))
    message(sprintf("Boot C-index:    %.4f (SD: %.4f)",
                    mean(boot_cindices, na.rm = TRUE),
                    sd(boot_cindices, na.rm = TRUE)))
    message(sprintf("\nAll plots saved to: %s", output_dir))
  }

  # Save summary CSV
  summary_out <- data.frame(
    metric = c("C-index", "C-index SE",
               "Log-rank p (median)", "Log-rank p (tertile)",
               "Permutation p-value",
               "Bootstrap C-index mean", "Bootstrap C-index SD",
               paste0("Sign consistency: ", names(sign_consistency))),
    value = c(
      round(cindex_continuous, 4), round(cindex_se, 4),
      format.pval(lr_pval, digits = 4),
      format.pval(lr_pval_tert, digits = 4),
      round(perm_pval, 4),
      round(mean(boot_cindices, na.rm = TRUE), 4),
      round(sd(boot_cindices, na.rm = TRUE), 4),
      round(sign_consistency, 3)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(summary_out,
            file.path(output_dir, "validation_summary.csv"),
            row.names = FALSE)

  return(results)
}


# ===========================================================================
# Example usage (uncomment and modify paths to run)
# ===========================================================================
#
# source("R/nested_cv.R")
# source("R/km_plotter.R")
#
# scl_common <- readRDS("Outputs/scl_common.rds")
#
# # --- After discovery: validate the stable genes ---
# # Suppose nested CV discovery found these genes stable:
# # (replace with your actual discovery results)
#
# # Option A: Use your original 4-gene signature with known coefficients
# val_results <- internal_validation(
#   cox_df = scl_common,
#   signature_genes = c("tfrc", "fam83f", "dlk1", "gng13"),
#   signature_coefs = c(tfrc = -0.9526, fam83f = -0.0855,
#                       dlk1 = -0.0006, gng13 = -0.2558),
#   n_permutations = 1000,
#   n_boot = 1000,
#   output_dir = "Outputs/internal_validation/four_gene"
# )
#
# # Option B: Let the function estimate coefficients via ridge
# # (better for a newly discovered signature)
# val_results <- internal_validation(
#   cox_df = scl_common,
#   signature_genes = c("gene1", "gene2", "gene3"),  # from discovery
#   signature_coefs = NULL,  # will be estimated
#   n_permutations = 1000,
#   n_boot = 1000,
#   output_dir = "Outputs/internal_validation/discovered_signature"
# )
#
#
# # ============================
# # INTERPRETING RESULTS
# # ============================
# #
# # What each analysis tells a reviewer:
# #
# # 1. KAPLAN-MEIER (median + tertile splits)
# #    Shows the risk score actually separates patients into distinct
# #    survival trajectories. Tertile split is stronger evidence than
# #    median split because it shows a dose-response: higher risk score
# #    = progressively worse survival. Reviewers expect this figure.
# #
# # 2. CALIBRATION
# #    Tests whether predicted survival probabilities match observed
# #    ones. Points near the diagonal = well-calibrated. Points above
# #    the line = model is too pessimistic. Below = too optimistic.
# #    This catches models that rank patients correctly (decent C-index)
# #    but give wrong absolute survival probabilities.
# #
# # 3. PERMUTATION TEST
# #    Answers: "could we have gotten this C-index by chance?" If
# #    the observed C-index is far to the right of the null distribution,
# #    the signature captures real signal. p < 0.01 is strong evidence.
# #
# # 4. BOOTSTRAP COEFFICIENT STABILITY
# #    Shows whether each gene's coefficient is consistent across
# #    resamples. If a gene's coefficient flips sign in 40% of
# #    bootstraps, it's not a reliable predictor. Sign consistency
# #    > 90% is strong; < 70% is concerning.
# #
# # 5. TIME-DEPENDENT AUC
# #    Shows discriminative ability at specific timepoints (1, 2, 3 years).
# #    More informative than a single C-index because it reveals whether
# #    the model predicts better for short-term or long-term survival.
# #    AUC > 0.7 at any timepoint is good for SCLC.
