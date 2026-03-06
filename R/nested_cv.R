# Name: nested_cv.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Perform nested cross-validation with varying master seeds to get
# an honest, unbiased estimate of model generalization performance. The outer
# loop holds out test folds that are never seen during gene selection or
# lambda tuning, eliminating the optimism bias present in standard CV.

suppressMessages({
  library(doParallel)
  library(glmnet)
  library(parallel)
  library(survival)
  library(tidyverse)
})

#' Run nested cross-validation for a penalized Cox model
#'
#' @param cox_df Data frame with columns: patient.id, vital.status, time,
#'   and gene expression columns.
#' @param candidate_genes Character vector of candidate gene names to consider.
#'   These are the full set of genes that the inner CV will select from.
#' @param n_outer_folds Number of outer CV folds (default 5).
#' @param n_inner_folds Number of inner CV folds (default 10).
#' @param master_seeds Integer vector of master seeds. The entire nested CV
#'   procedure is repeated for each seed. Use 50-100 seeds to characterize
#'   the distribution of performance.
#' @param my_alpha Elastic net mixing parameter (1 = lasso, 0 = ridge).
#' @param max_it Maximum iterations for glmnet (default 100000).
#' @param gene_num Maximum number of candidate genes to use (default 20).
#' @param progress_free Logical; if TRUE, use time2 column instead of time.
#' @param output_dir Directory to save results. Created if it doesn't exist.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with:
#'   - summary: data frame with mean, sd, median, and IQR of outer C-indices
#'   - per_seed: data frame with per-seed mean outer C-index
#'   - all_folds: data frame with every outer fold C-index
#'   - gene_stability: data frame showing how often each gene was selected
nested_cv <- function(cox_df,
                      candidate_genes,
                      n_outer_folds = 5,
                      n_inner_folds = 10,
                      master_seeds = 1:50,
                      my_alpha = 1,
                      max_it = 100000,
                      gene_num = 20,
                      progress_free = FALSE,
                      output_dir = "Outputs/nested_cv",
                      verbose = TRUE) {
  # Input validation
  if (is.null(cox_df) || !is.data.frame(cox_df)) {
    stop("cox_df must be a data.frame.")
  }
  if (is.null(candidate_genes) || length(candidate_genes) == 0) {
    stop("candidate_genes must be a non-empty character vector.")
  }

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set up parallel processing
  num_cores <- parallel::detectCores()
  registerDoParallel(cores = num_cores)

  # Prepare gene list: intersect candidates with available columns
  if (is.data.frame(candidate_genes)) {
    candidate_genes <- rownames(candidate_genes)
  }
  candidate_genes <- tolower(candidate_genes)
  candidate_genes <- intersect(candidate_genes, colnames(cox_df))
  candidate_genes <- head(candidate_genes, n = gene_num)

  if (length(candidate_genes) == 0) {
    stop("No candidate genes found in cox_df columns.")
  }

  # Determine time column
  time_col <- if (progress_free) "time2" else "time"
  if (!time_col %in% colnames(cox_df)) {
    stop(paste0("Column '", time_col, "' not found in cox_df."))
  }
  if (!"vital.status" %in% colnames(cox_df)) {
    stop("Column 'vital.status' not found in cox_df.")
  }

  # Storage for all results
  all_fold_results <- list()
  gene_selection_counts <- list()

  for (s_idx in seq_along(master_seeds)) {
    seed <- master_seeds[s_idx]
    set.seed(seed)

    if (verbose) {
      message(sprintf("Master seed %d/%d (seed = %d)",
                      s_idx, length(master_seeds), seed))
    }

    # Create outer fold assignments (stratified by vital.status)
    n <- nrow(cox_df)
    # Stratified sampling: balance censored/uncensored across folds
    deceased_idx <- which(cox_df$vital.status == 1)
    alive_idx <- which(cox_df$vital.status == 0)

    outer_folds <- integer(n)
    outer_folds[deceased_idx] <- sample(rep(1:n_outer_folds,
                                            length.out = length(deceased_idx)))
    outer_folds[alive_idx] <- sample(rep(1:n_outer_folds,
                                          length.out = length(alive_idx)))

    for (k in 1:n_outer_folds) {
      # Split into outer train and outer test
      test_idx <- which(outer_folds == k)
      train_idx <- which(outer_folds != k)

      train_data <- cox_df[train_idx, ]
      test_data <- cox_df[test_idx, ]

      # Check that test set has events
      if (sum(test_data$vital.status) == 0 || nrow(test_data) < 3) {
        if (verbose) {
          message(sprintf("  Fold %d: skipping (insufficient test events)", k))
        }
        next
      }

      # ----- Inner CV: feature selection + lambda tuning on train_data -----
      # Build model matrix from candidate genes
      available_genes <- intersect(candidate_genes, colnames(train_data))
      if (length(available_genes) == 0) next

      formula_str <- paste0("~", paste0(make.names(available_genes),
                                         collapse = "+"))
      colnames(train_data) <- make.names(colnames(train_data))
      colnames(test_data) <- make.names(colnames(test_data))

      train_x <- tryCatch(
        model.matrix(as.formula(formula_str), train_data),
        error = function(e) NULL
      )
      if (is.null(train_x)) next

      train_y <- Surv(
        time = train_data[[make.names(time_col)]],
        event = train_data$vital.status
      )

      # Inner CV to select lambda
      inner_folds <- sample(rep(1:n_inner_folds, length.out = nrow(train_data)))

      inner_cv <- tryCatch(
        cv.glmnet(
          x = train_x, y = train_y,
          nfolds = n_inner_folds,
          type.measure = "C",
          maxit = max_it,
          family = "cox",
          parallel = TRUE,
          alpha = my_alpha,
          foldid = inner_folds
        ),
        error = function(e) NULL
      )
      if (is.null(inner_cv)) next

      # Get selected genes from inner CV
      inner_coefs <- coef(inner_cv, s = "lambda.min")
      selected_idx <- which(as.vector(inner_coefs) != 0)
      selected_genes <- rownames(inner_coefs)[selected_idx]

      # Track gene selection frequency
      for (g in selected_genes) {
        gene_selection_counts[[g]] <- (gene_selection_counts[[g]] %||% 0) + 1
      }

      # ----- Evaluate on outer test fold -----
      # Refit on full training data with lambda.min from inner CV
      test_x <- tryCatch(
        model.matrix(as.formula(formula_str), test_data),
        error = function(e) NULL
      )
      if (is.null(test_x)) next

      test_y <- Surv(
        time = test_data[[make.names(time_col)]],
        event = test_data$vital.status
      )

      # Predict risk scores on test data
      test_pred <- tryCatch(
        as.vector(predict(inner_cv, newx = test_x, s = "lambda.min",
                          type = "link")),
        error = function(e) NULL
      )
      if (is.null(test_pred)) next

      # Calculate C-index on held-out test fold
      outer_cindex <- tryCatch({
        concordance_obj <- survival::concordance(test_y ~ test_pred)
        concordance_obj$concordance
      }, error = function(e) NA_real_)

      if (!is.na(outer_cindex)) {
        all_fold_results[[length(all_fold_results) + 1]] <- data.frame(
          master_seed = seed,
          outer_fold = k,
          outer_cindex = round(outer_cindex, 4),
          n_train = nrow(train_data),
          n_test = nrow(test_data),
          n_selected_genes = length(selected_genes),
          selected_genes = paste(selected_genes, collapse = ","),
          lambda_min = inner_cv$lambda.min,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Compile results
  if (length(all_fold_results) == 0) {
    stop("No valid results produced. Check your data and parameters.")
  }

  all_folds_df <- do.call(rbind, all_fold_results)

  # Per-seed summary (mean outer C-index per master seed)
  per_seed_df <- all_folds_df %>%
    group_by(master_seed) %>%
    summarise(
      mean_outer_cindex = mean(outer_cindex, na.rm = TRUE),
      sd_outer_cindex = sd(outer_cindex, na.rm = TRUE),
      n_folds_completed = n(),
      mean_n_genes = mean(n_selected_genes),
      .groups = "drop"
    )

  # Overall summary
  summary_df <- data.frame(
    mean_cindex = mean(per_seed_df$mean_outer_cindex),
    sd_cindex = sd(per_seed_df$mean_outer_cindex),
    median_cindex = median(per_seed_df$mean_outer_cindex),
    iqr_lower = quantile(per_seed_df$mean_outer_cindex, 0.25),
    iqr_upper = quantile(per_seed_df$mean_outer_cindex, 0.75),
    n_seeds = length(master_seeds),
    n_outer_folds = n_outer_folds,
    n_inner_folds = n_inner_folds,
    alpha = my_alpha,
    stringsAsFactors = FALSE
  )

  # Gene stability: how often each gene was selected across all outer folds
  total_fits <- nrow(all_folds_df)
  gene_stability_df <- data.frame(
    gene = names(gene_selection_counts),
    times_selected = unlist(gene_selection_counts),
    selection_frequency = round(unlist(gene_selection_counts) / total_fits, 3),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(selection_frequency))
  rownames(gene_stability_df) <- NULL

  # Save outputs
  write.csv(all_folds_df,
            file.path(output_dir, "nested_cv_all_folds.csv"),
            row.names = FALSE)
  write.csv(per_seed_df,
            file.path(output_dir, "nested_cv_per_seed.csv"),
            row.names = FALSE)
  write.csv(summary_df,
            file.path(output_dir, "nested_cv_summary.csv"),
            row.names = FALSE)
  write.csv(gene_stability_df,
            file.path(output_dir, "nested_cv_gene_stability.csv"),
            row.names = FALSE)

  if (verbose) {
    message("\n===== Nested CV Results =====")
    message(sprintf("Mean outer C-index:   %.4f", summary_df$mean_cindex))
    message(sprintf("SD across seeds:      %.4f", summary_df$sd_cindex))
    message(sprintf("Median outer C-index: %.4f", summary_df$median_cindex))
    message(sprintf("IQR:                  [%.4f, %.4f]",
                    summary_df$iqr_lower, summary_df$iqr_upper))
    message(sprintf("\nTop stable genes (selected in >50%% of fits):"))
    stable <- gene_stability_df %>% filter(selection_frequency > 0.5)
    if (nrow(stable) > 0) {
      for (i in 1:nrow(stable)) {
        message(sprintf("  %s: %.1f%%", stable$gene[i],
                        stable$selection_frequency[i] * 100))
      }
    } else {
      message("  (none above 50% threshold)")
    }
    message(sprintf("\nResults saved to: %s", output_dir))
  }

  return(list(
    summary = summary_df,
    per_seed = per_seed_df,
    all_folds = all_folds_df,
    gene_stability = gene_stability_df
  ))
}


# ===========================================================================
# Example usage (uncomment and modify paths to run)
# ===========================================================================
#
# source("R/cox_model.R")
#
# # Load your data (adjust paths as needed)
# scl_common <- readRDS("Outputs/scl_common.rds")
# # Or rebuild scl_common from main.Rmd data loading steps
#
# # --- Example 1: SDE genes ---
# sde_genes <- readRDS("Outputs/sde_genes_new_random_seed.rds")
# sde_results <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = sde_genes,
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:50,
#   my_alpha = 1,
#   gene_num = 20,
#   output_dir = "Outputs/nested_cv/sde"
# )
#
# # --- Example 2: Your 4-gene signature (tfrc, fam83f, dlk1, gng13) ---
# four_gene_results <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = c("tfrc", "fam83f", "dlk1", "gng13"),
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:100,
#   my_alpha = 0,  # ridge, since you have a fixed small gene set
#   gene_num = 4,
#   output_dir = "Outputs/nested_cv/four_gene"
# )
#
# # --- Example 3: MAD genes ---
# mad_genes <- readRDS("Outputs/mad_genes_new_random_seed.rds")
# mad_results <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = mad_genes,
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:50,
#   my_alpha = 1,
#   gene_num = 20,
#   output_dir = "Outputs/nested_cv/mad"
# )
#
# # --- Example 4: miRNA target genes ---
# mirna_genes <- readRDS("Outputs/mirna/rds/mirna_genes_200_mirnas_5_targets_up_.rds")
# mirna_results <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = mirna_genes,
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:50,
#   my_alpha = 1,
#   gene_num = 20,
#   output_dir = "Outputs/nested_cv/mirna"
# )
#
# # --- Compare methods ---
# comparison <- data.frame(
#   method = c("SDE", "4-gene", "MAD", "miRNA"),
#   mean_cindex = c(sde_results$summary$mean_cindex,
#                   four_gene_results$summary$mean_cindex,
#                   mad_results$summary$mean_cindex,
#                   mirna_results$summary$mean_cindex),
#   sd_cindex = c(sde_results$summary$sd_cindex,
#                 four_gene_results$summary$sd_cindex,
#                 mad_results$summary$sd_cindex,
#                 mirna_results$summary$sd_cindex)
# )
# print(comparison)
