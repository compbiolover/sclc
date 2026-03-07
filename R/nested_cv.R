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
#'   For discovery mode, set this high (e.g., 500-2000) to let the lasso
#'   select from a large pool.
#' @param lambda_rule Which lambda to use from inner CV: "lambda.1se"
#'   (default, more conservative/generalizable) or "lambda.min" (best
#'   in-sample performance but more prone to overfitting).
#' @param progress_free Logical; if TRUE, use time2 column instead of time.
#' @param output_dir Directory to save results. Created if it doesn't exist.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with:
#'   - summary: data frame with mean, sd, median, and IQR of outer C-indices
#'   - per_seed: data frame with per-seed mean outer C-index
#'   - all_folds: data frame with every outer fold C-index
#'   - gene_stability: data frame showing how often each gene was selected
#'     across all outer folds. In discovery mode, genes that appear with
#'     high frequency (>50%) across folds are strong signature candidates.
nested_cv <- function(cox_df,
                      candidate_genes,
                      n_outer_folds = 5,
                      n_inner_folds = 10,
                      master_seeds = 1:50,
                      my_alpha = 1,
                      max_it = 100000,
                      gene_num = 20,
                      lambda_rule = "lambda.1se",
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
  if (!lambda_rule %in% c("lambda.1se", "lambda.min")) {
    stop("lambda_rule must be 'lambda.1se' or 'lambda.min'.")
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
      inner_coefs <- coef(inner_cv, s = lambda_rule)
      selected_idx <- which(as.vector(inner_coefs) != 0)
      selected_genes <- rownames(inner_coefs)[selected_idx]

      # Track gene selection frequency
      for (g in selected_genes) {
        gene_selection_counts[[g]] <- (gene_selection_counts[[g]] %||% 0) + 1
      }

      # ----- Evaluate on outer test fold -----
      test_x <- tryCatch(
        model.matrix(as.formula(formula_str), test_data),
        error = function(e) NULL
      )
      if (is.null(test_x)) next

      test_y <- Surv(
        time = test_data[[make.names(time_col)]],
        event = test_data$vital.status
      )

      # Predict risk scores on test data using the chosen lambda rule
      test_pred <- tryCatch(
        as.vector(predict(inner_cv, newx = test_x, s = lambda_rule,
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
          lambda_selected = inner_cv[[lambda_rule]],
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
    lambda_rule = lambda_rule,
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
# # Load your data (adjust paths as needed)
# scl_common <- readRDS("Outputs/scl_common.rds")
# # Or rebuild scl_common from main.Rmd data loading steps
#
#
# # ============================
# # VALIDATION MODE
# # ============================
# # Test whether an existing gene signature generalizes honestly.
# # Use alpha = 0 (ridge) since genes are pre-selected — the inner CV
# # only tunes lambda, not which genes to include.
#
# # --- Your 4-gene signature (tfrc, fam83f, dlk1, gng13) ---
# four_gene_results <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = c("tfrc", "fam83f", "dlk1", "gng13"),
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:100,
#   my_alpha = 0,
#   gene_num = 4,
#   lambda_rule = "lambda.1se",
#   output_dir = "Outputs/nested_cv/four_gene_validation"
# )
#
#
# # ============================
# # DISCOVERY MODE
# # ============================
# # Let the lasso discover which genes matter from a large candidate pool.
# # Use alpha = 1 (lasso) so the inner CV performs gene selection AND
# # lambda tuning. The gene_stability output tells you which genes the
# # lasso consistently picks across all outer folds and seeds — those
# # are your new candidate signature genes.
#
# # --- Discovery from SDE genes ---
# sde_genes <- readRDS("Outputs/sde_genes_new_random_seed.rds")
# sde_discovery <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = sde_genes,
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:50,
#   my_alpha = 1,
#   gene_num = 500,       # let the lasso see a large pool
#   lambda_rule = "lambda.1se",
#   output_dir = "Outputs/nested_cv/sde_discovery"
# )
# # Check sde_discovery$gene_stability for consistently selected genes
#
# # --- Discovery from MAD genes ---
# mad_genes <- readRDS("Outputs/mad_genes_new_random_seed.rds")
# mad_discovery <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = mad_genes,
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:50,
#   my_alpha = 1,
#   gene_num = 500,
#   lambda_rule = "lambda.1se",
#   output_dir = "Outputs/nested_cv/mad_discovery"
# )
#
# # --- Discovery from miRNA target genes ---
# mirna_genes <- readRDS("Outputs/mirna/rds/mirna_genes_200_mirnas_5_targets_up_.rds")
# mirna_discovery <- nested_cv(
#   cox_df = scl_common,
#   candidate_genes = mirna_genes,
#   n_outer_folds = 5,
#   n_inner_folds = 10,
#   master_seeds = 1:50,
#   my_alpha = 1,
#   gene_num = 500,
#   lambda_rule = "lambda.1se",
#   output_dir = "Outputs/nested_cv/mirna_discovery"
# )
#
#
# # ============================
# # COMPARE METHODS
# # ============================
# comparison <- data.frame(
#   method = c("4-gene (validation)", "SDE (discovery)",
#              "MAD (discovery)", "miRNA (discovery)"),
#   mean_cindex = c(four_gene_results$summary$mean_cindex,
#                   sde_discovery$summary$mean_cindex,
#                   mad_discovery$summary$mean_cindex,
#                   mirna_discovery$summary$mean_cindex),
#   sd_cindex = c(four_gene_results$summary$sd_cindex,
#                 sde_discovery$summary$sd_cindex,
#                 mad_discovery$summary$sd_cindex,
#                 mirna_discovery$summary$sd_cindex)
# )
# print(comparison)
#
#
# # ============================
# # DISCOVERY WORKFLOW
# # ============================
# # After running discovery mode:
# #
# # 1. Look at gene_stability output:
# #    print(sde_discovery$gene_stability)
# #
# # 2. Genes with selection_frequency > 0.5 are strong candidates.
# #    These genes were independently selected by the lasso in >50%
# #    of outer folds across different seeds — meaning they're not
# #    artifacts of a particular train/test split.
# #
# # 3. Take those stable genes and run VALIDATION mode to get their
# #    honest C-index:
# #
# #    stable_genes <- sde_discovery$gene_stability %>%
# #      filter(selection_frequency > 0.5) %>%
# #      pull(gene)
# #
# #    new_sig_results <- nested_cv(
# #      cox_df = scl_common,
# #      candidate_genes = stable_genes,
# #      my_alpha = 0,
# #      gene_num = length(stable_genes),
# #      lambda_rule = "lambda.1se",
# #      master_seeds = 1:100,
# #      output_dir = "Outputs/nested_cv/new_signature_validation"
# #    )
# #
# # NOTE: This two-step workflow (discover then validate) still has
# # some optimism because both steps use the same dataset. The nested
# # CV C-index from the discovery step is already an honest estimate.
# # The validation step is mainly to confirm with ridge regression
# # that those specific genes carry signal on their own.
