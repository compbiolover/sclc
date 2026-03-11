# Name: nested_cv.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Perform nested cross-validation with varying master seeds to get
# an honest, unbiased estimate of model generalization performance. The outer
# loop holds out test folds that are never seen during gene selection or
# lambda tuning, eliminating the optimism bias present in standard CV.
#
# Seeds are parallelized via mclapply (fork-based) for M-series Macs.
# Each seed runs all outer folds sequentially; cv.glmnet runs single-threaded
# within each fork to avoid nested parallelism overhead.

suppressMessages({
  library(glmnet)
  library(parallel)
  library(survival)
  library(tidyverse)
})


# =========================================================================
# Helper: run nested CV for a single master seed
# =========================================================================
#
# This is the workhorse called by mclapply. It must be self-contained
# (no references to the parent environment beyond its arguments) so that
# forked processes work cleanly.

.run_one_seed <- function(seed, cox_df, candidate_genes, time_col,
                          n_outer_folds, n_inner_folds, n_inner_repeats,
                          my_alpha, max_it, lambda_rule, eval_times) {
  set.seed(seed)

  n <- nrow(cox_df)
  deceased_idx <- which(cox_df$vital.status == 1)
  alive_idx <- which(cox_df$vital.status == 0)

  # Stratified outer fold assignment

  outer_folds <- integer(n)
  outer_folds[deceased_idx] <- sample(rep(1:n_outer_folds,
                                          length.out = length(deceased_idx)))
  outer_folds[alive_idx] <- sample(rep(1:n_outer_folds,
                                        length.out = length(alive_idx)))

  fold_results <- list()
  gene_counts <- list()

  has_timeROC <- requireNamespace("timeROC", quietly = TRUE)

  for (k in 1:n_outer_folds) {
    test_idx <- which(outer_folds == k)
    train_idx <- which(outer_folds != k)

    train_data <- cox_df[train_idx, ]
    test_data <- cox_df[test_idx, ]

    # Need events in test set for concordance
    if (sum(test_data$vital.status) == 0 || nrow(test_data) < 3) next

    # Build model matrix
    available_genes <- intersect(candidate_genes, colnames(train_data))
    if (length(available_genes) == 0) next

    formula_str <- paste0("~", paste0(available_genes, collapse = "+"))

    train_x <- tryCatch(
      model.matrix(as.formula(formula_str), train_data),
      error = function(e) NULL
    )
    if (is.null(train_x)) next

    train_y <- Surv(
      time = train_data[[time_col]],
      event = train_data$vital.status
    )

    # Repeated inner CV: run n_inner_repeats times with different
    # stratified fold splits, then average the CV error curves to get
    # a more stable lambda estimate. This reduces noise in lambda
    # selection without biasing the outer fold evaluation.
    train_deceased <- which(train_data$vital.status == 1)
    train_alive <- which(train_data$vital.status == 0)

    inner_cv_fits <- list()
    for (rep_i in seq_len(n_inner_repeats)) {
      inner_folds <- integer(nrow(train_data))
      inner_folds[train_deceased] <- sample(rep(1:n_inner_folds,
                                                 length.out = length(train_deceased)))
      inner_folds[train_alive] <- sample(rep(1:n_inner_folds,
                                              length.out = length(train_alive)))

      fit <- tryCatch(
        cv.glmnet(
          x = train_x, y = train_y,
          nfolds = n_inner_folds,
          type.measure = "C",
          maxit = max_it,
          family = "cox",
          parallel = FALSE,
          alpha = my_alpha,
          foldid = inner_folds
        ),
        error = function(e) NULL
      )
      if (!is.null(fit)) inner_cv_fits[[length(inner_cv_fits) + 1]] <- fit
    }
    if (length(inner_cv_fits) == 0) next

    # Pick the best fit: average the CV error at each lambda across repeats,
    # then select the fit whose chosen lambda is closest to the averaged optimum.
    # For a single repeat, this just uses that fit directly.
    if (length(inner_cv_fits) == 1) {
      inner_cv <- inner_cv_fits[[1]]
    } else {
      # Find the common lambda grid (use the first fit's grid as reference)
      ref_lambdas <- inner_cv_fits[[1]]$lambda

      # Average the CV mean error (cvm) across repeats at each lambda
      cvm_matrix <- sapply(inner_cv_fits, function(fit) {
        # Interpolate each fit's cvm onto the reference lambda grid
        approx(x = log(fit$lambda), y = fit$cvm,
               xout = log(ref_lambdas), rule = 2)$y
      })
      avg_cvm <- rowMeans(cvm_matrix)

      # For type.measure = "C" (concordance), higher is better
      # cv.glmnet minimizes the loss, but for C-index it stores concordance
      # directly, so the best lambda maximizes avg_cvm
      best_idx <- which.max(avg_cvm)
      best_lambda <- ref_lambdas[best_idx]

      # Also compute the 1-SE rule on the averaged curve
      avg_cvsd <- if (ncol(cvm_matrix) > 1) {
        apply(cvm_matrix, 1, sd) / sqrt(ncol(cvm_matrix))
      } else {
        inner_cv_fits[[1]]$cvsd
      }
      best_cvm <- avg_cvm[best_idx]
      threshold_1se <- best_cvm - avg_cvsd[best_idx]
      # lambda.1se: largest lambda within 1 SE of the best
      candidates_1se <- which(avg_cvm >= threshold_1se)
      lambda_1se <- max(ref_lambdas[candidates_1se])

      # Use the first fit object but override its lambda selections
      inner_cv <- inner_cv_fits[[1]]
      inner_cv$lambda.min <- best_lambda
      inner_cv$lambda.1se <- lambda_1se
    }

    # Selected genes
    inner_coefs <- coef(inner_cv, s = lambda_rule)
    selected_idx <- which(as.vector(inner_coefs) != 0)
    selected_genes <- rownames(inner_coefs)[selected_idx]

    for (g in selected_genes) {
      gene_counts[[g]] <- (gene_counts[[g]] %||% 0L) + 1L
    }

    # Evaluate on outer test fold
    test_x <- tryCatch(
      model.matrix(as.formula(formula_str), test_data),
      error = function(e) NULL
    )
    if (is.null(test_x)) next

    test_y <- Surv(
      time = test_data[[time_col]],
      event = test_data$vital.status
    )

    test_pred <- tryCatch(
      as.vector(predict(inner_cv, newx = test_x, s = lambda_rule,
                        type = "link")),
      error = function(e) NULL
    )
    if (is.null(test_pred)) next

    # Global (Harrell's) C-index
    outer_cindex <- tryCatch({
      concordance_obj <- survival::concordance(test_y ~ test_pred)
      concordance_obj$concordance
    }, error = function(e) NA_real_)

    if (is.na(outer_cindex)) next

    # Time-dependent AUC via timeROC at each eval_time
    td_aucs <- rep(NA_real_, length(eval_times))
    names(td_aucs) <- paste0("auc_t", eval_times)

    if (has_timeROC && length(eval_times) > 0) {
      troc <- tryCatch(
        timeROC::timeROC(
          T = test_data[[time_col]],
          delta = test_data$vital.status,
          marker = test_pred,
          cause = 1,
          times = eval_times,
          iid = FALSE
        ),
        error = function(e) NULL
      )
      if (!is.null(troc)) {
        td_aucs <- troc$AUC[as.character(eval_times)]
      }
    }

    # Build result row
    row <- data.frame(
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

    # Append time-dependent AUC columns
    for (i in seq_along(eval_times)) {
      row[[paste0("auc_t", eval_times[i])]] <- round(td_aucs[i], 4)
    }

    fold_results[[length(fold_results) + 1]] <- row
  }

  list(fold_results = fold_results, gene_counts = gene_counts)
}


#' Run nested cross-validation for a penalized Cox model
#'
#' Seeds are parallelized via \code{parallel::mclapply} (fork-based),
#' which works well on M-series Macs. The inner \code{cv.glmnet} runs
#' single-threaded to avoid nested parallelism overhead.
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
#' @param n_inner_repeats Number of times to repeat the inner CV with
#'   different stratified fold splits (default 3). The CV error curves
#'   are averaged across repeats to produce a more stable lambda estimate.
#'   Only touches training data, so does not bias the outer fold evaluation.
#' @param eval_times Numeric vector of time points (in the same units as
#'   the time column) at which to compute time-dependent AUC. Requires
#'   the timeROC package. Default: c(12, 24, 36) for 1/2/3-year AUC.
#' @param progress_free Logical; if TRUE, use time2 column instead of time.
#' @param n_cores Number of cores for seed-level parallelism. Default uses
#'   all available cores. On M-series Macs, this uses fork-based mclapply.
#' @param output_dir Directory to save results. Created if it doesn't exist.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with:
#'   - summary: data frame with mean, sd, 95% CI, median, IQR of outer C-indices,
#'     plus mean time-dependent AUC at each eval_time
#'   - per_seed: data frame with per-seed mean outer C-index and AUCs
#'   - all_folds: data frame with every outer fold C-index and AUCs
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
                      n_inner_repeats = 3,
                      eval_times = c(12, 24, 36),
                      progress_free = FALSE,
                      n_cores = NULL,
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

  # Prepare gene list: intersect candidates with available columns
  if (is.data.frame(candidate_genes)) {
    candidate_genes <- rownames(candidate_genes)
  }
  candidate_genes <- tolower(candidate_genes)
  candidate_genes <- intersect(candidate_genes, colnames(cox_df))

  if (length(candidate_genes) > gene_num) {
    if (verbose) {
      message(sprintf("  Note: truncating %d candidates to gene_num = %d (first %d after intersection).",
                      length(candidate_genes), gene_num, gene_num))
    }
    candidate_genes <- head(candidate_genes, n = gene_num)
  }

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

  # Normalize column names once upfront to avoid repeated make.names() calls
  colnames(cox_df) <- make.names(colnames(cox_df))
  candidate_genes <- make.names(candidate_genes)
  time_col <- make.names(time_col)

  # Check for timeROC
  has_timeROC <- requireNamespace("timeROC", quietly = TRUE)
  if (!has_timeROC && length(eval_times) > 0) {
    if (verbose) {
      message("  timeROC package not installed — skipping time-dependent AUC.")
      message("  Install with: install.packages('timeROC')")
    }
    eval_times <- numeric(0)
  }

  # Filter eval_times to those within observed time range
  if (length(eval_times) > 0) {
    max_time <- max(cox_df[[time_col]], na.rm = TRUE)
    valid_times <- eval_times[eval_times < max_time]
    if (length(valid_times) < length(eval_times) && verbose) {
      dropped <- setdiff(eval_times, valid_times)
      message(sprintf("  Dropped eval_times > max observed time (%.1f): %s",
                      max_time, paste(dropped, collapse = ", ")))
    }
    eval_times <- valid_times
  }

  # Determine number of cores for seed-level parallelism
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores()
  }
  n_cores <- min(n_cores, length(master_seeds))

  if (verbose) {
    message(sprintf("Running nested CV: %d seeds x %d outer folds, %d cores",
                    length(master_seeds), n_outer_folds, n_cores))
    message(sprintf("  alpha = %s, lambda_rule = %s, %d candidate genes",
                    my_alpha, lambda_rule, length(candidate_genes)))
    message(sprintf("  Inner CV: %d folds x %d repeats (averaged for stable lambda)",
                    n_inner_folds, n_inner_repeats))
    if (length(eval_times) > 0) {
      message(sprintf("  Time-dependent AUC at: %s",
                      paste(eval_times, collapse = ", ")))
    }
  }

  # Run all seeds in parallel via mclapply (fork-based, ideal for macOS)
  seed_results <- parallel::mclapply(
    master_seeds,
    function(seed) {
      .run_one_seed(
        seed = seed,
        cox_df = cox_df,
        candidate_genes = candidate_genes,
        time_col = time_col,
        n_outer_folds = n_outer_folds,
        n_inner_folds = n_inner_folds,
        n_inner_repeats = n_inner_repeats,
        my_alpha = my_alpha,
        max_it = max_it,
        lambda_rule = lambda_rule,
        eval_times = eval_times
      )
    },
    mc.cores = n_cores,
    mc.set.seed = FALSE  # we set.seed() inside .run_one_seed
  )

  # Collect fold results and gene counts across all seeds
  all_fold_results <- list()
  gene_selection_counts <- list()

  for (res in seed_results) {
    if (is.null(res)) next

    # Fold results
    for (fr in res$fold_results) {
      all_fold_results[[length(all_fold_results) + 1]] <- fr
    }

    # Merge gene counts
    for (g in names(res$gene_counts)) {
      gene_selection_counts[[g]] <- (gene_selection_counts[[g]] %||% 0L) +
        res$gene_counts[[g]]
    }
  }

  # Compile results
  if (length(all_fold_results) == 0) {
    stop("No valid results produced. Check your data and parameters.")
  }

  all_folds_df <- do.call(rbind, all_fold_results)

  # Identify AUC columns for aggregation
  auc_cols <- grep("^auc_t", colnames(all_folds_df), value = TRUE)

  # Per-seed summary
  per_seed_df <- all_folds_df %>%
    group_by(master_seed) %>%
    summarise(
      mean_outer_cindex = mean(outer_cindex, na.rm = TRUE),
      sd_outer_cindex = sd(outer_cindex, na.rm = TRUE),
      n_folds_completed = n(),
      mean_n_genes = mean(n_selected_genes),
      across(all_of(auc_cols), ~ mean(.x, na.rm = TRUE),
             .names = "mean_{.col}"),
      .groups = "drop"
    )

  # Overall summary
  ci_vals <- per_seed_df$mean_outer_cindex
  ci_95_lower <- mean(ci_vals) - 1.96 * sd(ci_vals)
  ci_95_upper <- mean(ci_vals) + 1.96 * sd(ci_vals)

  summary_df <- data.frame(
    mean_cindex = mean(ci_vals),
    sd_cindex = sd(ci_vals),
    ci_95_lower = ci_95_lower,
    ci_95_upper = ci_95_upper,
    median_cindex = median(ci_vals),
    iqr_lower = quantile(ci_vals, 0.25),
    iqr_upper = quantile(ci_vals, 0.75),
    n_seeds = length(master_seeds),
    n_outer_folds = n_outer_folds,
    n_inner_folds = n_inner_folds,
    alpha = my_alpha,
    lambda_rule = lambda_rule,
    stringsAsFactors = FALSE
  )

  # Add time-dependent AUC summaries
  for (ac in auc_cols) {
    mean_col <- paste0("mean_", ac)
    if (mean_col %in% colnames(per_seed_df)) {
      vals <- per_seed_df[[mean_col]]
      summary_df[[paste0("mean_", ac)]] <- round(mean(vals, na.rm = TRUE), 4)
      summary_df[[paste0("sd_", ac)]] <- round(sd(vals, na.rm = TRUE), 4)
    }
  }

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
    message(sprintf("95%% CI:               [%.4f, %.4f]",
                    summary_df$ci_95_lower, summary_df$ci_95_upper))
    message(sprintf("SD across seeds:      %.4f", summary_df$sd_cindex))
    message(sprintf("Median outer C-index: %.4f", summary_df$median_cindex))
    message(sprintf("IQR:                  [%.4f, %.4f]",
                    summary_df$iqr_lower, summary_df$iqr_upper))

    # Time-dependent AUC summary
    if (length(auc_cols) > 0) {
      message("\nTime-dependent AUC (mean across seeds):")
      for (ac in auc_cols) {
        mean_col <- paste0("mean_", ac)
        sd_col <- paste0("sd_", ac)
        if (mean_col %in% colnames(summary_df)) {
          t_label <- gsub("auc_t", "t=", ac)
          message(sprintf("  %s: %.4f (SD %.4f)",
                          t_label, summary_df[[mean_col]], summary_df[[sd_col]]))
        }
      }
    }

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
#   eval_times = c(12, 24, 36),
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
#   eval_times = c(12, 24, 36),
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
#   eval_times = c(12, 24, 36),
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
#   eval_times = c(12, 24, 36),
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
