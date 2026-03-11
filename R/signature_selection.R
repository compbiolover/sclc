# Name: signature_selection.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Principled gene signature selection methods that go beyond
# arbitrary thresholds on selection frequency. Provides:
#   1. Stability selection with per-family error rate (PFER) control
#   2. Cross-pool consensus (convergent evidence across gene pools)
#   3. Univariate pre-screening with Benjamini-Hochberg correction
#   4. A combined selector that intersects all three

suppressMessages({
  library(glmnet)
  library(parallel)
  library(doParallel)
  library(survival)
  library(tidyverse)
})


# =========================================================================
# 1. STABILITY SELECTION (Meinshausen & Bühlmann, 2010)
# =========================================================================
#
# Instead of running lasso on the full dataset, we run it on many random
# subsamples of size n/2. For each subsample, we record which genes get
# nonzero coefficients. The selection probability for each gene is the
# fraction of subsamples in which it was selected.
#
# The key insight: Meinshausen & Bühlmann proved that for a threshold
# pi_thr and a regularization that selects at most q variables on average:
#
#   E[# falsely selected] <= q^2 / ((2*pi_thr - 1) * p)
#
# where p = total number of candidate genes. This gives you a principled
# threshold with error control, analogous to FDR for hypothesis testing.

#' Stability selection for Cox regression with PFER error control
#'
#' @param cox_df Data frame with vital.status, time, and gene columns.
#' @param candidate_genes Character vector of candidate gene names.
#' @param n_subsamples Number of subsamples to draw (default 200).
#' @param target_pfer Target per-family error rate — the expected number
#'   of falsely selected genes (default 1). Lower = more conservative.
#' @param my_alpha Elastic net alpha (default 1 = lasso).
#' @param gene_num Max genes to consider (default 500).
#' @param lambda_rule Lambda selection rule for glmnet (default "lambda.1se").
#' @param progress_free Logical; if TRUE, use time2 column.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A list with:
#'   - selected: character vector of genes passing the PFER-controlled threshold
#'   - selection_probs: named numeric vector of selection probabilities
#'   - threshold: the pi_thr used
#'   - pfer: the achieved PFER bound
#'   - all_probs: data frame of all genes and their selection probabilities
stability_selection <- function(cox_df,
                                candidate_genes,
                                n_subsamples = 200,
                                target_pfer = 1,
                                my_alpha = 1,
                                gene_num = 500,
                                lambda_rule = "lambda.1se",
                                progress_free = FALSE,
                                seed = 42,
                                verbose = TRUE) {
  set.seed(seed)

  # Prepare genes
  if (is.data.frame(candidate_genes)) {
    candidate_genes <- rownames(candidate_genes)
  }
  candidate_genes <- tolower(candidate_genes)
  candidate_genes <- intersect(candidate_genes, colnames(cox_df))
  candidate_genes <- head(candidate_genes, n = gene_num)
  p <- length(candidate_genes)

  if (p == 0) stop("No candidate genes found in cox_df.")

  time_col <- if (progress_free) "time2" else "time"
  n <- nrow(cox_df)
  subsample_size <- floor(n / 2)

  # Set up parallel
  num_cores <- parallel::detectCores()
  registerDoParallel(cores = num_cores)

  # Build formula
  formula_str <- paste0("~", paste0(make.names(candidate_genes), collapse = "+"))

  # Track selections
  selection_matrix <- matrix(0, nrow = n_subsamples, ncol = p)
  colnames(selection_matrix) <- make.names(candidate_genes)
  avg_selected <- numeric(n_subsamples)

  if (verbose) message("Running stability selection...")

  for (b in 1:n_subsamples) {
    if (verbose && b %% 50 == 0) {
      message(sprintf("  Subsample %d/%d", b, n_subsamples))
    }

    # Draw subsample of size n/2
    sub_idx <- sample(n, subsample_size, replace = FALSE)
    sub_data <- cox_df[sub_idx, ]
    colnames(sub_data) <- make.names(colnames(sub_data))

    sub_x <- tryCatch(
      model.matrix(as.formula(formula_str), sub_data),
      error = function(e) NULL
    )
    if (is.null(sub_x)) next

    sub_y <- Surv(
      time = sub_data[[make.names(time_col)]],
      event = sub_data$vital.status
    )

    cv_fit <- tryCatch(
      cv.glmnet(sub_x, sub_y, family = "cox", alpha = my_alpha,
                type.measure = "C", nfolds = 10, parallel = TRUE),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) next

    coefs <- coef(cv_fit, s = lambda_rule)
    selected <- rownames(coefs)[which(as.vector(coefs) != 0)]
    selection_matrix[b, colnames(selection_matrix) %in% selected] <- 1
    avg_selected[b] <- length(selected)
  }

  # Selection probabilities
  sel_probs <- colMeans(selection_matrix, na.rm = TRUE)
  names(sel_probs) <- candidate_genes

  # Average number of selected variables per subsample
  q <- mean(avg_selected[avg_selected > 0])

  # Compute PFER-controlled threshold
  # From Meinshausen & Bühlmann: E[V] <= q^2 / ((2*pi_thr - 1) * p)
  # Solving for pi_thr given target PFER:
  # target_pfer >= q^2 / ((2*pi_thr - 1) * p)
  # (2*pi_thr - 1) >= q^2 / (target_pfer * p)
  # pi_thr >= (1 + q^2 / (target_pfer * p)) / 2
  pi_thr <- (1 + q^2 / (target_pfer * p)) / 2
  pi_thr <- min(pi_thr, 0.99)  # Cap at 0.99
  pi_thr <- max(pi_thr, 0.6)   # Floor at 0.6 (minimum reasonable threshold)

  # Actual PFER at this threshold
  actual_pfer <- q^2 / ((2 * pi_thr - 1) * p)

  # Select genes above threshold
  selected_genes <- candidate_genes[sel_probs >= pi_thr]

  # Build results data frame
  all_probs_df <- data.frame(
    gene = candidate_genes,
    selection_prob = round(sel_probs, 3),
    above_threshold = sel_probs >= pi_thr,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(selection_prob))

  if (verbose) {
    message(sprintf("\n===== Stability Selection Results ====="))
    message(sprintf("Candidate genes (p):       %d", p))
    message(sprintf("Avg genes selected per sub: %.1f", q))
    message(sprintf("PFER-controlled threshold:  %.3f", pi_thr))
    message(sprintf("Achieved PFER bound:        %.2f (target: %.2f)",
                    actual_pfer, target_pfer))
    message(sprintf("Genes selected:            %d", length(selected_genes)))
    if (length(selected_genes) > 0) {
      message("Selected genes:")
      for (g in selected_genes) {
        message(sprintf("  %s (prob = %.3f)", g, sel_probs[g]))
      }
    } else {
      message("  No genes passed the threshold.")
      message("  Consider increasing target_pfer or gene_num.")
    }
  }

  return(list(
    selected = selected_genes,
    selection_probs = sort(sel_probs, decreasing = TRUE),
    threshold = pi_thr,
    pfer = actual_pfer,
    avg_q = q,
    all_probs = all_probs_df
  ))
}


# =========================================================================
# 2. CROSS-POOL CONSENSUS
# =========================================================================
#
# When you have multiple independent gene pools (SDE, MAD, miRNA targets),
# genes that appear as stable across multiple pools are more credible
# because they were identified through different biological rationales.

#' Find consensus genes across multiple stability selection results
#'
#' @param stability_results A named list of results from stability_selection()
#'   or nested_cv() (anything with a $selected or $gene_stability element).
#'   Names should describe the pool (e.g., "SDE", "MAD", "miRNA").
#' @param min_pools Minimum number of pools a gene must appear in to be
#'   considered consensus (default 2).
#' @param prob_threshold If using nested_cv results, the minimum selection
#'   frequency to consider a gene "selected" from that pool (default 0.5).
#' @param verbose Print results.
#'
#' @return A list with:
#'   - consensus_genes: character vector of genes appearing in >= min_pools
#'   - gene_pool_matrix: data frame showing which genes appear in which pools
#'   - pool_counts: how many pools each gene appears in
cross_pool_consensus <- function(stability_results,
                                 min_pools = 2,
                                 prob_threshold = 0.5,
                                 verbose = TRUE) {
  if (length(stability_results) < 2) {
    stop("Need at least 2 pools for cross-pool consensus.")
  }

  pool_names <- names(stability_results)
  if (is.null(pool_names)) {
    pool_names <- paste0("pool_", seq_along(stability_results))
  }

  # Extract selected genes from each pool
  pool_genes <- list()
  for (i in seq_along(stability_results)) {
    res <- stability_results[[i]]

    if (!is.null(res$selected)) {
      # From stability_selection()
      pool_genes[[pool_names[i]]] <- res$selected
    } else if (!is.null(res$gene_stability)) {
      # From nested_cv()
      stable <- res$gene_stability %>%
        filter(selection_frequency >= prob_threshold)
      pool_genes[[pool_names[i]]] <- stable$gene
    } else {
      warning(sprintf("Pool '%s' has no $selected or $gene_stability. Skipping.",
                      pool_names[i]))
    }
  }

  # Build gene-by-pool presence matrix
  all_genes <- unique(unlist(pool_genes))
  if (length(all_genes) == 0) {
    if (verbose) message("No genes selected in any pool.")
    return(list(consensus_genes = character(0),
                gene_pool_matrix = data.frame(),
                pool_counts = integer(0)))
  }

  pool_matrix <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  for (pool in names(pool_genes)) {
    pool_matrix[[pool]] <- pool_matrix$gene %in% pool_genes[[pool]]
  }
  pool_matrix$n_pools <- rowSums(pool_matrix[, -1, drop = FALSE])
  pool_matrix <- pool_matrix %>% arrange(desc(n_pools), gene)

  # Consensus genes
  consensus <- pool_matrix$gene[pool_matrix$n_pools >= min_pools]

  if (verbose) {
    message("\n===== Cross-Pool Consensus =====")
    message(sprintf("Pools: %s", paste(names(pool_genes), collapse = ", ")))
    message(sprintf("Genes per pool:"))
    for (pool in names(pool_genes)) {
      message(sprintf("  %s: %d genes", pool, length(pool_genes[[pool]])))
    }
    message(sprintf("\nConsensus (>= %d pools): %d genes", min_pools,
                    length(consensus)))
    if (length(consensus) > 0) {
      for (g in consensus) {
        pools_in <- names(pool_genes)[sapply(pool_genes, function(x) g %in% x)]
        message(sprintf("  %s (%s)", g, paste(pools_in, collapse = ", ")))
      }
    }
  }

  return(list(
    consensus_genes = consensus,
    gene_pool_matrix = pool_matrix,
    pool_counts = setNames(pool_matrix$n_pools, pool_matrix$gene)
  ))
}


# =========================================================================
# 3. UNIVARIATE PRE-SCREENING
# =========================================================================
#
# Before any multivariate model, test each gene individually with a
# univariate Cox model. Correct p-values with Benjamini-Hochberg.
# Only genes that are individually prognostic AND lasso-stable go
# into the final signature.

#' Univariate Cox screening with multiple testing correction
#'
#' @param cox_df Data frame with vital.status, time, and gene columns.
#' @param candidate_genes Character vector of gene names to screen.
#' @param fdr_threshold FDR threshold for significance (default 0.10).
#'   Use 0.10 rather than 0.05 because this is a pre-screen, not the
#'   final analysis, and you want to avoid being too aggressive.
#' @param gene_num Max genes to screen (default 500).
#' @param progress_free Logical; if TRUE, use time2 column.
#' @param verbose Print results.
#'
#' @return A list with:
#'   - significant_genes: genes passing the FDR threshold
#'   - results: data frame with gene, HR, p-value, BH-adjusted p-value
univariate_prescreen <- function(cox_df,
                                 candidate_genes,
                                 fdr_threshold = 0.10,
                                 gene_num = 500,
                                 progress_free = FALSE,
                                 verbose = TRUE) {
  # Prepare genes
  if (is.data.frame(candidate_genes)) {
    candidate_genes <- rownames(candidate_genes)
  }
  candidate_genes <- tolower(candidate_genes)
  candidate_genes <- intersect(candidate_genes, colnames(cox_df))
  candidate_genes <- head(candidate_genes, n = gene_num)

  if (length(candidate_genes) == 0) {
    stop("No candidate genes found in cox_df.")
  }

  time_col <- if (progress_free) "time2" else "time"

  if (verbose) {
    message(sprintf("Screening %d genes with univariate Cox models...",
                    length(candidate_genes)))
  }

  results_list <- list()

  for (gene in candidate_genes) {
    gene_clean <- make.names(gene)
    df_temp <- cox_df
    colnames(df_temp) <- make.names(colnames(df_temp))

    cox_result <- tryCatch({
      fit <- coxph(
        as.formula(sprintf("Surv(%s, vital.status) ~ %s",
                           make.names(time_col), gene_clean)),
        data = df_temp
      )
      s <- summary(fit)
      data.frame(
        gene = gene,
        hr = s$conf.int[1, 1],          # exp(coef)
        hr_lower = s$conf.int[1, 3],     # lower 95% CI
        hr_upper = s$conf.int[1, 4],     # upper 95% CI
        coef = s$coefficients[1, 1],     # raw coefficient
        se = s$coefficients[1, 3],       # SE
        z = s$coefficients[1, 4],        # z-score
        pval = s$coefficients[1, 5],     # p-value
        concordance = s$concordance[1],  # C-index
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(cox_result)) {
      results_list[[length(results_list) + 1]] <- cox_result
    }
  }

  results_df <- do.call(rbind, results_list)

  # BH correction
  results_df$pval_bh <- p.adjust(results_df$pval, method = "BH")
  results_df <- results_df %>% arrange(pval_bh)

  # Select significant genes
  sig_genes <- results_df$gene[results_df$pval_bh < fdr_threshold]

  if (verbose) {
    message(sprintf("\n===== Univariate Pre-Screening Results ====="))
    message(sprintf("Genes tested:         %d", nrow(results_df)))
    message(sprintf("FDR threshold:        %.2f", fdr_threshold))
    message(sprintf("Significant (BH):     %d", length(sig_genes)))
    if (length(sig_genes) > 0) {
      message("Top significant genes:")
      top <- head(results_df[results_df$pval_bh < fdr_threshold, ], 20)
      for (i in 1:nrow(top)) {
        message(sprintf("  %s: HR=%.3f, p_adj=%.4f, C=%.3f",
                        top$gene[i], top$hr[i], top$pval_bh[i],
                        top$concordance[i]))
      }
    } else {
      message("  No genes passed FDR threshold.")
      message("  Consider increasing fdr_threshold (e.g., 0.20).")
      # Show top genes anyway
      message("  Top 5 genes by raw p-value:")
      top5 <- head(results_df, 5)
      for (i in 1:nrow(top5)) {
        message(sprintf("    %s: HR=%.3f, p_raw=%.4f, p_adj=%.4f",
                        top5$gene[i], top5$hr[i], top5$pval[i], top5$pval_bh[i]))
      }
    }
  }

  return(list(
    significant_genes = sig_genes,
    results = results_df
  ))
}


# =========================================================================
# 4. COMBINED SIGNATURE SELECTOR
# =========================================================================
#
# Combines all three methods: a gene must be
#   (a) individually prognostic (univariate screen), AND
#   (b) lasso-stable (stability selection), AND
#   (c) optionally, selected across multiple gene pools (consensus)
#
# This triple filtering is conservative but produces signatures that
# are very difficult to dismiss as artifacts.

#' Select a gene signature using combined evidence
#'
#' @param cox_df Data frame with vital.status, time, and gene columns.
#' @param candidate_gene_pools Named list of character vectors, where each
#'   element is a pool of candidate genes (e.g., list(SDE = sde_genes,
#'   MAD = mad_genes, miRNA = mirna_genes)). If only one pool, pass it
#'   as list(pool = genes).
#' @param fdr_threshold FDR threshold for univariate screening (default 0.10).
#' @param target_pfer PFER target for stability selection (default 1).
#' @param n_subsamples Subsamples for stability selection (default 200).
#' @param require_consensus Logical; require genes to appear in >= 2 pools.
#'   Only applies when multiple pools are provided (default TRUE).
#' @param gene_num Max genes per pool (default 500).
#' @param lambda_rule Lambda rule for stability selection (default "lambda.1se").
#' @param progress_free Logical; if TRUE, use time2 column.
#' @param output_dir Directory to save results.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A list with:
#'   - final_signature: character vector of genes passing all filters
#'   - univariate_results: results from univariate pre-screening
#'   - stability_results: per-pool stability selection results
#'   - consensus_results: cross-pool consensus (if multiple pools)
#'   - filter_summary: data frame showing which genes passed which filters
select_signature <- function(cox_df,
                             candidate_gene_pools,
                             fdr_threshold = 0.10,
                             target_pfer = 1,
                             n_subsamples = 200,
                             require_consensus = TRUE,
                             gene_num = 500,
                             lambda_rule = "lambda.1se",
                             progress_free = FALSE,
                             output_dir = "Outputs/signature_selection",
                             seed = 42,
                             verbose = TRUE) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Handle single pool case
  if (is.character(candidate_gene_pools)) {
    candidate_gene_pools <- list(single_pool = candidate_gene_pools)
  }
  if (is.data.frame(candidate_gene_pools)) {
    candidate_gene_pools <- list(single_pool = rownames(candidate_gene_pools))
  }

  multi_pool <- length(candidate_gene_pools) > 1

  # ---- Step 1: Univariate pre-screen on all unique genes ----
  if (verbose) message("\n========== STEP 1: Univariate Pre-Screening ==========")

  all_candidates <- unique(tolower(unlist(candidate_gene_pools)))
  univar <- univariate_prescreen(
    cox_df = cox_df,
    candidate_genes = all_candidates,
    fdr_threshold = fdr_threshold,
    gene_num = length(all_candidates),
    progress_free = progress_free,
    verbose = verbose
  )
  univar_passing <- univar$significant_genes

  # ---- Step 2: Stability selection per pool ----
  if (verbose) message("\n========== STEP 2: Stability Selection ==========")

  stab_results <- list()
  for (pool_name in names(candidate_gene_pools)) {
    if (verbose) message(sprintf("\n--- Pool: %s ---", pool_name))
    stab_results[[pool_name]] <- stability_selection(
      cox_df = cox_df,
      candidate_genes = candidate_gene_pools[[pool_name]],
      n_subsamples = n_subsamples,
      target_pfer = target_pfer,
      gene_num = gene_num,
      lambda_rule = lambda_rule,
      progress_free = progress_free,
      seed = seed,
      verbose = verbose
    )
  }

  # All stability-selected genes
  all_stable <- unique(unlist(lapply(stab_results, function(x) x$selected)))

  # ---- Step 3: Cross-pool consensus (if multiple pools) ----
  consensus_result <- NULL
  consensus_genes <- all_stable  # default: all stable genes

  if (multi_pool && require_consensus) {
    if (verbose) message("\n========== STEP 3: Cross-Pool Consensus ==========")
    consensus_result <- cross_pool_consensus(
      stability_results = stab_results,
      min_pools = 2,
      verbose = verbose
    )
    consensus_genes <- consensus_result$consensus_genes
  }

  # ---- Step 4: Intersect all filters ----
  if (verbose) message("\n========== STEP 4: Final Intersection ==========")

  # Final signature = univariate significant AND stability-selected
  # AND (consensus if multi-pool)
  if (multi_pool && require_consensus) {
    final_sig <- intersect(univar_passing, consensus_genes)
  } else {
    final_sig <- intersect(univar_passing, all_stable)
  }

  # Build summary table showing all genes and which filters they passed
  all_genes_considered <- unique(c(univar_passing, all_stable, consensus_genes))
  filter_summary <- data.frame(
    gene = all_genes_considered,
    univariate_significant = all_genes_considered %in% univar_passing,
    stability_selected = all_genes_considered %in% all_stable,
    stringsAsFactors = FALSE
  )

  if (multi_pool && require_consensus) {
    filter_summary$cross_pool_consensus <- all_genes_considered %in% consensus_genes
  }

  filter_summary$in_final_signature <- all_genes_considered %in% final_sig

  # Add univariate stats
  filter_summary <- filter_summary %>%
    left_join(
      univar$results %>% select(gene, hr, pval_bh, concordance),
      by = "gene"
    ) %>%
    arrange(desc(in_final_signature), pval_bh)

  # Add stability probabilities (max across pools)
  filter_summary$max_stability_prob <- sapply(filter_summary$gene, function(g) {
    probs <- sapply(stab_results, function(s) {
      if (g %in% names(s$selection_probs)) s$selection_probs[g] else 0
    })
    max(probs)
  })

  if (verbose) {
    message(sprintf("\nUnivariate significant:   %d genes", length(univar_passing)))
    message(sprintf("Stability selected:       %d genes", length(all_stable)))
    if (multi_pool && require_consensus) {
      message(sprintf("Cross-pool consensus:     %d genes", length(consensus_genes)))
    }
    message(sprintf("Final signature:          %d genes", length(final_sig)))

    if (length(final_sig) > 0) {
      message("\nFinal signature genes:")
      for (g in final_sig) {
        row <- filter_summary[filter_summary$gene == g, ]
        message(sprintf("  %s: HR=%.3f, FDR=%.4f, stability=%.3f, C=%.3f",
                        g, row$hr, row$pval_bh, row$max_stability_prob,
                        row$concordance))
      }
    } else {
      message("\nNo genes passed all filters.")
      message("Consider relaxing fdr_threshold or target_pfer.")
      message("\nGenes passing at least one filter:")
      for (i in 1:min(10, nrow(filter_summary))) {
        row <- filter_summary[i, ]
        flags <- c()
        if (row$univariate_significant) flags <- c(flags, "U")
        if (row$stability_selected) flags <- c(flags, "S")
        if (multi_pool && require_consensus &&
            !is.null(row$cross_pool_consensus) && row$cross_pool_consensus) {
          flags <- c(flags, "C")
        }
        message(sprintf("  %s [%s]: HR=%.3f, FDR=%.4f, stab=%.3f",
                        row$gene, paste(flags, collapse = "+"),
                        row$hr, row$pval_bh, row$max_stability_prob))
      }
    }
  }

  # Save results
  write.csv(filter_summary,
            file.path(output_dir, "filter_summary.csv"),
            row.names = FALSE)
  write.csv(univar$results,
            file.path(output_dir, "univariate_results.csv"),
            row.names = FALSE)
  for (pool_name in names(stab_results)) {
    write.csv(stab_results[[pool_name]]$all_probs,
              file.path(output_dir,
                        sprintf("stability_%s.csv", pool_name)),
              row.names = FALSE)
  }

  if (verbose) {
    message(sprintf("\nResults saved to: %s", output_dir))
  }

  return(list(
    final_signature = final_sig,
    univariate_results = univar,
    stability_results = stab_results,
    consensus_results = consensus_result,
    filter_summary = filter_summary
  ))
}


# ===========================================================================
# Example usage (uncomment and modify paths to run)
# ===========================================================================
#
# scl_common <- readRDS("Outputs/scl_common.rds")
# sde_genes <- readRDS("Outputs/sde_genes_new_random_seed.rds")
# mad_genes <- readRDS("Outputs/mad_genes_new_random_seed.rds")
# mirna_genes <- readRDS("Outputs/mirna/rds/mirna_genes_200_mirnas_5_targets_up_.rds")
#
#
# # ============================
# # OPTION A: Full pipeline with all three pools
# # ============================
# sig <- select_signature(
#   cox_df = scl_common,
#   candidate_gene_pools = list(
#     SDE = sde_genes,
#     MAD = mad_genes,
#     miRNA = mirna_genes
#   ),
#   fdr_threshold = 0.10,
#   target_pfer = 1,
#   n_subsamples = 200,
#   require_consensus = TRUE,
#   gene_num = 500,
#   output_dir = "Outputs/signature_selection/all_pools"
# )
#
# # Your final signature:
# print(sig$final_signature)
#
# # See which genes passed which filters:
# print(sig$filter_summary)
#
#
# # ============================
# # OPTION B: Single pool (e.g., just SDE)
# # ============================
# sig_sde <- select_signature(
#   cox_df = scl_common,
#   candidate_gene_pools = list(SDE = sde_genes),
#   fdr_threshold = 0.10,
#   target_pfer = 1,
#   n_subsamples = 200,
#   require_consensus = FALSE,  # only 1 pool
#   gene_num = 500,
#   output_dir = "Outputs/signature_selection/sde_only"
# )
#
#
# # ============================
# # OPTION C: Run each step individually
# # ============================
#
# # Step 1: Univariate screen
# univar <- univariate_prescreen(
#   cox_df = scl_common,
#   candidate_genes = sde_genes,
#   fdr_threshold = 0.10
# )
#
# # Step 2: Stability selection
# stab <- stability_selection(
#   cox_df = scl_common,
#   candidate_genes = sde_genes,
#   n_subsamples = 200,
#   target_pfer = 1
# )
#
# # Step 3: Intersect manually
# robust_genes <- intersect(univar$significant_genes, stab$selected)
# print(robust_genes)
#
#
# # ============================
# # THEN: Validate with nested CV
# # ============================
# # Take the final signature and run nested_cv() to get the honest
# # C-index, then internal_validation() for the supporting figures:
# #
# # source("R/nested_cv.R")
# # source("R/internal_validation.R")
# #
# # ncv_results <- nested_cv(
# #   cox_df = scl_common,
# #   candidate_genes = sig$final_signature,
# #   my_alpha = 0,
# #   gene_num = length(sig$final_signature),
# #   master_seeds = 1:100,
# #   output_dir = "Outputs/nested_cv/final_signature"
# # )
# #
# # val_results <- internal_validation(
# #   cox_df = scl_common,
# #   signature_genes = sig$final_signature,
# #   output_dir = "Outputs/internal_validation/final_signature"
# # )
