# Name: run_discovery.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Self-contained runner script for the full discovery pipeline.
#   1. Load data
#   2. Run nested CV in discovery mode (lasso) across all gene pools
#   3. Identify stable genes from each pool
#   4. Run validation mode (ridge) on the stable genes
#   5. Run internal validation suite for publication figures
#   6. Compare all methods and print final results
#
# Usage: source("R/run_discovery.R") from the project root, or run
#   Rscript R/run_discovery.R from the command line.

# =========================================================================
# 0. CONFIGURATION — edit these paths and parameters as needed
# =========================================================================

config <- list(
  # --- Data paths ---
  scl_common_path = "Outputs/scl_common.rds",
  sde_genes_path  = "Outputs/sde_genes_new_random_seed.rds",
  mad_genes_path  = "Outputs/mad_genes_new_random_seed.rds",
  mirna_genes_path = "Outputs/mirna/rds/mirna_genes_200_mirnas_5_targets_up_.rds",

  # --- Nested CV parameters ---
  n_outer_folds = 5,
  n_inner_folds = 10,
  discovery_seeds = 1:50,       # 50 seeds for discovery (faster)
  validation_seeds = 1:100,     # 100 seeds for final validation (thorough)
  gene_num_discovery = 500,     # large pool for lasso to select from
  lambda_rule = "lambda.1se",   # conservative lambda (recommended)

  # --- Time-dependent AUC evaluation times ---
  eval_times = c(12, 24, 36),   # 1-year, 2-year, 3-year AUC

  # --- Gene stability threshold ---
  stability_threshold = 0.5,    # genes selected in >50% of fits

  # --- Output directories ---
  output_base = "Outputs/nested_cv"
)


# =========================================================================
# 1. SETUP
# =========================================================================

cat("============================================================\n")
cat("  SCLC Gene Signature Discovery Pipeline\n")
cat("============================================================\n\n")

# Source the pipeline functions
source("R/nested_cv.R")
source("R/internal_validation.R")
source("R/signature_selection.R")

# Load data
cat("Loading data...\n")

if (!file.exists(config$scl_common_path)) {
  stop(sprintf("Data file not found: %s\nRun main.Rmd first to generate scl_common.",
               config$scl_common_path))
}

scl_common <- readRDS(config$scl_common_path)
cat(sprintf("  Loaded scl_common: %d patients, %d columns\n",
            nrow(scl_common), ncol(scl_common)))

# Load gene pools (skip any that don't exist)
gene_pools <- list()

if (file.exists(config$sde_genes_path)) {
  gene_pools$SDE <- readRDS(config$sde_genes_path)
  cat(sprintf("  Loaded SDE genes: %d candidates\n", length(gene_pools$SDE)))
}
if (file.exists(config$mad_genes_path)) {
  gene_pools$MAD <- readRDS(config$mad_genes_path)
  cat(sprintf("  Loaded MAD genes: %d candidates\n", length(gene_pools$MAD)))
}
if (file.exists(config$mirna_genes_path)) {
  gene_pools$miRNA <- readRDS(config$mirna_genes_path)
  cat(sprintf("  Loaded miRNA genes: %d candidates\n", length(gene_pools$miRNA)))
}

if (length(gene_pools) == 0) {
  stop("No gene pool files found. Check paths in the config section.")
}


# =========================================================================
# 2. DISCOVERY MODE — lasso nested CV on each gene pool
# =========================================================================

cat("\n============================================================\n")
cat("  PHASE 1: Discovery (lasso nested CV per pool)\n")
cat("============================================================\n")

discovery_results <- list()

for (pool_name in names(gene_pools)) {
  cat(sprintf("\n--- Discovering from %s pool ---\n", pool_name))

  out_dir <- file.path(config$output_base, paste0(tolower(pool_name), "_discovery"))

  discovery_results[[pool_name]] <- nested_cv(
    cox_df           = scl_common,
    candidate_genes  = gene_pools[[pool_name]],
    n_outer_folds    = config$n_outer_folds,
    n_inner_folds    = config$n_inner_folds,
    master_seeds     = config$discovery_seeds,
    my_alpha         = 1,       # lasso for feature selection
    gene_num         = config$gene_num_discovery,
    lambda_rule      = config$lambda_rule,
    eval_times       = config$eval_times,
    output_dir       = out_dir
  )
}


# =========================================================================
# 3. IDENTIFY STABLE GENES
# =========================================================================

cat("\n============================================================\n")
cat("  PHASE 2: Gene Stability Analysis\n")
cat("============================================================\n")

stable_genes_per_pool <- list()

for (pool_name in names(discovery_results)) {
  res <- discovery_results[[pool_name]]
  stable <- res$gene_stability %>%
    dplyr::filter(selection_frequency >= config$stability_threshold)

  stable_genes_per_pool[[pool_name]] <- stable$gene

  cat(sprintf("\n%s pool: %d stable genes (>%.0f%% frequency)\n",
              pool_name, nrow(stable), config$stability_threshold * 100))
  if (nrow(stable) > 0) {
    for (i in 1:nrow(stable)) {
      cat(sprintf("  %s: %.1f%%\n", stable$gene[i],
                  stable$selection_frequency[i] * 100))
    }
  }
}

# Union of all stable genes across pools
all_stable <- unique(unlist(stable_genes_per_pool))
cat(sprintf("\nUnion of stable genes across all pools: %d\n", length(all_stable)))
cat(sprintf("  %s\n", paste(all_stable, collapse = ", ")))

# Cross-pool consensus: genes appearing in 2+ pools
if (length(gene_pools) >= 2) {
  consensus <- cross_pool_consensus(
    stability_results = discovery_results,
    min_pools = 2,
    prob_threshold = config$stability_threshold
  )
  consensus_genes <- consensus$consensus_genes
  cat(sprintf("\nCross-pool consensus genes (2+ pools): %d\n",
              length(consensus_genes)))
  if (length(consensus_genes) > 0) {
    cat(sprintf("  %s\n", paste(consensus_genes, collapse = ", ")))
  }
} else {
  consensus_genes <- all_stable
}


# =========================================================================
# 4. VALIDATION MODE — ridge nested CV on the discovered signature
# =========================================================================

cat("\n============================================================\n")
cat("  PHASE 3: Validation (ridge nested CV on stable genes)\n")
cat("============================================================\n")

# Validate the union of stable genes
if (length(all_stable) > 0) {
  cat(sprintf("\nValidating %d-gene signature: %s\n",
              length(all_stable), paste(all_stable, collapse = ", ")))

  validation_result <- nested_cv(
    cox_df           = scl_common,
    candidate_genes  = all_stable,
    n_outer_folds    = config$n_outer_folds,
    n_inner_folds    = config$n_inner_folds,
    master_seeds     = config$validation_seeds,
    my_alpha         = 0,       # ridge — genes are pre-selected
    gene_num         = length(all_stable),
    lambda_rule      = config$lambda_rule,
    eval_times       = config$eval_times,
    output_dir       = file.path(config$output_base, "stable_genes_validation")
  )
} else {
  cat("\nNo stable genes found. Consider:\n")
  cat("  - Lowering stability_threshold (e.g., 0.3)\n")
  cat("  - Increasing discovery_seeds (e.g., 1:100)\n")
  cat("  - Checking that gene names match between pools and scl_common\n")
  validation_result <- NULL
}

# Also validate existing 4-gene signature for comparison
cat("\nValidating existing 4-gene signature (TFRC, FAM83F, DLK1, GNG13)...\n")

four_gene_result <- nested_cv(
  cox_df           = scl_common,
  candidate_genes  = c("tfrc", "fam83f", "dlk1", "gng13"),
  n_outer_folds    = config$n_outer_folds,
  n_inner_folds    = config$n_inner_folds,
  master_seeds     = config$validation_seeds,
  my_alpha         = 0,
  gene_num         = 4,
  lambda_rule      = config$lambda_rule,
  eval_times       = config$eval_times,
  output_dir       = file.path(config$output_base, "four_gene_validation")
)


# =========================================================================
# 5. INTERNAL VALIDATION — publication figures for best signature
# =========================================================================

cat("\n============================================================\n")
cat("  PHASE 4: Internal Validation (figures & stats)\n")
cat("============================================================\n")

if (!is.null(validation_result) && length(all_stable) > 0) {
  cat(sprintf("\nRunning internal validation for %d-gene discovered signature...\n",
              length(all_stable)))

  iv_discovered <- internal_validation(
    cox_df           = scl_common,
    signature_genes  = all_stable,
    n_permutations   = 1000,
    n_boot           = 1000,
    output_dir       = "Outputs/internal_validation/discovered_signature"
  )
}

cat("\nRunning internal validation for 4-gene signature...\n")

iv_four_gene <- internal_validation(
  cox_df           = scl_common,
  signature_genes  = c("tfrc", "fam83f", "dlk1", "gng13"),
  n_permutations   = 1000,
  n_boot           = 1000,
  output_dir       = "Outputs/internal_validation/four_gene_signature"
)


# =========================================================================
# 6. COMPARISON TABLE
# =========================================================================

cat("\n============================================================\n")
cat("  FINAL RESULTS: Method Comparison\n")
cat("============================================================\n\n")

# Build comparison from all discovery runs + validation runs
comparison_rows <- list()

# Discovery C-indices (honest, from lasso nested CV)
for (pool_name in names(discovery_results)) {
  s <- discovery_results[[pool_name]]$summary
  comparison_rows[[length(comparison_rows) + 1]] <- data.frame(
    method    = sprintf("%s discovery (lasso)", pool_name),
    mean_cindex = round(s$mean_cindex, 4),
    sd_cindex   = round(s$sd_cindex, 4),
    median_cindex = round(s$median_cindex, 4),
    n_stable_genes = length(stable_genes_per_pool[[pool_name]]),
    stringsAsFactors = FALSE
  )
}

# Validation C-indices
if (!is.null(validation_result)) {
  s <- validation_result$summary
  comparison_rows[[length(comparison_rows) + 1]] <- data.frame(
    method    = sprintf("Discovered %d-gene (ridge)", length(all_stable)),
    mean_cindex = round(s$mean_cindex, 4),
    sd_cindex   = round(s$sd_cindex, 4),
    median_cindex = round(s$median_cindex, 4),
    n_stable_genes = length(all_stable),
    stringsAsFactors = FALSE
  )
}

s <- four_gene_result$summary
comparison_rows[[length(comparison_rows) + 1]] <- data.frame(
  method    = "4-gene TFRC/FAM83F/DLK1/GNG13 (ridge)",
  mean_cindex = round(s$mean_cindex, 4),
  sd_cindex   = round(s$sd_cindex, 4),
  median_cindex = round(s$median_cindex, 4),
  n_stable_genes = 4,
  stringsAsFactors = FALSE
)

comparison_df <- do.call(rbind, comparison_rows)

cat("Method Comparison (Nested CV C-index):\n")
cat("--------------------------------------------------------------\n")
print(comparison_df, row.names = FALSE, right = FALSE)

# Save comparison
write.csv(comparison_df,
          file.path(config$output_base, "method_comparison.csv"),
          row.names = FALSE)

cat(sprintf("\n\nAll results saved under: %s/\n", config$output_base))
cat("============================================================\n")
cat("  Discovery pipeline complete.\n")
cat("============================================================\n")
