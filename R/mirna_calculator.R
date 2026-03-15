# Name: mirna_calculator.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Multi-database miRNA target gene ranking with confidence weighting.
#
# Approach:
#   1. Identify cancer-dysregulated miRNAs via:
#      - PREFERRED: GEO differential expression (e.g., GSE19945 for SCLC)
#        with limma, using both up- and down-regulated miRNAs, weighted by
#        fold change magnitude
#      - FALLBACK: dbDEMC 3.0 + miRmap intersection (generic "lung cancer",
#        both up- and down-regulated)
#   2. Query three target prediction databases:
#      - TargetScan 8.0 (CNN-based biochemical model; McGeary et al., 2019)
#        Preferred: load bulk download files locally (no web scraping)
#        Fallback: query via hoardeR if bulk files unavailable
#      - miRDB (target scores; Wong & Wang, 2015)
#      - miRTarBase (experimentally validated; Huang et al., 2022)
#   3. Apply literature-established confidence thresholds:
#      - TargetScan: cumulative weighted context++ score <= -0.2
#      - miRDB: target score >= 80
#      - miRTarBase: strong experimental evidence (reporter assay, western blot, qPCR)
#   4. Require multi-database consensus: interaction must be supported by >= 2
#      databases, OR have strong experimental validation in miRTarBase alone
#   5. Weight gene scores by prediction confidence rather than binary counts
#
# Required local data files:
#   - Data/mirna_data/mirmap_mirnas.csv         (miRmap miRNA list)
#   - Data/mirna_data/dbdemc_3.0_high.txt       (dbDEMC upregulated miRNAs)
#   - Data/mirna_data/dbdemc_3.0_low.txt        (dbDEMC downregulated miRNAs)
#   - Data/mirna_data/miRDB_v6.0_prediction_result.txt (miRDB predictions)
#   - Data/mirna_data/miRTarBase.csv            (miRTarBase validated targets)
#
# TargetScan 8.0 bulk download (preferred over web scraping):
#   Download from: https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
#   - Data/mirna_data/Summary_Counts.default_predictions.txt
#       (default CNN-based predictions per miRNA family-gene pair)
#   - Data/mirna_data/miR_Family_Info.txt
#       (maps miRNA family IDs to individual mature miRNA names)
#   - Data/mirna_data/Predicted_Targets_Context_Scores.default_predictions.txt
#       (optional: per-site context++ scores for additional confidence weighting)
#
# References:
#   McGeary et al. (2019). Science. doi:10.1126/science.aav1741
#   Agarwal et al. (2015). eLife. doi:10.7554/eLife.05005
#   Wong & Wang (2015). Nucleic Acids Res. doi:10.1093/nar/gku1104
#   Huang et al. (2022). Nucleic Acids Res. doi:10.1093/nar/gkab1079


# =============================================================================
# Input validation
# =============================================================================

#' Validate mirna_calculator inputs
#'
#' @param mirmap_path Path to miRmap miRNAs CSV file.
#' @param dbdemc_high_path Path to dbDEMC upregulated miRNAs file.
#' @param dbdemc_low_path Path to dbDEMC downregulated miRNAs file (or NULL).
#' @param mirdb_path Path to miRDB predictions file (or NULL to skip).
#' @param mirtarbase_path Path to miRTarBase CSV file (or NULL to skip).
#' @param targetscan_path Path to TargetScan 8.0 bulk predictions file (or NULL
#'   to fall back to hoardeR web queries).
#' @param ts_family_path Path to TargetScan miR_Family_Info.txt file (or NULL).
#' @param cancer_type Cancer type string for dbDEMC filtering.
#' @param ts_context_threshold TargetScan context++ score threshold.
#' @param mirdb_score_threshold miRDB target score threshold.
#' @param min_databases Minimum number of databases for consensus.
#'
#' @return Invisible TRUE if all checks pass; stops with informative error otherwise.
validate_mirna_inputs <- function(mirmap_path,
                                  dbdemc_high_path,
                                  dbdemc_low_path = NULL,
                                  mirdb_path,
                                  mirtarbase_path,
                                  targetscan_path = NULL,
                                  ts_family_path = NULL,
                                  cancer_type,
                                  ts_context_threshold,
                                  mirdb_score_threshold,
                                  min_databases) {
  if (!is.character(mirmap_path) || length(mirmap_path) != 1) {
    stop("mirmap_path must be a single character string.")
  }
  if (!file.exists(mirmap_path)) {
    stop(paste0("miRmap file not found: ", mirmap_path))
  }

  if (!is.character(dbdemc_high_path) || length(dbdemc_high_path) != 1) {
    stop("dbdemc_high_path must be a single character string.")
  }
  if (!file.exists(dbdemc_high_path)) {
    stop(paste0("dbDEMC upregulated file not found: ", dbdemc_high_path))
  }

  if (!is.null(dbdemc_low_path)) {
    if (!is.character(dbdemc_low_path) || length(dbdemc_low_path) != 1) {
      stop("dbdemc_low_path must be a single character string or NULL.")
    }
    if (!file.exists(dbdemc_low_path)) {
      stop(paste0("dbDEMC downregulated file not found: ", dbdemc_low_path))
    }
  }

  if (!is.null(mirdb_path)) {
    if (!is.character(mirdb_path) || length(mirdb_path) != 1) {
      stop("mirdb_path must be a single character string or NULL.")
    }
    if (!file.exists(mirdb_path)) {
      stop(paste0("miRDB file not found: ", mirdb_path))
    }
  }

  if (!is.null(mirtarbase_path)) {
    if (!is.character(mirtarbase_path) || length(mirtarbase_path) != 1) {
      stop("mirtarbase_path must be a single character string or NULL.")
    }
    if (!file.exists(mirtarbase_path)) {
      stop(paste0("miRTarBase file not found: ", mirtarbase_path))
    }
  }

  if (!is.null(targetscan_path)) {
    if (!is.character(targetscan_path) || length(targetscan_path) != 1) {
      stop("targetscan_path must be a single character string or NULL.")
    }
    if (!file.exists(targetscan_path)) {
      stop(paste0("TargetScan bulk file not found: ", targetscan_path))
    }
  }

  if (!is.null(ts_family_path)) {
    if (!is.character(ts_family_path) || length(ts_family_path) != 1) {
      stop("ts_family_path must be a single character string or NULL.")
    }
    if (!file.exists(ts_family_path)) {
      stop(paste0("TargetScan miR family file not found: ", ts_family_path))
    }
  }

  if (!is.character(cancer_type) || length(cancer_type) != 1) {
    stop("cancer_type must be a single character string.")
  }

  if (!is.numeric(ts_context_threshold) || length(ts_context_threshold) != 1) {
    stop("ts_context_threshold must be a single numeric value.")
  }

  if (!is.numeric(mirdb_score_threshold) || length(mirdb_score_threshold) != 1 ||
      mirdb_score_threshold < 0 || mirdb_score_threshold > 100) {
    stop("mirdb_score_threshold must be a single numeric value between 0 and 100.")
  }

  if (!is.numeric(min_databases) || length(min_databases) != 1 ||
      min_databases < 1 || min_databases > 3) {
    stop("min_databases must be 1, 2, or 3.")
  }

  invisible(TRUE)
}


# =============================================================================
# miRNA filtering
# =============================================================================

#' Identify differentially expressed miRNAs from GEO expression data
#'
#' Downloads a miRNA expression dataset from GEO (default: GSE19945, which
#' contains 35 SCLC vs 8 normal lung samples) and performs differential
#' expression analysis with limma. This replaces the dbDEMC approach of using
#' a generic "lung cancer" label with actual SCLC-specific differential
#' expression, complete with fold changes and FDR-adjusted p-values.
#'
#' Both upregulated and downregulated miRNAs are returned, since:
#' - Upregulated miRNAs repress tumor suppressor targets
#' - Downregulated miRNAs de-repress oncogene targets
#' Both directions yield biologically relevant target genes.
#'
#' @param geo_accession GEO Series accession (default "GSE19945" for SCLC).
#' @param case_group Character string or regex pattern identifying case
#'   (disease) samples in the sample titles/characteristics.
#'   Default "SCLC|Small.cell lung cancer" to match GSE19945 metadata
#'   (which uses "tissue: Small-cell lung cancer").
#' @param control_group Character string or regex identifying control samples.
#'   Default "Normal|normal|Adjacent normal".
#' @param fdr_threshold FDR (Benjamini-Hochberg) threshold for significance.
#'   Default 0.05.
#' @param min_log2fc Minimum absolute log2 fold change. Default 1.0 (2-fold).
#' @param mirna_remove Character vector of miRNAs to exclude.
#' @param max_mirnas Maximum number of DE miRNAs to retain (ordered by
#'   significance, not fold change). Default 808.
#' @param de_results_path Optional path to a pre-computed DE results CSV file
#'   (columns: mirna, log2fc, adj_p_value). If provided, GEO download and
#'   limma analysis are skipped. This allows you to run the analysis once,
#'   save the results, and reuse them without re-downloading.
#' @param save_de_results Logical; if TRUE, save the full DE results table to
#'   output_path for future reuse via de_results_path.
#' @param output_path Path to save DE results CSV.
#' @param verbose Print diagnostic messages.
#'
#' @return A data frame with columns:
#'   - mirna: miRNA name (miRBase ID)
#'   - log2fc: log2 fold change (positive = upregulated in disease)
#'   - adj_p_value: BH-adjusted p-value
#'   - direction: "up" or "down"
#'   - abs_log2fc: absolute log2 fold change (for weighting)
#'   Sorted by adj_p_value (most significant first).
filter_mirnas_from_geo <- function(geo_accession = "GSE19945",
                                   case_group = "SCLC|Small.cell lung cancer",
                                   control_group = "Normal|normal|Adjacent normal",
                                   fdr_threshold = 0.05,
                                   min_log2fc = 1.0,
                                   mirna_remove = character(0),
                                   max_mirnas = 808,
                                   de_results_path = NULL,
                                   save_de_results = TRUE,
                                   output_path = "Outputs/mirna/sclc_de_mirnas.csv",
                                   verbose = TRUE) {
  # If pre-computed results exist, load and filter them
  if (!is.null(de_results_path) && file.exists(de_results_path)) {
    if (verbose) message(sprintf("Loading pre-computed DE results from %s",
                                  de_results_path))
    de_results <- read.csv(de_results_path, stringsAsFactors = FALSE)

    if (!all(c("mirna", "log2fc", "adj_p_value") %in% colnames(de_results))) {
      stop("DE results file must have columns: mirna, log2fc, adj_p_value")
    }

    de_sig <- de_results %>%
      dplyr::filter(adj_p_value < fdr_threshold & abs(log2fc) >= min_log2fc) %>%
      dplyr::mutate(
        direction = ifelse(log2fc > 0, "up", "down"),
        abs_log2fc = abs(log2fc)
      ) %>%
      dplyr::filter(!mirna %in% mirna_remove) %>%
      dplyr::arrange(adj_p_value)

    if (nrow(de_sig) > max_mirnas) {
      de_sig <- de_sig[1:max_mirnas, ]
    }

    if (verbose) {
      message(sprintf("DE miRNAs passing filters: %d (FDR < %.2f, |log2FC| >= %.1f)",
                      nrow(de_sig), fdr_threshold, min_log2fc))
      message(sprintf("  Upregulated: %d", sum(de_sig$direction == "up")))
      message(sprintf("  Downregulated: %d", sum(de_sig$direction == "down")))
    }

    return(de_sig)
  }

  # Otherwise, run fresh GEO download + limma analysis
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required. Install with: BiocManager::install('GEOquery')")
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required. Install with: BiocManager::install('limma')")
  }

  if (verbose) message(sprintf("Downloading %s from GEO...", geo_accession))
  gse <- GEOquery::getGEO(geo_accession, GSEMatrix = TRUE, getGPL = FALSE)

  if (is.list(gse)) {
    gse <- gse[[1]]
  }

  # Extract expression matrix and sample info
  expr_mat <- Biobase::exprs(gse)
  pheno <- Biobase::pData(gse)

  if (verbose) {
    message(sprintf("  Samples: %d, Features: %d", ncol(expr_mat), nrow(expr_mat)))
  }

  # Identify case and control samples from phenotype data
  # Search across all character columns for the group labels
  sample_labels <- rep(NA_character_, nrow(pheno))
  search_cols <- sapply(pheno, is.character)
  search_text <- apply(pheno[, search_cols, drop = FALSE], 1, paste, collapse = " ")

  sample_labels[grepl(case_group, search_text, ignore.case = TRUE)] <- "case"
  sample_labels[grepl(control_group, search_text, ignore.case = TRUE)] <- "control"

  # Remove ambiguous samples (matched both or neither)
  case_and_control <- sample_labels == "case" &
    grepl(control_group, search_text, ignore.case = TRUE)
  sample_labels[case_and_control] <- NA

  if (sum(sample_labels == "case", na.rm = TRUE) == 0) {
    stop(sprintf("No case samples matched pattern '%s'. Check case_group.",
                 case_group))
  }
  if (sum(sample_labels == "control", na.rm = TRUE) == 0) {
    stop(sprintf("No control samples matched pattern '%s'. Check control_group.",
                 control_group))
  }

  # Filter to labeled samples only
  keep <- !is.na(sample_labels)
  expr_mat <- expr_mat[, keep]
  sample_labels <- sample_labels[keep]

  if (verbose) {
    message(sprintf("  Case samples: %d, Control samples: %d",
                    sum(sample_labels == "case"),
                    sum(sample_labels == "control")))
  }

  # limma differential expression
  group <- factor(sample_labels, levels = c("control", "case"))
  design <- model.matrix(~ group)

  fit <- limma::lmFit(expr_mat, design)
  fit <- limma::eBayes(fit)
  results <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")

  # Extract miRNA names from feature annotation
  feat_data <- Biobase::fData(gse)
  if (!is.null(feat_data) && nrow(feat_data) > 0) {
    # Look for miRNA name column
    mirna_name_col <- grep("^miRNA|^ID$|^NAME$|^Gene\\.?Symbol|^Symbol",
                            colnames(feat_data), value = TRUE, ignore.case = TRUE)
    if (length(mirna_name_col) > 0) {
      # Use annotation to get miRNA names, matching by rownames
      common_features <- intersect(rownames(results), rownames(feat_data))
      results <- results[common_features, ]
      results$mirna <- feat_data[common_features, mirna_name_col[1]]
    } else {
      results$mirna <- rownames(results)
    }
  } else {
    results$mirna <- rownames(results)
  }

  # Build clean DE results table
  de_all <- data.frame(
    mirna = results$mirna,
    log2fc = results$logFC,
    adj_p_value = results$adj.P.Val,
    avg_expr = results$AveExpr,
    t_stat = results$t,
    stringsAsFactors = FALSE
  )

  # Remove rows with missing miRNA names
  de_all <- de_all[!is.na(de_all$mirna) & de_all$mirna != "", ]

  # Save full results for reuse
  if (save_de_results && !is.null(output_path)) {
    out_dir <- dirname(output_path)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    write.csv(de_all, output_path, row.names = FALSE)
    if (verbose) message(sprintf("  Full DE results saved to: %s", output_path))
  }

  # Filter to significant DE miRNAs
  de_sig <- de_all %>%
    dplyr::filter(adj_p_value < fdr_threshold & abs(log2fc) >= min_log2fc) %>%
    dplyr::mutate(
      direction = ifelse(log2fc > 0, "up", "down"),
      abs_log2fc = abs(log2fc)
    ) %>%
    dplyr::filter(!mirna %in% mirna_remove) %>%
    dplyr::arrange(adj_p_value)

  if (nrow(de_sig) > max_mirnas) {
    de_sig <- de_sig[1:max_mirnas, ]
  }

  if (verbose) {
    message(sprintf("\n===== SCLC Differential Expression Results ====="))
    message(sprintf("Total miRNAs tested: %d", nrow(de_all)))
    message(sprintf("DE miRNAs (FDR < %.2f, |log2FC| >= %.1f): %d",
                    fdr_threshold, min_log2fc, nrow(de_sig)))
    message(sprintf("  Upregulated in %s: %d", case_group,
                    sum(de_sig$direction == "up")))
    message(sprintf("  Downregulated in %s: %d", case_group,
                    sum(de_sig$direction == "down")))
    if (nrow(de_sig) > 0) {
      message("Top 10 DE miRNAs:")
      top <- head(de_sig, 10)
      for (i in seq_len(nrow(top))) {
        message(sprintf("  %s: log2FC=%.2f, FDR=%.1e (%s)",
                        top$mirna[i], top$log2fc[i], top$adj_p_value[i],
                        top$direction[i]))
      }
    }
  }

  de_sig
}


#' Filter miRNAs to cancer-dysregulated set (dbDEMC fallback)
#'
#' Intersects miRmap miRNAs with dbDEMC cancer-associated miRNAs. This is the
#' legacy approach — prefer filter_mirnas_from_geo() for SCLC-specific analysis
#' with proper differential expression statistics.
#'
#' Supports both dbDEMC 2.0 and 3.0 file formats. dbDEMC 3.0 files may include
#' additional columns (fold change, Q-value) which are used to prioritize
#' miRNAs when present, rather than treating all entries as a flat binary list.
#'
#' Integrates both up- and down-regulated miRNAs simultaneously by loading
#' both files when provided. Returns a data frame with direction information
#' (mirroring the GEO-based path's output structure).
#'
#' Limitations of this approach:
#' - dbDEMC groups all lung cancers together (NSCLC + SCLC)
#'
#' @param mirmap_path Path to miRmap CSV.
#' @param dbdemc_high_path Path to dbDEMC file for upregulated miRNAs.
#' @param dbdemc_low_path Path to dbDEMC file for downregulated miRNAs.
#'   If NULL, only upregulated miRNAs are used.
#' @param cancer_type Cancer type for dbDEMC filtering.
#' @param mirna_remove Character vector of miRNAs to exclude.
#' @param max_mirnas Maximum number of miRNAs to retain.
#' @param verbose Print diagnostic messages.
#'
#' @return Data frame with columns: mirna, direction, abs_log2fc (if available
#'   from dbDEMC 3.0, otherwise NA). When fold change data is available,
#'   miRNAs are sorted by absolute fold change (descending) so that the most
#'   dysregulated are retained first.
filter_cancer_mirnas <- function(mirmap_path,
                                 dbdemc_high_path,
                                 dbdemc_low_path = NULL,
                                 cancer_type,
                                 mirna_remove = character(0),
                                 max_mirnas = 808,
                                 verbose = TRUE) {
  mirmap_mirnas <- read.csv(mirmap_path, sep = ",")

  # Internal helper: parse a single dbDEMC file and return standardized df
  parse_dbdemc_file <- function(path, direction_label) {
    raw <- read.csv(path, sep = "\t")
    col_names <- colnames(raw)

    # Cancer type column
    cancer_col <- grep("^Cancer\\.?Type$|^cancer\\.?type$", col_names,
                        value = TRUE, ignore.case = TRUE)
    if (length(cancer_col) == 0) cancer_col <- col_names[1]
    cancer_col <- cancer_col[1]

    # Status column (up/down regulation)
    status_col <- grep("^Status$|^Expression$|^Regulation$", col_names,
                        value = TRUE, ignore.case = TRUE)
    if (length(status_col) == 0) status_col <- col_names[2]
    status_col <- status_col[1]

    # miRNA ID column
    mirna_col <- grep("miRBase|miRNA\\.?ID|miRNA\\.?Name", col_names,
                       value = TRUE, ignore.case = TRUE)
    if (length(mirna_col) == 0) mirna_col <- col_names[3]
    mirna_col <- mirna_col[1]

    # Optional fold change column (dbDEMC 3.0)
    fc_col <- grep("fold\\.?change|log\\.?fc|logFC|R\\.?fold", col_names,
                    value = TRUE, ignore.case = TRUE)
    has_fc <- length(fc_col) > 0
    if (has_fc) fc_col <- fc_col[1]

    filtered <- raw %>%
      dplyr::filter(.data[[cancer_col]] == cancer_type)

    # If file has a status column with mixed directions, filter to expected
    unique_statuses <- unique(toupper(filtered[[status_col]]))
    if (length(unique_statuses) > 1) {
      filtered <- filtered %>%
        dplyr::filter(toupper(.data[[status_col]]) == toupper(direction_label))
    }

    mirnas <- unique(filtered[[mirna_col]])
    mirnas <- mirnas[mirnas != "unknown"]

    result <- data.frame(mirna = mirnas, direction = direction_label,
                         stringsAsFactors = FALSE)

    if (has_fc && nrow(filtered) > 0) {
      fc_summary <- filtered %>%
        dplyr::filter(.data[[mirna_col]] %in% mirnas) %>%
        dplyr::group_by(.data[[mirna_col]]) %>%
        dplyr::summarise(
          abs_log2fc = mean(abs(as.numeric(.data[[fc_col]])), na.rm = TRUE),
          .groups = "drop"
        )
      result <- result %>%
        dplyr::left_join(
          fc_summary %>% dplyr::rename(mirna = 1),
          by = "mirna"
        )
    } else {
      result$abs_log2fc <- NA_real_
    }

    result
  }

  # Parse upregulated file
  up_df <- parse_dbdemc_file(dbdemc_high_path, "up")

  # Parse downregulated file if provided
  if (!is.null(dbdemc_low_path) && file.exists(dbdemc_low_path)) {
    down_df <- parse_dbdemc_file(dbdemc_low_path, "down")
    cancer_df <- dplyr::bind_rows(up_df, down_df)
  } else {
    cancer_df <- up_df
  }

  has_fc <- any(!is.na(cancer_df$abs_log2fc))

  # Intersect with miRmap miRNAs
  cancer_df <- cancer_df %>%
    dplyr::filter(mirna %in% mirmap_mirnas$mature_name) %>%
    dplyr::filter(!mirna %in% mirna_remove)

  # Deduplicate: if a miRNA appears in both up and down, keep the entry with
  # higher fold change (or first occurrence if no FC data)
  if (any(duplicated(cancer_df$mirna))) {
    if (has_fc) {
      cancer_df <- cancer_df %>%
        dplyr::arrange(dplyr::desc(abs_log2fc)) %>%
        dplyr::distinct(mirna, .keep_all = TRUE)
    } else {
      cancer_df <- cancer_df %>%
        dplyr::distinct(mirna, .keep_all = TRUE)
    }
  }

  # Sort by abs fold change if available
  if (has_fc) {
    cancer_df <- cancer_df %>%
      dplyr::arrange(dplyr::desc(abs_log2fc))
    if (verbose) message("  dbDEMC fold change data available; miRNAs ranked by dysregulation magnitude")
  }

  if (nrow(cancer_df) > max_mirnas) {
    cancer_df <- cancer_df[1:max_mirnas, ]
  }

  n_up <- sum(cancer_df$direction == "up")
  n_down <- sum(cancer_df$direction == "down")

  if (verbose) {
    message(sprintf("miRmap miRNAs: %d", nrow(mirmap_mirnas)))
    message(sprintf("dbDEMC %s miRNAs: %d up, %d down, %d total",
                    cancer_type, n_up, n_down, nrow(cancer_df)))
    if (has_fc) message("  Fold change data used for ranking")
  }

  cancer_df
}


# =============================================================================
# Database query functions
# =============================================================================

#' Load TargetScan 8.0 predictions from bulk download files
#'
#' This is the preferred method for accessing TargetScan data. It uses locally
#' downloaded prediction files from TargetScan 8.0, which include CNN-based
#' biochemical model scores (McGeary et al., 2019) that outperform the older
#' context++ scores by ~50%.
#'
#' Download files from:
#'   https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
#'
#' @param targetscan_path Path to TargetScan predictions file. This can be:
#'   - Summary_Counts.default_predictions.txt (CNN-based, recommended)
#'   - Predicted_Targets_Context_Scores.default_predictions.txt (context++)
#'   The file should be tab-delimited with columns including miRNA family/ID,
#'   gene symbol/ID, and a score column.
#' @param ts_family_path Path to miR_Family_Info.txt (maps miRNA family to
#'   individual mature miRNA IDs). If NULL, will attempt to match miRNAs by
#'   family name, which may reduce matches.
#' @param mirnas Character vector of miRNA names to filter to.
#' @param context_threshold Context++ score threshold for filtering. Interactions
#'   with cumulative weighted context++ scores <= this value are kept.
#'   Default -0.2 per Agarwal et al. (2015). Only applied when the file contains
#'   context++ scores; CNN-ranked files use their own ranking and all predictions
#'   from the file are kept.
#' @param species_id Taxonomy ID for species filtering (default 9606 for human).
#' @param verbose Print progress.
#'
#' @return Data frame with columns: mirna, gene, context_score, source.
load_targetscan_bulk <- function(targetscan_path,
                                 ts_family_path = NULL,
                                 mirnas,
                                 context_threshold = -0.2,
                                 species_id = 9606,
                                 verbose = TRUE) {
  if (verbose) message("Loading TargetScan 8.0 bulk predictions...")

  ts_data <- read.delim(targetscan_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "")

  # Identify columns flexibly across different TargetScan file formats
  col_names <- colnames(ts_data)

  # Gene column: "Gene Symbol", "Target gene", or "Gene.Symbol"
  gene_col <- grep("^Gene\\.?Symbol$|^Target\\.?gene$|^Gene\\.?ID$",
                    col_names, value = TRUE, ignore.case = TRUE)
  if (length(gene_col) == 0) {
    # Try positional fallback for Summary_Counts format
    if (ncol(ts_data) >= 2) {
      gene_col <- col_names[2]
      if (verbose) message(sprintf("  Using column '%s' as gene identifier", gene_col))
    } else {
      stop("Cannot identify gene column in TargetScan file.")
    }
  } else {
    gene_col <- gene_col[1]
  }

  # miRNA family column
  mirna_fam_col <- grep("^miR\\.?[Ff]amily$|^MiRNA\\.?family$|^miRNA$",
                         col_names, value = TRUE, ignore.case = TRUE)
  if (length(mirna_fam_col) == 0) {
    mirna_fam_col <- col_names[1]
    if (verbose) message(sprintf("  Using column '%s' as miRNA family", mirna_fam_col))
  } else {
    mirna_fam_col <- mirna_fam_col[1]
  }

  # Context++ score column (may not exist in CNN-ranked files)
  score_col <- grep("context\\+\\+|weighted\\.context|Cumulative\\.weighted",
                     col_names, value = TRUE, ignore.case = TRUE)

  # Species ID column for filtering
  species_col <- grep("^Species\\.?ID$|^Taxon", col_names, value = TRUE,
                       ignore.case = TRUE)
  if (length(species_col) > 0) {
    ts_data <- ts_data[ts_data[[species_col[1]]] == species_id, ]
  }

  # Build miRNA family -> mature miRNA name mapping
  family_to_mirna <- NULL
  if (!is.null(ts_family_path)) {
    fam_info <- read.delim(ts_family_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, comment.char = "")

    fam_mirna_col <- grep("^MiRBase\\.?ID$|^miRBase\\.?ID$|^Mature\\.?ID",
                           colnames(fam_info), value = TRUE, ignore.case = TRUE)
    fam_family_col <- grep("^miR\\.?[Ff]amily$|^MiRNA\\.?family",
                            colnames(fam_info), value = TRUE, ignore.case = TRUE)
    fam_species_col <- grep("^Species\\.?ID$|^Taxon",
                             colnames(fam_info), value = TRUE, ignore.case = TRUE)

    if (length(fam_mirna_col) > 0 && length(fam_family_col) > 0) {
      fam_subset <- fam_info
      if (length(fam_species_col) > 0) {
        fam_subset <- fam_subset[fam_subset[[fam_species_col[1]]] == species_id, ]
      }
      family_to_mirna <- data.frame(
        family = fam_subset[[fam_family_col[1]]],
        mirna = fam_subset[[fam_mirna_col[1]]],
        stringsAsFactors = FALSE
      )
      # Filter to requested miRNAs
      family_to_mirna <- family_to_mirna[family_to_mirna$mirna %in% mirnas, ]
    }
  }

  # Map families to individual miRNAs
  if (!is.null(family_to_mirna) && nrow(family_to_mirna) > 0) {
    # Join: expand each family row to its member miRNAs
    relevant_families <- unique(family_to_mirna$family)
    ts_data <- ts_data[ts_data[[mirna_fam_col]] %in% relevant_families, ]

    # Expand to individual miRNAs
    ts_expanded <- merge(
      ts_data,
      family_to_mirna,
      by.x = mirna_fam_col,
      by.y = "family",
      allow.cartesian = TRUE
    )
    mirna_col_name <- "mirna"
  } else {
    # No family mapping: try direct miRNA name matching
    ts_data <- ts_data[ts_data[[mirna_fam_col]] %in% mirnas, ]
    ts_expanded <- ts_data
    ts_expanded$mirna <- ts_expanded[[mirna_fam_col]]
    mirna_col_name <- "mirna"
  }

  if (nrow(ts_expanded) == 0) {
    if (verbose) message("  TargetScan bulk: 0 interactions matched your miRNAs")
    return(data.frame(mirna = character(), gene = character(),
                      context_score = numeric(), source = character(),
                      stringsAsFactors = FALSE))
  }

  # Build output
  result <- data.frame(
    mirna = ts_expanded[[mirna_col_name]],
    gene = ts_expanded[[gene_col]],
    stringsAsFactors = FALSE
  )

  # Add context++ score if available and apply threshold
  if (length(score_col) > 0) {
    result$context_score <- as.numeric(ts_expanded[[score_col[1]]])
    result <- result[!is.na(result$context_score) &
                       result$context_score <= context_threshold, ]
  } else {
    result$context_score <- NA_real_
  }

  result$source <- "targetscan"
  result <- unique(result)

  if (verbose) {
    message(sprintf("  TargetScan bulk: %d interactions from %d miRNAs",
                    nrow(result), length(unique(result$mirna))))
    if (length(score_col) > 0) {
      message(sprintf("  Score column used: %s (threshold <= %.2f)",
                      score_col[1], context_threshold))
    } else {
      message("  No context++ score column found; using CNN-ranked predictions as-is")
    }
  }

  result
}


#' Query TargetScan for miRNA targets with context++ score filtering (web API)
#'
#' Fallback method that queries TargetScan via the hoardeR package's web
#' scraping interface. Prefer load_targetscan_bulk() with locally downloaded
#' TargetScan 8.0 files for better performance and access to CNN-based scores
#' (McGeary et al., 2019).
#'
#' @param mirnas Character vector of miRNA names.
#' @param ts_org Species ("Human").
#' @param ts_version TargetScan release version.
#' @param max_targets Maximum targets per miRNA to retrieve.
#' @param context_threshold Context++ score threshold (interactions with
#'   scores <= this value are kept). Default -0.2 per Agarwal et al. (2015).
#' @param verbose Print progress.
#'
#' @return Data frame with columns: mirna, gene, context_score, source.
query_targetscan <- function(mirnas,
                             ts_org = "Human",
                             ts_version = "8.0",
                             max_targets = 500,
                             context_threshold = -0.2,
                             verbose = TRUE) {
  if (!requireNamespace("hoardeR", quietly = TRUE)) {
    stop("Package 'hoardeR' is required for TargetScan queries.")
  }
  if (!requireNamespace("furrr", quietly = TRUE)) {
    stop("Package 'furrr' is required for parallel TargetScan queries.")
  }

  future::plan(future::multisession, workers = min(11, future::availableCores()))

  if (verbose) message(sprintf("Querying TargetScan for %d miRNAs...",
                                length(mirnas)))

  results <- furrr::future_map(mirnas, function(m) {
    tryCatch({
      ts_result <- hoardeR::targetScan(
        mirna = m,
        species = ts_org,
        release = ts_version,
        maxOut = max_targets
      )

      if (!"Ortholog" %in% names(ts_result)) return(NULL)

      df <- ts_result %>%
        dplyr::select(dplyr::any_of(c("Ortholog", "csScore"))) %>%
        dplyr::rename(gene = Ortholog)

      if ("csScore" %in% names(df)) {
        df <- df %>%
          dplyr::mutate(
            context_score = as.numeric(gsub("[^0-9.\\-]", "", csScore))
          ) %>%
          dplyr::select(-csScore) %>%
          dplyr::filter(!is.na(context_score) & context_score <= context_threshold)
      } else {
        df$context_score <- NA_real_
      }

      df$mirna <- m
      df$source <- "targetscan"
      df
    }, error = function(e) {
      if (verbose) message(sprintf("  TargetScan error for %s: %s", m,
                                    e$message))
      NULL
    })
  })

  out <- dplyr::bind_rows(results[!sapply(results, is.null)])
  if (verbose) message(sprintf("  TargetScan: %d interactions from %d miRNAs",
                                nrow(out), length(unique(out$mirna))))
  out
}


#' Load miRDB predictions with score filtering
#'
#' @param mirdb_path Path to miRDB prediction file (tab-separated:
#'   miRNA, gene symbol, target score).
#' @param mirnas Character vector of miRNAs to filter to.
#' @param score_threshold Minimum target score (default 80, per Wong & Wang 2015).
#' @param verbose Print progress.
#'
#' @return Data frame with columns: mirna, gene, mirdb_score, source.
load_mirdb <- function(mirdb_path,
                       mirnas,
                       score_threshold = 80,
                       verbose = TRUE) {
  if (verbose) message("Loading miRDB predictions...")

  mirdb <- read.delim(mirdb_path, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE)

  # miRDB format: V1 = miRNA, V2 = target (gene symbol or RefSeq), V3 = score
  # Some versions use RefSeq IDs; adjust column names accordingly
  if (ncol(mirdb) >= 3) {
    colnames(mirdb)[1:3] <- c("mirna", "gene", "mirdb_score")
  } else {
    stop("miRDB file must have at least 3 columns (miRNA, gene/target, score).")
  }

  mirdb$mirdb_score <- as.numeric(mirdb$mirdb_score)

  out <- mirdb %>%
    dplyr::filter(mirna %in% mirnas) %>%
    dplyr::filter(mirdb_score >= score_threshold) %>%
    dplyr::mutate(source = "mirdb") %>%
    dplyr::select(mirna, gene, mirdb_score, source)

  if (verbose) message(sprintf("  miRDB: %d high-confidence interactions (score >= %d)",
                                nrow(out), score_threshold))
  out
}


#' Load miRTarBase experimentally validated interactions
#'
#' @param mirtarbase_path Path to miRTarBase CSV file.
#' @param mirnas Character vector of miRNAs to filter to.
#' @param strong_only Logical; if TRUE (default), only keep interactions with
#'   strong experimental evidence (reporter assay, western blot, qPCR).
#' @param verbose Print progress.
#'
#' @return Data frame with columns: mirna, gene, evidence_type, source.
load_mirtarbase <- function(mirtarbase_path,
                            mirnas,
                            strong_only = TRUE,
                            verbose = TRUE) {
  if (verbose) message("Loading miRTarBase validated interactions...")

  mtb <- read.csv(mirtarbase_path, stringsAsFactors = FALSE)

  # miRTarBase standard columns (may vary slightly by version)
  mirna_col <- grep("^miRNA$|^miRTarBase.ID|^miRNA\\.ID", names(mtb),
                     value = TRUE, ignore.case = TRUE)
  gene_col <- grep("^Target.Gene$|^Target\\.Gene|^Gene\\.Symbol",
                    names(mtb), value = TRUE, ignore.case = TRUE)
  evidence_col <- grep("^Experiments$|^Support\\.Type|^Evidence",
                        names(mtb), value = TRUE, ignore.case = TRUE)

  if (length(mirna_col) == 0 || length(gene_col) == 0) {
    stop(paste0("Cannot identify miRNA/gene columns in miRTarBase file. ",
                "Found columns: ", paste(names(mtb), collapse = ", ")))
  }

  mtb_clean <- data.frame(
    mirna = mtb[[mirna_col[1]]],
    gene = mtb[[gene_col[1]]],
    stringsAsFactors = FALSE
  )

  if (length(evidence_col) > 0) {
    mtb_clean$evidence_type <- mtb[[evidence_col[1]]]
  } else {
    mtb_clean$evidence_type <- "unknown"
  }

  strong_methods <- c("Luciferase reporter assay", "Reporter assay",
                       "Western blot", "qRT-PCR", "Northern blot",
                       "Proteomics", "In-situ hybridization")

  if (strong_only && "evidence_type" %in% names(mtb_clean)) {
    mtb_clean <- mtb_clean %>%
      dplyr::filter(grepl(paste(strong_methods, collapse = "|"),
                          evidence_type, ignore.case = TRUE))
  }

  out <- mtb_clean %>%
    dplyr::filter(mirna %in% mirnas) %>%
    dplyr::mutate(source = "mirtarbase") %>%
    dplyr::select(mirna, gene, evidence_type, source) %>%
    dplyr::distinct(mirna, gene, .keep_all = TRUE)

  if (verbose) {
    evidence_label <- if (strong_only) "strong-evidence" else "all"
    message(sprintf("  miRTarBase: %d %s validated interactions",
                    nrow(out), evidence_label))
  }
  out
}


# =============================================================================
# Multi-database consensus
# =============================================================================

#' Build consensus miRNA-gene interactions from multiple databases
#'
#' @param targetscan_df Data frame from query_targetscan().
#' @param mirdb_df Data frame from load_mirdb() (or NULL).
#' @param mirtarbase_df Data frame from load_mirtarbase() (or NULL).
#' @param min_databases Minimum databases supporting an interaction (default 2).
#'   Exception: miRTarBase strong-evidence interactions always pass regardless.
#' @param verbose Print summary.
#'
#' @return Data frame with columns: mirna, gene, n_databases, databases,
#'   confidence_score, has_experimental_validation.
build_consensus <- function(targetscan_df,
                            mirdb_df = NULL,
                            mirtarbase_df = NULL,
                            min_databases = 2,
                            verbose = TRUE) {
  # Standardize all to mirna + gene + source
  all_interactions <- list()

  if (!is.null(targetscan_df) && nrow(targetscan_df) > 0) {
    ts_pairs <- targetscan_df %>%
      dplyr::select(mirna, gene) %>%
      dplyr::mutate(source = "targetscan") %>%
      dplyr::distinct()
    all_interactions[["targetscan"]] <- ts_pairs
  }

  if (!is.null(mirdb_df) && nrow(mirdb_df) > 0) {
    mdb_pairs <- mirdb_df %>%
      dplyr::select(mirna, gene) %>%
      dplyr::mutate(source = "mirdb") %>%
      dplyr::distinct()
    all_interactions[["mirdb"]] <- mdb_pairs
  }

  if (!is.null(mirtarbase_df) && nrow(mirtarbase_df) > 0) {
    mtb_pairs <- mirtarbase_df %>%
      dplyr::select(mirna, gene) %>%
      dplyr::mutate(source = "mirtarbase") %>%
      dplyr::distinct()
    all_interactions[["mirtarbase"]] <- mtb_pairs
  }

  n_dbs <- length(all_interactions)
  if (n_dbs == 0) {
    stop("No interactions found in any database.")
  }
  if (n_dbs == 1 && min_databases > 1) {
    warning(paste0("Only 1 database provided but min_databases = ", min_databases,
                   ". Setting min_databases = 1."))
    min_databases <- 1
  }

  combined <- dplyr::bind_rows(all_interactions)

  # Count databases per mirna-gene pair
  consensus <- combined %>%
    dplyr::group_by(mirna, gene) %>%
    dplyr::summarise(
      n_databases = dplyr::n_distinct(source),
      databases = paste(sort(unique(source)), collapse = ","),
      .groups = "drop"
    )

  # Check for experimental validation
  experimentally_validated <- character(0)
  if (!is.null(mirtarbase_df) && nrow(mirtarbase_df) > 0) {
    experimentally_validated <- paste(mirtarbase_df$mirna,
                                      mirtarbase_df$gene, sep = ":")
  }

  consensus <- consensus %>%
    dplyr::mutate(
      has_experimental_validation = paste(mirna, gene, sep = ":") %in%
        experimentally_validated
    )

  # Apply consensus filter: >= min_databases OR experimentally validated
  consensus_filtered <- consensus %>%
    dplyr::filter(n_databases >= min_databases | has_experimental_validation)

  # Build confidence score: weighted sum of database support
  # TargetScan context++ score (normalized to 0-1 where lower raw = higher confidence)
  ts_scores <- NULL
  if (!is.null(targetscan_df) && nrow(targetscan_df) > 0 &&
      "context_score" %in% names(targetscan_df)) {
    ts_scores <- targetscan_df %>%
      dplyr::select(mirna, gene, context_score) %>%
      dplyr::group_by(mirna, gene) %>%
      dplyr::summarise(ts_confidence = mean(abs(context_score), na.rm = TRUE),
                        .groups = "drop")
  }

  # miRDB score (normalized to 0-1)
  mdb_scores <- NULL
  if (!is.null(mirdb_df) && nrow(mirdb_df) > 0) {
    mdb_scores <- mirdb_df %>%
      dplyr::select(mirna, gene, mirdb_score) %>%
      dplyr::group_by(mirna, gene) %>%
      dplyr::summarise(mdb_confidence = mean(mirdb_score, na.rm = TRUE) / 100,
                        .groups = "drop")
  }

  # Join confidence scores
  if (!is.null(ts_scores)) {
    consensus_filtered <- consensus_filtered %>%
      dplyr::left_join(ts_scores, by = c("mirna", "gene"))
  } else {
    consensus_filtered$ts_confidence <- NA_real_
  }

  if (!is.null(mdb_scores)) {
    consensus_filtered <- consensus_filtered %>%
      dplyr::left_join(mdb_scores, by = c("mirna", "gene"))
  } else {
    consensus_filtered$mdb_confidence <- NA_real_
  }

  # Composite confidence: mean of available scores + database count bonus
  consensus_filtered <- consensus_filtered %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      confidence_score = {
        scores <- c(ts_confidence, mdb_confidence)
        scores <- scores[!is.na(scores)]
        base <- if (length(scores) > 0) mean(scores) else 0.5
        db_bonus <- (n_databases - 1) * 0.1
        exp_bonus <- ifelse(has_experimental_validation, 0.2, 0)
        min(base + db_bonus + exp_bonus, 1.0)
      }
    ) %>%
    dplyr::ungroup()

  if (verbose) {
    message(sprintf("\n===== Multi-Database Consensus ====="))
    message(sprintf("Total unique interactions across all databases: %d",
                    nrow(consensus)))
    message(sprintf("Passing consensus (>= %d DBs or experimentally validated): %d",
                    min_databases, nrow(consensus_filtered)))
    message(sprintf("  Supported by 1 database: %d",
                    sum(consensus_filtered$n_databases == 1)))
    message(sprintf("  Supported by 2 databases: %d",
                    sum(consensus_filtered$n_databases == 2)))
    if (n_dbs >= 3) {
      message(sprintf("  Supported by 3 databases: %d",
                      sum(consensus_filtered$n_databases == 3)))
    }
    message(sprintf("  Experimentally validated: %d",
                    sum(consensus_filtered$has_experimental_validation)))
    message(sprintf("Unique genes in consensus: %d",
                    length(unique(consensus_filtered$gene))))
    message(sprintf("Unique miRNAs in consensus: %d",
                    length(unique(consensus_filtered$mirna))))
  }

  consensus_filtered
}


# =============================================================================
# Gene ranking
# =============================================================================

#' Compute confidence-weighted miRNA gene rankings
#'
#' Instead of counting how many miRNAs target each gene (biased toward
#' well-studied hub genes), this weights each interaction by its prediction
#' confidence and normalizes to a probability distribution.
#'
#' @param consensus_df Data frame from build_consensus().
#' @param verbose Print summary.
#'
#' @return Named numeric vector of gene scores (sums to 1), sorted descending.
compute_gene_ranking <- function(consensus_df, verbose = TRUE) {
  gene_scores <- consensus_df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      raw_score = sum(confidence_score, na.rm = TRUE),
      n_mirnas = dplyr::n_distinct(mirna),
      n_databases_mean = mean(n_databases),
      has_validation = any(has_experimental_validation),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(raw_score))

  # Normalize to sum to 1
  total <- sum(gene_scores$raw_score)
  ranking <- gene_scores$raw_score / total
  names(ranking) <- gene_scores$gene
  ranking <- sort(ranking, decreasing = TRUE)

  if (verbose) {
    message(sprintf("\n===== Gene Ranking Summary ====="))
    message(sprintf("Total genes ranked: %d", length(ranking)))
    message(sprintf("Top 10 genes:"))
    for (i in seq_len(min(10, length(ranking)))) {
      g <- names(ranking)[i]
      row <- gene_scores[gene_scores$gene == g, ]
      message(sprintf("  %2d. %s: score=%.4f, miRNAs=%d, mean_dbs=%.1f, validated=%s",
                      i, g, ranking[i], row$n_mirnas, row$n_databases_mean,
                      ifelse(row$has_validation, "yes", "no")))
    }
  }

  ranking
}


# =============================================================================
# Main function
# =============================================================================

#' Multi-database miRNA target gene calculator
#'
#' Identifies cancer-dysregulated miRNAs, queries their gene targets across
#' TargetScan, miRDB, and miRTarBase, applies literature-established confidence
#' thresholds, requires multi-database consensus, and produces confidence-
#' weighted gene rankings.
#'
#' miRNA filtering (Step 1) has two modes:
#'
#'   **GEO-based (preferred for SCLC):** Downloads SCLC miRNA expression data
#'   from GEO (default: GSE19945, 35 SCLC vs 8 normal) and runs limma
#'   differential expression. This provides SCLC-specific miRNAs with actual
#'   fold changes and FDR-adjusted p-values, and uses BOTH directions (up and
#'   down) simultaneously. Set use_geo = TRUE (default).
#'
#'   **dbDEMC fallback:** Uses pre-compiled cancer-miRNA associations from
#'   dbDEMC 3.0 (Yang et al., 2022). Integrates both up- and down-regulated
#'   miRNAs simultaneously. Groups all lung cancers together (not SCLC-specific).
#'   When fold change data is available in the dbDEMC 3.0 file, miRNAs are
#'   ranked by dysregulation magnitude. Set use_geo = FALSE.
#'
#' For TargetScan (Step 2), the preferred approach is to download the TargetScan
#' 8.0 bulk prediction files locally. These include CNN-based biochemical model
#' scores (McGeary et al., 2019, Science) that outperform the older context++
#' scores (Agarwal et al., 2015) by ~50%.
#'
#' @param use_geo Logical; if TRUE (default), use GEO-based differential
#'   expression for miRNA filtering. If FALSE, fall back to dbDEMC.
#' @param geo_accession GEO accession for miRNA expression data.
#'   Default "GSE19945" (35 SCLC vs 8 normal lung, Agilent miRNA array).
#' @param case_group Regex pattern to identify disease samples in GEO metadata.
#'   Default "SCLC".
#' @param control_group Regex pattern to identify control samples.
#'   Default "Normal|normal|Adjacent normal".
#' @param de_fdr_threshold FDR threshold for GEO differential expression.
#'   Default 0.05.
#' @param de_min_log2fc Minimum absolute log2 fold change for DE miRNAs.
#'   Default 1.0 (2-fold change).
#' @param de_results_path Path to pre-computed DE results CSV. If provided and
#'   the file exists, GEO download is skipped. Allows one-time analysis with
#'   reuse across runs.
#' @param targetscan_path Path to TargetScan 8.0 bulk predictions file, or NULL
#'   to fall back to hoardeR web queries. Download from:
#'   https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
#'   Recommended file: Summary_Counts.default_predictions.txt (CNN-based).
#' @param ts_family_path Path to TargetScan miR_Family_Info.txt for mapping
#'   miRNA families to individual mature miRNA names. Download from same page.
#' @param ts_org TargetScan species (default "Human"). Only used when falling
#'   back to hoardeR web queries.
#' @param ts_version TargetScan release (default "8.0"). Only used with hoardeR.
#' @param ts_max_targets Max targets per miRNA from hoardeR web queries.
#' @param ts_context_threshold TargetScan context++ score threshold. Interactions
#'   with scores <= this value are kept. Default -0.2 per Agarwal et al. (2015).
#' @param mirdb_path Path to miRDB predictions file, or NULL to skip miRDB.
#' @param mirdb_score_threshold miRDB minimum target score. Default 80
#'   per Wong & Wang (2015).
#' @param mirtarbase_path Path to miRTarBase CSV, or NULL to skip miRTarBase.
#' @param mirtarbase_strong_only Logical; only use strong experimental evidence.
#' @param min_databases Minimum databases for consensus (default 2).
#'   Interactions with strong miRTarBase evidence always pass regardless.
#' @param cancer_type Cancer type for dbDEMC filtering (only used when
#'   use_geo = FALSE).
#' @param mirna_remove Character vector of miRNAs to exclude.
#' @param max_mirnas Maximum miRNAs to use.
#' @param mirmap_path Path to miRmap miRNAs CSV (only used when use_geo = FALSE).
#' @param dbdemc_high_path Path to dbDEMC 3.0 upregulated miRNAs file
#'   (only used when use_geo = FALSE). Also accepts dbDEMC 2.0 format.
#' @param dbdemc_low_path Path to dbDEMC 3.0 downregulated miRNAs file
#'   (only used when use_geo = FALSE). If NULL, only upregulated are used.
#' @param save_outputs Logical; save ranking and consensus data to files.
#' @param output_dir Directory for output files.
#' @param output_prefix Prefix for output filenames.
#' @param verbose Print progress messages.
#'
#' @return Named numeric vector of gene scores (sums to 1), sorted descending.
#'   The consensus data frame is attached as attribute "consensus_details".
#'   The DE miRNA table is attached as attribute "de_mirnas" (includes
#'   direction and, when available, abs_log2fc for each miRNA).
mirna_calculator <- function(use_geo = TRUE,
                             geo_accession = "GSE19945",
                             case_group = "SCLC|Small.cell lung cancer",
                             control_group = "Normal|normal|Adjacent normal",
                             de_fdr_threshold = 0.05,
                             de_min_log2fc = 1.0,
                             de_results_path = NULL,
                             targetscan_path = "Data/mirna_data/Summary_Counts.default_predictions.txt",
                             ts_family_path = "Data/mirna_data/miR_Family_Info.txt",
                             ts_org = "Human",
                             ts_version = "8.0",
                             ts_max_targets = 500,
                             ts_context_threshold = -0.2,
                             mirdb_path = "Data/mirna_data/miRDB_v6.0_prediction_result.txt",
                             mirdb_score_threshold = 80,
                             mirtarbase_path = "Data/mirna_data/miRTarBase.csv",
                             mirtarbase_strong_only = TRUE,
                             min_databases = 2,
                             cancer_type = "lung cancer",
                             mirna_remove = "hsa-miR-129-1-3p",
                             max_mirnas = 808,
                             mirmap_path = "Data/mirna_data/mirmap_mirnas.csv",
                             dbdemc_high_path = "Data/mirna_data/dbdemc_3.0_high.txt",
                             dbdemc_low_path = "Data/mirna_data/dbdemc_3.0_low.txt",
                             save_outputs = TRUE,
                             output_dir = "Outputs/mirna",
                             output_prefix = "mirna_consensus",
                             verbose = TRUE) {
  suppressMessages({
    library(dplyr)
  })

  if (save_outputs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # =========================================================================
  # Step 1: Filter miRNAs to cancer-associated set
  # =========================================================================
  if (verbose) message("\n========== Step 1: Identify cancer-associated miRNAs ==========")

  de_mirna_df <- NULL

  if (use_geo) {
    # GEO-based: SCLC-specific limma differential expression (preferred)
    if (verbose) message("Using GEO-based differential expression (SCLC-specific)")

    de_mirna_df <- filter_mirnas_from_geo(
      geo_accession = geo_accession,
      case_group = case_group,
      control_group = control_group,
      fdr_threshold = de_fdr_threshold,
      min_log2fc = de_min_log2fc,
      mirna_remove = mirna_remove,
      max_mirnas = max_mirnas,
      de_results_path = de_results_path,
      save_de_results = save_outputs,
      output_path = file.path(output_dir, "sclc_de_mirnas.csv"),
      verbose = verbose
    )

    common_mirnas <- de_mirna_df$mirna

    if (length(common_mirnas) == 0) {
      stop("No miRNAs passed differential expression filters. ",
           "Try relaxing de_fdr_threshold or de_min_log2fc.")
    }
  } else {
    # dbDEMC 3.0 fallback (legacy approach)
    if (verbose) {
      message("Using dbDEMC 3.0 filtering (legacy approach)")
      message("  Note: dbDEMC groups all lung cancers together (not SCLC-specific).")
      message("  Consider use_geo = TRUE for SCLC-specific filtering.")
    }

    # Validate dbDEMC-specific inputs
    validate_mirna_inputs(
      mirmap_path = mirmap_path,
      dbdemc_high_path = dbdemc_high_path,
      dbdemc_low_path = dbdemc_low_path,
      mirdb_path = mirdb_path,
      mirtarbase_path = mirtarbase_path,
      targetscan_path = targetscan_path,
      ts_family_path = ts_family_path,
      cancer_type = cancer_type,
      ts_context_threshold = ts_context_threshold,
      mirdb_score_threshold = mirdb_score_threshold,
      min_databases = min_databases
    )

    de_mirna_df <- filter_cancer_mirnas(
      mirmap_path = mirmap_path,
      dbdemc_high_path = dbdemc_high_path,
      dbdemc_low_path = dbdemc_low_path,
      cancer_type = cancer_type,
      mirna_remove = mirna_remove,
      max_mirnas = max_mirnas,
      verbose = verbose
    )

    common_mirnas <- de_mirna_df$mirna

    if (length(common_mirnas) == 0) {
      stop("No miRNAs passed cancer filtering. Check cancer_type and dbDEMC files.")
    }
  }

  # =========================================================================
  # Step 2: Query target prediction databases
  # =========================================================================
  if (verbose) message("\n========== Step 2: Query target prediction databases ==========")

  # Prefer TargetScan 8.0 bulk files (CNN-based model, McGeary et al. 2019)
  # Fall back to hoardeR web queries if bulk files not available
  if (!is.null(targetscan_path)) {
    if (verbose) message("Using TargetScan 8.0 bulk predictions (CNN-based model)")
    ts_df <- load_targetscan_bulk(
      targetscan_path = targetscan_path,
      ts_family_path = ts_family_path,
      mirnas = common_mirnas,
      context_threshold = ts_context_threshold,
      verbose = verbose
    )
  } else {
    if (verbose) message("No bulk TargetScan file provided; falling back to hoardeR web queries")
    ts_df <- query_targetscan(
      mirnas = common_mirnas,
      ts_org = ts_org,
      ts_version = ts_version,
      max_targets = ts_max_targets,
      context_threshold = ts_context_threshold,
      verbose = verbose
    )
  }

  mirdb_df <- NULL
  if (!is.null(mirdb_path)) {
    mirdb_df <- load_mirdb(
      mirdb_path = mirdb_path,
      mirnas = common_mirnas,
      score_threshold = mirdb_score_threshold,
      verbose = verbose
    )
  }

  mirtarbase_df <- NULL
  if (!is.null(mirtarbase_path)) {
    mirtarbase_df <- load_mirtarbase(
      mirtarbase_path = mirtarbase_path,
      mirnas = common_mirnas,
      strong_only = mirtarbase_strong_only,
      verbose = verbose
    )
  }

  # =========================================================================
  # Step 3: Build consensus
  # =========================================================================
  if (verbose) message("\n========== Step 3: Build multi-database consensus ==========")
  consensus_df <- build_consensus(
    targetscan_df = ts_df,
    mirdb_df = mirdb_df,
    mirtarbase_df = mirtarbase_df,
    min_databases = min_databases,
    verbose = verbose
  )

  # =========================================================================
  # Step 4: Weight by miRNA dysregulation magnitude
  # =========================================================================
  if (!is.null(de_mirna_df) && "abs_log2fc" %in% names(de_mirna_df)) {
    if (verbose) message("\n========== Step 4: Weight by miRNA dysregulation ==========")

    # Join fold change information to consensus interactions
    fc_lookup <- de_mirna_df %>%
      dplyr::select(mirna, abs_log2fc, direction) %>%
      dplyr::distinct(mirna, .keep_all = TRUE)

    consensus_df <- consensus_df %>%
      dplyr::left_join(fc_lookup, by = "mirna") %>%
      dplyr::mutate(
        # Scale confidence by how strongly the miRNA is dysregulated
        # Normalize abs_log2fc to [0, 1] range relative to max observed
        mirna_dysregulation_weight = ifelse(
          !is.na(abs_log2fc),
          abs_log2fc / max(abs_log2fc, na.rm = TRUE),
          0.5  # default for miRNAs without FC data
        ),
        # Multiply confidence by dysregulation weight
        confidence_score = confidence_score * (0.5 + 0.5 * mirna_dysregulation_weight)
      )

    if (verbose) {
      message(sprintf("  Interactions with FC weighting: %d", nrow(consensus_df)))
      message(sprintf("  miRNAs with FC data: %d / %d",
                      sum(!is.na(consensus_df$abs_log2fc)),
                      nrow(consensus_df)))
    }
  }

  # =========================================================================
  # Step 5: Compute confidence-weighted gene ranking
  # =========================================================================
  step_num <- if (!is.null(de_mirna_df) && "abs_log2fc" %in% names(de_mirna_df)) 5 else 4
  if (verbose) message(sprintf("\n========== Step %d: Compute gene rankings ==========",
                                step_num))
  ranking <- compute_gene_ranking(consensus_df, verbose = verbose)

  # =========================================================================
  # Step 6: Save outputs
  # =========================================================================
  if (save_outputs) {
    ranking_path <- file.path(output_dir,
                               paste0(output_prefix, "_ranking.csv"))
    ranking_rds <- file.path(output_dir,
                              paste0(output_prefix, "_ranking.rds"))
    consensus_path <- file.path(output_dir,
                                 paste0(output_prefix, "_interactions.csv"))

    write.csv(data.frame(gene = names(ranking), score = unname(ranking)),
              file = ranking_path, row.names = FALSE)
    saveRDS(ranking, file = ranking_rds)
    write.csv(consensus_df, file = consensus_path, row.names = FALSE)

    if (verbose) {
      message(sprintf("\nOutputs saved to %s:", output_dir))
      message(sprintf("  Ranking:      %s", basename(ranking_path)))
      message(sprintf("  Ranking RDS:  %s", basename(ranking_rds)))
      message(sprintf("  Interactions: %s", basename(consensus_path)))
    }
  }

  attr(ranking, "consensus_details") <- consensus_df
  if (!is.null(de_mirna_df)) {
    attr(ranking, "de_mirnas") <- de_mirna_df
  }
  return(ranking)
}
