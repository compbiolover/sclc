# Name: mirna_calculator.R
# Author: Andrew Willems <awillems@vols.utk.edu>
# Purpose: Multi-database miRNA target gene ranking with confidence weighting.
#
# Approach:
#   1. Filter miRNAs to those dysregulated in the cancer of interest (dbDEMC)
#      AND confirmed in miRmap
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
#   - Data/mirna_data/dbdemc_2.0_high.txt       (dbDEMC upregulated miRNAs)
#   - Data/mirna_data/dbdemc_2.0_low.txt        (dbDEMC downregulated miRNAs)
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
#' @param dbdemc_path Path to dbDEMC data file.
#' @param mirdb_path Path to miRDB predictions file (or NULL to skip).
#' @param mirtarbase_path Path to miRTarBase CSV file (or NULL to skip).
#' @param targetscan_path Path to TargetScan 8.0 bulk predictions file (or NULL
#'   to fall back to hoardeR web queries).
#' @param ts_family_path Path to TargetScan miR_Family_Info.txt file (or NULL).
#' @param cancer_type Cancer type string for dbDEMC filtering.
#' @param status Direction of dysregulation ("up" or "down").
#' @param ts_context_threshold TargetScan context++ score threshold.
#' @param mirdb_score_threshold miRDB target score threshold.
#' @param min_databases Minimum number of databases for consensus.
#'
#' @return Invisible TRUE if all checks pass; stops with informative error otherwise.
validate_mirna_inputs <- function(mirmap_path,
                                  dbdemc_path,
                                  mirdb_path,
                                  mirtarbase_path,
                                  targetscan_path = NULL,
                                  ts_family_path = NULL,
                                  cancer_type,
                                  status,
                                  ts_context_threshold,
                                  mirdb_score_threshold,
                                  min_databases) {
  if (!is.character(mirmap_path) || length(mirmap_path) != 1) {
    stop("mirmap_path must be a single character string.")
  }
  if (!file.exists(mirmap_path)) {
    stop(paste0("miRmap file not found: ", mirmap_path))
  }

  if (!is.character(dbdemc_path) || length(dbdemc_path) != 1) {
    stop("dbdemc_path must be a single character string.")
  }
  if (!file.exists(dbdemc_path)) {
    stop(paste0("dbDEMC file not found: ", dbdemc_path))
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

  if (!status %in% c("up", "down")) {
    stop("status must be 'up' or 'down'.")
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

#' Filter miRNAs to cancer-dysregulated set
#'
#' Intersects miRmap miRNAs with dbDEMC cancer-associated miRNAs.
#'
#' @param mirmap_path Path to miRmap CSV.
#' @param dbdemc_path Path to dbDEMC file.
#' @param cancer_type Cancer type for dbDEMC filtering.
#' @param status "up" or "down" regulation.
#' @param mirna_remove Character vector of miRNAs to exclude.
#' @param max_mirnas Maximum number of miRNAs to retain.
#' @param verbose Print diagnostic messages.
#'
#' @return Character vector of filtered miRNA names.
filter_cancer_mirnas <- function(mirmap_path,
                                 dbdemc_path,
                                 cancer_type,
                                 status,
                                 mirna_remove = character(0),
                                 max_mirnas = 808,
                                 verbose = TRUE) {
  mirmap_mirnas <- read.csv(mirmap_path, sep = ",")
  cancer_data <- read.csv(dbdemc_path, sep = "\t")

  cancer_data <- cancer_data %>%
    dplyr::filter(Cancer.Type == cancer_type) %>%
    dplyr::filter(Status == toupper(status))

  cancer_mirnas <- unique(cancer_data$miRBase.Update.ID)
  cancer_mirnas <- cancer_mirnas[cancer_mirnas != "unknown"]

  common_mirnas <- intersect(mirmap_mirnas$mature_name, cancer_mirnas)
  common_mirnas <- common_mirnas[!common_mirnas %in% mirna_remove]

  if (length(common_mirnas) > max_mirnas) {
    common_mirnas <- common_mirnas[1:max_mirnas]
  }

  if (verbose) {
    message(sprintf("miRmap miRNAs: %d", nrow(mirmap_mirnas)))
    message(sprintf("dbDEMC %s-%s miRNAs: %d", cancer_type, status,
                    length(cancer_mirnas)))
    message(sprintf("Common miRNAs (after filtering): %d",
                    length(common_mirnas)))
  }

  common_mirnas
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
#' Queries TargetScan, miRDB, and miRTarBase for miRNA-gene interactions,
#' applies literature-established confidence thresholds, requires multi-database
#' consensus, and produces confidence-weighted gene rankings.
#'
#' For TargetScan, the preferred approach is to download the TargetScan 8.0
#' bulk prediction files locally. These include CNN-based biochemical model
#' scores (McGeary et al., 2019, Science) that outperform the older context++
#' scores (Agarwal et al., 2015) by ~50%. If bulk files are provided via
#' targetscan_path, they are used instead of web scraping via hoardeR.
#'
#' @param targetscan_path Path to TargetScan 8.0 bulk predictions file, or NULL
#'   to fall back to hoardeR web queries. Download from:
#'   https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi
#'   Recommended file: Summary_Counts.default_predictions.txt (CNN-based).
#' @param ts_family_path Path to TargetScan miR_Family_Info.txt for mapping
#'   miRNA families to individual mature miRNA names. Only needed when using
#'   targetscan_path. Download from the same page.
#' @param ts_org TargetScan species (default "Human"). Only used when falling
#'   back to hoardeR web queries.
#' @param ts_version TargetScan release (default "8.0"). Only used when falling
#'   back to hoardeR web queries.
#' @param ts_max_targets Max targets to retrieve per miRNA from TargetScan.
#'   Only used when falling back to hoardeR web queries.
#' @param ts_context_threshold TargetScan context++ score threshold. Interactions
#'   with scores <= this value are kept. Default -0.2 per Agarwal et al. (2015).
#' @param mirdb_path Path to miRDB predictions file, or NULL to skip miRDB.
#' @param mirdb_score_threshold miRDB minimum target score. Default 80
#'   per Wong & Wang (2015).
#' @param mirtarbase_path Path to miRTarBase CSV, or NULL to skip miRTarBase.
#' @param mirtarbase_strong_only Logical; only use strong experimental evidence.
#' @param min_databases Minimum databases for consensus (default 2).
#'   Interactions with strong miRTarBase evidence always pass regardless.
#' @param cancer_type Cancer type for dbDEMC filtering.
#' @param status Dysregulation direction ("up" or "down").
#' @param mirna_remove Character vector of miRNAs to exclude.
#' @param max_mirnas Maximum miRNAs to use.
#' @param mirmap_path Path to miRmap miRNAs CSV.
#' @param dbdemc_path Path to dbDEMC file. If NULL, auto-selects based on
#'   cancer_up parameter.
#' @param cancer_up Logical; TRUE for upregulated, FALSE for downregulated.
#'   Used to select dbDEMC file when dbdemc_path is NULL.
#' @param save_outputs Logical; save ranking and consensus data to files.
#' @param output_dir Directory for output files.
#' @param output_prefix Prefix for output filenames.
#' @param verbose Print progress messages.
#'
#' @return Named numeric vector of gene scores (sums to 1), sorted descending.
#'   Also invisibly returns a list with the full consensus data frame
#'   as an attribute "consensus_details".
mirna_calculator <- function(targetscan_path = "Data/mirna_data/Summary_Counts.default_predictions.txt",
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
                             status = "up",
                             mirna_remove = "hsa-miR-129-1-3p",
                             max_mirnas = 808,
                             mirmap_path = "Data/mirna_data/mirmap_mirnas.csv",
                             dbdemc_path = NULL,
                             cancer_up = TRUE,
                             save_outputs = TRUE,
                             output_dir = "Outputs/mirna",
                             output_prefix = "mirna_consensus",
                             verbose = TRUE) {
  suppressMessages({
    library(dplyr)
  })

  # Auto-select dbDEMC file if not provided
  if (is.null(dbdemc_path)) {
    dbdemc_path <- if (cancer_up) {
      "Data/mirna_data/dbdemc_2.0_high.txt"
    } else {
      "Data/mirna_data/dbdemc_2.0_low.txt"
    }
  }

  # Validate inputs
  validate_mirna_inputs(
    mirmap_path = mirmap_path,
    dbdemc_path = dbdemc_path,
    mirdb_path = mirdb_path,
    mirtarbase_path = mirtarbase_path,
    targetscan_path = targetscan_path,
    ts_family_path = ts_family_path,
    cancer_type = cancer_type,
    status = status,
    ts_context_threshold = ts_context_threshold,
    mirdb_score_threshold = mirdb_score_threshold,
    min_databases = min_databases
  )

  if (save_outputs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Step 1: Filter miRNAs to cancer-associated set
  if (verbose) message("\n========== Step 1: Filter cancer-associated miRNAs ==========")
  common_mirnas <- filter_cancer_mirnas(
    mirmap_path = mirmap_path,
    dbdemc_path = dbdemc_path,
    cancer_type = cancer_type,
    status = status,
    mirna_remove = mirna_remove,
    max_mirnas = max_mirnas,
    verbose = verbose
  )

  if (length(common_mirnas) == 0) {
    stop("No miRNAs passed cancer filtering. Check cancer_type and status.")
  }

  # Step 2: Query databases
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

  # Step 3: Build consensus
  if (verbose) message("\n========== Step 3: Build multi-database consensus ==========")
  consensus_df <- build_consensus(
    targetscan_df = ts_df,
    mirdb_df = mirdb_df,
    mirtarbase_df = mirtarbase_df,
    min_databases = min_databases,
    verbose = verbose
  )

  # Step 4: Compute confidence-weighted gene ranking
  if (verbose) message("\n========== Step 4: Compute gene rankings ==========")
  ranking <- compute_gene_ranking(consensus_df, verbose = verbose)

  # Step 5: Save outputs
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
  return(ranking)
}
