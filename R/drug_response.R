# ============================================================================
# drug_response.R
#
# WS3 (therapy-response arm): test whether the EMT cell-state axis predicts
# drug sensitivity / gene dependency across SCLC cell lines. This is where the
# project's thesis predicts signal -- EMT in SCLC drives CHEMORESISTANCE
# (Chen et al. 2024, PMID 38834595) -- and where the WS0 power analysis said we
# are adequately powered (|r| >= 0.35 at ~60 lines).
#
#   emt_drug_correlation()  per-feature correlation of the EMT axis with a
#                           cell-line x feature matrix (drug AUC/IC50, or CRISPR
#                           gene-effect), with BH FDR across features.
#   top_emt_associations()  convenience summary of the strongest hits.
#
# ORIENTATION (state it in your data!): for a sensitivity matrix where HIGHER =
# MORE RESISTANT (e.g. drug AUC, or log IC50), a POSITIVE correlation means
# EMT-high lines are more resistant -- the chemoresistance hypothesis. Flip the
# sign if your metric is "higher = more sensitive".
#
# Conventions match the other modules: pure functions, {cli} errors, tidy out.
# ============================================================================

#' tibble if available else data.frame
#' @keywords internal
.dr_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Coerce an EMT input to a named numeric vector (cell line -> score)
#' @keywords internal
.dr_emt_vector <- function(emt, score = "consensus") {
  if (is.numeric(emt) && !is.null(names(emt))) return(emt)
  if (is.data.frame(emt)) {
    id <- if ("sample" %in% names(emt)) "sample" else names(emt)[1]
    if (!score %in% names(emt)) {
      cli::cli_abort("EMT data.frame must contain a {.val {score}} column (or pass a named numeric vector).")
    }
    return(stats::setNames(emt[[score]], as.character(emt[[id]])))
  }
  cli::cli_abort("{.arg emt} must be a named numeric vector or a data.frame with sample + {score}.")
}

# ============================================================================
# 1. EMT vs per-feature correlation (drugs or dependencies)
# ============================================================================

#' Correlate the EMT axis with a cell-line x feature matrix
#'
#' For each feature (drug or gene), correlates the per-line EMT score with the
#' feature across the cell lines they share, then BH-adjusts the p-values across
#' all features.
#'
#' @param emt Named numeric vector (cell line -> EMT score) or a data.frame with
#'   `sample` + a score column.
#' @param features A numeric matrix/data.frame with cell lines in ROWS (rownames
#'   = cell line IDs) and features in columns (e.g. drug AUC, gene-effect).
#' @param method Correlation method: "spearman" (default) or "pearson".
#' @param min_n Minimum number of shared, non-missing cell lines required to
#'   test a feature. Default 10.
#' @param score Score column name when `emt` is a data.frame. Default "consensus".
#' @return A tibble: `feature`, `n`, `rho`, `p`, `q` (BH FDR), sorted by `q`.
#'   Positive `rho` = EMT-high lines have higher feature value (e.g. more
#'   resistant, if the matrix is oriented higher = more resistant).
#' @export
emt_drug_correlation <- function(emt, features, method = c("spearman", "pearson"),
                                 min_n = 10, score = "consensus") {
  method <- match.arg(method)
  emt_v <- .dr_emt_vector(emt, score)
  if (is.data.frame(features)) features <- as.matrix(features)
  if (!is.matrix(features) || is.null(rownames(features))) {
    cli::cli_abort("{.arg features} must be a matrix/data.frame with cell-line IDs as rownames.")
  }
  if (is.null(colnames(features))) colnames(features) <- paste0("feature_", seq_len(ncol(features)))
  storage.mode(features) <- "numeric"

  shared <- intersect(names(emt_v), rownames(features))
  if (length(shared) < min_n) {
    cli::cli_abort(c(
      "x" = "Only {length(shared)} cell lines shared between EMT scores and features (need >= {min_n}).",
      "i" = "Check that the ID schemes match (e.g. {names(emt_v)[1]} vs {rownames(features)[1]})."
    ))
  }
  ev <- emt_v[shared]
  fm <- features[shared, , drop = FALSE]

  res <- lapply(seq_len(ncol(fm)), function(j) {
    y <- fm[, j]
    ok <- is.finite(ev) & is.finite(y)
    n <- sum(ok)
    if (n < min_n) return(NULL)
    ct <- suppressWarnings(stats::cor.test(ev[ok], y[ok], method = method))
    data.frame(feature = colnames(fm)[j], n = n,
               rho = unname(ct$estimate), p = ct$p.value,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res) || nrow(res) == 0) {
    cli::cli_abort("No feature had >= {min_n} shared non-missing cell lines.")
  }
  res$q <- stats::p.adjust(res$p, method = "BH")
  res <- res[order(res$q, res$p), , drop = FALSE]
  rownames(res) <- NULL
  .dr_as_tbl(res)
}

# ============================================================================
# 2. Convenience summary
# ============================================================================

#' Summarize the strongest EMT-feature associations
#'
#' @param cor_table Output of [emt_drug_correlation()].
#' @param q_max Keep features with `q` at or below this FDR. Default 0.05.
#' @param highlight Optional character vector of feature names to always report
#'   (e.g. c("cisplatin", "etoposide")), regardless of `q_max`.
#' @return A tibble of the significant (and highlighted) associations.
#' @export
top_emt_associations <- function(cor_table, q_max = 0.05, highlight = NULL) {
  if (!all(c("feature", "rho", "q") %in% names(cor_table))) {
    cli::cli_abort("{.arg cor_table} must come from emt_drug_correlation().")
  }
  keep <- cor_table$q <= q_max
  if (!is.null(highlight)) {
    keep <- keep | tolower(cor_table$feature) %in% tolower(highlight)
  }
  out <- cor_table[keep, , drop = FALSE]
  out[order(out$q, out$p), , drop = FALSE]
}
