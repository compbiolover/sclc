#' Calculate Median Absolute Deviation (MAD) Rankings
#'
#' @param denoised_sc A numeric matrix of denoised single-cell expression data
#' @param normalize Logical indicating whether to normalize MAD scores
#' @param min_expression Minimum MAD score threshold
#' @param genes_are_rows Logical indicating if genes are in rows
#'
#' @return A named numeric vector of MAD rankings
#' @export
mad_calculator <- function(denoised_sc,
                           normalize = TRUE,
                           min_expression = 0,
                           genes_are_rows = TRUE) {
  # Input validation
  validate_mad_inputs(denoised_sc, normalize, min_expression, genes_are_rows)

  # Calculate MAD scores
  margin <- if (genes_are_rows) 1 else 2
  mads <- apply(denoised_sc, margin, mad)

  # Filter by minimum expression - changed to >= for inclusive threshold
  mads <- mads[mads >= min_expression]

  # Return early if all genes are filtered out
  if (length(mads) == 0) {
    return(numeric(0))
  }

  # Sort by decreasing MAD score
  mad_ranking <- mads[order(mads, decreasing = TRUE)]

  # Normalize if requested
  if (normalize) {
    mad_ranking <- normalize_mad_scores(mad_ranking)
  }

  mad_ranking
}

#' Validate inputs for MAD calculation
#'
#' @param denoised_sc Expression matrix to validate
#' @param normalize Normalization flag to validate
#' @param min_expression Minimum expression threshold to validate
#' @param genes_are_rows Gene orientation flag to validate
#'
#' @noRd
validate_mad_inputs <- function(denoised_sc, normalize, min_expression,
                                genes_are_rows) {
  if (!is.matrix(denoised_sc)) {
    cli::cli_abort(c(
      "x" = "denoised_sc must be a matrix",
      "i" = "You supplied a {.cls {class(denoised_sc)}}"
    ))
  }

  if (!is.numeric(denoised_sc)) {
    cli::cli_abort(c(
      "x" = "denoised_sc must contain numeric values",
      "i" = "Current matrix contains {.cls {storage.mode(denoised_sc)}}"
    ))
  }

  if (!is.logical(normalize)) {
    cli::cli_abort(c(
      "x" = "normalize must be TRUE or FALSE",
      "i" = "You supplied {.val {normalize}}"
    ))
  }

  if (!is.logical(genes_are_rows)) {
    cli::cli_abort(c(
      "x" = "genes_are_rows must be TRUE or FALSE",
      "i" = "You supplied {.val {genes_are_rows}}"
    ))
  }

  if (!is.numeric(min_expression) || length(min_expression) != 1) {
    cli::cli_abort(c(
      "x" = "min_expression must be a single numeric value",
      "i" = "You supplied {.val {min_expression}}"
    ))
  }
}

#' Normalize MAD scores
#'
#' @param mad_scores Vector of MAD scores to normalize
#'
#' @noRd
normalize_mad_scores <- function(mad_scores) {
  abs(mad_scores) / sum(abs(mad_scores))
}