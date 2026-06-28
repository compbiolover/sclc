# ============================================================================
# emt_subtype_map.R
#
# Map the EMT cell-state axis (from R/emt_scoring.R) onto SCLC biology:
#   - call_sclc_subtype()   robust A/N/P transcription-factor subtype calls
#   - ne_score()            neuroendocrine (NE) vs non-NE score (Zhang 2018)
#   - map_emt_to_subtype()  cross-tabulate the EMT axis vs subtype / NE
#
# Why this exists
#   The legacy subtype caller (R/main.Rmd) used raw-expression argmax over
#   ASCL1/NEUROD1/POU2F3/YAP1. That is biased (raw magnitudes differ per gene,
#   so the most highly-expressed marker tends to win) and YAP1 is no longer a
#   lineage subtype -- YAP1-expressing SCLC lines were reclassified as
#   SMARCA4-deficient undifferentiated tumors (Meder et al. 2024, Clin Cancer
#   Res, PMID 38328215). Here subtype is called on PER-GENE Z-SCORES across
#   samples among A/N/P only, with a separate SMARCA4-UT QC flag for high-YAP1 /
#   low-NE samples.
#
# Biological expectation (positive control; Groves et al. 2023, PMID 36900269):
#   SCLC-A leans epithelial (low EMT), non-NE samples lean mesenchymal (high
#   EMT); EMT axis should correlate negatively with the NE score.
#
# Conventions match R/emt_scoring.R: pure functions, {cli} errors, tidy output.
# Self-contained (own small internal helpers) so it can be sourced alone.
# ============================================================================

# ---- internal helpers ------------------------------------------------------

#' Coerce expression input to a genes-in-rows numeric matrix
#' @keywords internal
.sm_gene_matrix <- function(expr, genes_are_rows = TRUE) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) {
    cli::cli_abort(c("x" = "{.arg expr} must be a numeric matrix or data.frame.",
                     "i" = "You supplied a {.cls {class(expr)}}."))
  }
  if (!genes_are_rows) expr <- t(expr)
  if (is.null(rownames(expr))) {
    cli::cli_abort("{.arg expr} must have gene symbols as rownames.")
  }
  if (is.null(colnames(expr))) colnames(expr) <- paste0("sample_", seq_len(ncol(expr)))
  expr
}

#' Row z-score (across samples) for one gene, matched case-insensitively
#' @keywords internal
.sm_zrow <- function(expr, gene) {
  hit <- which(toupper(rownames(expr)) == toupper(gene))
  if (length(hit) == 0) return(NULL)
  x <- as.numeric(expr[hit[1], ])
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(stats::setNames(rep(0, length(x)), colnames(expr)))
  stats::setNames((x - mean(x, na.rm = TRUE)) / s, colnames(expr))
}

#' tibble if available, else data.frame
#' @keywords internal
.sm_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

# ============================================================================
# 1. SCLC subtype caller (A / N / P) + SMARCA4-UT flag
# ============================================================================

#' Call SCLC transcription-factor subtypes (A / N / P)
#'
#' Assigns each sample to SCLC-A (ASCL1), SCLC-N (NEUROD1) or SCLC-P (POU2F3)
#' by the largest PER-GENE z-score (computed across samples) among the three
#' lineage transcription factors. YAP1 is deliberately NOT a subtype: instead a
#' `smarca4_ut_flag` marks samples with high YAP1 and uniformly low A/N/P, which
#' are candidate SMARCA4-deficient undifferentiated tumors rather than SCLC.
#'
#' @param expr Numeric expression matrix (log-scale recommended). Gene matching
#'   is case-insensitive, so lowercase symbols (as in `scl_common`) are fine.
#' @param markers Named list mapping subtype letters to marker symbols.
#'   Default `list(A = "ASCL1", N = "NEUROD1", P = "POU2F3")`.
#' @param yap1_gene YAP1 symbol used only for the SMARCA4-UT flag. Default "YAP1".
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @param confidence_min Minimum z-score gap between the top and second marker
#'   for a confident call; below this the call is kept but `low_confidence` is
#'   TRUE. Default 0.25.
#' @param yap1_flag_z YAP1 z-score above which (with all A/N/P z below
#'   `anp_flag_z`) a sample is flagged SMARCA4-UT. Default 1.
#' @param anp_flag_z Upper bound on the max A/N/P z-score for the SMARCA4-UT
#'   flag. Default 0.
#' @return A tibble: `sample`, `subtype` (A/N/P), `z_A`,`z_N`,`z_P`, `margin`,
#'   `low_confidence`, `yap1_z`, `smarca4_ut_flag`.
#' @export
call_sclc_subtype <- function(expr,
                              markers = list(A = "ASCL1", N = "NEUROD1", P = "POU2F3"),
                              yap1_gene = "YAP1", genes_are_rows = TRUE,
                              confidence_min = 0.25, yap1_flag_z = 1,
                              anp_flag_z = 0) {
  expr <- .sm_gene_matrix(expr, genes_are_rows)
  if (ncol(expr) < 3) {
    cli::cli_abort("Subtype calling needs >= 3 samples for across-sample z-scores (got {ncol(expr)}).")
  }
  z <- lapply(markers, function(g) .sm_zrow(expr, g))
  missing <- names(z)[vapply(z, is.null, logical(1))]
  if (length(missing) > 0) {
    cli::cli_abort(c("x" = "Subtype marker(s) not found in the matrix: {.val {unlist(markers[missing])}}.",
                     "i" = "Check gene symbols / ID type."))
  }
  zmat <- do.call(cbind, z)                       # samples x {A,N,P}
  colnames(zmat) <- names(markers)
  best <- max.col(zmat, ties.method = "first")
  subtype <- names(markers)[best]
  sorted <- t(apply(zmat, 1, sort, decreasing = TRUE))
  margin <- sorted[, 1] - sorted[, 2]
  max_anp <- sorted[, 1]

  yz <- .sm_zrow(expr, yap1_gene)
  yap1_z <- if (is.null(yz)) rep(NA_real_, ncol(expr)) else as.numeric(yz)
  smarca4_ut <- !is.na(yap1_z) & yap1_z >= yap1_flag_z & max_anp <= anp_flag_z

  out <- data.frame(
    sample = colnames(expr), subtype = subtype,
    z_A = zmat[, "A"], z_N = zmat[, "N"], z_P = zmat[, "P"],
    margin = margin, low_confidence = margin < confidence_min,
    yap1_z = yap1_z, smarca4_ut_flag = smarca4_ut,
    row.names = NULL, stringsAsFactors = FALSE
  )
  if (any(smarca4_ut)) {
    cli::cli_warn("{sum(smarca4_ut)} sample(s) flagged as candidate SMARCA4-UT (high YAP1, low A/N/P) -- review before treating as SCLC.")
  }
  .sm_as_tbl(out)
}

# ============================================================================
# 2. Neuroendocrine (NE) score -- Zhang et al. 2018
# ============================================================================

#' Load the vendored Zhang 2018 NE template
#'
#' @param path TSV with columns `gene`, `ne_ref`, `nonne_ref` (reference mean
#'   expression of each of the 50 genes in NE and non-NE samples). See
#'   Data/sclc_signatures/PROVENANCE.md for the source (Zhang et al. 2018,
#'   Transl Lung Cancer Res, PMID 29535911).
#' @return A data.frame, or NULL if the file is absent.
#' @export
load_ne_template <- function(path = "Data/sclc_signatures/zhang_ne_50.tsv") {
  if (!file.exists(path)) return(NULL)
  utils::read.delim(path, stringsAsFactors = FALSE)
}

#' Neuroendocrine (NE) score (Zhang et al. 2018)
#'
#' Per sample, `NE = (cor_NE - cor_nonNE) / 2`, where `cor_NE`/`cor_nonNE` are
#' the Pearson correlations between the sample's expression of the 50 template
#' genes and the NE / non-NE reference vectors. Ranges roughly [-1, 1]; positive
#' = neuroendocrine, negative = non-NE.
#'
#' @param expr Numeric expression matrix.
#' @param ne_template Output of [load_ne_template()] (`gene`, `ne_ref`,
#'   `nonne_ref`). Required -- the 50 genes + references are not hand-typed.
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @return Named numeric vector; higher = more neuroendocrine.
#' @export
ne_score <- function(expr, ne_template = load_ne_template(), genes_are_rows = TRUE) {
  if (is.null(ne_template) || !is.data.frame(ne_template) ||
      !all(c("gene", "ne_ref", "nonne_ref") %in% names(ne_template))) {
    cli::cli_abort(c(
      "x" = "ne_score() needs the Zhang 2018 NE template (gene, ne_ref, nonne_ref).",
      "i" = "Vendor it at Data/sclc_signatures/zhang_ne_50.tsv -- see PROVENANCE.md."
    ))
  }
  expr <- .sm_gene_matrix(expr, genes_are_rows)
  ridx <- match(toupper(ne_template$gene), toupper(rownames(expr)))
  keep <- !is.na(ridx)
  if (sum(keep) < 10) {
    cli::cli_abort("Only {sum(keep)}/{nrow(ne_template)} NE template genes found; too few to score.")
  }
  if (sum(keep) < nrow(ne_template)) {
    cli::cli_warn("NE score: matched {sum(keep)}/{nrow(ne_template)} template genes.")
  }
  sub <- expr[ridx[keep], , drop = FALSE]
  ne_ref <- ne_template$ne_ref[keep]
  nonne_ref <- ne_template$nonne_ref[keep]
  scores <- vapply(seq_len(ncol(sub)), function(j) {
    v <- sub[, j]
    (stats::cor(v, ne_ref, use = "pairwise.complete.obs") -
     stats::cor(v, nonne_ref, use = "pairwise.complete.obs")) / 2
  }, numeric(1))
  stats::setNames(scores, colnames(expr))
}

# ============================================================================
# 3. Map the EMT axis onto subtype / NE
# ============================================================================

#' Cross-tabulate the EMT axis against SCLC subtype and NE score
#'
#' Joins per-sample EMT consensus scores with subtype calls (and optionally NE
#' scores) and summarizes how the EMT axis distributes across SCLC biology --
#' the WS2 deliverable / positive control.
#'
#' @param emt_scores Output of `emt_scoring::compute_emt_scores()` (needs
#'   `sample` and `consensus`; `emt_state` used if present).
#' @param subtypes Output of [call_sclc_subtype()] (needs `sample`, `subtype`).
#' @param ne Optional named numeric vector or data.frame (`sample`, `ne`) of NE
#'   scores from [ne_score()].
#' @return A list with `per_sample` (merged tidy table) and `by_subtype`
#'   (summary: n, mean/median EMT consensus, and EMT~NE correlation if NE given).
#' @export
map_emt_to_subtype <- function(emt_scores, subtypes, ne = NULL) {
  if (!all(c("sample", "consensus") %in% names(emt_scores))) {
    cli::cli_abort("{.arg emt_scores} must contain `sample` and `consensus` columns.")
  }
  if (!all(c("sample", "subtype") %in% names(subtypes))) {
    cli::cli_abort("{.arg subtypes} must contain `sample` and `subtype` columns.")
  }
  per <- merge(
    data.frame(sample = emt_scores$sample, emt_consensus = emt_scores$consensus,
               emt_state = if ("emt_state" %in% names(emt_scores)) as.character(emt_scores$emt_state) else NA,
               stringsAsFactors = FALSE),
    data.frame(sample = subtypes$sample, subtype = subtypes$subtype,
               smarca4_ut_flag = if ("smarca4_ut_flag" %in% names(subtypes)) subtypes$smarca4_ut_flag else NA,
               stringsAsFactors = FALSE),
    by = "sample"
  )
  if (!is.null(ne)) {
    if (is.data.frame(ne)) {
      if (!all(c("sample", "ne") %in% names(ne))) {
        cli::cli_abort("{.arg ne} data.frame must contain `sample` and `ne` columns.")
      }
      ne_df <- ne[, c("sample", "ne")]
    } else {
      if (is.null(names(ne))) {
        cli::cli_abort("{.arg ne} vector must be named by sample.")
      }
      ne_df <- data.frame(sample = names(ne), ne = as.numeric(ne), stringsAsFactors = FALSE)
    }
    per <- merge(per, ne_df, by = "sample", all.x = TRUE)
  }
  if (nrow(per) == 0) {
    cli::cli_abort("No samples in common between `emt_scores` and `subtypes`.")
  }

  subs <- split(per, per$subtype)
  by_sub <- do.call(rbind, lapply(names(subs), function(s) {
    d <- subs[[s]]
    row <- data.frame(
      subtype = s, n = nrow(d),
      mean_emt = mean(d$emt_consensus, na.rm = TRUE),
      median_emt = stats::median(d$emt_consensus, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (!is.null(ne) && sum(!is.na(d$ne)) > 2) {
      row$cor_emt_ne <- stats::cor(d$emt_consensus, d$ne,
                                   method = "spearman", use = "pairwise.complete.obs")
    }
    row
  }))
  list(per_sample = .sm_as_tbl(per), by_subtype = .sm_as_tbl(by_sub))
}
