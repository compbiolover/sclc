# ============================================================================
# emt_scoring.R
#
# Quantify EMT as a continuous cell-state axis from an expression matrix
# (bulk or single-cell). Several validated, complementary scores are computed
# so conclusions are not method-dependent, then combined into a single
# consensus EMT axis.
#
# COMMON ORIENTATION (important): every scorer in this file returns values
# oriented so that HIGHER = MORE MESENCHYMAL (more EMT). Methods whose native
# orientation is epithelial-high (e.g. 76GS) are negated internally so they
# agree with the rest. This makes the scores directly comparable and lets us
# rank-average them into the consensus axis.
#
# Methods implemented
#   - score_76gs()          Byers et al. 2013, Clin Cancer Res (PMID 23172883).
#                           Weighted sum; weights = within-dataset Pearson
#                           correlation of each signature gene with CDH1.
#   - score_ks()            Tan et al. 2014, EMBO Mol Med (PMID 24711451).
#                           Signed two-sample Kolmogorov-Smirnov statistic of
#                           mesenchymal vs epithelial gene CDFs, per sample.
#   - score_hallmark_emt()  GSVA/ssGSEA enrichment of MSigDB HALLMARK_EMT.
#   - score_mlr()           George et al. 2017, Cancer Res (PMID 29038175).
#                           Multinomial hybrid-E/M score in [0, 2]. EXPERIMENTAL:
#                           requires the published coefficient table.
#   - score_emt_singlecell() Sparse-robust per-cell scoring via UCell/AUCell.
#
# Orchestration
#   - load_emt_signatures()    read the vendored gene sets (Data/emt_signatures).
#   - compute_emt_scores()     run available methods -> tidy table + consensus.
#   - emt_method_concordance() Spearman agreement across methods (QC).
#
# Conventions follow R/power_analysis.R: pure functions, roxygen2 docs, {cli}
# errors, tidy tibble outputs, Bioconductor packages as soft (suggested) deps
# with an actionable error if missing.
# ============================================================================

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

#' Coerce expression input to a genes-in-rows numeric matrix
#' @keywords internal
.as_gene_matrix <- function(expr, genes_are_rows = TRUE) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) {
    cli::cli_abort(c(
      "x" = "{.arg expr} must be a numeric matrix or data.frame.",
      "i" = "You supplied a {.cls {class(expr)}}."
    ))
  }
  if (!genes_are_rows) expr <- t(expr)
  if (is.null(rownames(expr))) {
    cli::cli_abort("{.arg expr} must have gene names as rownames (genes_are_rows = {genes_are_rows}).")
  }
  if (is.null(colnames(expr))) colnames(expr) <- paste0("sample_", seq_len(ncol(expr)))
  expr
}

#' Intersect a gene set with the matrix, warning about coverage
#' @keywords internal
.match_genes <- function(genes, expr, set_name = "signature", min_genes = 3) {
  present <- intersect(unique(genes), rownames(expr))
  n_in <- length(present)
  n_set <- length(unique(genes))
  if (n_in < min_genes) {
    cli::cli_abort(c(
      "x" = "Only {n_in}/{n_set} {set_name} genes found in the expression matrix.",
      "i" = "Need at least {min_genes}. Check gene-ID type (symbols vs Ensembl) and orientation."
    ))
  }
  if (n_in < n_set) {
    cli::cli_warn("{set_name}: matched {n_in}/{n_set} genes ({round(100 * n_in / n_set)}% coverage).")
  }
  present
}

#' tibble if available, else data.frame
#' @keywords internal
.emt_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Require a soft dependency with an actionable error
#' @keywords internal
.emt_require_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    how <- if (bioc) {
      sprintf("BiocManager::install('%s')", pkg)
    } else {
      sprintf("install.packages('%s')", pkg)
    }
    cli::cli_abort(c("x" = "Package {.pkg {pkg}} is required for this function.",
                     "i" = "Install it with {.code {how}}."))
  }
}

#' Rank-normalize a numeric vector to (0, 1), preserving orientation
#' @keywords internal
.rank01 <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  r / (sum(!is.na(x)) + 1)
}

# ---------------------------------------------------------------------------
# 1. 76-gene signature (Byers 2013)
# ---------------------------------------------------------------------------

#' 76GS EMT score (Byers et al. 2013)
#'
#' Weighted sum of signature-gene expression, where each gene's weight is its
#' Pearson correlation (across samples) with CDH1 (E-cadherin). Natively the
#' score is epithelial-high; it is mean-centered and **negated** here so that
#' higher = more mesenchymal, matching the rest of this module.
#'
#' @param expr Numeric expression matrix (log-scale recommended).
#' @param genes_76gs Character vector of the 76GS genes.
#' @param cdh1_gene Symbol used for E-cadherin. Default "CDH1".
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @return Named numeric vector (one per sample); higher = more mesenchymal.
#' @export
score_76gs <- function(expr, genes_76gs, cdh1_gene = "CDH1",
                       genes_are_rows = TRUE) {
  expr <- .as_gene_matrix(expr, genes_are_rows)
  if (ncol(expr) < 3) {
    cli::cli_abort("76GS needs >= 3 samples to compute within-dataset CDH1 correlations (got {ncol(expr)}).")
  }
  if (!cdh1_gene %in% rownames(expr)) {
    cli::cli_abort(c("x" = "76GS requires {.val {cdh1_gene}} in the matrix to weight genes.",
                     "i" = "Supply a matrix containing E-cadherin or set {.arg cdh1_gene}."))
  }
  present <- .match_genes(genes_76gs, expr, "76GS")
  cdh1 <- expr[cdh1_gene, ]
  sub <- expr[present, , drop = FALSE]
  weights <- apply(sub, 1, function(g)
    suppressWarnings(stats::cor(g, cdh1, use = "pairwise.complete.obs")))
  weights[is.na(weights)] <- 0           # constant / all-NA genes contribute nothing
  # na.rm so a few missing expression values contribute 0 rather than poisoning
  # the whole sample's score.
  epithelial_score <- colSums(weights * sub, na.rm = TRUE)
  centered <- epithelial_score - mean(epithelial_score)
  stats::setNames(-centered, colnames(expr))   # negate -> higher = mesenchymal
}

# ---------------------------------------------------------------------------
# 2. KS score (Tan 2014)
# ---------------------------------------------------------------------------

#' Signed two-sample KS statistic for one sample
#' @keywords internal
.ks_signed <- function(epi_vals, mes_vals) {
  epi_vals <- epi_vals[is.finite(epi_vals)]
  mes_vals <- mes_vals[is.finite(mes_vals)]
  if (length(epi_vals) < 1 || length(mes_vals) < 1) return(NA_real_)
  grid <- sort(c(epi_vals, mes_vals))
  d <- stats::ecdf(epi_vals)(grid) - stats::ecdf(mes_vals)(grid)
  d_plus <- max(d)     # epithelial CDF above mesenchymal -> mesenchymal higher
  d_minus <- max(-d)   # mesenchymal CDF above epithelial -> epithelial higher
  if (d_plus >= d_minus) d_plus else -d_minus
}

#' KS EMT score (Tan et al. 2014)
#'
#' Per sample, compares the empirical CDFs of mesenchymal vs epithelial
#' signature-gene expression with a signed two-sample Kolmogorov-Smirnov
#' statistic. Returns values in [-1, 1] with higher = more mesenchymal.
#'
#' @param expr Numeric expression matrix.
#' @param epithelial_genes,mesenchymal_genes Character vectors.
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @return Named numeric vector in [-1, 1]; higher = more mesenchymal.
#' @export
score_ks <- function(expr, epithelial_genes, mesenchymal_genes,
                     genes_are_rows = TRUE) {
  expr <- .as_gene_matrix(expr, genes_are_rows)
  epi <- .match_genes(epithelial_genes, expr, "KS epithelial")
  mes <- .match_genes(mesenchymal_genes, expr, "KS mesenchymal")
  scores <- vapply(seq_len(ncol(expr)), function(j) {
    .ks_signed(expr[epi, j], expr[mes, j])
  }, numeric(1))
  stats::setNames(scores, colnames(expr))
}

# ---------------------------------------------------------------------------
# 3. Hallmark EMT via GSVA / ssGSEA
# ---------------------------------------------------------------------------

#' Hallmark EMT enrichment score (GSVA / ssGSEA)
#'
#' Single-sample enrichment of the MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_
#' TRANSITION gene set. The set is mesenchymal-oriented, so higher = more
#' mesenchymal natively. Supports both the modern GSVA parameter-object API
#' (>= 1.50) and the legacy `gsva(method=)` signature.
#'
#' @param expr Numeric expression matrix.
#' @param hallmark_genes Character vector of the Hallmark EMT genes.
#' @param method "gsva" or "ssgsea". Default "gsva".
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @return Named numeric vector; higher = more mesenchymal.
#' @export
score_hallmark_emt <- function(expr, hallmark_genes,
                               method = c("gsva", "ssgsea"),
                               genes_are_rows = TRUE) {
  method <- match.arg(method)
  .emt_require_pkg("GSVA", bioc = TRUE)
  expr <- .as_gene_matrix(expr, genes_are_rows)
  present <- .match_genes(hallmark_genes, expr, "Hallmark EMT")
  gsl <- list(HALLMARK_EMT = present)

  res <- if (method == "gsva" && exists("gsvaParam", where = asNamespace("GSVA"))) {
    GSVA::gsva(GSVA::gsvaParam(expr, gsl))
  } else if (method == "ssgsea" && exists("ssgseaParam", where = asNamespace("GSVA"))) {
    GSVA::gsva(GSVA::ssgseaParam(expr, gsl))
  } else {
    # legacy API
    GSVA::gsva(expr, gsl, method = method, verbose = FALSE)
  }
  stats::setNames(as.numeric(res["HALLMARK_EMT", ]), colnames(expr))
}

# ---------------------------------------------------------------------------
# 4. MLR hybrid-E/M score (George 2017) -- EXPERIMENTAL
# ---------------------------------------------------------------------------

#' MLR hybrid-E/M score (George et al. 2017) -- experimental
#'
#' Computes the multinomial logistic-regression EMT score in [0, 2]
#' (0 = epithelial, 1 = hybrid E/M, 2 = mesenchymal). This requires the
#' published model: predictor genes, normalizer genes, and class coefficients,
#' supplied via [load_emt_signatures()] (the `mlr` element). Flagged
#' experimental because the score is sensitive to the exact published
#' coefficients and to cross-platform normalization.
#'
#' @param expr Numeric expression matrix.
#' @param mlr_model A data.frame with columns: `gene`, `role`
#'   ("predictor"/"normalizer"), and one coefficient column per non-reference
#'   class (e.g. `coef_hybrid`, `coef_M`). Per-class intercepts may be supplied
#'   as an optional row with `role == "intercept"`; if that row is absent,
#'   intercepts default to 0.
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @return Named numeric vector in [0, 2]; higher = more mesenchymal.
#' @export
score_mlr <- function(expr, mlr_model, genes_are_rows = TRUE) {
  expr <- .as_gene_matrix(expr, genes_are_rows)
  if (is.null(mlr_model) || !is.data.frame(mlr_model) ||
      !all(c("gene", "role") %in% names(mlr_model))) {
    cli::cli_abort(c(
      "x" = "score_mlr() needs the published George 2017 coefficient table.",
      "i" = "Provide it via load_emt_signatures()$mlr (gene, role, coef_* columns).",
      "i" = "See Data/emt_signatures/PROVENANCE.md for the source."
    ))
  }
  coef_cols <- grep("^coef_", names(mlr_model), value = TRUE)
  if (length(coef_cols) < 1) {
    cli::cli_abort("MLR model table has no `coef_*` columns; cannot score.")
  }
  normalizers <- mlr_model$gene[mlr_model$role == "normalizer"]
  predictors  <- mlr_model$gene[mlr_model$role == "predictor"]
  norm_present <- .match_genes(normalizers, expr, "MLR normalizers", min_genes = 1)
  pred_present <- .match_genes(predictors, expr, "MLR predictors", min_genes = 1)

  baseline <- colMeans(expr[norm_present, , drop = FALSE])
  x <- sweep(expr[pred_present, , drop = FALSE], 2, baseline, "-")

  # linear predictors per non-reference class
  cls <- mlr_model[mlr_model$role == "predictor", , drop = FALSE]
  rownames(cls) <- cls$gene
  cls <- cls[pred_present, , drop = FALSE]
  intercepts <- mlr_model[mlr_model$role == "intercept", coef_cols, drop = FALSE]

  eta <- vapply(coef_cols, function(cc) {
    b0 <- if (nrow(intercepts) == 1) as.numeric(intercepts[[cc]]) else 0
    b0 + colSums(cls[[cc]] * x)
  }, numeric(ncol(expr)))
  eta <- cbind(reference = 0, eta)                 # reference class = epithelial
  eta <- eta - apply(eta, 1, max)                  # stabilize softmax (avoid overflow)
  ex <- exp(eta)
  probs <- ex / rowSums(ex)                         # softmax
  class_levels <- 0:(ncol(eta) - 1)                # 0 = E, 1 = hybrid, 2 = M, ...
  score <- as.numeric(probs %*% class_levels)
  stats::setNames(score, colnames(expr))
}

# ---------------------------------------------------------------------------
# 5. Single-cell scoring (UCell / AUCell)
# ---------------------------------------------------------------------------

#' Per-cell EMT score for single-cell data (UCell / AUCell)
#'
#' Scores each cell for epithelial and mesenchymal programs with a
#' rank-based method robust to single-cell sparsity, then returns the
#' mesenchymal-minus-epithelial contrast (higher = more mesenchymal).
#'
#' @param mat Genes-in-rows expression matrix (raw or normalized counts).
#'   A sparse `dgCMatrix` is accepted and passed through.
#' @param epithelial_genes,mesenchymal_genes Character vectors.
#' @param method "UCell" or "AUCell". Default "UCell".
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @param ncores Cores for UCell. Default 1.
#' @return Named numeric vector (one per cell); higher = more mesenchymal.
#' @export
score_emt_singlecell <- function(mat, epithelial_genes, mesenchymal_genes,
                                 method = c("UCell", "AUCell"),
                                 genes_are_rows = TRUE, ncores = 1) {
  method <- match.arg(method)
  if (!genes_are_rows) {
    # Only reach for Matrix on sparse/Matrix inputs; base t() handles dense
    # matrices without hard-requiring the Matrix package.
    mat <- if (inherits(mat, "Matrix")) Matrix::t(mat) else t(mat)
  }
  if (is.null(rownames(mat))) cli::cli_abort("{.arg mat} must have gene names as rownames.")
  cells <- colnames(mat)
  if (is.null(cells)) cells <- paste0("cell_", seq_len(ncol(mat)))
  sigs <- list(
    EMT_epithelial  = .match_genes(epithelial_genes, mat, "single-cell epithelial", min_genes = 1),
    EMT_mesenchymal = .match_genes(mesenchymal_genes, mat, "single-cell mesenchymal", min_genes = 1)
  )
  if (method == "UCell") {
    .emt_require_pkg("UCell", bioc = TRUE)
    sc <- UCell::ScoreSignatures_UCell(mat, features = sigs, ncores = ncores)
    epi <- sc[, "EMT_epithelial_UCell"]
    mes <- sc[, "EMT_mesenchymal_UCell"]
  } else {
    .emt_require_pkg("AUCell", bioc = TRUE)
    rk <- AUCell::AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE)
    au <- AUCell::AUCell_calcAUC(sigs, rk, verbose = FALSE)
    am <- AUCell::getAUC(au)
    epi <- am["EMT_epithelial", ]
    mes <- am["EMT_mesenchymal", ]
  }
  stats::setNames(as.numeric(mes - epi), cells)
}

# ---------------------------------------------------------------------------
# 6. Signature loading
# ---------------------------------------------------------------------------

#' Load the vendored EMT gene sets
#'
#' Reads the canonical, versioned gene sets from `dir`. Missing files yield
#' `NULL` for that element; the corresponding scorer then errors with an
#' actionable message rather than silently using a wrong list.
#'
#' @param dir Directory holding the vendored TSVs. Default
#'   "Data/emt_signatures". See its PROVENANCE.md for sources/citations.
#' @return A named list: `gs_76`, `ks_epithelial`, `ks_mesenchymal`,
#'   `hallmark`, `mlr` (data.frame or NULL).
#' @export
load_emt_signatures <- function(dir = "Data/emt_signatures") {
  read_genes <- function(file) {
    path <- file.path(dir, file)
    if (!file.exists(path)) return(NULL)
    df <- utils::read.delim(path, stringsAsFactors = FALSE)
    if ("gene" %in% names(df)) df$gene else df[[1]]
  }
  read_table <- function(file) {
    path <- file.path(dir, file)
    if (!file.exists(path)) return(NULL)
    utils::read.delim(path, stringsAsFactors = FALSE)
  }
  list(
    gs_76          = read_genes("byers_76gs.tsv"),
    ks_epithelial  = read_genes("tan_ks_epithelial.tsv"),
    ks_mesenchymal = read_genes("tan_ks_mesenchymal.tsv"),
    hallmark       = read_genes("hallmark_emt.tsv"),
    mlr            = read_table("george_mlr.tsv")
  )
}

# ---------------------------------------------------------------------------
# 7. Orchestrator + concordance QC
# ---------------------------------------------------------------------------

#' Compute multiple EMT scores and a consensus EMT axis
#'
#' Runs the requested methods, assembles a tidy per-sample table, rank-
#' normalizes each method to (0, 1) (orientation preserved: higher =
#' mesenchymal), averages them into a `consensus` EMT axis, and assigns
#' tertile states (E / hybrid / M).
#'
#' @param expr Numeric expression matrix.
#' @param signatures Output of [load_emt_signatures()].
#' @param methods Which scorers to run. Any of "76gs", "ks", "hallmark", "mlr".
#'   By default runs every method whose gene set is available.
#' @param genes_are_rows Whether genes are in rows. Default TRUE.
#' @param hallmark_method Passed to [score_hallmark_emt()]. Default "gsva".
#' @return A tibble: `sample`, one column per method, `consensus`, `emt_state`.
#' @export
compute_emt_scores <- function(expr, signatures = load_emt_signatures(),
                               methods = c("76gs", "ks", "hallmark", "mlr"),
                               genes_are_rows = TRUE,
                               hallmark_method = "gsva") {
  expr <- .as_gene_matrix(expr, genes_are_rows)
  samples <- colnames(expr)
  out <- list(sample = samples)

  run <- function(flag, have, fun) {
    if (flag %in% methods && isTRUE(have)) {
      res <- tryCatch(fun(), error = function(e) {
        cli::cli_warn("Skipping {.val {flag}}: {conditionMessage(e)}")
        NULL
      })
      if (!is.null(res)) out[[flag]] <<- unname(res[samples])
    }
  }

  run("76gs", !is.null(signatures$gs_76),
      function() score_76gs(expr, signatures$gs_76, genes_are_rows = TRUE))
  run("ks", !is.null(signatures$ks_epithelial) && !is.null(signatures$ks_mesenchymal),
      function() score_ks(expr, signatures$ks_epithelial,
                          signatures$ks_mesenchymal, genes_are_rows = TRUE))
  run("hallmark", !is.null(signatures$hallmark),
      function() score_hallmark_emt(expr, signatures$hallmark,
                                    method = hallmark_method, genes_are_rows = TRUE))
  run("mlr", !is.null(signatures$mlr),
      function() score_mlr(expr, signatures$mlr, genes_are_rows = TRUE))

  method_cols <- setdiff(names(out), "sample")
  if (length(method_cols) == 0) {
    cli::cli_abort(c("x" = "No EMT methods could be computed.",
                     "i" = "Check that Data/emt_signatures/ is populated (see PROVENANCE.md)."))
  }
  df <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)

  norm <- vapply(method_cols, function(m) .rank01(df[[m]]), numeric(nrow(df)))
  df$consensus <- rowMeans(norm, na.rm = TRUE)
  # Tertile states from the rank of the consensus axis. Cutting the
  # rank-normalized score on FIXED breaks (not data quantiles) keeps cut()
  # valid even when consensus is constant or has few distinct values, where
  # quantile() breaks would collide and error.
  df$emt_state <- cut(.rank01(df$consensus), breaks = c(0, 1/3, 2/3, 1),
                      labels = c("E", "hybrid", "M"), include.lowest = TRUE)
  if (length(method_cols) == 1) {
    cli::cli_warn("Only one EMT method available ({method_cols}); consensus equals it. Cross-method QC is not meaningful.")
  }
  .emt_as_tbl(df)
}

#' Cross-method concordance (Spearman) -- WS1 QC positive control
#'
#' @param scores Output of [compute_emt_scores()].
#' @return A Spearman correlation matrix among the per-method columns. Strong
#'   positive off-diagonal values indicate the methods agree on the EMT axis.
#' @export
emt_method_concordance <- function(scores) {
  method_cols <- setdiff(names(scores), c("sample", "consensus", "emt_state"))
  if (length(method_cols) < 2) {
    cli::cli_abort("Need >= 2 method columns to assess concordance; found {length(method_cols)}.")
  }
  m <- as.matrix(scores[, method_cols, drop = FALSE])
  stats::cor(m, method = "spearman", use = "pairwise.complete.obs")
}
