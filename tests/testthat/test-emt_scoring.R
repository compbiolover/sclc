# tests/testthat/test-emt_scoring.R
#
# Tests for the EMT scoring engine. Uses a synthetic expression matrix with a
# known epithelial->mesenchymal gradient and synthetic gene sets, so the tests
# do not depend on the vendored canonical lists. Every scorer is checked for
# the common orientation (higher = more mesenchymal).

if (!exists("score_76gs", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "emt_scoring.R"))
}

# ---- synthetic data: 6 samples on an E -> M gradient -----------------------
make_gradient <- function() {
  epi_genes <- c("CDH1", "EPCAM", "CLDN4", "CLDN7")
  mes_genes <- c("VIM", "ZEB1", "FN1", "CDH2")
  set.seed(42)
  m <- 1:6                                   # mesenchymal level per sample
  epi_mat <- t(sapply(epi_genes, function(g) (7 - m) + stats::rnorm(6, 0, 0.05)))
  mes_mat <- t(sapply(mes_genes, function(g) m + stats::rnorm(6, 0, 0.05)))
  expr <- rbind(epi_mat, mes_mat)
  colnames(expr) <- paste0("s", m)
  list(expr = expr, epi = epi_genes, mes = mes_genes, all = c(epi_genes, mes_genes))
}

test_that("score_76gs is oriented mesenchymal-high", {
  d <- make_gradient()
  s <- score_76gs(d$expr, genes_76gs = d$all)
  expect_length(s, 6)
  expect_equal(names(s), colnames(d$expr))
  # most mesenchymal sample (s6) scores higher than most epithelial (s1)
  expect_gt(s[["s6"]], s[["s1"]])
  # monotone increasing along the gradient
  expect_true(cor(s, 1:6, method = "spearman") > 0.9)
})

test_that("score_76gs tolerates a few missing values (no NA poisoning)", {
  d <- make_gradient()
  expr <- d$expr
  expr["VIM", "s3"] <- NA          # a few scattered NAs
  expr["CDH1", "s5"] <- NA
  s <- score_76gs(expr, genes_76gs = d$all)
  expect_true(all(is.finite(s)))   # not all-NA despite missing values
  expect_gt(s[["s6"]], s[["s1"]])  # orientation preserved
})

# Larger matrix (n=40) where CDH1 is deterministically DECOUPLED from the E/M
# gradient (alternating pattern), mimicking neuroendocrine SCLC where CDH1 does
# not track epithelial state.
make_uninformative_cdh1 <- function(n = 40) {
  epi <- c("CDH1", "EPCAM", "CLDN4", "CLDN7"); mes <- c("VIM", "ZEB1", "FN1", "CDH2")
  set.seed(13); lvl <- seq(0, 1, length.out = n)
  e <- t(sapply(epi, function(g) (1 - lvl) * 5 + stats::rnorm(n, 0, 0.05)))  # epithelial down
  m <- t(sapply(mes, function(g) lvl * 5 + stats::rnorm(n, 0, 0.05)))        # mesenchymal up
  expr <- rbind(e, m); colnames(expr) <- paste0("c", seq_len(n))
  expr["CDH1", ] <- rep_len(c(-1, 1), n)                                     # decouple CDH1
  list(expr = expr, epi = epi, mes = mes, all = c(epi, mes))
}

test_that("cdh1_informativeness detects when CDH1 tracks epithelial state", {
  d <- make_gradient()
  expect_gt(cdh1_informativeness(d$expr, d$epi, d$mes), 0.5)   # CDH1 covaries w/ epithelial
  u <- make_uninformative_cdh1()
  expect_lt(cdh1_informativeness(u$expr, u$epi, u$mes), 0.2)   # CDH1 decoupled
})

test_that("compute_emt_scores drops 76GS from consensus when CDH1 uninformative", {
  u <- make_uninformative_cdh1()
  sigs <- list(gs_76 = u$all, ks_epithelial = u$epi, ks_mesenchymal = u$mes,
               hallmark = NULL, mlr = NULL)
  res <- expect_warning(
    compute_emt_scores(u$expr, signatures = sigs, methods = c("76gs", "ks")),
    "mesenchymal markers")
  expect_true("76gs" %in% names(res))                            # column kept for transparency
  expect_false("76gs" %in% attr(res, "consensus_methods"))       # but excluded from consensus
  # informative CDH1 -> 76GS retained in the consensus
  d <- make_gradient()
  sigs2 <- list(gs_76 = d$all, ks_epithelial = d$epi, ks_mesenchymal = d$mes,
                hallmark = NULL, mlr = NULL)
  res2 <- compute_emt_scores(d$expr, signatures = sigs2, methods = c("76gs", "ks"))
  expect_true("76gs" %in% attr(res2, "consensus_methods"))
})

test_that("consensus guard warns and skips when too few mesenchymal markers present", {
  d <- make_gradient()
  sigs <- list(gs_76 = d$all, ks_epithelial = d$epi, ks_mesenchymal = d$mes,
               hallmark = NULL, mlr = NULL)
  res <- expect_warning(
    compute_emt_scores(d$expr, signatures = sigs, methods = c("76gs", "ks"),
                       mes_markers = c("NOPE1", "NOPE2")),
    "Skipping the consensus validity check")
  expect_setequal(attr(res, "consensus_methods"), c("76gs", "ks"))  # nothing dropped
})

test_that("score_ks is in [-1, 1] and oriented mesenchymal-high", {
  d <- make_gradient()
  s <- score_ks(d$expr, epithelial_genes = d$epi, mesenchymal_genes = d$mes)
  expect_true(all(s >= -1 & s <= 1))
  expect_gt(s[["s6"]], s[["s1"]])
  expect_gt(s[["s6"]], 0)   # clearly mesenchymal
  expect_lt(s[["s1"]], 0)   # clearly epithelial
})

test_that("compute_emt_scores builds a consensus axis and tertile states", {
  d <- make_gradient()
  sigs <- list(gs_76 = d$all, ks_epithelial = d$epi, ks_mesenchymal = d$mes,
               hallmark = NULL, mlr = NULL)
  res <- compute_emt_scores(d$expr, signatures = sigs)
  expect_true(all(c("sample", "76gs", "ks", "consensus", "emt_state") %in% names(res)))
  expect_equal(nrow(res), 6)
  # consensus increases along the gradient
  expect_true(cor(res$consensus, 1:6, method = "spearman") > 0.9)
  # tertile extremes are correctly called
  expect_equal(as.character(res$emt_state[res$sample == "s1"]), "E")
  expect_equal(as.character(res$emt_state[res$sample == "s6"]), "M")
})

test_that("emt_method_concordance returns positive cross-method agreement", {
  d <- make_gradient()
  sigs <- list(gs_76 = d$all, ks_epithelial = d$epi, ks_mesenchymal = d$mes)
  res <- compute_emt_scores(d$expr, signatures = sigs)
  cc <- emt_method_concordance(res)
  expect_true(is.matrix(cc))
  expect_gt(cc["76gs", "ks"], 0.5)
})

test_that("input validation and missing-gene handling are strict", {
  d <- make_gradient()
  expect_error(score_76gs("not a matrix", d$all), "numeric matrix")
  # CDH1 required for 76GS
  no_cdh1 <- d$expr[setdiff(rownames(d$expr), "CDH1"), ]
  expect_error(score_76gs(no_cdh1, d$all), "CDH1")
  # 76GS needs >= 3 samples (correlations across samples)
  expect_error(score_76gs(d$expr[, 1:2], d$all), ">= 3 samples")
  # too few matched genes
  expect_error(score_ks(d$expr, c("NOPE1", "NOPE2"), d$mes), "epithelial")
})

test_that("compute_emt_scores errors when no signatures are available", {
  d <- make_gradient()
  empty <- list(gs_76 = NULL, ks_epithelial = NULL, ks_mesenchymal = NULL,
                hallmark = NULL, mlr = NULL)
  expect_error(compute_emt_scores(d$expr, signatures = empty), "No EMT methods")
})

test_that("score_mlr refuses to run without published coefficients", {
  d <- make_gradient()
  bad <- data.frame(gene = d$all, role = "predictor")  # no coef_* columns
  expect_error(score_mlr(d$expr, bad), "coef")
})

test_that("score_hallmark_emt (GSVA) runs and is oriented mesenchymal-high", {
  d <- make_gradient()
  skip_if_not_installed("GSVA")
  # mesenchymal set as a proxy "EMT/Hallmark" set: enrichment should rise along
  # the E -> M gradient.
  s <- score_hallmark_emt(d$expr, hallmark_genes = d$mes)
  expect_length(s, 6)
  expect_gt(s[["s6"]], s[["s1"]])
})

# ---- single-cell: synthetic sparse E/M cell mixture -------------------------
make_sc_gradient <- function() {
  set.seed(7)
  epi <- c("CDH1", "EPCAM", "CLDN4", "CLDN7")
  mes <- c("VIM", "ZEB1", "FN1", "CDH2")
  bg  <- paste0("BG", seq_len(200))                  # background genes for ranking
  genes <- c(epi, mes, bg)
  n_epi <- 15L; n_mes <- 15L; ncell <- n_epi + n_mes
  m <- matrix(stats::rpois(length(genes) * ncell, 1), nrow = length(genes),
              dimnames = list(genes, paste0("c", seq_len(ncell))))
  m[epi, seq_len(n_epi)] <- m[epi, seq_len(n_epi)] +
    stats::rpois(length(epi) * n_epi, 20)
  m[mes, (n_epi + 1):ncell] <- m[mes, (n_epi + 1):ncell] +
    stats::rpois(length(mes) * n_mes, 20)
  list(mat = Matrix::Matrix(m, sparse = TRUE), epi = epi, mes = mes,
       epi_cells = seq_len(n_epi), mes_cells = (n_epi + 1):ncell)
}

test_that("score_emt_singlecell (UCell) scores mesenchymal cells higher", {
  skip_if_not_installed("UCell")
  d <- make_sc_gradient()
  s <- score_emt_singlecell(d$mat, d$epi, d$mes, method = "UCell")
  expect_length(s, ncol(d$mat))
  expect_gt(mean(s[d$mes_cells]), mean(s[d$epi_cells]))
})

test_that("score_emt_singlecell (AUCell) scores mesenchymal cells higher", {
  skip_if_not_installed("AUCell")
  d <- make_sc_gradient()
  s <- score_emt_singlecell(d$mat, d$epi, d$mes, method = "AUCell")
  expect_length(s, ncol(d$mat))
  expect_gt(mean(s[d$mes_cells]), mean(s[d$epi_cells]))
})

test_that("compute_emt_scores integrates 3 methods (76GS + KS + Hallmark)", {
  skip_if_not_installed("GSVA")
  d <- make_gradient()
  sigs <- list(gs_76 = d$all, ks_epithelial = d$epi, ks_mesenchymal = d$mes,
               hallmark = d$mes, mlr = NULL)
  res <- compute_emt_scores(d$expr, signatures = sigs)
  expect_true(all(c("76gs", "ks", "hallmark", "consensus") %in% names(res)))
  expect_true(cor(res$consensus, 1:6, method = "spearman") > 0.9)
  cc <- emt_method_concordance(res)
  expect_equal(dim(cc), c(3, 3))
})
