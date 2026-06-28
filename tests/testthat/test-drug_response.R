# tests/testthat/test-drug_response.R
#
# Tests for the WS3 drug-response module on synthetic cell-line data where one
# drug's resistance rises with EMT, one falls, and others are noise.

if (!exists("emt_drug_correlation", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "drug_response.R"))
}

make_drug_data <- function(n = 80) {
  set.seed(11)
  lines <- paste0("CL", seq_len(n))
  emt <- stats::setNames(stats::rnorm(n), lines)
  feats <- cbind(
    resistant = 0.6 * emt + stats::rnorm(n, 0, 0.8),   # resistance up with EMT (+)
    sensitive = -0.6 * emt + stats::rnorm(n, 0, 0.8),  # resistance down with EMT (-)
    noise     = stats::rnorm(n),                       # no association
    sparse    = NA_real_                               # too few values -> excluded
  )
  feats[1:5, "sparse"] <- stats::rnorm(5)              # only 5 non-NA
  rownames(feats) <- lines
  list(emt = emt, feats = feats, lines = lines)
}

test_that("emt_drug_correlation returns a tidy FDR-adjusted table", {
  d <- make_drug_data()
  res <- emt_drug_correlation(d$emt, d$feats)
  expect_true(all(c("feature", "n", "rho", "p", "q") %in% names(res)))
  expect_true(all(res$q >= res$p - 1e-9))                 # BH q >= p
  expect_false("sparse" %in% res$feature)                 # < min_n excluded
})

test_that("it recovers the EMT->resistance signal with correct sign", {
  d <- make_drug_data()
  res <- emt_drug_correlation(d$emt, d$feats)
  r <- res[res$feature == "resistant", ]
  s <- res[res$feature == "sensitive", ]
  expect_gt(r$rho, 0); expect_lt(r$q, 0.05)               # EMT-high = more resistant
  expect_lt(s$rho, 0); expect_lt(s$q, 0.05)               # EMT-high = more sensitive
  expect_gt(res$q[res$feature == "noise"], 0.05)          # noise not significant
})

test_that("emt_drug_correlation accepts a data.frame EMT input", {
  d <- make_drug_data()
  emt_df <- data.frame(sample = names(d$emt), consensus = unname(d$emt),
                       stringsAsFactors = FALSE)
  res <- emt_drug_correlation(emt_df, d$feats)
  expect_equal(nrow(res), 3)                              # resistant, sensitive, noise
})

test_that("emt_drug_correlation errors on insufficient overlap", {
  d <- make_drug_data()
  bad <- d$feats; rownames(bad) <- paste0("X", rownames(bad))
  expect_error(emt_drug_correlation(d$emt, bad), "shared")
})

test_that("emt_drug_correlation rejects non-numeric features and EMT df without sample", {
  d <- make_drug_data()
  ch <- d$feats; ch[] <- as.character(ch)                 # character matrix, rownames kept
  expect_error(emt_drug_correlation(d$emt, ch), "numeric")
  df_nosample <- data.frame(id = names(d$emt), consensus = unname(d$emt),
                            stringsAsFactors = FALSE)
  expect_error(emt_drug_correlation(df_nosample, d$feats), "sample")
})

test_that("top_emt_associations validates required columns", {
  d <- make_drug_data()
  res <- emt_drug_correlation(d$emt, d$feats)
  expect_error(top_emt_associations(res[, c("feature", "rho", "q")]),
               "must come from emt_drug_correlation")
})

test_that("top_emt_associations filters by FDR and honors highlight", {
  d <- make_drug_data()
  res <- emt_drug_correlation(d$emt, d$feats)
  top <- top_emt_associations(res, q_max = 0.05)
  expect_true(all(top$q <= 0.05))
  expect_true("noise" %in% top_emt_associations(res, q_max = 0.05,
                                                highlight = "noise")$feature)
})
