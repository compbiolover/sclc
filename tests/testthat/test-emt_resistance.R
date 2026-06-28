# tests/testthat/test-emt_resistance.R
#
# Tests for the WS4 dynamic-resistance module on synthetic paired-CDX data where
# the resistant state is, by construction, more mesenchymal (higher EMT
# location), more heterogeneous (higher dispersion), and richer in mesenchymal
# cells. No external data required.

if (!exists("emt_resistance_shift", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "emt_resistance.R"))
}

# n_models paired CDX models; each has `n_per` sensitive and `n_per` resistant
# cells. Resistant cells get a positive EMT offset (location), wider spread
# (dispersion), and per-model baselines so pairing is meaningful.
make_cdx_data <- function(n_models = 6, n_per = 60, seed = 7,
                          loc = 0.9, sd_sens = 1.0, sd_res = 1.7) {
  set.seed(seed)
  models <- paste0("SC", seq_len(n_models))
  rows <- lapply(seq_along(models), function(i) {
    base <- stats::rnorm(1, 0, 0.5)                 # per-model baseline EMT level
    sens <- stats::rnorm(n_per, base, sd_sens)
    res  <- stats::rnorm(n_per, base + loc, sd_res) # up in mean AND spread
    data.frame(
      cell = c(paste0(models[i], "_S", seq_len(n_per)),
               paste0(models[i], "_R", seq_len(n_per))),
      model = models[i],
      condition = rep(c("sensitive", "resistant"), each = n_per),
      emt = c(sens, res),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Build the (named EMT vector, cell_meta) pair prepare_resistance_emt() consumes.
as_inputs <- function(d, condition = d$condition) {
  emt <- stats::setNames(d$emt, d$cell)
  meta <- data.frame(model = d$model, condition = condition,
                     row.names = d$cell, stringsAsFactors = FALSE)
  list(emt = emt, meta = meta)
}

# ---- condition parsing -----------------------------------------------------

test_that(".er_parse_condition reads keywords, numerics, and explicit labels", {
  kw <- .er_parse_condition(c("chemo-naive", "Relapsed", "treatment-naive", "CR"))
  expect_equal(as.character(kw), c("sensitive", "resistant", "sensitive", "resistant"))
  expect_identical(levels(kw), c("sensitive", "resistant"))

  num <- .er_parse_condition(c(0, 1, 0, 1))
  expect_equal(as.character(num), c("sensitive", "resistant", "sensitive", "resistant"))

  lab <- .er_parse_condition(c("baselineX", "LY2606368", "baselineX"),
                             sensitive = "baselineX", resistant = "LY2606368")
  expect_equal(as.character(lab), c("sensitive", "resistant", "sensitive"))
})

test_that(".er_parse_condition handles underscore-delimited tokens and character 0/1", {
  # Underscores must delimit tokens (\\b would treat "_" as a word char and miss these).
  us <- .er_parse_condition(c("SC4_CR", "SC4_CS", "SC16_cr", "SC16_naive"))
  expect_equal(as.character(us), c("resistant", "sensitive", "resistant", "sensitive"))
  # Character "0"/"1" (as TSV imports often produce) read like numeric 0/1.
  ch <- .er_parse_condition(c("0", "1", "0"))
  expect_equal(as.character(ch), c("sensitive", "resistant", "sensitive"))
  # A drug-named resistant arm without keywords stays unparseable (needs labels).
  expect_error(.er_parse_condition(c("SC4_baseline", "SC4_Talazoparib")), "Could not parse")
})

test_that(".er_parse_condition errors on unparseable values and half-specified labels", {
  expect_error(.er_parse_condition(c("naive", "weird-state")), "Could not parse")
  expect_error(.er_parse_condition(c(0, 2)), "0 .sensitive. / 1 .resistant.")
  expect_error(.er_parse_condition(c("a", "b"), sensitive = "a"), "together")
  # Overlapping explicit labels are rejected, not silently overwritten.
  expect_error(
    .er_parse_condition(c("x", "y"), sensitive = c("x", "y"), resistant = "y"),
    "share value"
  )
})

test_that(".er_parse_condition keeps NA conditions missing (all branches)", {
  # `%in%` returns FALSE for NA, so explicit-label parsing leaves NA as NA
  # (not silently coerced to a state) -- it is dropped later by prepare_*.
  expl <- .er_parse_condition(c("base", "drug", NA),
                              sensitive = "base", resistant = "drug")
  expect_equal(as.character(expl), c("sensitive", "resistant", NA))
  kw <- .er_parse_condition(c("untreated", NA, "relapsed"))
  expect_equal(as.character(kw), c("sensitive", NA, "resistant"))
})

# ---- prepare ---------------------------------------------------------------

test_that("prepare_resistance_emt assembles a tidy paired frame", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  expect_setequal(names(prepared), c("cell", "model", "condition", "emt"))
  expect_identical(levels(factor(prepared$condition)), c("sensitive", "resistant"))
  expect_equal(nrow(prepared), nrow(d))
  expect_length(attr(prepared, "unpaired_models"), 0)   # every model is paired
})

test_that("prepare_resistance_emt drops tiny groups and flags unpaired models", {
  d <- make_cdx_data(n_models = 3, n_per = 60)
  # Make model SC3 sensitive-only by relabelling its resistant cells, and shrink
  # SC2's resistant group below min_cells.
  d$condition[d$model == "SC3" & d$condition == "resistant"] <- "sensitive"
  drop_cells <- d$cell[d$model == "SC2" & d$condition == "resistant"][-(1:5)]
  d <- d[!d$cell %in% drop_cells, ]                      # SC2 resistant -> 5 cells
  inp <- as_inputs(d)
  prepared <- suppressWarnings(prepare_resistance_emt(inp$emt, inp$meta, min_cells = 20))
  expect_false(any(prepared$model == "SC2" & prepared$condition == "resistant"))
  expect_true("SC3" %in% attr(prepared, "unpaired_models"))
  expect_true("SC2" %in% attr(prepared, "unpaired_models"))  # lost its resistant arm
})

test_that("prepare_resistance_emt validates ids", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  expect_error(prepare_resistance_emt(unname(inp$emt), inp$meta), "named")
  expect_error(
    prepare_resistance_emt(inp$emt, inp$meta[, "model", drop = FALSE]),
    "Missing .* column"
  )
  expect_error(
    prepare_resistance_emt(stats::setNames(1:3, paste0("z", 1:3)), inp$meta),
    "No cell ids shared"
  )
})

# ---- location: EMT shifts up with resistance --------------------------------

test_that("emt_resistance_shift detects the upward EMT shift", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  res <- emt_resistance_shift(prepared)
  expect_equal(nrow(res$per_model), 6)
  expect_true(all(res$per_model$emt_delta > 0))         # every model goes up
  expect_equal(res$test$n_up, 6)
  expect_equal(res$test$method, "paired Wilcoxon signed-rank")
  expect_lt(res$test$p_value, 0.05)                     # 6/6 up, one-sided
  expect_gt(res$test$estimate, 0)
})

test_that("emt_resistance_shift omits the test but returns deltas with too few pairs", {
  d <- make_cdx_data(n_models = 2)
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  res <- suppressWarnings(emt_resistance_shift(prepared, min_models = 3))
  expect_equal(nrow(res$per_model), 2)
  expect_true(is.na(res$test$p_value))
  expect_true(is.na(res$test$method))
  expect_equal(res$test$n_models, 2)
})

# ---- dispersion: heterogeneity widens with resistance -----------------------

test_that("emt_heterogeneity_shift detects wider EMT spread under resistance", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  res <- emt_heterogeneity_shift(prepared, dispersion = "sd")
  expect_true(all(res$per_model$sd_delta > 0))
  expect_lt(res$test$p_value, 0.05)
  expect_equal(res$test$quantity, "sd")
  # robust dispersion option works and is labelled
  res_mad <- emt_heterogeneity_shift(prepared, dispersion = "mad")
  expect_equal(res_mad$test$quantity, "mad")
  expect_true("mad_delta" %in% names(res_mad$per_model))
})

# ---- composition: mesenchymal fraction rises with resistance ----------------

test_that("emt_state_composition detects the rise in mesenchymal-cell fraction", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  res <- emt_state_composition(prepared)
  expect_true(all(res$per_model$mes_frac_delta >= 0))
  # both fractions must lie within [0, 1] (check both bounds for both columns)
  fr <- c(res$per_model$mes_frac_sensitive, res$per_model$mes_frac_resistant)
  expect_true(all(fr >= 0 & fr <= 1))
  expect_lt(res$test$p_value, 0.05)
  expect_false(is.null(attr(res$test, "threshold")))
})

test_that("emt_dispersion_permutation flags a real dispersion shift, not a null one", {
  inp <- as_inputs(make_cdx_data())                      # resistant sd 1.7 > sens 1.0
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  res <- emt_dispersion_permutation(prepared, dispersion = "sd", n_perm = 200, seed = 1)
  expect_equal(res$observed$n_pos, 6)                    # all 6 models widen
  expect_lt(res$p_pos, 0.05)
  expect_lt(res$p_mean, 0.05)
  # reproducible under a fixed seed
  res2 <- emt_dispersion_permutation(prepared, dispersion = "sd", n_perm = 200, seed = 1)
  expect_equal(res$p_mean, res2$p_mean)

  # null: equal dispersion (and no location shift) -> not significant
  null_inp <- as_inputs(make_cdx_data(loc = 0, sd_sens = 1, sd_res = 1, seed = 3))
  null_prep <- prepare_resistance_emt(null_inp$emt, null_inp$meta)
  null_res <- emt_dispersion_permutation(null_prep, dispersion = "sd", n_perm = 200, seed = 1)
  expect_gt(null_res$p_mean, 0.05)
})

test_that("emt_dispersion_downsample holds under equal-n subsampling", {
  inp <- as_inputs(make_cdx_data())
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  res <- emt_dispersion_downsample(prepared, dispersion = "sd", n_cells = 30,
                                   n_rep = 50, seed = 1)
  expect_equal(res$n_cells, 30)
  expect_equal(res$n_models, 6)
  expect_gt(res$frac_all_up, 0.8)                        # widening survives downsampling
  expect_gt(res$mean_delta, 0)
  # default n_cells is the smallest group; can't exceed it
  expect_error(emt_dispersion_downsample(prepared, n_cells = 10^6), "exceeds the smallest")
})

test_that("emt_state_composition validates threshold and prepared", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  expect_error(emt_state_composition(prepared, threshold = c(1, 2)), "single finite")
  expect_error(emt_resistance_shift(data.frame(a = 1)), "prepare_resistance_emt")
})

test_that("analyses accept a condition factor carrying unused levels", {
  d <- make_cdx_data()
  inp <- as_inputs(d)
  prepared <- prepare_resistance_emt(inp$emt, inp$meta)
  # Re-level with a spurious unused level, as can happen after subsetting.
  prepared$condition <- factor(as.character(prepared$condition),
                               levels = c("sensitive", "resistant", "unused"))
  expect_true(any(levels(prepared$condition) == "unused"))
  res <- emt_resistance_shift(prepared)                 # must not reject on the unused level
  expect_equal(nrow(res$per_model), 6)
})
