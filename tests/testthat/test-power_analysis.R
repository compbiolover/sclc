# tests/testthat/test-power_analysis.R
#
# Tests for the a-priori power / feasibility module. Validates the statistical
# formulas against closed-form expectations, internal consistency (round-trips),
# expected monotonicity, and -- when available -- the {pwr} reference package.

# Make the test robust whether the project is loaded as a package or sourced.
if (!exists("power_cox", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "power_analysis.R"))
}

test_that("mdes_cox matches the closed-form Schoenfeld value", {
  # exp((z_{0.975} + z_{0.80}) / sqrt(events)) with sigma = 1
  expected <- exp((qnorm(0.975) + qnorm(0.80)) / sqrt(50))
  got <- mdes_cox(events = 50)$mdes_hr
  expect_equal(got, expected, tolerance = 1e-8)
  # sanity: ~1.49 per SD at 50 events
  expect_equal(round(got, 2), 1.49)
})

test_that("events_required_cox inverts power_cox / mdes_cox (round-trip)", {
  hr <- 1.5
  d <- events_required_cox(hr = hr)$events_required
  # With exactly that many events the power should be >= target (rounding up).
  expect_gte(power_cox(events = d, hr = hr)$power, 0.80)
  # And the MDES at d events should be at or below the HR we solved for.
  expect_lte(mdes_cox(events = d)$mdes_hr, hr + 1e-6)
})

test_that("power_cox is monotonically increasing in events and in effect size", {
  p_events <- power_cox(events = c(20, 40, 80, 160), hr = 1.5)$power
  expect_true(all(diff(p_events) > 0))
  p_hr <- power_cox(events = 50, hr = c(1.2, 1.5, 2.0))$power
  expect_true(all(diff(p_hr) > 0))
  expect_true(all(p_events >= 0 & p_events <= 1))
})

test_that("correlation power and MDES are internally consistent", {
  # n = 60, r = 0.35 -> ~0.79 power (Fisher z)
  expect_equal(round(power_correlation(n = 60, r = 0.35)$power, 2), 0.79)
  # MDES round-trips through n_required_correlation
  r_min <- mdes_correlation(n = 60)$mdes_r
  n_back <- n_required_correlation(r = r_min)$n_required
  expect_lte(abs(n_back - 60), 1)
})

test_that("interaction MDES is larger (harder) than the main-effect MDES", {
  main <- mdes_cox(events = 200)$mdes_hr
  inter <- mdes_interaction_cox(events = 200)$mdes_hr_interaction
  expect_gt(inter, main)
  # 4x variance inflation -> ~2x the log-HR requirement
  expect_equal(log(inter) / log(main), 2, tolerance = 0.05)
})

test_that("run_power_analysis returns one row per cohort with effect strings", {
  rep <- run_power_analysis()
  expect_equal(nrow(rep), nrow(emt_cohort_assumptions()))
  expect_true(all(c("cohort", "claim", "min_detectable_effect",
                    "detectable_value") %in% names(rep)))
  expect_false(any(is.na(rep$detectable_value)))
})

test_that("input validation rejects nonsense", {
  expect_error(power_cox(events = -1, hr = 1.5))
  expect_error(power_cox(events = 50, hr = -2))
  expect_error(power_correlation(n = 2, r = 0.3))   # n must exceed 3
  expect_error(mdes_correlation(n = 60, power = 1.5))
})

test_that("base-R formulas agree with {pwr} when installed", {
  skip_if_not_installed("pwr")
  cmp <- validate_against_packages()
  expect_true(all(cmp$abs_diff < 1e-3))
})
