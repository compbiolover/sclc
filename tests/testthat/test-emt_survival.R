# tests/testthat/test-emt_survival.R
#
# Tests for the WS3 survival module on synthetic data where higher EMT causes
# higher hazard (HR per SD ~ exp(0.6) = 1.82). No real cohort data needed.

if (!exists("emt_cox", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "emt_survival.R"))
}

make_surv <- function(n = 200) {
  set.seed(7)
  emt <- stats::rnorm(n)
  t_event <- stats::rexp(n, rate = exp(0.6 * emt))   # hazard rises with EMT
  t_cens <- stats::rexp(n, rate = 0.25)
  time <- pmin(t_event, t_cens); event <- as.integer(t_event <= t_cens)
  state <- cut(emt, stats::quantile(emt, c(0, 1/3, 2/3, 1)),
               labels = c("E", "hybrid", "M"), include.lowest = TRUE)
  list(
    emt = data.frame(sample = paste0("s", seq_len(n)), consensus = emt,
                     emt_state = as.character(state), stringsAsFactors = FALSE),
    clin = data.frame(sample = paste0("s", seq_len(n)), os_time = time,
                      os_status = ifelse(event == 1, "1:DECEASED", "0:LIVING"),
                      stage = sample(c("I", "II", "III"), n, replace = TRUE),
                      cohort = rep(c("A", "B"), length.out = n),
                      stringsAsFactors = FALSE)
  )
}

test_that(".sv_parse_event handles numeric, logical, and cBioPortal strings", {
  expect_equal(.sv_parse_event(c("1:DECEASED", "0:LIVING", "DECEASED", "LIVING")),
               c(1L, 0L, 1L, 0L))
  expect_equal(.sv_parse_event(c(1, 0, 2)), c(1L, 0L, 1L))
  expect_equal(.sv_parse_event(c(TRUE, FALSE)), c(1L, 0L))
})

test_that(".sv_parse_event prefers explicit prefix and matches whole words only", {
  expect_equal(.sv_parse_event(c("1:LIVING", "0:DECEASED")), c(1L, 0L))  # numeric prefix wins
  expect_true(is.na(.sv_parse_event("DEADLINE")))                        # 'DEAD' substring != whole word
  expect_equal(.sv_parse_event(c("Censored", "Relapse")), c(0L, 1L))
})

test_that("emt_cox / emt_added_value error on a constant EMT score", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status",
                        covariates = "stage")
  s$emt_consensus <- 1
  expect_error(emt_cox(s), "constant")
  expect_error(emt_added_value(s, base_covariates = "stage"), "constant")
})

test_that("emt_km errors without >= 2 groups or without events", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status")
  s1 <- s; s1$emt_state <- "M"
  expect_error(emt_km(s1), ">= 2 non-empty groups")
  s2 <- s; s2$.event <- 0L
  expect_error(emt_km(s2), "No events")
})

test_that("pool_emt_survival warns on a flat cohort but still fits", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status",
                        covariates = "cohort")
  s$emt_consensus[s$cohort == "B"] <- 0.5          # cohort B has no EMT variation
  res <- expect_warning(pool_emt_survival(s, cohort = "cohort"),
                        "no within-cohort EMT variation")
  expect_true(is.finite(res$coef$hr[res$coef$term == "emt_consensus"]))
})

test_that("prepare_survival joins and parses event/time", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status",
                        covariates = c("stage", "cohort"))
  expect_true(all(c("emt_consensus", ".time", ".event", "stage", "cohort") %in% names(s)))
  expect_true(all(s$.event %in% c(0L, 1L)))
  expect_equal(nrow(s), nrow(d$emt))
  expect_error(prepare_survival(d$emt, d$clin, time_col = "nope", event_col = "os_status"),
               "Missing clinical column")
})

test_that("prepare_survival guards against duplicate IDs and zero overlap", {
  d <- make_surv()
  dup <- d$emt; dup <- rbind(dup, dup[1, ])               # duplicate a sample id
  expect_error(prepare_survival(dup, d$clin, time_col = "os_time", event_col = "os_status"),
               "Duplicate sample IDs")
  nomatch <- d$clin; nomatch$sample <- paste0("X", nomatch$sample)  # no overlap
  expect_error(prepare_survival(d$emt, nomatch, time_col = "os_time", event_col = "os_status"),
               "No overlapping sample IDs")
})

test_that("emt_cox recovers the positive EMT-hazard association (HR>1)", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status")
  res <- emt_cox(s)
  row <- res$coef[res$coef$term == "emt_consensus", ]
  expect_gt(row$hr, 1)
  expect_lt(row$p, 0.05)
  expect_gt(res$c_index, 0.5)
  expect_gt(res$events, 5)
})

test_that("emt_km separates the EMT tertiles (log-rank)", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status")
  km <- emt_km(s)
  expect_length(km$n_by_group, 3)
  expect_lt(km$logrank_p, 0.05)
})

test_that("emt_added_value shows EMT adds over a base model", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status",
                        covariates = "stage")
  av <- emt_added_value(s, base_covariates = "stage")
  expect_lt(av$lrt_p, 0.05)
  expect_gt(av$delta_c_index, 0)
})

test_that("pool_emt_survival fits a cohort-stratified model", {
  d <- make_surv()
  s <- prepare_survival(d$emt, d$clin, time_col = "os_time", event_col = "os_status",
                        covariates = "cohort")
  res <- pool_emt_survival(s, cohort = "cohort")
  expect_equal(res$n_cohorts, 2)
  expect_gt(res$coef$hr[res$coef$term == "emt_consensus"], 1)
})
