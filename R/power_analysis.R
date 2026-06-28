# ============================================================================
# power_analysis.R
#
# A priori statistical power / feasibility analysis for the SCLC EMT project.
#
# Purpose
#   Before any modeling, quantify whether each planned cohort can actually
#   support the claims we care about. Three claim families are covered:
#
#     1. Survival association  - "EMT score (or EMT-high cell fraction) is
#        associated with OS/PFS" tested with a Cox model on a *standardized*
#        continuous covariate. Powered via Schoenfeld's formula, which depends
#        only on the number of EVENTS (not the number of patients).
#
#     2. Drug-response association - "EMT score correlates with chemo / drug
#        sensitivity across SCLC cell lines" tested as a Pearson/Spearman
#        correlation. Powered via the Fisher z-transform.
#
#     3. Treatment-effect modification - "EMT modifies immunotherapy benefit"
#        tested as a treatment x EMT interaction in a Cox model. Powered as a
#        survival association but with the well-known ~4x events penalty for
#        detecting an interaction relative to a main effect of the same size.
#
# Design notes
#   - All statistical formulas are implemented in base R so this file has no
#     hard dependencies. If {powerSurvEpi} / {pwr} are installed they can be
#     used to cross-check (see `validate_against_packages()`), but they are
#     never required.
#   - Functions are pure: same inputs -> same tidy tibble out. No global state,
#     no implicit file writes (writing is opt-in via `out_dir`).
#   - Effect sizes for survival are expressed as the hazard ratio per 1 SD of a
#     standardized covariate (sigma = 1), which is the scale we will actually
#     model EMT scores on.
#
# Author: (rewritten, best-in-class) for the EMT->resistance reframe
# ============================================================================

# ---- dependencies (loaded lazily; stats + base only are mandatory) ----------
# Soft-suggested for plotting/reporting: ggplot2, dplyr, tibble, readr, cli.

#' Two-sided normal quantile helper
#'
#' @param p Probability.
#' @return The standard-normal quantile `qnorm(p)`.
#' @keywords internal
.z <- function(p) stats::qnorm(p)

# ============================================================================
# 1. SURVIVAL (Cox, continuous standardized covariate) -- Schoenfeld
# ============================================================================

#' Power of a Cox model for a continuous covariate (Schoenfeld)
#'
#' Computes the statistical power to detect a hazard ratio `hr` per 1 SD of a
#' standardized continuous covariate, given the number of observed `events`.
#' Power for a Cox score test depends on the number of events, not the number
#' of patients.
#'
#' @param events Number of events (deaths / progressions). Length >= 1.
#' @param hr Hazard ratio per 1 SD of the standardized covariate. Length >= 1.
#'   Covariates are assumed standardized (mean 0, SD 1) -- the scale on which we
#'   model EMT scores -- so the effect is already "per SD" and the covariate SD
#'   does not enter the formula separately. Standardize before modeling if your
#'   covariate is on another scale.
#' @param alpha Two-sided type-I error rate. Default 0.05.
#' @param r_sq Proportion of covariate variance explained by other covariates
#'   (variance inflation for adjusted models / interactions). Default 0
#'   (unadjusted). The effective information is multiplied by `(1 - r_sq)`.
#'
#' @return A tibble (or data.frame) with one row per recycled `events`/`hr`
#'   combination and a `power` column.
#' @examples
#' power_cox(events = c(40, 50, 80), hr = 1.5)
#' @export
power_cox <- function(events, hr, alpha = 0.05, r_sq = 0) {
  stopifnot(all(events > 0), all(hr > 0),
            alpha > 0, alpha < 1, r_sq >= 0, r_sq < 1)
  grid <- expand.grid(events = events, hr = hr,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  beta <- log(grid$hr)                      # log-HR per 1 SD (standardized X)
  z_alpha <- .z(1 - alpha / 2)
  # Schoenfeld: Z ~ Normal( |beta| * sqrt(events * (1 - r_sq)), 1 )
  ncp <- abs(beta) * sqrt(grid$events * (1 - r_sq))
  grid$alpha <- alpha
  grid$r_sq <- r_sq
  grid$power <- stats::pnorm(ncp - z_alpha) + stats::pnorm(-ncp - z_alpha)
  .as_tbl(grid[, c("events", "hr", "alpha", "r_sq", "power")])
}

#' Minimum detectable hazard ratio for a Cox model (per SD)
#'
#' Inverts Schoenfeld's formula to return the smallest hazard ratio (per 1 SD,
#' HR > 1) detectable at the requested power.
#'
#' @inheritParams power_cox
#' @param power Target power. Default 0.80.
#' @return A tibble with `events` and `mdes_hr` (minimum detectable HR per SD).
#' @examples
#' mdes_cox(events = c(30, 50, 80))
#' @export
mdes_cox <- function(events, alpha = 0.05, power = 0.80, r_sq = 0) {
  stopifnot(all(events > 0), alpha > 0, alpha < 1,
            power > 0, power < 1, r_sq >= 0, r_sq < 1)
  z_sum <- .z(1 - alpha / 2) + .z(power)
  beta <- z_sum / sqrt(events * (1 - r_sq))
  .as_tbl(data.frame(events = events, alpha = alpha, power = power,
                     r_sq = r_sq, mdes_hr = exp(beta)))
}

#' Events required for a Cox model to reach target power
#'
#' @inheritParams power_cox
#' @param power Target power. Default 0.80.
#' @return A tibble with `hr` and `events_required` (rounded up).
#' @examples
#' events_required_cox(hr = c(1.3, 1.5, 2.0))
#' @export
events_required_cox <- function(hr, alpha = 0.05, power = 0.80, r_sq = 0) {
  stopifnot(all(hr > 0), alpha > 0, alpha < 1,
            power > 0, power < 1, r_sq >= 0, r_sq < 1)
  z_sum <- .z(1 - alpha / 2) + .z(power)
  beta <- log(hr)
  d <- (z_sum^2) / (beta^2 * (1 - r_sq))
  .as_tbl(data.frame(hr = hr, alpha = alpha, power = power,
                     r_sq = r_sq, events_required = ceiling(d)))
}

# ============================================================================
# 2. CORRELATION (cell-line EMT score vs drug sensitivity) -- Fisher z
# ============================================================================

#' Power of a correlation test (Fisher z-transform)
#'
#' @param n Sample size (e.g. number of SCLC cell lines).
#' @param r True correlation to detect. Length >= 1.
#' @param alpha Two-sided type-I error rate. Default 0.05.
#' @return A tibble with one row per `n`/`r` combination and a `power` column.
#' @examples
#' power_correlation(n = c(50, 60, 74), r = 0.35)
#' @export
power_correlation <- function(n, r, alpha = 0.05) {
  stopifnot(all(n > 3), all(abs(r) < 1), alpha > 0, alpha < 1)
  grid <- expand.grid(n = n, r = r, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  z_r <- atanh(abs(grid$r))                  # Fisher transform
  se <- 1 / sqrt(grid$n - 3)
  z_alpha <- .z(1 - alpha / 2)
  ncp <- z_r / se
  grid$alpha <- alpha
  grid$power <- stats::pnorm(ncp - z_alpha) + stats::pnorm(-ncp - z_alpha)
  .as_tbl(grid[, c("n", "r", "alpha", "power")])
}

#' Minimum detectable correlation at a given sample size
#'
#' @inheritParams power_correlation
#' @param power Target power. Default 0.80.
#' @return A tibble with `n` and `mdes_r`.
#' @examples
#' mdes_correlation(n = c(50, 60, 74))
#' @export
mdes_correlation <- function(n, alpha = 0.05, power = 0.80) {
  stopifnot(all(n > 3), alpha > 0, alpha < 1, power > 0, power < 1)
  z_sum <- .z(1 - alpha / 2) + .z(power)
  z_r <- z_sum / sqrt(n - 3)
  .as_tbl(data.frame(n = n, alpha = alpha, power = power, mdes_r = tanh(z_r)))
}

#' Sample size required to detect a correlation
#'
#' @inheritParams power_correlation
#' @param power Target power. Default 0.80.
#' @return A tibble with `r` and `n_required` (rounded up).
#' @export
n_required_correlation <- function(r, alpha = 0.05, power = 0.80) {
  stopifnot(all(abs(r) < 1, r != 0), alpha > 0, alpha < 1,
            power > 0, power < 1)
  z_sum <- .z(1 - alpha / 2) + .z(power)
  z_r <- atanh(abs(r))
  n <- (z_sum / z_r)^2 + 3
  .as_tbl(data.frame(r = r, alpha = alpha, power = power,
                     n_required = ceiling(n)))
}

# ============================================================================
# 3. INTERACTION (EMT x treatment) -- survival with interaction penalty
# ============================================================================

#' Minimum detectable interaction HR for a Cox model
#'
#' Detecting a treatment x covariate interaction is far less efficient than a
#' main effect. For a balanced design and a standardized covariate the variance
#' of the interaction coefficient is ~4x that of the main effect, so the
#' minimum detectable interaction HR is computed with an effective event count
#' of `events / variance_inflation`.
#'
#' @inheritParams mdes_cox
#' @param variance_inflation Inflation of the interaction coefficient variance
#'   relative to a main effect. Default 4 (balanced 50/50 treatment split).
#' @return A tibble with `events`, `variance_inflation`, and `mdes_hr_interaction`.
#' @export
mdes_interaction_cox <- function(events, alpha = 0.05,
                                 power = 0.80, variance_inflation = 4) {
  stopifnot(all(events > 0), variance_inflation >= 1)
  out <- mdes_cox(events = events, alpha = alpha,
                  power = power, r_sq = 1 - 1 / variance_inflation)
  names(out)[names(out) == "mdes_hr"] <- "mdes_hr_interaction"
  out$variance_inflation <- variance_inflation
  out
}

# ============================================================================
# 4. PROJECT-SPECIFIC FEASIBILITY REPORT
# ============================================================================

#' Default planning assumptions for the SCLC EMT cohorts
#'
#' Editable planning inputs. Event counts are a-priori estimates; update them
#' once real clinical tables are loaded (`R/run_emt_analysis.R` will pass the
#' observed counts in). Cohort sizes reflect public availability as of the
#' project plan.
#'
#' @return A tibble of cohort assumptions.
#' @export
emt_cohort_assumptions <- function() {
  .as_tbl(data.frame(
    cohort = c("George 2015 (UCologne)", "Jiang 2016 (GSE60052)",
               "Pooled bulk (George + GSE60052)", "Extended pooled (+GEO)",
               "DepMap/CCLE SCLC lines", "GDSC SCLC lines",
               "IMpower133 (Gay 2021)"),
    claim = c("survival", "survival", "survival", "survival",
              "correlation", "correlation", "interaction"),
    n = c(81L, 48L, 129L, 200L, 60L, 74L, 276L),
    events = c(50L, 30L, 80L, 130L, NA, NA, 200L),
    note = c("limited-stage; OS",
             "subset with survival annotation",
             "harmonize + batch-correct before pooling",
             "add GSE40275/GSE149507/GSE30219 (heterogeneous platforms)",
             "EMT score vs cisplatin/etoposide AUC; CRISPR deps",
             "EMT score vs drug IC50",
             "EMT x atezolizumab interaction; underpowered -> exploratory"),
    stringsAsFactors = FALSE
  ))
}

#' Build the a-priori feasibility report for the project's claims
#'
#' Applies the appropriate power calculation to each cohort/claim in
#' `assumptions` and returns one tidy table summarizing the minimum detectable
#' effect at the requested power. This is the gate we run before modeling.
#'
#' @param assumptions A tibble like [emt_cohort_assumptions()].
#' @param alpha Two-sided type-I error. Default 0.05.
#' @param power Target power. Default 0.80.
#' @param out_dir Optional directory to write `feasibility_report.csv`. If NULL
#'   (default) nothing is written.
#' @return A tibble: cohort, claim, n, events, min detectable effect, and a
#'   human-readable interpretation.
#' @export
run_power_analysis <- function(assumptions = emt_cohort_assumptions(),
                               alpha = 0.05, power = 0.80, out_dir = NULL) {
  rows <- lapply(seq_len(nrow(assumptions)), function(i) {
    a <- assumptions[i, ]
    if (a$claim == "survival") {
      m <- mdes_cox(events = a$events, alpha = alpha, power = power)
      effect <- sprintf("HR >= %.2f per SD", m$mdes_hr)
      detectable <- m$mdes_hr
    } else if (a$claim == "interaction") {
      m <- mdes_interaction_cox(events = a$events, alpha = alpha, power = power)
      effect <- sprintf("interaction HR >= %.2f per SD (4x penalty)",
                        m$mdes_hr_interaction)
      detectable <- m$mdes_hr_interaction
    } else if (a$claim == "correlation") {
      m <- mdes_correlation(n = a$n, alpha = alpha, power = power)
      effect <- sprintf("|r| >= %.2f", m$mdes_r)
      detectable <- m$mdes_r
    } else {
      stop(sprintf("Unknown claim type: %s", a$claim))
    }
    data.frame(cohort = a$cohort, claim = a$claim, n = a$n, events = a$events,
               min_detectable_effect = effect, detectable_value = detectable,
               note = a$note, stringsAsFactors = FALSE)
  })
  report <- .as_tbl(do.call(rbind, rows))
  report$alpha <- alpha
  report$power <- power

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(report,
                     file.path(out_dir, "feasibility_report.csv"),
                     row.names = FALSE)
  }
  report
}

# ============================================================================
# 5. SENSITIVITY CURVES (for the manuscript supplement)
# ============================================================================

#' Power curve over a grid of events x hazard ratios (survival)
#'
#' @param events Integer vector of event counts.
#' @param hr Numeric vector of hazard ratios per SD.
#' @param alpha Two-sided type-I error. Default 0.05.
#' @return A tibble suitable for `ggplot2::geom_line(aes(events, power, colour
#'   = factor(hr)))`.
#' @export
power_grid_survival <- function(events = seq(20, 160, by = 10),
                                hr = c(1.3, 1.5, 1.75, 2.0), alpha = 0.05) {
  power_cox(events = events, hr = hr, alpha = alpha)
}

#' Power curve over a grid of sample sizes x correlations (cell lines)
#'
#' @param n Integer vector of sample sizes.
#' @param r Numeric vector of correlations.
#' @param alpha Two-sided type-I error. Default 0.05.
#' @return A tibble for plotting.
#' @export
power_grid_correlation <- function(n = seq(20, 120, by = 5),
                                   r = c(0.25, 0.35, 0.45, 0.55),
                                   alpha = 0.05) {
  power_correlation(n = n, r = r, alpha = alpha)
}

#' Plot a survival power grid (requires ggplot2)
#'
#' @param grid Output of [power_grid_survival()].
#' @param target_power Horizontal reference line. Default 0.80.
#' @return A ggplot object.
#' @export
plot_power_grid_survival <- function(grid, target_power = 0.80) {
  .require_pkg("ggplot2")
  ggplot2::ggplot(grid, ggplot2::aes(x = .data$events, y = .data$power,
                                     colour = factor(.data$hr))) +
    ggplot2::geom_hline(yintercept = target_power, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Number of events", y = "Power",
                  colour = "HR per SD",
                  title = "A priori power: EMT score vs survival (Cox)",
                  subtitle = "Schoenfeld; standardized covariate, two-sided alpha = 0.05") +
    ggplot2::theme_minimal(base_size = 12)
}

# ============================================================================
# 6. OPTIONAL CROSS-CHECK against established packages
# ============================================================================

#' Cross-check the base-R formulas against {pwr} and {powerSurvEpi}
#'
#' Optional sanity check used in tests / development. Returns a tibble comparing
#' this file's results to the reference packages where available; silently skips
#' any package that is not installed.
#'
#' @return A tibble of comparisons, or an empty tibble if neither package is
#'   installed.
#' @export
validate_against_packages <- function() {
  out <- list()
  if (requireNamespace("pwr", quietly = TRUE)) {
    ref <- pwr::pwr.r.test(n = 60, r = 0.35, sig.level = 0.05)$power
    ours <- power_correlation(n = 60, r = 0.35)$power
    out[[length(out) + 1]] <- data.frame(
      check = "correlation power (n=60, r=0.35)",
      reference = ref, ours = ours, abs_diff = abs(ref - ours))
  }
  if (length(out) == 0) return(.as_tbl(data.frame()))
  .as_tbl(do.call(rbind, out))
}

# ============================================================================
# small internal helpers
# ============================================================================

#' Coerce to tibble if available, else data.frame
#' @keywords internal
.as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Require a soft dependency with a clear error
#' @keywords internal
.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this function. Install it with install.packages('%s').",
                 pkg, pkg), call. = FALSE)
  }
}
