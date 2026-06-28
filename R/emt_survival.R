# ============================================================================
# emt_survival.R
#
# WS3 (survival arm): connect the validated EMT cell-state axis (from
# R/emt_scoring.R -> compute_emt_scores()$consensus) to clinical outcome.
#
#   prepare_survival()  join EMT scores + clinical, parse the event indicator
#   emt_cox()           Cox model for the EMT axis (HR per SD, CI, p, C-index)
#   emt_km()            Kaplan-Meier by EMT tertile state + log-rank test
#   emt_added_value()   does EMT add to a base model? (LRT + delta C-index)
#   pool_emt_survival() pooled multi-cohort Cox, baseline hazard stratified by
#                       cohort (the WS0 power analysis said to pool for power)
#
# The EMT axis is standardized (mean 0, SD 1) before modeling, so hazard ratios
# are "per 1 SD of the EMT consensus" -- the scale the WS0 feasibility analysis
# was powered on. Uses base {survival}; plotting helpers are optional.
#
# Conventions match the other modules: pure functions, {cli} errors, tidy out.
# ============================================================================

#' tibble if available else data.frame
#' @keywords internal
.sv_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Parse a survival event indicator to 0/1
#'
#' Accepts numeric 0/1, logical, or cBioPortal-style strings like
#' "1:DECEASED" / "0:LIVING" / "DECEASED" / "LIVING".
#' @keywords internal
.sv_parse_event <- function(x) {
  if (is.numeric(x)) return(as.integer(x != 0))
  if (is.logical(x)) return(as.integer(x))
  s <- toupper(trimws(as.character(x)))
  ev <- rep(NA_integer_, length(s))
  ev[grepl("^1|DECEASED|DEAD|EVENT|PROGRESS|RELAPSE|YES|TRUE", s)] <- 1L
  ev[grepl("^0|LIVING|ALIVE|CENSOR|NO|FALSE", s)] <- 0L
  ev
}

# ============================================================================
# 1. Assemble the survival frame
# ============================================================================

#' Join EMT scores with clinical data into a survival-analysis frame
#'
#' @param emt_scores Output of `compute_emt_scores()` (needs `sample`,
#'   `consensus`; `emt_state` carried through if present).
#' @param clinical A data.frame of clinical data with a sample id, a time
#'   column, and an event column.
#' @param sample_col,time_col,event_col Column names in `clinical`.
#' @param covariates Optional character vector of extra `clinical` columns to
#'   retain as model covariates (e.g. stage, subtype).
#' @return A tibble with `sample`, `emt_consensus`, `emt_state` (if present),
#'   `.time`, `.event` (0/1), and any requested covariates; rows with missing
#'   time/event dropped.
#' @export
prepare_survival <- function(emt_scores, clinical, sample_col = "sample",
                             time_col, event_col, covariates = NULL) {
  if (!all(c("sample", "consensus") %in% names(emt_scores))) {
    cli::cli_abort("{.arg emt_scores} must contain `sample` and `consensus`.")
  }
  need <- c(sample_col, time_col, event_col, covariates)
  miss <- setdiff(need, names(clinical))
  if (length(miss) > 0) cli::cli_abort("Missing clinical column(s): {.val {miss}}.")

  e <- data.frame(sample = as.character(emt_scores$sample),
                  emt_consensus = emt_scores$consensus,
                  stringsAsFactors = FALSE)
  if ("emt_state" %in% names(emt_scores)) e$emt_state <- as.character(emt_scores$emt_state)
  cl <- clinical[, need, drop = FALSE]
  names(cl)[1:3] <- c("sample", ".time", ".event")
  cl$sample <- as.character(cl$sample)
  cl$.time <- suppressWarnings(as.numeric(cl$.time))
  cl$.event <- .sv_parse_event(cl$.event)

  df <- merge(e, cl, by = "sample")
  ok <- is.finite(df$.time) & df$.time >= 0 & !is.na(df$.event)
  dropped <- sum(!ok)
  if (dropped > 0) cli::cli_warn("Dropping {dropped} sample(s) with missing/invalid time or event.")
  df <- df[ok, , drop = FALSE]
  if (nrow(df) == 0) cli::cli_abort("No samples with usable survival data after the join.")
  .sv_as_tbl(df)
}

# ============================================================================
# 2. Cox model for the EMT axis
# ============================================================================

#' Cox model for the EMT consensus axis
#'
#' Fits `Surv(time, event) ~ EMT [+ covariates]` with the EMT axis standardized,
#' so the hazard ratio is per 1 SD of the EMT consensus.
#'
#' @param data A frame from [prepare_survival()] (`.time`, `.event`,
#'   `emt_consensus`).
#' @param score Column to model. Default "emt_consensus".
#' @param covariates Optional character vector of covariate columns.
#' @param standardize Standardize `score` to per-SD. Default TRUE.
#' @return A list: `coef` (tidy: term, hr, ci_low, ci_high, p), `n`, `events`,
#'   `c_index`, and the `fit` (coxph object).
#' @export
emt_cox <- function(data, score = "emt_consensus", covariates = NULL,
                    standardize = TRUE) {
  .sv_require("survival")
  for (col in c(".time", ".event", score, covariates)) {
    if (!col %in% names(data)) cli::cli_abort("Column {.val {col}} not found; run prepare_survival() first.")
  }
  d <- as.data.frame(data)
  d <- d[stats::complete.cases(d[, c(".time", ".event", score, covariates)]), , drop = FALSE]
  if (sum(d$.event) < 5) {
    cli::cli_warn("Only {sum(d$.event)} events; Cox estimates will be unstable (see WS0 power analysis).")
  }
  if (standardize) d[[score]] <- as.numeric(scale(d[[score]]))
  rhs <- paste(c(score, covariates), collapse = " + ")
  fit <- survival::coxph(stats::as.formula(sprintf("survival::Surv(.time, .event) ~ %s", rhs)),
                         data = d)
  s <- summary(fit)
  coef <- data.frame(
    term = rownames(s$coefficients),
    hr = s$coefficients[, "exp(coef)"],
    ci_low = s$conf.int[, "lower .95"],
    ci_high = s$conf.int[, "upper .95"],
    p = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL, stringsAsFactors = FALSE
  )
  list(coef = .sv_as_tbl(coef), n = fit$n, events = fit$nevent,
       c_index = unname(survival::concordance(fit)$concordance), fit = fit)
}

# ============================================================================
# 3. Kaplan-Meier by EMT tertile state
# ============================================================================

#' Kaplan-Meier survival by EMT state (E / hybrid / M) + log-rank
#'
#' @param data A frame from [prepare_survival()].
#' @param group Grouping column. Default "emt_state".
#' @return A list: `fit` (survfit), `logrank_p`, `n_by_group`,
#'   `median_survival` (per group).
#' @export
emt_km <- function(data, group = "emt_state") {
  .sv_require("survival")
  if (!group %in% names(data)) cli::cli_abort("Group column {.val {group}} not found.")
  d <- as.data.frame(data)
  d <- d[!is.na(d[[group]]) & stats::complete.cases(d[, c(".time", ".event")]), , drop = FALSE]
  d[[group]] <- factor(d[[group]])
  f <- stats::as.formula(sprintf("survival::Surv(.time, .event) ~ %s", group))
  fit <- survival::survfit(f, data = d)
  sd_fit <- survival::survdiff(f, data = d)
  logrank_p <- stats::pchisq(sd_fit$chisq, df = length(sd_fit$n) - 1, lower.tail = FALSE)
  med <- summary(fit)$table
  med_surv <- if (is.matrix(med)) med[, "median"] else stats::setNames(med["median"], levels(d[[group]]))
  list(fit = fit, logrank_p = logrank_p,
       n_by_group = table(d[[group]]), median_survival = med_surv)
}

# ============================================================================
# 4. Added value over a base model
# ============================================================================

#' Does the EMT axis add prognostic value beyond a base model?
#'
#' Compares a base Cox model (covariates only) with base + EMT via a likelihood-
#' ratio test and the change in C-index.
#'
#' @param data A frame from [prepare_survival()].
#' @param score EMT column. Default "emt_consensus".
#' @param base_covariates Character vector of base-model covariates (>= 1).
#' @param standardize Standardize `score`. Default TRUE.
#' @return A list: `lrt_p`, `c_index_base`, `c_index_full`, `delta_c_index`.
#' @export
emt_added_value <- function(data, score = "emt_consensus", base_covariates,
                            standardize = TRUE) {
  .sv_require("survival")
  if (length(base_covariates) < 1) cli::cli_abort("Provide >= 1 base covariate.")
  d <- as.data.frame(data)
  d <- d[stats::complete.cases(d[, c(".time", ".event", score, base_covariates)]), , drop = FALSE]
  if (standardize) d[[score]] <- as.numeric(scale(d[[score]]))
  surv <- "survival::Surv(.time, .event)"
  base <- survival::coxph(stats::as.formula(paste(surv, "~", paste(base_covariates, collapse = " + "))), data = d)
  full <- survival::coxph(stats::as.formula(paste(surv, "~", paste(c(base_covariates, score), collapse = " + "))), data = d)
  lrt <- stats::anova(base, full, test = "LRT")
  list(lrt_p = lrt$`Pr(>|Chi|)`[2],
       c_index_base = unname(survival::concordance(base)$concordance),
       c_index_full = unname(survival::concordance(full)$concordance),
       delta_c_index = unname(survival::concordance(full)$concordance -
                              survival::concordance(base)$concordance))
}

# ============================================================================
# 5. Pooled multi-cohort Cox (stratified baseline hazard)
# ============================================================================

#' Pooled EMT survival across cohorts (cohort-stratified baseline hazard)
#'
#' The WS0 power analysis showed single SCLC cohorts only detect HR >= ~1.5/SD;
#' pooling is required for adequate power. This standardizes the EMT axis WITHIN
#' each cohort (so scores are comparable) and fits a Cox model with a separate
#' baseline hazard per cohort via `strata(cohort)`.
#'
#' @param data A frame from [prepare_survival()] with an extra cohort column.
#' @param cohort Cohort column name.
#' @param score EMT column. Default "emt_consensus".
#' @param covariates Optional extra covariates.
#' @return Same shape as [emt_cox()] (HR is per within-cohort SD of EMT).
#' @export
pool_emt_survival <- function(data, cohort, score = "emt_consensus",
                              covariates = NULL) {
  .sv_require("survival")
  if (!cohort %in% names(data)) cli::cli_abort("Cohort column {.val {cohort}} not found.")
  d <- as.data.frame(data)
  d <- d[stats::complete.cases(d[, c(".time", ".event", score, cohort, covariates)]), , drop = FALSE]
  # standardize EMT within each cohort
  d[[score]] <- stats::ave(d[[score]], d[[cohort]],
                           FUN = function(x) as.numeric(scale(x)))
  rhs <- paste(c(score, covariates, sprintf("survival::strata(%s)", cohort)), collapse = " + ")
  fit <- survival::coxph(stats::as.formula(sprintf("survival::Surv(.time, .event) ~ %s", rhs)), data = d)
  s <- summary(fit)
  coef <- data.frame(
    term = rownames(s$coefficients), hr = s$coefficients[, "exp(coef)"],
    ci_low = s$conf.int[, "lower .95"], ci_high = s$conf.int[, "upper .95"],
    p = s$coefficients[, "Pr(>|z|)"], row.names = NULL, stringsAsFactors = FALSE
  )
  list(coef = .sv_as_tbl(coef), n = fit$n, events = fit$nevent,
       c_index = unname(survival::concordance(fit)$concordance), fit = fit,
       n_cohorts = length(unique(d[[cohort]])))
}

#' Require a package with an actionable error
#' @keywords internal
.sv_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cli::cli_abort(c("x" = "Package {.pkg {pkg}} is required.",
                     "i" = "Install it with {.code install.packages('{pkg}')}."))
  }
}
