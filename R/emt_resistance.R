# ============================================================================
# emt_resistance.R
#
# WS4 (dynamic-resistance arm): does EMT *arise in SCLC as a function of
# treatment*? Tests whether the validated per-cell EMT axis (from
# R/emt_scoring.R -> score_emt_singlecell(); higher = more mesenchymal) rises
# when an SCLC tumor evolves from chemo-sensitive to chemo-resistant.
#
# Motivating data: Stewart et al. 2020, Nat Cancer (PMID 33521652), CTC-derived
# xenografts (CDX) profiled by scRNA-seq as PAIRED chemo-sensitive vs
# chemo-resistant states of the same model. Their headline finding is that
# intratumoral HETEROGENEITY increases with resistance, with EMT named as one of
# the heterogeneous resistance programs. So we ask three complementary questions,
# not just one:
#
#   emt_resistance_shift()      LOCATION: does the EMT axis shift up (resistant
#                               more mesenchymal than sensitive)?
#   emt_heterogeneity_shift()   DISPERSION: does per-cell EMT heterogeneity widen
#                               with resistance (the paper's central claim)?
#   emt_state_composition()     COMPOSITION: does the mesenchymal-cell fraction
#                               increase with resistance?
#   prepare_resistance_emt()    assemble per-cell EMT + (model, condition) input.
#
# STATISTICAL STANCE (important). Cells within one model x condition are
# pseudoreplicates, not independent samples -- a per-cell test would report
# absurd significance off a handful of tumors. Every test here therefore
# AGGREGATES each model x condition to one number (a robust pseudobulk summary)
# and tests PAIRED ACROSS MODELS (resistant - sensitive). The model, not the
# cell, is the unit of replication. With the few CDX models such studies carry,
# we report the per-model deltas and sign-consistency alongside the paired test,
# because a Wilcoxon p-value on ~5 pairs is weak on its own.
#
# Conventions match the other modules: pure functions, {cli} errors, tidy out,
# self-contained (own small helpers) so it can be sourced alone.
# ============================================================================

# ---- internal helpers ------------------------------------------------------

#' tibble if available else data.frame
#' @keywords internal
.er_as_tbl <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Coerce a per-cell EMT input to a named numeric vector (cell -> score)
#' @keywords internal
.er_emt_vector <- function(emt, score = "consensus") {
  if (is.numeric(emt)) {
    if (is.null(names(emt))) {
      cli::cli_abort("{.arg emt_cells} numeric vector must be named by cell id.")
    }
    if (anyDuplicated(names(emt))) {
      cli::cli_abort("{.arg emt_cells} has duplicate cell ids: {.val {unique(names(emt)[duplicated(names(emt))])}}.")
    }
    return(emt)
  }
  if (is.data.frame(emt)) {
    if (!all(c("cell", score) %in% names(emt))) {
      cli::cli_abort("EMT data.frame must contain {.val cell} and {.val {score}} columns (or pass a named numeric vector).")
    }
    if (!is.numeric(emt[[score]])) cli::cli_abort("EMT {.val {score}} column must be numeric.")
    if (anyDuplicated(emt$cell)) {
      cli::cli_abort("EMT {.val cell} column has duplicate ids: {.val {unique(emt$cell[duplicated(emt$cell)])}}.")
    }
    return(stats::setNames(emt[[score]], as.character(emt$cell)))
  }
  cli::cli_abort("{.arg emt_cells} must be a named numeric vector or a data.frame with `cell` + `{score}`.")
}

#' Parse a treatment-condition vector to an ordered factor sensitive < resistant
#'
#' Accepts the labels SCLC resistance studies actually use. Numeric/logical are
#' read as 0 = sensitive, 1 = resistant. Explicit `sensitive`/`resistant` label
#' vectors override the keyword parsing entirely.
#' @keywords internal
.er_parse_condition <- function(x, sensitive = NULL, resistant = NULL) {
  s_raw <- as.character(x)
  nonNA <- s_raw[!is.na(s_raw)]
  if (!is.null(sensitive) || !is.null(resistant)) {
    if (is.null(sensitive) || is.null(resistant)) {
      cli::cli_abort("Pass {.arg sensitive_labels} and {.arg resistant_labels} together, or neither.")
    }
    overlap <- intersect(as.character(sensitive), as.character(resistant))
    if (length(overlap) > 0) {
      cli::cli_abort(c(
        "x" = "{.arg sensitive_labels} and {.arg resistant_labels} share value(s): {.val {overlap}}.",
        "i" = "A label can denote only one state; remove the overlap."
      ))
    }
    out <- rep(NA_character_, length(s_raw))
    out[s_raw %in% as.character(sensitive)] <- "sensitive"
    out[s_raw %in% as.character(resistant)] <- "resistant"
  } else if (is.logical(x) ||
             (length(nonNA) > 0 && all(grepl("^\\s*[01]\\s*$", nonNA)))) {
    # Logical, numeric 0/1, OR *character* "0"/"1" (TSV imports often yield
    # character columns) -> 1 = resistant, 0 = sensitive.
    xi <- if (is.logical(x)) as.integer(x) else suppressWarnings(as.integer(trimws(s_raw)))
    out <- ifelse(xi == 1L, "resistant", "sensitive")
  } else if (is.numeric(x)) {
    cli::cli_abort(c("x" = "Numeric condition must be 0 (sensitive) / 1 (resistant).",
                     "i" = "Saw values: {.val {sort(unique(x))}}. Pass labels via {.arg sensitive_labels}/{.arg resistant_labels} instead."))
  } else {
    # Keyword match on TOKENS split at any non-alphanumeric separator, so "_",
    # "-", "/", and spaces all delimit (e.g. "SC4_CR" -> tokens "sc4","cr").
    # `\\b` would NOT do this, since "_" is a regex word character.
    res_kw <- c("resistant", "resistance", "chemoresistant", "relapse",
                "relapsed", "refractory", "cr")
    sen_kw <- c("sensitive", "chemosensitive", "naive", "chemonaive",
                "untreated", "baseline", "pretreatment", "cs")
    classify <- function(v) {
      if (is.na(v) || !nzchar(v)) return(NA_character_)
      toks <- strsplit(v, "[^a-z0-9]+")[[1]]
      toks <- toks[nzchar(toks)]
      if (any(toks %in% res_kw)) return("resistant")
      if (any(toks %in% sen_kw)) return("sensitive")
      # Substring fallback for glued compound words (only the long, unambiguous
      # stems -- never the 2-letter cr/cs, which must stay token-exact).
      if (grepl("resist|relaps|refractor", v)) return("resistant")
      if (grepl("sensitiv|naive|untreat|baseline|pretreat", v)) return("sensitive")
      NA_character_
    }
    out <- vapply(tolower(trimws(s_raw)), classify, character(1), USE.NAMES = FALSE)
  }
  unparsed <- unique(s_raw[is.na(out) & !is.na(s_raw)])
  if (length(unparsed) > 0) {
    cli::cli_abort(c(
      "x" = "Could not parse {length(unparsed)} condition value(s) as sensitive/resistant: {.val {unparsed}}.",
      "i" = "Supply explicit {.arg sensitive_labels} / {.arg resistant_labels}."
    ))
  }
  factor(out, levels = c("sensitive", "resistant"))
}

#' Aggregate per-cell EMT to one value per model x condition with `fun`
#' @keywords internal
.er_aggregate <- function(prepared, fun) {
  ag <- stats::aggregate(emt ~ model + condition, data = prepared, FUN = fun)
  n  <- stats::aggregate(emt ~ model + condition, data = prepared,
                         FUN = function(z) sum(is.finite(z)))
  names(n)[3] <- "n"
  merge(ag, n, by = c("model", "condition"))
}

#' Paired model-level test of resistant vs sensitive, returned as a tidy row +
#' a per-model delta table. `value_name` labels the aggregated quantity.
#' @keywords internal
.er_paired_test <- function(agg, value_name, min_models = 3,
                            alternative = "greater") {
  agg$condition <- as.character(agg$condition)
  sens <- agg[agg$condition == "sensitive", c("model", "emt", "n"), drop = FALSE]
  res  <- agg[agg$condition == "resistant", c("model", "emt", "n"), drop = FALSE]
  if (nrow(sens) == 0 || nrow(res) == 0) {
    cli::cli_abort(c("x" = "Need both sensitive and resistant states to compare {value_name}.",
                     "i" = "Found {nrow(sens)} sensitive and {nrow(res)} resistant model-group(s)."))
  }
  names(sens) <- c("model", "v_sensitive", "n_sensitive")
  names(res)  <- c("model", "v_resistant", "n_resistant")
  # Inner-join keeps only models carrying BOTH states; an unpaired model cannot
  # contribute to a paired test and would bias a between-group comparison by
  # model composition.
  paired <- merge(sens, res, by = "model")
  paired <- paired[is.finite(paired$v_sensitive) & is.finite(paired$v_resistant), , drop = FALSE]
  n_models <- nrow(paired)
  per_model <- data.frame(
    model       = paired$model,
    n_sensitive = paired$n_sensitive,
    n_resistant = paired$n_resistant,
    sensitive   = paired$v_sensitive,
    resistant   = paired$v_resistant,
    delta       = paired$v_resistant - paired$v_sensitive,
    row.names = NULL, stringsAsFactors = FALSE
  )
  # Label the aggregated columns by the quantity tested (emt / sd / mes_frac ...).
  names(per_model)[4] <- paste0(value_name, "_sensitive")
  names(per_model)[5] <- paste0(value_name, "_resistant")
  names(per_model)[6] <- paste0(value_name, "_delta")

  n_up   <- sum(per_model[[paste0(value_name, "_delta")]] > 0)
  n_down <- sum(per_model[[paste0(value_name, "_delta")]] < 0)

  if (n_models < min_models) {
    cli::cli_warn(c(
      "!" = "Only {n_models} paired model(s) for {value_name} (< {min_models}); the paired test is omitted.",
      "i" = "Per-model deltas and sign counts are still returned."
    ))
    test <- data.frame(
      quantity = value_name, method = NA_character_, n_models = n_models,
      n_up = n_up, n_down = n_down, statistic = NA_real_,
      estimate = if (n_models > 0) stats::median(per_model[[paste0(value_name, "_delta")]]) else NA_real_,
      p_value = NA_real_, alternative = alternative,
      stringsAsFactors = FALSE
    )
    return(list(per_model = .er_as_tbl(per_model), test = .er_as_tbl(test)))
  }

  wt <- suppressWarnings(stats::wilcox.test(
    per_model[[paste0(value_name, "_resistant")]],
    per_model[[paste0(value_name, "_sensitive")]],
    paired = TRUE, alternative = alternative, conf.int = TRUE, exact = FALSE
  ))
  # Hodges-Lehmann median of paired deltas; fall back to the plain median if the
  # CI/estimate could not be computed (e.g. heavy ties).
  est <- if (!is.null(wt$estimate)) unname(wt$estimate) else
    stats::median(per_model[[paste0(value_name, "_delta")]])
  test <- data.frame(
    quantity = value_name, method = "paired Wilcoxon signed-rank",
    n_models = n_models, n_up = n_up, n_down = n_down,
    statistic = unname(wt$statistic),
    estimate = est,
    p_value = wt$p.value, alternative = alternative,
    stringsAsFactors = FALSE
  )
  list(per_model = .er_as_tbl(per_model), test = .er_as_tbl(test))
}

# ============================================================================
# 1. Assemble the per-cell resistance frame
# ============================================================================

#' Assemble a per-cell EMT-vs-resistance analysis frame
#'
#' Joins per-cell EMT scores (from `score_emt_singlecell()`) with per-cell
#' metadata carrying the model id and the treatment condition, returning the
#' tidy frame the WS4 analyses consume.
#'
#' @param emt_cells Named numeric vector (cell id -> EMT score; higher = more
#'   mesenchymal) or a data.frame with `cell` + a score column.
#' @param cell_meta A data.frame of per-cell metadata.
#' @param model_col,condition_col,cell_col Column names in `cell_meta` for the
#'   model id, the treatment condition, and the cell id. `cell_col` defaults to
#'   rownames if `NULL`.
#' @param sensitive_labels,resistant_labels Optional character vectors giving the
#'   exact `condition_col` values that denote each state, overriding keyword
#'   parsing (use when your labels are study-specific, e.g. drug names for the
#'   resistant arm).
#' @param score Score column name when `emt_cells` is a data.frame. Default
#'   "consensus".
#' @param min_cells Minimum finite cells required per model x condition group;
#'   smaller groups are dropped with a warning (a handful of cells gives an
#'   unstable per-group summary). Default 20.
#' @return A tibble: `cell`, `model`, `condition` (ordered sensitive < resistant),
#'   `emt`. Models lacking both conditions are kept but flagged in
#'   `attr(., "unpaired_models")`.
#' @export
prepare_resistance_emt <- function(emt_cells, cell_meta,
                                   model_col = "model",
                                   condition_col = "condition",
                                   cell_col = NULL,
                                   sensitive_labels = NULL,
                                   resistant_labels = NULL,
                                   score = "consensus", min_cells = 20) {
  emt_v <- .er_emt_vector(emt_cells, score)
  if (!is.data.frame(cell_meta)) cli::cli_abort("{.arg cell_meta} must be a data.frame.")
  # When cell_col is given, require it too so a typo'd/missing id column fails
  # fast here rather than as a confusing "No cell ids shared" downstream.
  need <- c(model_col, condition_col, if (!is.null(cell_col)) cell_col)
  miss <- setdiff(need, names(cell_meta))
  if (length(miss) > 0) cli::cli_abort("Missing {.arg cell_meta} column(s): {.val {miss}}.")

  cells <- if (is.null(cell_col)) rownames(cell_meta) else as.character(cell_meta[[cell_col]])
  if (is.null(cells)) {
    cli::cli_abort(c("x" = "No cell ids found in {.arg cell_meta}.",
                     "i" = "Set {.arg cell_col} or give the data.frame rownames matching {.arg emt_cells}."))
  }
  if (anyDuplicated(cells)) {
    cli::cli_abort("{.arg cell_meta} has duplicate cell ids: {.val {unique(cells[duplicated(cells)])}}.")
  }

  shared <- intersect(names(emt_v), cells)
  if (length(shared) == 0) {
    cli::cli_abort(c(
      "x" = "No cell ids shared between {.arg emt_cells} and {.arg cell_meta}.",
      "i" = "Check the id schemes match (e.g. {.val {names(emt_v)[1]}} vs {.val {cells[1]}})."
    ))
  }
  if (length(shared) < length(emt_v)) {
    cli::cli_warn("Matched {length(shared)}/{length(emt_v)} EMT-scored cells to metadata.")
  }
  idx <- match(shared, cells)
  cond <- .er_parse_condition(cell_meta[[condition_col]][idx],
                              sensitive = sensitive_labels, resistant = resistant_labels)

  df <- data.frame(
    cell = shared,
    model = as.character(cell_meta[[model_col]][idx]),
    condition = cond,
    emt = as.numeric(emt_v[shared]),
    row.names = NULL, stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$emt), , drop = FALSE]
  if (nrow(df) == 0) cli::cli_abort("No cells with finite EMT scores after the join.")

  # Drop cells with a missing model or condition BEFORE grouping: an NA in either
  # makes interaction() emit an NA group that table()/`%in%` then silently discard,
  # so without this the cell count would shrink with no explanation.
  bad_meta <- is.na(df$model) | is.na(df$condition)
  if (any(bad_meta)) {
    cli::cli_warn("Dropping {sum(bad_meta)} cell(s) with a missing model or condition.")
    df <- df[!bad_meta, , drop = FALSE]
  }
  if (nrow(df) == 0) cli::cli_abort("No cells left after dropping missing model/condition.")

  # Drop tiny model x condition groups (unstable per-group summaries).
  grp <- interaction(df$model, df$condition, drop = TRUE)
  keep_grp <- names(which(table(grp) >= min_cells))
  dropped <- sum(!(as.character(grp) %in% keep_grp))
  if (dropped > 0) {
    cli::cli_warn("Dropping {dropped} cell(s) in model x condition groups with < {min_cells} cells.")
  }
  df <- df[as.character(grp) %in% keep_grp, , drop = FALSE]
  if (nrow(df) == 0) {
    cli::cli_abort("No model x condition group reached {min_cells} cells; lower {.arg min_cells} or check inputs.")
  }

  # Flag models that lack one of the two states (cannot contribute to pairing).
  tab <- table(df$model, df$condition)
  unpaired <- rownames(tab)[tab[, "sensitive"] == 0 | tab[, "resistant"] == 0]
  if (length(unpaired) > 0) {
    cli::cli_warn(c(
      "!" = "{length(unpaired)} model(s) have only one treatment state and cannot be paired: {.val {unpaired}}.",
      "i" = "Paired tests use the remaining {length(setdiff(rownames(tab), unpaired))} model(s)."
    ))
  }
  out <- .er_as_tbl(df)
  attr(out, "unpaired_models") <- unpaired
  out
}

# ============================================================================
# 2. LOCATION: does the EMT axis shift up with resistance?
# ============================================================================

#' Does the EMT axis rise from chemo-sensitive to chemo-resistant?
#'
#' Summarizes each model x condition by a robust location of its per-cell EMT
#' (median by default -- a pseudobulk of the axis), then tests resistant vs
#' sensitive PAIRED across models. Positive deltas mean the resistant state is
#' more mesenchymal -- i.e. EMT arose with treatment.
#'
#' @param prepared Output of [prepare_resistance_emt()].
#' @param summary_fun Per-group location summary. Default `stats::median`.
#' @param alternative Paired-test alternative. Default "greater" (the directional
#'   hypothesis that resistance raises EMT); use "two.sided" to stay agnostic.
#' @param min_models Minimum paired models for the Wilcoxon test. Default 3.
#' @return A list: `per_model` (per-model sensitive/resistant location + delta)
#'   and `test` (one-row paired-test summary with sign counts).
#' @export
emt_resistance_shift <- function(prepared, summary_fun = stats::median,
                                 alternative = c("greater", "two.sided", "less"),
                                 min_models = 3) {
  alternative <- match.arg(alternative)
  .er_check_prepared(prepared)
  agg <- .er_aggregate(prepared, function(z) summary_fun(z[is.finite(z)]))
  .er_paired_test(agg, value_name = "emt", min_models = min_models,
                  alternative = alternative)
}

# ============================================================================
# 3. DISPERSION: does EMT heterogeneity widen with resistance?
# ============================================================================

#' Does per-cell EMT heterogeneity increase with resistance?
#'
#' This is the dispersion analogue of [emt_resistance_shift()], and the question
#' Stewart et al. 2020 actually emphasize: resistance coexists with WIDER
#' intratumoral spread of cell states. Summarizes each model x condition by a
#' dispersion of its per-cell EMT, then tests resistant vs sensitive paired
#' across models. Positive deltas mean the resistant state is more heterogeneous.
#'
#' @param prepared Output of [prepare_resistance_emt()].
#' @param dispersion Which spread to use: "sd" (default), "var", "iqr", or "mad"
#'   (median absolute deviation). MAD/IQR are robust to outlier cells.
#' @param alternative,min_models As in [emt_resistance_shift()].
#' @return A list: `per_model` and `test`, labelled by the chosen dispersion.
#' @export
emt_heterogeneity_shift <- function(prepared,
                                    dispersion = c("sd", "var", "iqr", "mad"),
                                    alternative = c("greater", "two.sided", "less"),
                                    min_models = 3) {
  dispersion <- match.arg(dispersion)
  alternative <- match.arg(alternative)
  .er_check_prepared(prepared)
  disp_fun <- .er_disp_fun(dispersion)
  agg <- .er_aggregate(prepared, function(z) {
    z <- z[is.finite(z)]
    if (length(z) < 2) return(NA_real_)
    disp_fun(z)
  })
  .er_paired_test(agg, value_name = dispersion, min_models = min_models,
                  alternative = alternative)
}

#' Dispersion function by name (sd / var / iqr / mad)
#' @keywords internal
.er_disp_fun <- function(dispersion) {
  switch(dispersion,
    sd  = function(z) stats::sd(z),
    var = function(z) stats::var(z),
    iqr = function(z) stats::IQR(z),
    mad = function(z) stats::mad(z)
  )
}

#' Per-model dispersion deltas (resistant - sensitive), paired models only
#' @keywords internal
.er_model_dispersion_deltas <- function(prepared, disp_fun) {
  ag <- stats::aggregate(emt ~ model + condition, data = prepared, FUN = function(z) {
    z <- z[is.finite(z)]; if (length(z) < 2) NA_real_ else disp_fun(z)
  })
  ag$condition <- as.character(ag$condition)
  s <- ag[ag$condition == "sensitive", c("model", "emt")]
  r <- ag[ag$condition == "resistant", c("model", "emt")]
  m <- merge(s, r, by = "model", suffixes = c("_s", "_r"))
  m <- m[is.finite(m$emt_s) & is.finite(m$emt_r), , drop = FALSE]
  stats::setNames(m$emt_r - m$emt_s, m$model)
}

# ============================================================================
# 3b. ROBUSTNESS: is the dispersion shift real, or a cell-count artifact?
# ============================================================================

#' Permutation null for the EMT dispersion shift
#'
#' Asks whether the observed widening of per-cell EMT under resistance exceeds
#' what random relabeling of the same cells produces. Within EACH model the
#' sensitive/resistant labels are shuffled (group sizes preserved, so this also
#' controls for the unequal cell counts between conditions), the per-model
#' dispersion deltas are recomputed, and two aggregate statistics are compared
#' to their null: the MEAN per-model delta and the NUMBER of models with a
#' positive delta (sign-consistency).
#'
#' @param prepared Output of [prepare_resistance_emt()].
#' @param dispersion "sd" (default), "var", "iqr", or "mad".
#' @param n_perm Number of permutations. Default 1000.
#' @param seed RNG seed (results are otherwise non-reproducible). Default 1.
#' @return A list: `observed` (mean_delta, n_pos, n_models, per_model deltas),
#'   `p_mean` and `p_pos` (one-sided empirical p-values, `(1+#ge)/(1+n_perm)`),
#'   `dispersion`, `n_perm`.
#' @export
emt_dispersion_permutation <- function(prepared,
                                       dispersion = c("sd", "var", "iqr", "mad"),
                                       n_perm = 1000, seed = 1) {
  dispersion <- match.arg(dispersion)
  .er_check_prepared(prepared)
  disp_fun <- .er_disp_fun(dispersion)
  obs <- .er_model_dispersion_deltas(prepared, disp_fun)
  if (length(obs) < 2) {
    cli::cli_abort("Need >= 2 paired models for a permutation test; found {length(obs)}.")
  }
  obs_mean <- mean(obs); obs_pos <- sum(obs > 0); n_models <- length(obs)
  models <- unique(prepared$model)
  idx_by_model <- lapply(models, function(m) which(prepared$model == m))
  cond <- prepared$condition
  set.seed(seed)
  perm_mean <- numeric(n_perm); perm_pos <- integer(n_perm)
  for (p in seq_len(n_perm)) {
    pp <- prepared
    cc <- cond
    for (ix in idx_by_model) cc[ix] <- sample(cc[ix])   # shuffle labels within model
    pp$condition <- cc
    d <- .er_model_dispersion_deltas(pp, disp_fun)
    perm_mean[p] <- mean(d); perm_pos[p] <- sum(d > 0)
  }
  list(
    observed = list(mean_delta = obs_mean, n_pos = obs_pos, n_models = n_models,
                    per_model = obs),
    p_mean = (1 + sum(perm_mean >= obs_mean)) / (1 + n_perm),
    p_pos  = (1 + sum(perm_pos  >= obs_pos )) / (1 + n_perm),
    dispersion = dispersion, n_perm = n_perm
  )
}

#' Equal-n downsampling stability for the EMT dispersion shift
#'
#' Recomputes the per-model dispersion deltas on repeated subsamples that take
#' the SAME number of cells from every model x condition group, removing the
#' cell-count imbalance as an explanation. Reports how often all models still go
#' up and the mean delta across repeats.
#'
#' @param prepared Output of [prepare_resistance_emt()].
#' @param dispersion "sd" (default), "var", "iqr", or "mad".
#' @param n_cells Cells to sample per group. Default: the smallest group size.
#' @param n_rep Number of subsample repeats. Default 200.
#' @param seed RNG seed. Default 1.
#' @return A list: `n_cells`, `n_models`, `mean_delta` (averaged over repeats),
#'   `frac_all_up` (fraction of repeats where every model's delta > 0),
#'   `mean_n_up` (average count of up models), `dispersion`, `n_rep`.
#' @export
emt_dispersion_downsample <- function(prepared,
                                      dispersion = c("sd", "var", "iqr", "mad"),
                                      n_cells = NULL, n_rep = 200, seed = 1) {
  dispersion <- match.arg(dispersion)
  .er_check_prepared(prepared)
  disp_fun <- .er_disp_fun(dispersion)
  tab <- table(prepared$model, prepared$condition)
  paired_models <- rownames(tab)[tab[, "sensitive"] > 0 & tab[, "resistant"] > 0]
  if (length(paired_models) < 2) {
    cli::cli_abort("Need >= 2 paired models for downsampling; found {length(paired_models)}.")
  }
  sub <- prepared[prepared$model %in% paired_models, , drop = FALSE]
  min_grp <- min(tab[paired_models, c("sensitive", "resistant")])
  if (is.null(n_cells)) n_cells <- min_grp
  if (n_cells > min_grp) {
    cli::cli_abort("{.arg n_cells} ({n_cells}) exceeds the smallest group ({min_grp}); sampling is without replacement.")
  }
  key <- paste(sub$model, sub$condition, sep = "\r")
  rows_by_grp <- split(seq_len(nrow(sub)), key)
  set.seed(seed)
  rep_mean <- numeric(n_rep); rep_up <- integer(n_rep)
  for (i in seq_len(n_rep)) {
    take <- unlist(lapply(rows_by_grp, function(r) sample(r, n_cells)), use.names = FALSE)
    d <- .er_model_dispersion_deltas(sub[take, , drop = FALSE], disp_fun)
    rep_mean[i] <- mean(d); rep_up[i] <- sum(d > 0)
  }
  list(
    n_cells = n_cells, n_models = length(paired_models),
    mean_delta = mean(rep_mean),
    frac_all_up = mean(rep_up == length(paired_models)),
    mean_n_up = mean(rep_up), dispersion = dispersion, n_rep = n_rep
  )
}

# ============================================================================
# 4. COMPOSITION: does the mesenchymal-cell fraction increase with resistance?
# ============================================================================

#' Does the fraction of mesenchymal cells increase with resistance?
#'
#' Calls each cell mesenchymal if its EMT score exceeds `threshold`, computes the
#' mesenchymal fraction per model x condition, then tests resistant vs sensitive
#' paired across models. By default `threshold` is the cohort-wide upper-tertile
#' cut of the per-cell EMT axis, so "mesenchymal" means "in the top third of the
#' EMT axis across all cells" -- a fixed reference shared by both states.
#'
#' @param prepared Output of [prepare_resistance_emt()].
#' @param threshold EMT cut above which a cell is mesenchymal. Default: the
#'   `2/3` quantile of `prepared$emt` across all cells.
#' @param alternative,min_models As in [emt_resistance_shift()].
#' @return A list: `per_model` (mesenchymal fraction per state + delta) and
#'   `test`. `attr(test, "threshold")` records the cut used.
#' @export
emt_state_composition <- function(prepared, threshold = NULL,
                                  alternative = c("greater", "two.sided", "less"),
                                  min_models = 3) {
  alternative <- match.arg(alternative)
  .er_check_prepared(prepared)
  if (is.null(threshold)) {
    threshold <- stats::quantile(prepared$emt[is.finite(prepared$emt)], 2/3, names = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold)) {
    cli::cli_abort("{.arg threshold} must be a single finite number.")
  }
  agg <- .er_aggregate(prepared, function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0) return(NA_real_)
    mean(z > threshold)
  })
  res <- .er_paired_test(agg, value_name = "mes_frac", min_models = min_models,
                         alternative = alternative)
  attr(res$test, "threshold") <- threshold
  res
}

#' Validate a prepared resistance frame
#' @keywords internal
.er_check_prepared <- function(prepared) {
  if (!is.data.frame(prepared) ||
      !all(c("cell", "model", "condition", "emt") %in% names(prepared))) {
    cli::cli_abort(c("x" = "{.arg prepared} must come from prepare_resistance_emt().",
                     "i" = "Expected columns: cell, model, condition, emt."))
  }
  # Check the actual values present, not factor levels: a factor carried over
  # from subsetting may keep unused levels that are not in the data.
  vals <- setdiff(unique(as.character(prepared$condition)), NA)
  extra <- setdiff(vals, c("sensitive", "resistant"))
  if (length(extra) > 0) {
    cli::cli_abort("{.arg prepared}$condition has unexpected value(s): {.val {extra}}. Expected sensitive / resistant.")
  }
  invisible(TRUE)
}
