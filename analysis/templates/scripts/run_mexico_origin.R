#!/usr/bin/env Rscript

# Mexico admission-origin analysis from grouped, identifier-free source counts.
# PROCED=2 is an ED-origin proxy and PROCED=1 is outpatient-consultation
# origin. Other and missing codes remain in the extraction but are not assigned
# to either model arm. No population offset is used.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE) else
  normalizePath("work/repair_20261112/scripts/run_mexico_origin.R", winslash = "/", mustWork = FALSE)
repair_root <- dirname(dirname(script_path))
repair_rlib <- file.path(repair_root, "Rlib")
if (dir.exists(repair_rlib)) .libPaths(c(repair_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

engine_path <- file.path(repair_root, "scripts", "repair_engine.R")
if (!file.exists(engine_path)) stop(sprintf("repair engine not found: %s", engine_path), call. = FALSE)
source(engine_path, local = FALSE)

canonical_state <- function(x) {
  z <- trimws(as.character(x))
  numeric <- grepl("^[0-9]+$", z)
  z[numeric] <- sprintf("%02d", as.integer(z[numeric]))
  z
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (nzchar(Sys.which("sha256sum"))) {
    output <- system2(Sys.which("sha256sum"), path, stdout = TRUE, stderr = FALSE)
    if (!length(output)) stop(sprintf("sha256sum failed for %s", path), call. = FALSE)
    return(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]])
  }
  if (nzchar(Sys.which("certutil"))) {
    output <- system2(Sys.which("certutil"), c("-hashfile", shQuote(path), "SHA256"), stdout = TRUE, stderr = FALSE)
    value <- gsub("[[:space:]]", "", output[grepl("^[[:xdigit:]]{64}$", output)])
    if (length(value)) return(tolower(value[[1L]]))
  }
  NA_character_
}

fread_gzip <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE) && .Platform$OS.type == "unix") {
    return(data.table::fread(cmd = paste("gzip -dc --", shQuote(path)), ...))
  }
  data.table::fread(path, ...)
}

read_model_inputs <- function(root, grouped_path = NULL, calendar_path = NULL) {
  if (is.null(grouped_path)) grouped_path <- file.path(root, "data", "mx_origin_grouped_2015_2023.csv.gz")
  if (is.null(calendar_path)) {
    updated_calendar <- file.path(root, "data", "mx_updated_calendar.csv.gz")
    calendar_path <- if (file.exists(updated_calendar)) updated_calendar else file.path(root, "data", "mx_calendar.csv")
  }
  if (!file.exists(grouped_path)) stop(sprintf("grouped origin input not found: %s", grouped_path), call. = FALSE)
  if (!file.exists(calendar_path)) stop(sprintf("Mexico calendar not found: %s", calendar_path), call. = FALSE)
  clinical <- fread_gzip(
    grouped_path,
    select = c("state", "date", "dx3", "origin_raw", "origin_group", "los_group", "type_service", "n"),
    colClasses = list(character = c("state", "date", "dx3", "origin_raw", "origin_group", "los_group", "type_service")),
    showProgress = TRUE
  )
  calendar <- fread_gzip(
    calendar_path,
    select = c("date", "unit", "E3"),
    colClasses = list(character = c("date", "unit")),
    showProgress = TRUE
  )
  if (!all(c("date", "unit", "E3") %in% names(calendar))) stop("mx_calendar.csv must include date, unit, and precomputed E3", call. = FALSE)
  clinical[, state := canonical_state(state)]
  clinical[, date := as.Date(date)]
  clinical[, dx3 := toupper(trimws(as.character(dx3)))]
  clinical[, origin_group := trimws(as.character(origin_group))]
  clinical[, los_group := trimws(as.character(los_group))]
  clinical[, n := suppressWarnings(as.numeric(n))]
  if (anyNA(clinical$date) || any(!is.finite(clinical$n) | clinical$n < 0)) stop("origin grouped input has invalid dates or counts", call. = FALSE)
  calendar[, state := canonical_state(unit)]
  calendar[, date := as.Date(date)]
  calendar[, E3 := suppressWarnings(as.numeric(E3))]
  if (anyNA(calendar$date)) stop("Mexico calendar contains invalid dates", call. = FALSE)
  duplicate <- calendar[, .N, by = .(state, date)][N > 1L]
  if (nrow(duplicate)) stop(sprintf("calendar contains %d duplicate state-day rows", nrow(duplicate)), call. = FALSE)
  list(clinical = clinical, calendar = calendar, grouped_path = grouped_path, calendar_path = calendar_path)
}

daily_outcomes <- function(clinical, calendar, event_window = c("2015-01-01", "2022-12-31"), los_filter = NULL) {
  date_lo <- as.Date(event_window[[1L]])
  date_hi <- as.Date(event_window[[2L]])
  d <- clinical[date >= date_lo & date <= date_hi]
  if (!is.null(los_filter)) d <- d[los_group == los_filter]
  relevant_origins <- c("ED_origin", "outpatient_consultation_origin")
  d <- d[origin_group %in% relevant_origins & dx3 %in% c("K80", "N20", "N21", "N22", "N23")]
  daily <- d[, .(
    y_K80_ED = sum(n[dx3 == "K80" & origin_group == "ED_origin"]),
    y_N20_N23_ED = sum(n[dx3 %in% c("N20", "N21", "N22", "N23") & origin_group == "ED_origin"]),
    y_N20_only_ED = sum(n[dx3 == "N20" & origin_group == "ED_origin"]),
    y_K80_OP = sum(n[dx3 == "K80" & origin_group == "outpatient_consultation_origin"]),
    y_N20_N23_OP = sum(n[dx3 %in% c("N20", "N21", "N22", "N23") & origin_group == "outpatient_consultation_origin"])
  ), by = .(state, date)]
  cal <- calendar[date >= date_lo & date <= date_hi, .(state, date, E3)]
  panel <- merge(cal, daily, by = c("state", "date"), all.x = TRUE, sort = FALSE)
  model_cols <- c("y_K80_ED", "y_N20_N23_ED", "y_N20_only_ED", "y_K80_OP", "y_N20_N23_OP")
  for (nm in model_cols) {
    set(panel, which(is.na(panel[[nm]])), nm, 0)
    set(panel, j = nm, value = as.numeric(panel[[nm]]))
  }
  setorder(panel, state, date)
  panel
}

unit_target <- function(d) {
  eligible <- d$E3[is.finite(d$E3)]
  if (length(eligible) < 400L || length(unique(eligible)) < 2L) stop("fewer than 400 E3-eligible days or fewer than two E3 values", call. = FALSE)
  list(
    basis = new_temperature_basis(eligible, df = 3L),
    p95 = unname(stats::quantile(eligible, 0.95, type = 7)),
    median = stats::median(eligible),
    exposure_days = length(eligible)
  )
}

linear_effect <- function(L, beta, V) {
  L <- setNames(as.numeric(L), names(beta))
  estimate <- as.numeric(L %*% beta)
  variance <- as.numeric(t(L) %*% V %*% L)
  if (!is.finite(variance) || variance <= 0) stop("primary HAC contrast variance is non-positive or non-finite", call. = FALSE)
  se <- sqrt(variance)
  list(
    estimate = estimate,
    variance = variance,
    se = se,
    rr = exp(estimate),
    lo = exp(estimate - stats::qnorm(0.975) * se),
    hi = exp(estimate + stats::qnorm(0.975) * se),
    L = L
  )
}

basis_endpoint_vector <- function(fit_result, endpoint, basis_name) {
  fit <- fit_result$fit
  template <- as.data.frame(fit_result$data[1L, , drop = FALSE])
  template$disease <- factor(endpoint, levels = levels(fit_result$data$disease))
  template$B1 <- 0; template$B2 <- 0; template$B3 <- 0
  tt <- stats::delete.response(stats::terms(fit))
  x0 <- stats::model.matrix(tt, data = template, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  template[[basis_name]] <- 1
  x1 <- stats::model.matrix(tt, data = template, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  beta_names <- names(stats::coef(fit))
  as.numeric(x1[, beta_names, drop = FALSE] - x0[, beta_names, drop = FALSE])
}

formal_basis_test <- function(fit_result, endpoint_names) {
  # Rows are the three E3 spline-basis components of
  # (N20-N23 vs K80)_ED minus (N20-N23 vs K80)_outpatient.
  pairs <- stats::setNames(vector("list", 3L), c("B1", "B2", "B3"))
  for (basis_name in c("B1", "B2", "B3")) {
    vec <- function(endpoint) basis_endpoint_vector(fit_result, endpoint, basis_name)
    pairs[[basis_name]] <- (vec("y_N20_N23_ED") - vec("y_K80_ED")) -
      (vec("y_N20_N23_OP") - vec("y_K80_OP"))
  }
  Lmat <- do.call(rbind, pairs)
  beta <- stats::coef(fit_result$fit)
  V <- fit_result$hac_vcov$hac14
  colnames(Lmat) <- names(beta)
  basis_estimate <- as.numeric(Lmat %*% beta)
  covariance <- Lmat %*% V %*% t(Lmat)
  chi <- as.numeric(tryCatch(t(basis_estimate) %*% solve(covariance, basis_estimate), error = function(e) NA_real_))
  list(
    estimate = setNames(basis_estimate, rownames(Lmat)),
    covariance = covariance,
    chi_square_3df = chi,
    p_value = if (is.finite(chi)) stats::pchisq(chi, df = 3L, lower.tail = FALSE) else NA_real_,
    L = Lmat
  )
}

summarize_primary_fit <- function(d, unit, target) {
  outcomes <- c("y_K80_ED", "y_N20_N23_ED", "y_K80_OP", "y_N20_N23_OP")
  eligible <- is.finite(d$E3)
  events <- vapply(outcomes, function(nm) sum(d[[nm]][eligible]), numeric(1))
  if (target$exposure_days < 400L || any(events < 30L)) {
    return(list(status = "ineligible_threshold", unit = unit, n_days = target$exposure_days,
                events = events, failure = "Requires at least 30 events in each of four cells and at least 400 E3-eligible days."))
  }
  fit <- fit_stacked_contrast(
    data = d, outcomes = outcomes, exposure = "E3", date = "date", offset = NULL,
    basis = target$basis, p95 = target$p95, median = target$median,
    nuisance = "separate", hac_lags = 14L, full_covariance = TRUE
  )
  beta <- stats::coef(fit$fit)
  V <- fit$hac_vcov$hac14
  names_effects <- lapply(fit$exposure_effects, function(value) setNames(as.numeric(value), names(beta)))
  ed_L <- names_effects$y_N20_N23_ED - names_effects$y_K80_ED
  op_L <- names_effects$y_N20_N23_OP - names_effects$y_K80_OP
  primary_L <- ed_L - op_L
  primary <- linear_effect(primary_L, beta, V)
  broad_ed <- linear_effect(ed_L, beta, V)
  broad_op <- linear_effect(op_L, beta, V)
  pair_rows <- list()
  pairs <- utils::combn(outcomes, 2L, simplify = FALSE)
  for (pair in pairs) {
    denominator <- pair[[1L]]
    numerator <- pair[[2L]]
    L <- names_effects[[numerator]] - names_effects[[denominator]]
    label <- paste0(numerator, "_vs_", denominator)
    effect <- linear_effect(L, beta, V)
    engine_value <- fit$contrasts[[label]]$coefficient_contrast$estimate
    if (!isTRUE(all.equal(unname(effect$estimate), unname(engine_value), tolerance = 1e-9))) {
      stop(sprintf("full-L and engine pairwise contrast disagree for %s", label), call. = FALSE)
    }
    pair_rows[[label]] <- list(
      numerator = numerator, denominator = denominator,
      logRRR = effect$estimate, variance_hac14 = effect$variance, se_hac14 = effect$se,
      RRR = effect$rr, lo = effect$lo, hi = effect$hi,
      coefficient_contrast_L = as.list(effect$L)
    )
  }
  formal <- formal_basis_test(fit, names_effects)
  list(
    status = "ok", unit = unit, n_days = target$exposure_days,
    events = events, p95 = target$p95, median = target$median,
    primary = primary, broad_ed = broad_ed, broad_op = broad_op,
    formal = formal, fit = fit, beta = beta, V = V,
    exposure_effects = names_effects, all_pairs = pair_rows,
    audit = fit$audit
  )
}

summarize_strict_ed_fit <- function(d, unit, target) {
  outcomes <- c("y_K80_ED", "y_N20_only_ED")
  eligible <- is.finite(d$E3)
  events <- vapply(outcomes, function(nm) sum(d[[nm]][eligible]), numeric(1))
  if (target$exposure_days < 400L || any(events < 30L)) {
    return(list(status = "ineligible_threshold", unit = unit, n_days = target$exposure_days,
                events = events, failure = "Strict ED N20-only fit requires at least 30 events in each cell and 400 E3-eligible days."))
  }
  fit <- fit_stacked_contrast(
    data = d, outcomes = outcomes, exposure = "E3", date = "date", offset = NULL,
    basis = target$basis, p95 = target$p95, median = target$median,
    nuisance = "separate", hac_lags = 14L, full_covariance = TRUE
  )
  beta <- stats::coef(fit$fit)
  V <- fit$hac_vcov$hac14
  effects <- lapply(fit$exposure_effects, function(value) setNames(as.numeric(value), names(beta)))
  L <- effects$y_N20_only_ED - effects$y_K80_ED
  effect <- linear_effect(L, beta, V)
  list(status = "ok", unit = unit, n_days = target$exposure_days, events = events,
       p95 = target$p95, median = target$median, effect = effect,
       fit = fit, beta = beta, V = V, exposure_effects = effects, audit = fit$audit)
}

unit_model_worker <- function(unit) {
  d <- panels_by_state[[unit]]
  target <- tryCatch(unit_target(d), error = function(e) e)
  if (inherits(target, "error")) {
    skipped <- list(status = "failed", unit = unit, n_days = sum(is.finite(d$E3)), failure = conditionMessage(target))
    return(list(unit = unit, primary = skipped, strict = skipped, shortstay = skipped))
  }
  primary <- tryCatch(summarize_primary_fit(d, unit, target), error = function(e) list(
    status = "failed", unit = unit, n_days = target$exposure_days, failure = conditionMessage(e)
  ))
  strict <- tryCatch(summarize_strict_ed_fit(d, unit, target), error = function(e) list(
    status = "failed", unit = unit, n_days = target$exposure_days, failure = conditionMessage(e)
  ))
  short_panel <- shortstay_panels[[unit]]
  shortstay <- tryCatch(summarize_primary_fit(short_panel, unit, unit_target(short_panel)), error = function(e) list(
    status = "failed", unit = unit, n_days = sum(is.finite(short_panel$E3)), failure = conditionMessage(e)
  ))
  list(unit = unit, primary = primary, strict = strict, shortstay = shortstay)
}

result_to_row <- function(result, model_name) {
  base <- data.table(
    model = model_name,
    state = result$unit,
    status = result$status,
    n_eligible_days = result$n_days %||% NA_real_,
    failure = result$failure %||% NA_character_
  )
  if (!identical(result$status, "ok")) return(base)
  if (model_name %in% c("primary_broad_four_endpoint", "sensitivity_recorded_los_1to30")) {
    base[, `:=`(
      events_K80_ED = unname(result$events[["y_K80_ED"]]),
      events_N20_N23_ED = unname(result$events[["y_N20_N23_ED"]]),
      events_K80_outpatient = unname(result$events[["y_K80_OP"]]),
      events_N20_N23_outpatient = unname(result$events[["y_N20_N23_OP"]]),
      p95_E3 = result$p95, median_E3 = result$median,
      logRRR_N20_N23_vs_K80_ED = result$broad_ed$estimate,
      var_N20_N23_vs_K80_ED_hac14 = result$broad_ed$variance,
      RRR_N20_N23_vs_K80_ED = result$broad_ed$rr,
      lo_N20_N23_vs_K80_ED = result$broad_ed$lo,
      hi_N20_N23_vs_K80_ED = result$broad_ed$hi,
      logRRR_N20_N23_vs_K80_outpatient = result$broad_op$estimate,
      var_N20_N23_vs_K80_outpatient_hac14 = result$broad_op$variance,
      RRR_N20_N23_vs_K80_outpatient = result$broad_op$rr,
      lo_N20_N23_vs_K80_outpatient = result$broad_op$lo,
      hi_N20_N23_vs_K80_outpatient = result$broad_op$hi,
      logRRR_ED_minus_outpatient = result$primary$estimate,
      var_ED_minus_outpatient_hac14 = result$primary$variance,
      se_ED_minus_outpatient_hac14 = result$primary$se,
      RRR_ratio_ED_to_outpatient = result$primary$rr,
      lo_RRR_ratio_ED_to_outpatient = result$primary$lo,
      hi_RRR_ratio_ED_to_outpatient = result$primary$hi,
      joint_wald_chisq_3df = result$formal$chi_square_3df,
      joint_wald_p_3df = result$formal$p_value,
      joint_basis_estimates_json = as.character(toJSON(result$formal$estimate, auto_unbox = TRUE, digits = NA)),
      joint_basis_covariance_json = as.character(toJSON(unclass(result$formal$covariance), auto_unbox = TRUE, digits = NA)),
      all_pair_effects_json = as.character(toJSON(result$all_pairs, auto_unbox = TRUE, digits = NA)),
      exposure_effects_json = as.character(toJSON(result$exposure_effects, auto_unbox = TRUE, digits = NA)),
      formula = result$audit$formula,
      convergence_warnings_json = as.character(toJSON(result$audit$warnings, auto_unbox = TRUE, digits = NA))
    )]
  } else {
    base[, `:=`(
      events_K80_ED = unname(result$events[["y_K80_ED"]]),
      events_N20_only_ED = unname(result$events[["y_N20_only_ED"]]),
      p95_E3 = result$p95, median_E3 = result$median,
      logRRR_N20_only_vs_K80_ED = result$effect$estimate,
      var_N20_only_vs_K80_ED_hac14 = result$effect$variance,
      se_N20_only_vs_K80_ED_hac14 = result$effect$se,
      RRR_N20_only_vs_K80_ED = result$effect$rr,
      lo_N20_only_vs_K80_ED = result$effect$lo,
      hi_N20_only_vs_K80_ED = result$effect$hi,
      formula = result$audit$formula,
      convergence_warnings_json = as.character(toJSON(result$audit$warnings, auto_unbox = TRUE, digits = NA))
    )]
  }
  base
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

write_covariance_audit <- function(result, model_name, outdir) {
  if (!identical(result$status, "ok")) return(NA_character_)
  path <- file.path(outdir, "covariance", paste0(model_name, "_state_", result$unit, "_hac14.json"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  payload <- list(
    model = model_name,
    state = result$unit,
    covariance = "Bartlett calendar HAC, lag 14",
    coefficient_names = names(result$beta),
    covariance_hac14 = unclass(result$V),
    exposure_effects_L = result$exposure_effects,
    all_pair_effects = result$all_pairs %||% NULL,
    primary_difference_L = result$primary$L %||% NULL,
    primary_variance = result$primary$variance %||% NULL,
    joint_basis_L = result$formal$L %||% NULL,
    joint_basis_estimate = result$formal$estimate %||% NULL,
    joint_basis_covariance = unclass(result$formal$covariance %||% matrix(numeric(), 0L, 0L)),
    joint_wald_chisq_3df = result$formal$chi_square_3df %||% NULL
  )
  write_json(payload, path, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA)
  path
}

run_synthetic_test <- function(outdir = NULL) {
  set.seed(118)
  n <- 1100L
  d <- data.table(
    date = seq(as.Date("2015-01-01"), by = "day", length.out = n),
    E3 = stats::rnorm(n, 20, 4)
  )
  d[, `:=`(
    y_K80_ED = stats::rpois(.N, exp(1.2 + 0.012 * E3)),
    y_N20_N23_ED = stats::rpois(.N, exp(1.1 + 0.018 * E3)),
    y_K80_OP = stats::rpois(.N, exp(1.3 + 0.010 * E3)),
    y_N20_N23_OP = stats::rpois(.N, exp(1.15 + 0.014 * E3))
  )]
  target <- unit_target(d)
  fit <- fit_stacked_contrast(
    d, c("y_K80_ED", "y_N20_N23_ED", "y_K80_OP", "y_N20_N23_OP"),
    exposure = "E3", date = "date", offset = NULL,
    basis = target$basis, p95 = target$p95, median = target$median,
    nuisance = "separate", hac_lags = 14L, full_covariance = TRUE
  )
  beta <- stats::coef(fit$fit); V <- fit$hac_vcov$hac14
  exposure_effects <- lapply(fit$exposure_effects, function(value) setNames(as.numeric(value), names(beta)))
  L <- (exposure_effects$y_N20_N23_ED - exposure_effects$y_K80_ED) -
    (exposure_effects$y_N20_N23_OP - exposure_effects$y_K80_OP)
  direct <- as.numeric(L %*% beta)
  endpoint_direct <- as.numeric((exposure_effects$y_N20_N23_ED - exposure_effects$y_K80_ED) %*% beta) -
    as.numeric((exposure_effects$y_N20_N23_OP - exposure_effects$y_K80_OP) %*% beta)
  stopifnot(isTRUE(all.equal(direct, endpoint_direct, tolerance = 1e-10)))
  projected <- as.numeric(fit$score$daily_scores %*% (fit$score$Hinv %*% L))
  score_meat <- .bartlett_hac_meat(matrix(projected, ncol = 1L), fit$score$daily_dates, lag = 14L)
  projected_variance <- as.numeric(score_meat[1L, 1L])
  full_l_variance <- as.numeric(t(L) %*% V %*% L)
  stopifnot(isTRUE(all.equal(full_l_variance, projected_variance, tolerance = 1e-9)))
  interaction <- formal_basis_test(fit, exposure_effects)
  direct_basis <- do.call(rbind, lapply(c("B1", "B2", "B3"), function(bn) {
    v <- function(ep) basis_endpoint_vector(fit, ep, bn)
    (v("y_N20_N23_ED") - v("y_K80_ED")) - (v("y_N20_N23_OP") - v("y_K80_OP"))
  }))
  stopifnot(isTRUE(all.equal(unname(interaction$L), unname(direct_basis), tolerance = 1e-10)))
  stopifnot(is.finite(interaction$chi_square_3df), is.finite(interaction$p_value))
  result <- list(
    status = "passed",
    fixture = "synthetic four-endpoint daily count panel",
    n_days = n,
    checks = c(
      "primary ED-minus-outpatient full-L estimate equals difference of within-origin exposure effects",
      "primary HAC14 variance evaluated from the full contrast vector and joint covariance",
      "3-df E3 spline interaction basis matches independently assembled endpoint-basis vectors",
      "all six engine pairwise contrasts are returned"
    ),
    primary_logRRR_difference = direct,
    primary_variance_hac14 = full_l_variance,
    primary_variance_from_projected_daily_scores = projected_variance,
    joint_wald_chisq_3df = interaction$chi_square_3df,
    joint_wald_p_3df = interaction$p_value,
    pairwise_contrast_labels = names(fit$contrasts)
  )
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write_json(result, file.path(outdir, "mx_origin_synthetic_fixture_test.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  }
  result
}

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(root = repair_root, grouped = NULL, calendar = NULL,
              outdir = file.path(repair_root, "results", "mx_origin"),
              workers = 8L, max_units = Inf, fixture_only = FALSE)
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (token %in% c("--fixture-only", "--synthetic-only")) {
      out$fixture_only <- TRUE; i <- i + 1L; next
    }
    if (!startsWith(token, "--") || i == length(args)) stop(sprintf("unexpected argument: %s", token), call. = FALSE)
    key <- sub("^--", "", token); value <- args[[i + 1L]]
    if (key == "root") out$root <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "grouped") out$grouped <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "calendar") out$calendar <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "outdir") out$outdir <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "workers") out$workers <- as.integer(value)
    else if (key == "max-units") out$max_units <- as.integer(value)
    else stop(sprintf("unknown option: %s", token), call. = FALSE)
    i <- i + 2L
  }
  if (out$workers < 1L || out$workers > 8L) stop("--workers must be 1-8", call. = FALSE)
  if (!is.infinite(out$max_units) && out$max_units < 1L) stop("--max-units must be positive", call. = FALSE)
  out
}

pool_effect <- function(per_state, estimate_col, variance_col, analysis_name) {
  attempted <- per_state[model == analysis_name]
  if (anyDuplicated(attempted$state)) {
    stop(sprintf("multiple rows per state reached pooled contrast '%s'", analysis_name), call. = FALSE)
  }
  fit <- attempted[status == "ok" & is.finite(get(estimate_col)) & is.finite(get(variance_col)) & get(variance_col) > 0]
  failed <- attempted[status != "ok" | !is.finite(get(estimate_col)) | !is.finite(get(variance_col)) | get(variance_col) <= 0]
  rows <- list()
  for (method in c("DL", "REML_HK")) {
    pooled <- pool_meta(fit[[estimate_col]], fit[[variance_col]], method = method)
    rows[[length(rows) + 1L]] <- data.table(
      analysis = analysis_name,
      pooling = method,
      status = if (is.null(pooled$failure)) "ok" else "partial_or_failed",
      units_attempted = nrow(attempted),
      units_fit = nrow(fit),
      units_failed_or_ineligible = nrow(failed),
      state_coverage_fraction = if (nrow(attempted)) nrow(fit) / nrow(attempted) else NA_real_,
      state_codes_fit_json = as.character(toJSON(fit$state, auto_unbox = TRUE)),
      state_codes_failed_or_ineligible_json = as.character(toJSON(failed$state, auto_unbox = TRUE)),
      k = pooled$k,
      logRRR = pooled$estimate,
      se = pooled$se,
      RRR = pooled$RR,
      lo = pooled$RR95CI[[1L]],
      hi = pooled$RR95CI[[2L]],
      I2 = pooled$I2,
      tau2 = pooled$tau2,
      concordance_positive = pooled$concordance,
      failure = pooled$failure %||% NA_character_
    )
  }
  rbindlist(rows, fill = TRUE)
}

main <- function() {
  args <- parse_cli()
  if (args$fixture_only) {
    result <- run_synthetic_test(args$outdir)
    cat("synthetic four-endpoint fixture:", result$status, "\n")
    return(invisible(result))
  }
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
  inputs <- read_model_inputs(args$root, args$grouped, args$calendar)
  primary_panel <- daily_outcomes(inputs$clinical, inputs$calendar)
  short_panel <- daily_outcomes(inputs$clinical, inputs$calendar, los_filter = "1to30")
  setnames(primary_panel, "state", "unit")
  setnames(short_panel, "state", "unit")
  panels_by_state <<- split(primary_panel, by = "unit", keep.by = TRUE)
  shortstay_panels <<- split(short_panel, by = "unit", keep.by = TRUE)
  state_codes <- sort(names(panels_by_state))
  expected_codes <- sprintf("%02d", 1:32)
  if (!identical(state_codes, expected_codes)) {
    stop(sprintf("Mexico calendar state coverage differs from 01-32: %s", paste(state_codes, collapse = ",")), call. = FALSE)
  }
  if (!is.infinite(args$max_units)) state_codes <- head(state_codes, as.integer(args$max_units))
  cat(sprintf("Running %d Mexico resident states with %d worker(s)\n", length(state_codes), min(args$workers, length(state_codes))))

  results <- if (min(args$workers, length(state_codes)) == 1L) {
    lapply(state_codes, unit_model_worker)
  } else {
    cl <- parallel::makeCluster(min(args$workers, length(state_codes)), outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, "engine_path", envir = .GlobalEnv)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({library(data.table); library(jsonlite)})
      source(engine_path, local = FALSE)
      NULL
    })
    parallel::clusterExport(cl, c(
      "panels_by_state", "shortstay_panels", "unit_model_worker", "unit_target",
      "summarize_primary_fit", "summarize_strict_ed_fit", "formal_basis_test",
      "basis_endpoint_vector", "linear_effect"
    ), envir = .GlobalEnv)
    parallel::parLapply(cl, state_codes, unit_model_worker)
  }

  primary_rows <- rbindlist(lapply(results, function(z) result_to_row(z$primary, "primary_broad_four_endpoint")), fill = TRUE)
  strict_rows <- rbindlist(lapply(results, function(z) result_to_row(z$strict, "strict_ED_N20_only_vs_K80")), fill = TRUE)
  short_rows <- rbindlist(lapply(results, function(z) result_to_row(z$shortstay, "sensitivity_recorded_los_1to30")), fill = TRUE)
  if (anyDuplicated(primary_rows$state) || anyDuplicated(strict_rows$state) || anyDuplicated(short_rows$state)) {
    stop("unit-level model export contains duplicate state rows", call. = FALSE)
  }
  if (nrow(primary_rows) != length(state_codes) || nrow(strict_rows) != length(state_codes) || nrow(short_rows) != length(state_codes)) {
    stop("unit-level model export row count differs from attempted state count", call. = FALSE)
  }
  support_rows <- rbindlist(list(
    primary_rows[status == "ok", .(state, analysis = "broad_N20_N23_vs_K80_within_ED_origin", status,
                                      n_eligible_days, events_numerator = events_N20_N23_ED, events_denominator = events_K80_ED,
                                      logRRR = logRRR_N20_N23_vs_K80_ED, variance_hac14 = var_N20_N23_vs_K80_ED_hac14,
                                      RRR = RRR_N20_N23_vs_K80_ED, lo = lo_N20_N23_vs_K80_ED, hi = hi_N20_N23_vs_K80_ED)],
    primary_rows[status == "ok", .(state, analysis = "broad_N20_N23_vs_K80_within_outpatient_origin", status,
                                      n_eligible_days, events_numerator = events_N20_N23_outpatient, events_denominator = events_K80_outpatient,
                                      logRRR = logRRR_N20_N23_vs_K80_outpatient, variance_hac14 = var_N20_N23_vs_K80_outpatient_hac14,
                                      RRR = RRR_N20_N23_vs_K80_outpatient, lo = lo_N20_N23_vs_K80_outpatient, hi = hi_N20_N23_vs_K80_outpatient)],
    strict_rows[, .(state, analysis = "N20_only_vs_K80_within_ED_origin", status, n_eligible_days,
                      events_numerator = events_N20_only_ED, events_denominator = events_K80_ED,
                      logRRR = logRRR_N20_only_vs_K80_ED, variance_hac14 = var_N20_only_vs_K80_ED_hac14,
                      RRR = RRR_N20_only_vs_K80_ED, lo = lo_N20_only_vs_K80_ED, hi = hi_N20_only_vs_K80_ED,
                      failure)]
  ), fill = TRUE)
  support_duplicates <- support_rows[, .N, by = .(analysis, state)][N != 1L]
  if (nrow(support_duplicates)) stop("supportive contrast export contains duplicate state rows", call. = FALSE)
  primary_rows[, `:=`(model = "primary_broad_four_endpoint", effect = "ED_minus_outpatient_N20_N23_vs_K80")]
  short_rows[, `:=`(model = "sensitivity_recorded_los_1to30", effect = "ED_minus_outpatient_N20_N23_vs_K80")]

  covariance_paths <- list()
  primary_covariance_paths <- vapply(results, function(z) write_covariance_audit(z$primary, "primary_broad_four_endpoint", args$outdir), character(1))
  strict_covariance_paths <- vapply(results, function(z) write_covariance_audit(z$strict, "strict_ED_N20_only_vs_K80", args$outdir), character(1))
  short_covariance_paths <- vapply(results, function(z) write_covariance_audit(z$shortstay, "sensitivity_recorded_los_1to30", args$outdir), character(1))
  covariance_paths$primary <- setNames(primary_covariance_paths, state_codes)
  covariance_paths$strict_ed <- setNames(strict_covariance_paths, state_codes)
  covariance_paths$los1to30 <- setNames(short_covariance_paths, state_codes)
  primary_rows[, covariance_hac14_file := primary_covariance_paths]
  strict_rows[, covariance_hac14_file := strict_covariance_paths]
  short_rows[, covariance_hac14_file := short_covariance_paths]

  fwrite(primary_rows, file.path(args$outdir, "mx_origin_primary_per_state.csv"), na = "")
  fwrite(support_rows, file.path(args$outdir, "mx_origin_supportive_per_state.csv"), na = "")
  fwrite(strict_rows, file.path(args$outdir, "mx_origin_strict_ed_per_state.csv"), na = "")
  fwrite(short_rows, file.path(args$outdir, "mx_origin_los1to30_sensitivity_per_state.csv"), na = "")
  failure_rows <- rbindlist(list(
    primary_rows[status != "ok", .(model, state, status, n_eligible_days, failure)],
    strict_rows[status != "ok", .(model, state, status, n_eligible_days, failure)],
    short_rows[status != "ok", .(model, state, status, n_eligible_days, failure)]
  ), fill = TRUE)
  fwrite(failure_rows, file.path(args$outdir, "mx_origin_failures.csv"), na = "")

  main_meta <- primary_rows[, .(
    model = "primary_ED_minus_outpatient",
    state = state,
    status = status,
    estimate = logRRR_ED_minus_outpatient,
    variance = var_ED_minus_outpatient_hac14
  )]
  sensitivity_meta <- short_rows[, .(
    model = "sensitivity_recorded_los_1to30",
    state = state,
    status = status,
    estimate = logRRR_ED_minus_outpatient,
    variance = var_ED_minus_outpatient_hac14
  )]
  supportive_meta <- support_rows[, .(model = analysis, state, status, estimate = logRRR, variance = variance_hac14)]
  meta_input <- rbindlist(list(main_meta, sensitivity_meta, supportive_meta), fill = TRUE)
  meta_duplicates <- meta_input[, .N, by = .(model, state)][N != 1L]
  if (nrow(meta_duplicates)) stop("duplicate state contrast rows detected before meta-analysis", call. = FALSE)
  meta_names <- unique(meta_input$model)
  pooled <- rbindlist(lapply(meta_names, function(nm) pool_effect(meta_input, "estimate", "variance", nm)), fill = TRUE)
  fwrite(pooled, file.path(args$outdir, "mx_origin_pooled.csv"), na = "")

  support_counts <- inputs$clinical[, .(n = sum(n)), by = .(origin_group, origin_raw, type_service)]
  support_counts[, origin_raw := as.character(origin_raw)]
  audit <- list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    specification_path = file.path(args$root, "admin", "mx_origin_specification.json"),
    specification_sha256 = sha256_file(file.path(args$root, "admin", "mx_origin_specification.json")),
    input_grouped_path = inputs$grouped_path,
    input_grouped_sha256 = sha256_file(inputs$grouped_path),
    input_calendar_path = inputs$calendar_path,
    input_calendar_sha256 = sha256_file(inputs$calendar_path),
    input_grouped_rows = nrow(inputs$clinical),
    source_event_year_dx3_totals_2015_2022 = inputs$clinical[date >= as.Date("2015-01-01") & date <= as.Date("2022-12-31"), .(n = sum(n)), by = .(date_year = format(date, "%Y"), dx3)],
    origin_and_service_code_coverage = support_counts,
    model_years = 2015:2022,
    state_codes_attempted = state_codes,
    unit_thresholds = list(primary = ">=30 events in each 4 cells and >=400 E3-eligible days", strict_N20_only_ED = ">=30 in each of 2 cells and >=400 E3-eligible days", los_sensitivity = "same as primary"),
    exposure = list(field = "precomputed E3 from mx_calendar.csv", basis_df = 3L, reference = "state-specific P95 versus median from all-origin eligible calendar days in 2015-2022; same values shared across all outcomes and sensitivity"),
    model = list(family = "quasi-Poisson", nuisance = "separate endpoint-specific day-of-week, month, and 12-df natural spline of time", offset = "none", covariance = "calendar-date Bartlett HAC lag 14", pooling = c("DerSimonian-Laird", "REML with Hartung-Knapp")),
    interpretation_limits = c("PROCED=2 is an ED-origin proxy, not proof of disease onset, acuity, triage severity, or mechanism.", "PROCED 1 and all non-ED codes are not elective.", "TIPSERV is not an urgency classifier.", "LOS sensitivity denotes recorded DIAS_ESTA 1-30 days, not hospital certification."),
    origin_data_totals_by_group = support_counts,
    units_fit = list(
      primary = primary_rows[status == "ok", state],
      primary_failed_or_ineligible = primary_rows[status != "ok", .(state, status, failure)],
      strict_ed = strict_rows[status == "ok", state],
      strict_ed_failed_or_ineligible = strict_rows[status != "ok", .(state, status, failure)],
      los1to30 = short_rows[status == "ok", state],
      los1to30_failed_or_ineligible = short_rows[status != "ok", .(state, status, failure)]
    ),
    pooled_summary = pooled,
    covariance_hac14_files = covariance_paths,
    patient_level_rows_used_or_emitted = FALSE
  )
  write_json(audit, file.path(args$outdir, "mx_origin_audit.json"), auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA)
  cat("Mexico origin analysis finished; results:", args$outdir, "\n")
  cat("Primary states fit:", sum(primary_rows$status == "ok"), "/", nrow(primary_rows), "\n")
  cat("Primary model failures/ineligible:", sum(primary_rows$status != "ok"), "\n")
  cat("Primary cohort: 2015-2022; common full-period state-specific E3 references; no offset; HAC14\n")
}

if (sys.nframe() == 0L) main()
