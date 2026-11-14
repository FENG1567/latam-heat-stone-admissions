#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(splines)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1L]]),
  winslash = "/", mustWork = FALSE) else normalizePath("run_temporal_m04.R", winslash = "/", mustWork = FALSE)
temporal_root_default <- dirname(dirname(dirname(script_path)))
base_root_default <- dirname(temporal_root_default)

parse_options <- function() {
  out <- list(stage = "verify", base_root = base_root_default, temporal_root = temporal_root_default,
              outdir = NULL, verify_manifest = NULL, primary_evidence_log = NULL,
              holiday_csv = NULL, workers = 4L,
              production_authorized = FALSE)
  raw <- commandArgs(trailingOnly = TRUE)
  i <- 1L
  while (i <= length(raw)) {
    token <- raw[[i]]
    if (token == "--production-authorized") {
      out$production_authorized <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(token, "--") || i == length(raw)) stop(sprintf("Invalid argument: %s", token), call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
    if (!key %in% names(out)) stop(sprintf("Unknown option --%s", gsub("_", "-", key, fixed = TRUE)), call. = FALSE)
    value <- raw[[i + 1L]]
    if (key == "workers") out[[key]] <- as.integer(value) else out[[key]] <- value
    i <- i + 2L
  }
  if (!out$stage %in% c("verify", "main", "clinical", "all")) stop("--stage must be verify, main, clinical, or all", call. = FALSE)
  if (!out$workers %in% 1:4) stop("--workers must be from 1 to 4", call. = FALSE)
  if (out$stage != "verify" && !out$production_authorized) stop("production stages require --production-authorized", call. = FALSE)
  out$base_root <- normalizePath(out$base_root, winslash = "/", mustWork = TRUE)
  out$temporal_root <- normalizePath(out$temporal_root, winslash = "/", mustWork = TRUE)
  if (is.null(out$holiday_csv)) out$holiday_csv <- file.path(out$temporal_root, "results/temporal/inputs/national_holidays.csv")
  out$holiday_manifest <- file.path(dirname(out$holiday_csv), "national_holidays_manifest.json")
  if (is.null(out$verify_manifest)) out$verify_manifest <- file.path(out$temporal_root, "results/temporal/verification_manifest.json")
  if (is.null(out$primary_evidence_log)) {
    out$primary_evidence_log <- file.path(out$temporal_root, "logs/temporal/m04_verify_20261113_rev4.log")
  }
  if (is.null(out$outdir)) {
    run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
    out$outdir <- file.path(out$temporal_root, "results/temporal", paste0("m04_", out$stage, "_", run_id))
  }
  out$holiday_csv <- normalizePath(out$holiday_csv, winslash = "/", mustWork = TRUE)
  out$holiday_manifest <- normalizePath(out$holiday_manifest, winslash = "/", mustWork = TRUE)
  out$verify_manifest <- normalizePath(out$verify_manifest, winslash = "/", mustWork = FALSE)
  out$primary_evidence_log <- normalizePath(out$primary_evidence_log, winslash = "/", mustWork = out$stage == "verify")
  out$outdir <- normalizePath(out$outdir, winslash = "/", mustWork = FALSE)
  if (out$stage != "verify") {
    if (!file.exists(out$verify_manifest)) stop(sprintf("verified baseline manifest not found: %s", out$verify_manifest), call. = FALSE)
    verified <- jsonlite::fromJSON(out$verify_manifest, simplifyVector = FALSE)
    engine_now <- sha256_file(file.path(out$temporal_root, "scripts/temporal/repair_engine_temporal.R"))
    if (!isTRUE(verified$passed) || !identical(verified$engine_sha256, engine_now)) {
      stop("baseline verification is not passed for the current temporal engine", call. = FALSE)
    }
  }
  out
}

sha256_file <- function(path) {
  cmd <- Sys.which("sha256sum")
  if (!nzchar(cmd)) stop("sha256sum is required", call. = FALSE)
  ans <- suppressWarnings(system2(cmd, shQuote(path), stdout = TRUE, stderr = TRUE))
  status <- attr(ans, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("sha256sum failed for %s", path), call. = FALSE)
  sub("[[:space:]].*$", "", ans[[1L]])
}

write_json_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".m04-json-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writeLines(as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
    digits = NA, pretty = TRUE)), tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop(sprintf("failed to atomically write %s", path), call. = FALSE)
}

log_event <- function(event, ...) {
  cat(as.character(jsonlite::toJSON(c(list(event = event,
    time_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")), list(...)),
    auto_unbox = TRUE, null = "null", na = "null", digits = NA)), "\n", sep = "")
  flush.console()
}

load_validated_definitions <- function(path) {
  wanted <- c("parse_date_safe", "canonical_unit", "normalize_dx4", "map_sex", "scenario_identity",
    "fread_path", "read_calendar", "read_clinical", "build_daily_panel", "unit_weather_target",
    "design_at_endpoint", "pure_basis_coef", "linear_effect", "interaction_results")
  exprs <- parse(path, keep.source = FALSE)
  found <- character()
  for (expr in exprs) {
    if (is.call(expr) && as.character(expr[[1L]]) %in% c("<-", "=") && is.symbol(expr[[2L]]) &&
        as.character(expr[[2L]]) %in% wanted) {
      eval(expr, envir = .GlobalEnv)
      found <- c(found, as.character(expr[[2L]]))
    }
  }
  missing <- setdiff(wanted, found)
  if (length(missing)) stop(sprintf("validated clinical definitions missing: %s", paste(missing, collapse = ", ")))
  invisible(found)
}

read_holidays <- function(path, country, dates) {
  country_id <- as.character(country)
  h <- fread(path, select = c("country", "date", "national_holiday"),
             colClasses = list(character = c("country", "date")))
  country_keep <- h[["country"]] == country_id
  h <- h[which(country_keep)]
  if (!nrow(h)) stop(sprintf("holiday calendar has no country coverage for %s", country_id), call. = FALSE)
  h[, date := as.Date(date)]
  h[, national_holiday := as.integer(national_holiday)]
  if (anyNA(h$date) || anyNA(h$national_holiday) || any(!h$national_holiday %in% c(0L, 1L)) ||
      anyDuplicated(h$date)) stop(sprintf("invalid holiday calendar for %s", country), call. = FALSE)
  idx <- match(as.Date(dates), h$date)
  if (anyNA(idx)) stop(sprintf("holiday calendar does not cover all %s calendar dates", country), call. = FALSE)
  h$national_holiday[idx]
}

fread_safe <- function(path, ...) {
  is_gzip <- grepl("\\.gz$", path, ignore.case = TRUE)
  gzip <- if (is_gzip && .Platform$OS.type == "unix") Sys.which("gzip") else ""
  if (nzchar(gzip) && !requireNamespace("R.utils", quietly = TRUE)) {
    return(data.table::fread(cmd = paste(shQuote(gzip), "-dc --", shQuote(path)), ...))
  }
  data.table::fread(path, ...)
}

read_main_calendar <- function(path, country) {
  required <- c("unit", "date", "E3", "y_K80", "y_K81", "y_urinary")
  cal <- fread_safe(path, select = required, colClasses = list(character = c("unit", "date")))
  cal[, unit := as.character(unit)]
  cal[, date := as.Date(date)]
  for (nm in setdiff(required, c("unit", "date"))) set(cal, j = nm, value = as.numeric(cal[[nm]]))
  if (anyNA(cal$date) || anyDuplicated(cal[, .(unit, date)])) stop(sprintf("invalid %s unit-date calendar", country), call. = FALSE)
  setorder(cal, unit, date)
  cal[, national_holiday := read_holidays(opt$holiday_csv, country, date)]
  cal
}

read_target_manifest <- function(country) {
  path <- if (country == "BR") file.path(opt$base_root, "results/fits/primary_cache_manifest.csv") else
    file.path(opt$base_root, "results/mx_updated/fits/fit_cache_manifest.csv")
  if (!file.exists(path)) stop(sprintf("primary target-cache manifest missing: %s", path), call. = FALSE)
  z <- fread(path, colClasses = list(character = c("country", "unit", "scenario", "cache_path")))
  if (country == "BR") z <- z[country == "BR" & scenario == "primary_threeoutcome_E3_separate"] else
    z <- z[country == "MX" & scenario == "primary_complete_source_2015_2022"]
  if (anyDuplicated(z$unit)) stop(sprintf("duplicate %s primary cache units", country), call. = FALSE)
  attr(z, "path") <- path
  z
}

target_from_cache <- function(row) {
  cache <- readRDS(row$cache_path)
  fit <- cache$engine_fit %||% cache$fit
  if (is.null(fit)) stop(sprintf("cache has no fit object: %s", row$cache_path), call. = FALSE)
  basis <- if (!is.null(cache$engine_fit)) fit$basis else cache$target_basis %||% fit$basis
  p95 <- as.numeric(cache$p95 %||% fit$p95)
  median <- as.numeric(cache$median %||% fit$median)
  if (length(p95) != 1L || length(median) != 1L || !is.finite(p95) || !is.finite(median) || is.null(basis)) {
    stop(sprintf("primary cache target is incomplete: %s", row$cache_path), call. = FALSE)
  }
  if (abs(p95 - as.numeric(row$p95)) > 1e-10 || abs(median - as.numeric(row$median)) > 1e-10) {
    stop(sprintf("primary-cache manifest target mismatch: %s", row$cache_path), call. = FALSE)
  }
  list(cache = cache, fit = fit, basis = basis, p95 = p95, median = median)
}

primary_eligible <- function(d) {
  outcomes <- c("y_K80", "y_K81", "y_urinary")
  keep <- is.finite(d$E3) & !is.na(d$date)
  for (nm in outcomes) keep <- keep & is.finite(d[[nm]]) & d[[nm]] >= 0
  events <- vapply(outcomes, function(nm) sum(d[[nm]][keep]), numeric(1))
  list(keep = keep, n = sum(keep), events = events)
}

primary_cache_fit <- function(target, unit_data) {
  fit <- target$fit
  contrasts <- fit$contrasts[["y_urinary_vs_y_K80"]]
  if (is.null(contrasts)) stop("primary cached fit lacks y_urinary_vs_y_K80")
  list(estimate = as.numeric(contrasts$coefficient_contrast$by_covariance$hac14["estimate"]),
       variance = as.numeric(contrasts$coefficient_contrast$by_covariance$hac14["variance"]))
}

clinical_margins <- function(fit, endpoints) {
  out <- list()
  for (g in c("elective", "urgent")) {
    k <- endpoints[[paste0("K80_", g)]]
    u <- endpoints[[paste0("N20_", g)]]
    lk <- as.numeric(design_at_endpoint(fit, k, fit$p95) - design_at_endpoint(fit, k, fit$median))
    lu <- as.numeric(design_at_endpoint(fit, u, fit$p95) - design_at_endpoint(fit, u, fit$median))
    q <- linear_effect(lu - lk, fit$fit, fit$hac_vcov$hac14)
    out[[g]] <- q
  }
  out
}

read_old_clinical_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing old clinical result: %s", path), call. = FALSE)
  z <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  data.table::setDT(z)
  if ("unit" %in% names(z)) data.table::set(z, j = "unit", value = as.character(z[["unit"]]))
  z
}

read_primary_preflight_evidence <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hit <- grep('"event":"primary_preflight_metrics"', lines, fixed = TRUE)
  if (length(hit) != 1L) stop(sprintf("expected one primary-preflight event in %s", path), call. = FALSE)
  record <- jsonlite::fromJSON(lines[[hit]], simplifyVector = FALSE)
  if (!identical(record$event, "primary_preflight_metrics")) stop("primary evidence event label is invalid", call. = FALSE)
  payload_names <- setdiff(names(record), c("event", "time_utc"))
  if (length(payload_names) != 1L || !is.list(record[[payload_names[[1L]]]])) {
    stop("primary evidence event has an unexpected payload", call. = FALSE)
  }
  metrics <- record[[payload_names[[1L]]]]
  if (!identical(as.character(metrics$unit), "11001")) stop("primary evidence is not for BR unit 11001", call. = FALSE)
  numeric_names <- c("p95_calendar", "p95_cache", "median_calendar", "median_cache",
    "p95_abs_difference", "median_abs_difference", "basis_probe_max_abs_difference",
    "two_outcome_logRRR_hac14", "cached_three_outcome_logRRR_hac14", "estimate_abs_difference",
    "two_outcome_variance_hac14", "cached_three_outcome_variance_hac14", "variance_abs_difference",
    "refit_three_outcome_logRRR_hac14", "refit_three_outcome_variance_hac14",
    "refit_three_outcome_estimate_abs_difference", "refit_three_outcome_variance_abs_difference",
    "n_complete")
  missing <- setdiff(c("unit", numeric_names), names(metrics))
  if (length(missing)) stop(sprintf("primary evidence lacks fields: %s", paste(missing, collapse = ", ")))
  values <- vapply(numeric_names, function(nm) as.numeric(metrics[[nm]]), numeric(1))
  if (any(!is.finite(values))) stop("primary evidence contains non-finite metrics")
  derived_differences <- c(
    abs(values[["two_outcome_logRRR_hac14"]] - values[["cached_three_outcome_logRRR_hac14"]]),
    abs(values[["two_outcome_variance_hac14"]] - values[["cached_three_outcome_variance_hac14"]]),
    abs(values[["refit_three_outcome_logRRR_hac14"]] - values[["cached_three_outcome_logRRR_hac14"]]),
    abs(values[["refit_three_outcome_variance_hac14"]] - values[["cached_three_outcome_variance_hac14"]])
  )
  recorded_differences <- c(values[["estimate_abs_difference"]], values[["variance_abs_difference"]],
    values[["refit_three_outcome_estimate_abs_difference"]], values[["refit_three_outcome_variance_abs_difference"]])
  if (max(abs(derived_differences - recorded_differences)) > 1e-12) {
    stop("primary evidence estimate/variance deltas are internally inconsistent")
  }
  pass_differences <- c(values[["p95_abs_difference"]], values[["median_abs_difference"]],
    values[["basis_probe_max_abs_difference"]], values[["refit_three_outcome_estimate_abs_difference"]],
    values[["refit_three_outcome_variance_abs_difference"]])
  if (values[["n_complete"]] < 250 || max(pass_differences) > 1e-7) {
    stop("saved three-outcome primary preflight metrics do not pass the original tolerance")
  }
  metrics$unit <- as.character(metrics$unit)
  for (nm in numeric_names) metrics[[nm]] <- as.numeric(metrics[[nm]])
  list(metrics = metrics, diagnostics = data.table::as.data.table(metrics))
}

verify_real_baselines <- function() {
  validated_path <- file.path(opt$base_root, "scripts/run_clinical_validated.R")
  engine_path <- file.path(opt$temporal_root, "scripts/temporal/repair_engine_temporal.R")
  load_validated_definitions(validated_path)
  source(engine_path, local = .GlobalEnv)
  real_unit <- "11001"
  primary_evidence <- read_primary_preflight_evidence(opt$primary_evidence_log)
  main_preflight <- primary_evidence$diagnostics
  main_metrics <- primary_evidence$metrics
  br_manifest <- read_target_manifest("BR")
  br_cache_row <- br_manifest[unit == real_unit]
  if (nrow(br_cache_row) != 1L) stop("verification BR unit is absent from primary target manifest")
  if (abs(as.numeric(br_cache_row$p95) - main_metrics$p95_cache) > 1e-10 ||
      abs(as.numeric(br_cache_row$median) - main_metrics$median_cache) > 1e-10) {
    stop("saved primary evidence targets differ from the current BR cache manifest")
  }
  primary_check <- data.table(
    check = "BR_11001_three_outcome_refit_matches_cached_primary_and_two_outcome_delta_recorded",
    status = "passed", unit = real_unit, target_p95 = main_metrics$p95_cache, target_median = main_metrics$median_cache,
    target_basis_max_probe_difference = main_metrics$basis_probe_max_abs_difference,
    target_p95_difference = main_metrics$p95_abs_difference, target_median_difference = main_metrics$median_abs_difference,
    two_vs_cached_three_estimate_abs_difference = main_metrics$estimate_abs_difference,
    two_vs_cached_three_variance_abs_difference = main_metrics$variance_abs_difference,
    refit_three_vs_cached_estimate_abs_difference = main_metrics$refit_three_outcome_estimate_abs_difference,
    refit_three_vs_cached_variance_abs_difference = main_metrics$refit_three_outcome_variance_abs_difference,
    primary_n_complete = main_metrics$n_complete, formula = NA_character_
  )

  clinical_calendar <- read_calendar(file.path(opt$base_root, "data/br_calendar.csv.gz"))$data
  clinical <- read_clinical(file.path(opt$base_root, "results/brazil_clinical_aggregate_all_admission_years_by_source_year.csv.gz"))
  scenario <- list(country = "brazil", name = "matched_source_2014_2020",
    label = "matched-source replication; source years 2014-2020", years = 2014:2020,
    analysis_years = 2014:2020, offset_col = "population", modes = "admission_type")
  panel_object <- build_daily_panel(clinical, clinical_calendar, scenario)
  panel <- panel_object$panel
  panel[, national_holiday := read_holidays(opt$holiday_csv, "BR", date)]
  cd <- panel[unit == real_unit]
  old_path <- file.path(opt$base_root, "results/clinical_full_BR/clinical_brazil_interaction_per_unit.csv")
  old <- read_old_clinical_csv(old_path)
  old <- old[scenario == "matched_source_2014_2020" & nuisance == "separate" & mode == "admission_type" & unit == real_unit]
  if (nrow(old) != 1L || old$status != "ok") stop("verification clinical unit is not a successful old baseline unit")
  old_group_path <- file.path(opt$base_root, "results/clinical_full_BR/clinical_brazil_interaction_group_per_unit.csv")
  old_group <- read_old_clinical_csv(old_group_path)
  old_group <- old_group[scenario == "matched_source_2014_2020" & nuisance == "separate" &
    mode == "admission_type" & unit == real_unit & status == "ok"]
  setorder(old_group, group_order)
  expected_groups <- c("elective", "urgent")
  if (nrow(old_group) != 2L || anyDuplicated(old_group$group) ||
      !identical(as.character(old_group$group), expected_groups) || any(old_group$group_order != seq_along(expected_groups))) {
    stop("verification flat clinical margins do not contain one ordered elective/urgent row each")
  }
  target <- unit_weather_target(cd)
  if (abs(target$p95 - old$p95) > 1e-10 || abs(target$median - old$median) > 1e-10) {
    stop("clinical full-period temperature target does not match the old baseline")
  }
  endpoints <- list(K80_elective = "y_K80_elective", N20_elective = "y_N20_elective",
                    K80_urgent = "y_K80_urgent", N20_urgent = "y_N20_urgent")
  outcomes <- unlist(endpoints, use.names = FALSE)
  counts <- vapply(outcomes, function(nm) sum(cd[[nm]], na.rm = TRUE), numeric(1))
  if (any(counts <= 0)) stop("verification clinical unit violates positive-count baseline eligibility")
  clinical_fit <- fit_stacked_contrast(
    data = cd, outcomes = outcomes, exposure = "E3", date = "date", offset = "population",
    basis = target$basis, p95 = target$p95, median = target$median, nuisance = "separate",
    hac_lags = c(14L, 28L), full_covariance = TRUE
  )
  clinical_fit$p95 <- target$p95; clinical_fit$median <- target$median
  direct <- interaction_results(clinical_fit, endpoints, "urgent", "elective", "admission_type")
  margins <- clinical_margins(clinical_fit, endpoints)
  direct_est_diff <- abs(direct$logRRR_difference - as.numeric(old$logRRR_difference))
  direct_var_diff <- abs(direct$variance_hac14 - as.numeric(old$variance_hac14))
  margin_diffs <- unlist(lapply(expected_groups, function(g) {
    old_margin <- old_group[group == g]
    z <- margins[[g]]
    c(setNames(abs(as.numeric(old_margin$logRRR) - as.numeric(z["estimate"])), paste0(g, "_estimate")),
      setNames(abs(as.numeric(old_margin$variance_hac14) - as.numeric(z["variance"])), paste0(g, "_variance")))
  }))
  if (any(!is.finite(c(direct_est_diff, direct_var_diff, margin_diffs))) ||
      max(c(direct_est_diff, direct_var_diff, margin_diffs)) > 1e-7) {
    stop("real-unit clinical direct/marginal baseline verification failed")
  }
  clinical_check <- data.table(
    check = "BR_11001_validated_four_outcome_baseline_wrapper_reproduces_direct_and_margins",
    status = "passed", unit = real_unit, target_p95 = target$p95, target_median = target$median,
    direct_estimate_abs_difference = direct_est_diff, direct_variance_abs_difference = direct_var_diff,
    margin_max_estimate_or_variance_abs_difference = max(margin_diffs), n_days = nrow(cd),
    events_K80_elective = counts[["y_K80_elective"]], events_N20_elective = counts[["y_N20_elective"]],
    events_K80_urgent = counts[["y_K80_urgent"]], events_N20_urgent = counts[["y_N20_urgent"]],
    formula = clinical_fit$audit$formula
  )
  check_dt <- rbindlist(list(primary_check, clinical_check), fill = TRUE)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  fwrite(main_preflight, file.path(opt$outdir, "primary_preflight_diagnostics.csv"), na = "")
  fwrite(check_dt, file.path(opt$outdir, "real_unit_baseline_checks.csv"), na = "")
  inputs <- c(file.path(opt$base_root, "data/br_calendar.csv.gz"),
    file.path(opt$base_root, "results/brazil_clinical_aggregate_all_admission_years_by_source_year.csv.gz"),
    attr(br_manifest, "path"), br_cache_row$cache_path, opt$primary_evidence_log,
    file.path(opt$base_root, "results/clinical_full_BR/clinical_brazil_interaction_per_unit.csv"),
    old_group_path,
    file.path(opt$temporal_root, "scripts/temporal/repair_engine_temporal.R"), validated_path,
    opt$holiday_csv, opt$holiday_manifest)
  input_hashes <- lapply(inputs, sha256_file); names(input_hashes) <- basename(inputs)
  manifest <- list(passed = TRUE, verified_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    engine_sha256 = sha256_file(engine_path), runner_sha256 = sha256_file(script_path),
    R = R.version.string, workers = 1L, native_threads = 1L, unit = real_unit,
    primary_evidence_log = opt$primary_evidence_log, primary_evidence_reused_without_refit = TRUE,
    clinical_unit_verified_in_old_411_success_set = TRUE, checks = check_dt, input_hashes = input_hashes)
  write_json_atomic(manifest, opt$verify_manifest)
  log_event("real_unit_verification_passed", manifest = opt$verify_manifest,
    primary_two_vs_three_estimate_abs_difference = main_metrics$estimate_abs_difference,
    primary_two_vs_three_variance_abs_difference = main_metrics$variance_abs_difference,
    clinical_direct_estimate_abs_difference = direct_est_diff,
    clinical_direct_variance_abs_difference = direct_var_diff,
    clinical_margin_max_abs_difference = max(margin_diffs))
  invisible(manifest)
}

main_attempt <- function(country, nuisance, unit, d, cache_row) {
  fail <- function(reason, n = NA_integer_, events = c(y_K80 = NA_real_, y_K81 = NA_real_, y_urinary = NA_real_)) data.table(
    country = country, unit = unit, nuisance = nuisance, contrast = "y_urinary_vs_y_K80",
    covariance = NA_character_, status = "failed", failure = reason, n_complete = n,
    events_K80 = events[["y_K80"]], events_K81 = events[["y_K81"]], events_urinary = events[["y_urinary"]])
  if (is.null(cache_row) || nrow(cache_row) != 1L) return(fail("primary_target_cache_missing"))
  tryCatch({
    d <- copy(d)
    eligible <- primary_eligible(d)
    event_named <- eligible$events
    if (eligible$n < 250L || any(event_named < 40)) return(fail("primary_paired_eligibility_failed", eligible$n, event_named))
    target <- target_from_cache(cache_row)
    fit <- fit_stacked_contrast(
      data = d[eligible$keep], outcomes = c("y_K80", "y_K81", "y_urinary"), exposure = "E3", date = "date",
      basis = target$basis, p95 = target$p95, median = target$median, nuisance = nuisance,
      hac_lags = c(14L, 28L), full_covariance = TRUE
    )
    contrast <- fit$contrasts[["y_urinary_vs_y_K80"]]$coefficient_contrast
    rows <- lapply(names(contrast$by_covariance), function(cov_name) {
      q <- contrast$by_covariance[[cov_name]]
      data.table(country = country, unit = unit, nuisance = nuisance, contrast = "y_urinary_vs_y_K80",
        covariance = cov_name, status = "ok", failure = NA_character_, estimate = as.numeric(q["estimate"]),
        variance = as.numeric(q["variance"]), se = as.numeric(q["se"]), RR = as.numeric(q["RR"]),
        lo = as.numeric(q["lo"]), hi = as.numeric(q["hi"]), n_complete = eligible$n,
        events_K80 = event_named[["y_K80"]], events_K81 = event_named[["y_K81"]],
        events_urinary = event_named[["y_urinary"]], p95 = target$p95, median = target$median,
        flexible_time_df = fit$flexible_time_df, national_holiday_days = fit$national_holiday_days,
        target_cache_path = cache_row$cache_path, target_n_complete = as.integer(cache_row$n_complete),
        formula = fit$audit$formula, rank = fit$audit$rank, columns = fit$audit$columns,
        warnings = paste(fit$audit$warnings, collapse = " | "))
    })
    rbindlist(rows, fill = TRUE)
  }, error = function(e) fail(conditionMessage(e)))
}

run_main_models <- function() {
  source(file.path(opt$temporal_root, "scripts/temporal/repair_engine_temporal.R"), local = .GlobalEnv)
  results <- list(); memberships <- list()
  countries <- if (opt$stage == "main" || opt$stage == "all") c("BR", "MX") else character()
  for (country in countries) {
    log_event("main_calendar_read_start", country = country)
    cal_path <- file.path(opt$base_root, if (country == "BR") "data/br_calendar.csv.gz" else "data/mx_updated_calendar.csv.gz")
    cal <- read_main_calendar(cal_path, country)
    target_manifest <- read_target_manifest(country)
    units <- sort(intersect(unique(cal$unit), target_manifest$unit))
    if (!length(units)) stop(sprintf("no smoke unit has a %s primary cache", country))
    units <- units[[1L]]
    nuisances <- if (country == "MX") c("yearmonth", "flexible") else "flexible"
    for (nuisance in nuisances) {
      log_event("main_model_start", country = country, nuisance = nuisance, units = length(units), workers = opt$workers)
      worker <- function(unit) {
        unit_id <- as.character(unit)
        target_row <- target_manifest[which(target_manifest[["unit"]] == unit_id)]
        d <- cal[which(cal[["unit"]] == unit_id)]
        out <- main_attempt(country, nuisance, unit, d, target_row)
        list(rows = out)
      }
      groups <- split(units, ceiling(seq_along(units) / 20L))
      chunks <- list()
      completed <- 0L
      for (block in groups) {
        block_result <- parallel::mclapply(block, worker, mc.cores = min(opt$workers, length(block)),
          mc.preschedule = FALSE, mc.set.seed = FALSE)
        chunks[[length(chunks) + 1L]] <- rbindlist(lapply(block_result, `[[`, "rows"), fill = TRUE)
        completed <- completed + length(block)
        latest <- rbindlist(chunks, fill = TRUE)
        fwrite(latest, file.path(opt$outdir, sprintf("checkpoint_main_%s_%s.csv", country, nuisance)), na = "")
        log_event("main_model_progress", country = country, nuisance = nuisance, completed = completed,
          total = length(units), successful_units = uniqueN(latest[status == "ok", unit]),
          failed_units = uniqueN(latest[status == "failed", unit]))
      }
      result <- rbindlist(chunks, fill = TRUE)
      results[[paste(country, nuisance, sep = "_")]] <- result
    }
    target_units <- unique(target_manifest$unit)
    membership <- data.table(country = country, unit = units, primary_target_cache_member = units %in% target_units)
    for (nuisance in nuisances) {
      z <- results[[paste(country, nuisance, sep = "_")]]
      membership[, (paste0(nuisance, "_success")) := unit %in% unique(z[status == "ok", unit])]
    }
    memberships[[country]] <- membership
    rm(cal); gc(verbose = FALSE)
  }
  per_unit <- rbindlist(results, fill = TRUE)
  failures <- per_unit[status != "ok"]
  fwrite(per_unit, file.path(opt$outdir, "main_per_unit.csv"), na = "")
  fwrite(failures, file.path(opt$outdir, "main_failures.csv"), na = "")
  membership <- rbindlist(memberships, fill = TRUE)
  fwrite(membership, file.path(opt$outdir, "main_unit_membership.csv"), na = "")
  pooled <- list()
  groups <- unique(per_unit[status == "ok", .(country, nuisance, contrast, covariance)])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- per_unit[country == g$country & nuisance == g$nuisance & contrast == g$contrast & covariance == g$covariance & status == "ok"]
    if (anyDuplicated(z$unit)) stop("duplicate successful main rows within a pooling cell", call. = FALSE)
    for (method in c("DL", "REML_HK")) {
      q <- pool_meta(z$estimate, z$variance, method = method)
      pooled[[length(pooled) + 1L]] <- data.table(country = g$country, nuisance = g$nuisance,
        contrast = g$contrast, covariance = g$covariance, pooling = method, units_attempted = uniqueN(per_unit[country == g$country & nuisance == g$nuisance, unit]),
        units_fit = uniqueN(z$unit), k = q$k, estimate = q$estimate, se = q$se, RR = q$RR,
        lo = q$RR95CI[1L], hi = q$RR95CI[2L], I2 = q$I2, tau2 = q$tau2,
        concordance = q$concordance, failure = if (is.null(q$failure)) NA_character_ else q$failure)
    }
  }
  pooled_dt <- rbindlist(pooled, fill = TRUE)
  fwrite(pooled_dt, file.path(opt$outdir, "main_pooled_DL_REML.csv"), na = "")
  log_event("main_models_complete", unit_attempt_rows = nrow(membership), fits_ok = uniqueN(per_unit[status == "ok", paste(country, nuisance, unit)]), failures = nrow(failures))
  list(per_unit = per_unit, pooled = pooled_dt, membership = membership, failures = failures)
}

all_attempts_from_old <- function(success_path, failure_path, scenario, expected_units) {
  scenario_id <- as.character(scenario)
  success <- fread(success_path, colClasses = list(character = c("unit", "scenario")))
  success <- success[which(success[["scenario"]] == scenario_id & success[["nuisance"]] == "separate" &
    success[["mode"]] == "admission_type")]
  failure <- fread(failure_path, colClasses = list(character = c("unit", "scenario")))
  failure <- failure[which(failure[["scenario"]] == scenario_id & failure[["nuisance"]] == "separate" &
    failure[["mode"]] == "admission_type")]
  if (anyDuplicated(success$unit) || anyDuplicated(failure$unit) || length(intersect(success$unit, failure$unit))) {
    stop(sprintf("old bridge %s has duplicate or overlapping success/failure units", scenario), call. = FALSE)
  }
  if (!setequal(c(success$unit, failure$unit), expected_units)) {
    stop(sprintf("old bridge %s does not cover the complete expected unit set", scenario), call. = FALSE)
  }
  attempt <- data.table(unit = expected_units)
  attempt[, status := "not_recorded"]
  add <- function(x) {
    for (nm in setdiff(names(x), "unit")) {
      if (!nm %in% names(attempt)) set(attempt, j = nm, value = x[[nm]][rep(NA_integer_, nrow(attempt))])
      attempt[match(x$unit, unit), (nm) := x[[nm]]]
    }
  }
  if (nrow(success)) add(success)
  if (nrow(failure)) add(failure)
  setorder(attempt, unit)
  attempt
}

clinical_attempt <- function(unit, d, model_name, nuisance, drop_2020 = FALSE, expected_target = NULL) {
  target <- NULL
  fail <- function(reason) data.table(country = "brazil", scenario = "latebilling_reconstruction_2014_2021",
    model = model_name, nuisance = nuisance, unit = unit, status = "failed", failure = reason,
    n_days = if (is.null(d)) NA_integer_ else nrow(d))
  tryCatch({
    target <- unit_weather_target(d)
    if (!is.null(expected_target) && all(is.finite(expected_target))) {
      if (abs(target$p95 - expected_target[["p95"]]) > 1e-10 ||
          abs(target$median - expected_target[["median"]]) > 1e-10) {
        stop("full-clinical-period temperature target differs from the old latebilling bridge")
      }
    }
    dfit <- if (drop_2020) d[as.integer(format(date, "%Y")) != 2020L] else d
    endpoints <- list(K80_elective = "y_K80_elective", N20_elective = "y_N20_elective",
                      K80_urgent = "y_K80_urgent", N20_urgent = "y_N20_urgent")
    outcomes <- unlist(endpoints, use.names = FALSE)
    counts <- vapply(outcomes, function(nm) sum(dfit[[nm]], na.rm = TRUE), numeric(1))
    if (any(!is.finite(counts) | counts <= 0)) stop(sprintf("interaction endpoint has zero/invalid counts: %s",
      paste(names(counts)[!is.finite(counts) | counts <= 0], collapse = ",")))
    fit <- fit_stacked_contrast(
      data = dfit, outcomes = outcomes, exposure = "E3", date = "date", offset = "population",
      basis = target$basis, p95 = target$p95, median = target$median, nuisance = nuisance,
      hac_lags = c(14L, 28L), full_covariance = TRUE
    )
    fit$p95 <- target$p95; fit$median <- target$median
    effect <- interaction_results(fit, endpoints, "urgent", "elective", "admission_type")
    margin <- clinical_margins(fit, endpoints)
    time_df <- fit$flexible_time_df
    row <- data.table(country = "brazil", scenario = "latebilling_reconstruction_2014_2021", model = model_name,
      nuisance = nuisance, unit = unit, status = "ok", failure = NA_character_, n_days = nrow(dfit),
      events_K80_elective = sum(dfit$y_K80_elective), events_N20_elective = sum(dfit$y_N20_elective),
      events_K80_urgent = sum(dfit$y_K80_urgent), events_N20_urgent = sum(dfit$y_N20_urgent),
      p95 = target$p95, median = target$median, target_exposure_days = target$exposure_days,
      flexible_time_df = time_df, national_holiday_days = if (nuisance == "flexible") sum(dfit$national_holiday == 1L) else NA_integer_,
      logRRR_difference = effect$logRRR_difference, variance_hac14 = effect$variance_hac14,
      se_hac14 = effect$se_hac14, RR_ratio = effect$RR_ratio_hac14,
      RR_lo = effect$RR_ratio_lo, RR_hi = effect$RR_ratio_hi,
      joint_wald_chisq3 = effect$joint_wald_chisq3, joint_wald_p = effect$joint_wald_p,
      joint_basis_coefficients_json = as.character(toJSON(effect$joint_basis_coefficients, auto_unbox = TRUE)),
      joint_basis_covariance_json = as.character(toJSON(effect$joint_basis_covariance, auto_unbox = TRUE)),
      formula = fit$audit$formula, rank = fit$audit$rank, columns = fit$audit$columns,
      warnings = paste(fit$audit$warnings, collapse = " | "))
    margins <- rbindlist(lapply(names(margin), function(g) {
      q <- margin[[g]]
      data.table(country = "brazil", scenario = "latebilling_reconstruction_2014_2021", model = model_name,
        nuisance = nuisance, unit = unit, group = g, estimate = as.numeric(q["estimate"]),
        variance_hac14 = as.numeric(q["variance"]), se_hac14 = as.numeric(q["se"]),
        RRR = as.numeric(q["RR"]), lo = as.numeric(q["lo"]), hi = as.numeric(q["hi"]),
        p95 = target$p95, median = target$median)
    }), fill = TRUE)
    list(row = row, margins = margins)
  }, error = function(e) list(row = fail(conditionMessage(e)), margins = data.table()))
}

pool_clinical <- function(per_unit, margins) {
  out <- list()
  groups <- unique(per_unit[status == "ok", .(model, nuisance)])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- per_unit[model == g$model & nuisance == g$nuisance & status == "ok"]
    if (anyDuplicated(z$unit)) stop("duplicate clinical interaction successes", call. = FALSE)
    for (method in c("DL", "REML_HK")) {
      q <- pool_meta(z$logRRR_difference, z$variance_hac14, method = method)
      out[[length(out) + 1L]] <- data.table(estimand = "urgent_minus_elective_logRRR_difference",
        model = g$model, nuisance = g$nuisance, pooling = method, units_fit = uniqueN(z$unit), k = q$k,
        estimate = q$estimate, se = q$se, RR_ratio = q$RR, lo = q$RR95CI[1L], hi = q$RR95CI[2L],
        I2 = q$I2, tau2 = q$tau2, concordance = q$concordance,
        failure = if (is.null(q$failure)) NA_character_ else q$failure)
    }
  }
  margin_groups <- if (nrow(margins)) unique(margins[, .(model, nuisance, group)]) else data.table()
  for (i in seq_len(nrow(margin_groups))) {
    g <- margin_groups[i]
    z <- margins[model == g$model & nuisance == g$nuisance & group == g$group]
    for (method in c("DL", "REML_HK")) {
      q <- pool_meta(z$estimate, z$variance_hac14, method = method)
      out[[length(out) + 1L]] <- data.table(estimand = paste0("within_", g$group, "_urinary_vs_K80"),
        model = g$model, nuisance = g$nuisance, pooling = method, units_fit = uniqueN(z$unit), k = q$k,
        estimate = q$estimate, se = q$se, RR_ratio = q$RR, lo = q$RR95CI[1L], hi = q$RR95CI[2L],
        I2 = q$I2, tau2 = q$tau2, concordance = q$concordance,
        failure = if (is.null(q$failure)) NA_character_ else q$failure)
    }
  }
  rbindlist(out, fill = TRUE)
}

bridge_margins <- function(rows, old_group_rows, scenario_name) {
  wanted_units <- as.character(rows$unit)
  wanted_groups <- c("elective", "urgent")
  old <- old_group_rows[scenario == scenario_name & nuisance == "separate" & mode == "admission_type" &
    status == "ok" & unit %in% wanted_units & group %in% wanted_groups]
  if (anyDuplicated(old[, .(unit, group)]) || uniqueN(old$unit) != length(wanted_units) ||
      nrow(old) != length(wanted_units) * length(wanted_groups) ||
      any(old[, .N, by = unit]$N != length(wanted_groups))) {
    stop(sprintf("flat old clinical margins do not match successful %s units", scenario_name), call. = FALSE)
  }
  setorder(old, unit, group_order)
  old[, .(unit, group, estimate = as.numeric(logRRR), variance_hac14 = as.numeric(variance_hac14),
          se_hac14 = as.numeric(se_hac14), RRR = as.numeric(RRR_hac14),
          lo = as.numeric(RR_lo_hac14), hi = as.numeric(RR_hi_hac14))]
}

run_clinical_models <- function() {
  validated_path <- file.path(opt$base_root, "scripts/run_clinical_validated.R")
  load_validated_definitions(validated_path)
  source(file.path(opt$temporal_root, "scripts/temporal/repair_engine_temporal.R"), local = .GlobalEnv)
  calendar <- read_calendar(file.path(opt$base_root, "data/br_calendar.csv.gz"))$data
  clinical <- read_clinical(file.path(opt$base_root, "results/brazil_clinical_aggregate_all_admission_years_by_source_year.csv.gz"))
  scenario <- list(country = "brazil", name = "latebilling_reconstruction_2014_2021",
    label = "late-billing reconstruction; source years 2014-2021; source-year 2021 contributes late 2020 admissions",
    years = 2014:2021, analysis_years = 2014:2020, offset_col = "population", modes = "admission_type")
  built <- build_daily_panel(clinical, calendar, scenario)
  panel <- built$panel
  panel[, national_holiday := read_holidays(opt$holiday_csv, "BR", date)]
  units <- sort(unique(panel$unit))
  if (length(units) != 450L) stop(sprintf("expected 450 BR clinical units, found %d", length(units)), call. = FALSE)
  if (!is.finite(built$audit$source_2021_admission_2020_rows_included) || built$audit$source_2021_admission_2020_rows_included <= 0) {
    stop("preferred late-billing cohort does not include source-year 2021 admission-2020 rows", call. = FALSE)
  }
  source_dir <- file.path(opt$base_root, "results/clinical_full_BR")
  success_path <- file.path(source_dir, "clinical_brazil_interaction_per_unit.csv")
  failure_path <- file.path(source_dir, "clinical_brazil_interaction_failures.csv")
  old_group_path <- file.path(source_dir, "clinical_brazil_interaction_group_per_unit.csv")
  old_group_rows <- read_old_clinical_csv(old_group_path)
  standard <- all_attempts_from_old(success_path, failure_path, "matched_source_2014_2020", units)
  bridge <- all_attempts_from_old(success_path, failure_path, "latebilling_reconstruction_2014_2021", units)
  target_rows <- bridge[status == "ok" & is.finite(p95) & is.finite(median), .(unit, p95, median)]
  standard_ok <- sum(standard$status == "ok"); standard_failed <- sum(standard$status == "failed")
  if (standard_ok != 411L || standard_failed != 39L) {
    stop(sprintf("expected old clinical standard 411/39 success/failure, got %d/%d", standard_ok, standard_failed))
  }
  fwrite(standard, file.path(opt$outdir, "clinical_original_standard_all_attempts.csv"), na = "")
  fwrite(bridge, file.path(opt$outdir, "clinical_latebilling_separate_baseline_all_attempts.csv"), na = "")

  variants <- list(
    list(name = "flexible_holiday_full_2014_2020", nuisance = "flexible", drop_2020 = FALSE),
    list(name = "separate_exclude_2020_full_targets", nuisance = "separate", drop_2020 = TRUE)
  )
  result_rows <- list(); margin_rows <- list()
  for (variant in variants) {
    log_event("clinical_model_start", model = variant$name, nuisance = variant$nuisance,
      units = length(units), time_df_rule = if (variant$nuisance == "flexible") "6*unique analysis years" else "legacy ns(time,12)",
      holidays = variant$nuisance == "flexible", source_2021_latebill_rows = built$audit$source_2021_admission_2020_rows_included)
    worker <- function(unit) {
      unit_id <- as.character(unit)
      d <- panel[which(panel[["unit"]] == unit_id)]
      expected <- target_rows[which(target_rows[["unit"]] == unit_id)]
      expected_target <- if (nrow(expected) == 1L) c(p95 = expected$p95, median = expected$median) else NULL
      clinical_attempt(unit, d, variant$name, variant$nuisance, variant$drop_2020, expected_target)
    }
    groups <- split(units, ceiling(seq_along(units) / 20L))
    chunks <- list(); margin_chunks <- list(); completed <- 0L
    for (block in groups) {
      block_result <- parallel::mclapply(block, worker, mc.cores = min(opt$workers, length(block)),
        mc.preschedule = FALSE, mc.set.seed = FALSE)
      chunks[[length(chunks) + 1L]] <- rbindlist(lapply(block_result, `[[`, "row"), fill = TRUE)
      margin_chunks[[length(margin_chunks) + 1L]] <- rbindlist(lapply(block_result, `[[`, "margins"), fill = TRUE)
      completed <- completed + length(block)
      latest <- rbindlist(chunks, fill = TRUE)
      fwrite(latest, file.path(opt$outdir, paste0("checkpoint_clinical_", variant$name, ".csv")), na = "")
      log_event("clinical_model_progress", model = variant$name, completed = completed, total = length(units),
        successful_units = uniqueN(latest[status == "ok", unit]), failed_units = uniqueN(latest[status == "failed", unit]))
    }
    result_rows[[variant$name]] <- rbindlist(chunks, fill = TRUE)
    margin_rows[[variant$name]] <- rbindlist(margin_chunks, fill = TRUE)
  }
  per_unit <- rbindlist(result_rows, fill = TRUE)
  margins <- rbindlist(margin_rows, fill = TRUE)
  failures <- per_unit[status != "ok"]
  pooled <- pool_clinical(per_unit, margins)
  fwrite(per_unit, file.path(opt$outdir, "clinical_variant_per_unit.csv"), na = "")
  fwrite(margins, file.path(opt$outdir, "clinical_variant_margins_per_unit.csv"), na = "")
  fwrite(failures, file.path(opt$outdir, "clinical_variant_failures.csv"), na = "")
  fwrite(pooled, file.path(opt$outdir, "clinical_variant_pooled_DL_REML.csv"), na = "")

  standard_status <- standard$status[match(units, standard$unit)]
  bridge_status <- bridge$status[match(units, bridge$unit)]
  if (anyNA(standard_status) || anyNA(bridge_status)) stop("clinical old-baseline membership is incomplete")
  membership <- data.table(unit = units, original_standard_success = standard_status == "ok",
                           latebilling_separate_success = bridge_status == "ok")
  for (variant in variants) {
    z <- result_rows[[variant$name]]
    status <- z$status[match(units, z$unit)]
    if (anyNA(status)) stop(sprintf("new variant membership incomplete: %s", variant$name))
    membership[, (paste0(variant$name, "_success")) := status == "ok"]
  }
  membership[, common_standard_flex := original_standard_success & flexible_holiday_full_2014_2020_success]
  membership[, common_latebilling_flex := latebilling_separate_success & flexible_holiday_full_2014_2020_success]
  membership[, common_standard_exclude2020 := original_standard_success & separate_exclude_2020_full_targets_success]
  membership[, common_latebilling_exclude2020 := latebilling_separate_success & separate_exclude_2020_full_targets_success]
  membership[, common_all_four := original_standard_success & latebilling_separate_success &
    flexible_holiday_full_2014_2020_success & separate_exclude_2020_full_targets_success]
  fwrite(membership, file.path(opt$outdir, "clinical_exact_unit_membership.csv"), na = "")

  # Pool the old baseline and new variant separately on exactly shared successful units;
  # these rows describe matched membership and are not a paired test of their difference.
  bridge_rows <- list()
  reference_sets <- list(
    original_standard = list(data = standard, model = "matched_source_2014_2020"),
    latebilling_separate = list(data = bridge, model = "latebilling_reconstruction_2014_2021")
  )
  for (variant in variants) for (ref_name in names(reference_sets)) {
    ref <- reference_sets[[ref_name]]$data
    ref_success <- ref[status == "ok"]
    var_success <- result_rows[[variant$name]][status == "ok"]
    common <- intersect(ref_success$unit, var_success$unit)
    if (!length(common)) next
    ref_rows <- ref_success[unit %in% common]
    var_rows <- var_success[unit %in% common]
    ref_margin <- bridge_margins(ref_rows, old_group_rows, reference_sets[[ref_name]]$model)
    var_margin <- margin_rows[[variant$name]][unit %in% common]
    for (role in c("old_baseline", "new_variant")) {
      direct <- if (role == "old_baseline") ref_rows else var_rows
      for (method in c("DL", "REML_HK")) {
        q <- pool_meta(direct$logRRR_difference, direct$variance_hac14, method = method)
        bridge_rows[[length(bridge_rows) + 1L]] <- data.table(reference = ref_name, variant = variant$name,
          role = role, estimand = "urgent_minus_elective_logRRR_difference", pooling = method,
          common_success_units = length(common), k = q$k, estimate = q$estimate, se = q$se,
          RR_ratio = q$RR, lo = q$RR95CI[1L], hi = q$RR95CI[2L], I2 = q$I2, tau2 = q$tau2,
          concordance = q$concordance,
          failure = if (is.null(q$failure)) NA_character_ else q$failure)
      }
      mm <- if (role == "old_baseline") ref_margin else var_margin
      for (group_name in c("elective", "urgent")) {
        mg <- mm[group == group_name]
        for (method in c("DL", "REML_HK")) {
          q <- pool_meta(mg$estimate, mg$variance_hac14, method = method)
          bridge_rows[[length(bridge_rows) + 1L]] <- data.table(reference = ref_name, variant = variant$name,
            role = role, estimand = paste0("within_", group_name, "_urinary_vs_K80"), pooling = method,
            common_success_units = length(common), k = q$k, estimate = q$estimate, se = q$se,
            RR_ratio = q$RR, lo = q$RR95CI[1L], hi = q$RR95CI[2L], I2 = q$I2, tau2 = q$tau2,
            concordance = q$concordance,
            failure = if (is.null(q$failure)) NA_character_ else q$failure)
        }
      }
    }
  }
  bridge_dt <- rbindlist(bridge_rows, fill = TRUE)
  fwrite(bridge_dt, file.path(opt$outdir, "clinical_common_success_baseline_bridge_pooled.csv"), na = "")
  write_json_atomic(built$audit, file.path(opt$outdir, "clinical_panel_source_audit.json"))
  log_event("clinical_models_complete", units_expected = length(units), old_standard_success = standard_ok,
    old_standard_failed = standard_failed, flex_success = uniqueN(per_unit[model == variants[[1L]]$name & status == "ok", unit]),
    flex_failed = uniqueN(per_unit[model == variants[[1L]]$name & status == "failed", unit]),
    exclude2020_success = uniqueN(per_unit[model == variants[[2L]]$name & status == "ok", unit]),
    exclude2020_failed = uniqueN(per_unit[model == variants[[2L]]$name & status == "failed", unit]),
    common_all_four = sum(membership$common_all_four))
  list(per_unit = per_unit, margins = margins, pooled = pooled, failures = failures,
       membership = membership, bridge_pooled = bridge_dt, audit = built$audit)
}

opt <- parse_options()
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(opt$outdir, all.files = TRUE, no.. = TRUE))) stop(sprintf("refusing to overwrite nonempty output directory: %s", opt$outdir))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", BLAS_NUM_THREADS = "1")
data.table::setDTthreads(1L)
log_event("m04_temporal_start", stage = opt$stage, base_root = opt$base_root, temporal_root = opt$temporal_root,
  outdir = opt$outdir, workers = opt$workers, native_threads = 1L, runner = script_path)

if (opt$stage == "verify") {
  result <- verify_real_baselines()
} else {
  outputs <- list()
  if (opt$stage %in% c("main", "all")) outputs$main <- run_main_models()
  if (opt$stage %in% c("clinical", "all")) outputs$clinical <- run_clinical_models()
  manifest_files <- c(file.path(opt$base_root, "data/br_calendar.csv.gz"),
    file.path(opt$base_root, "data/mx_updated_calendar.csv.gz"),
    file.path(opt$base_root, "results/brazil_clinical_aggregate_all_admission_years_by_source_year.csv.gz"),
    file.path(opt$base_root, "results/fits/primary_cache_manifest.csv"),
    file.path(opt$base_root, "results/mx_updated/fits/fit_cache_manifest.csv"),
    file.path(opt$base_root, "results/clinical_full_BR/clinical_brazil_interaction_per_unit.csv"),
    file.path(opt$base_root, "results/clinical_full_BR/clinical_brazil_interaction_failures.csv"),
    file.path(opt$base_root, "results/clinical_full_BR/clinical_brazil_interaction_group_per_unit.csv"),
    file.path(opt$temporal_root, "scripts/temporal/run_temporal_m04.R"),
    file.path(opt$temporal_root, "scripts/temporal/repair_engine_temporal.R"), opt$holiday_csv, opt$holiday_manifest)
  file_hashes <- lapply(manifest_files, sha256_file); names(file_hashes) <- make.unique(basename(manifest_files))
  run_manifest <- list(run_id = basename(opt$outdir), completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stage = opt$stage, workers = opt$workers, native_threads = 1L, R = R.version.string,
    verified_baseline_manifest = opt$verify_manifest, input_file_sha256 = file_hashes,
    holiday_calendar_scope = "National public holidays only; no assumption all hospitals close; 2006-2023 coverage",
    holiday_indicator = "holidays 0.108, observed=TRUE, subdiv=NULL; national_holiday only",
    holiday_generator_manifest = jsonlite::fromJSON(opt$holiday_manifest, simplifyVector = FALSE),
    flexible_formula = "disease*(B1+B2+B3+dow+ns(time,6*unique retained analysis years)+national_holiday); offset(population) only for clinical interaction",
    eligibility = list(main = "same primary paired K80/K81/urinary complete dates; >=250 dates and >=40 events per endpoint",
      clinical = "all four strict N20/K80 admission-type totals positive plus existing GLM full-rank/convergence audit; no threshold relaxation"),
    analyses = list("MX complete-source 2015-2022 year-month fixed effects",
      "BR and MX primary urinary/K80 flexible-time/national-holiday sensitivity with cached primary targets",
      "BR strict N20/K80 latebilling clinical urgency interaction, full 2014-2020 flexible-time/holiday",
      "BR strict N20/K80 latebilling clinical urgency interaction excluding 2020, separate legacy calendar controls and full-period targets"))
  file.copy(opt$holiday_manifest, file.path(opt$outdir, "national_holidays_manifest.json"), overwrite = FALSE)
  write_json_atomic(run_manifest, file.path(opt$outdir, "run_manifest.json"))
  log_event("m04_temporal_complete", manifest = file.path(opt$outdir, "run_manifest.json"))
}
