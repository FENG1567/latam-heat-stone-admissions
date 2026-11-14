#!/usr/bin/env Rscript

# Dictionary-supported secondary Mexico admission-origin analysis.
# Reuses the locked model helpers in run_mexico_origin.R while keeping the
# 2015-2022 main cohort and its output provenance unchanged.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE) else
  normalizePath("work/repair_20261112/scripts/run_mexico_origin_dictionary_window.R", winslash = "/", mustWork = FALSE)
repair_root <- dirname(dirname(script_path))
repair_rlib <- file.path(repair_root, "Rlib")
if (dir.exists(repair_rlib)) .libPaths(c(repair_rlib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# source() keeps run_mexico_origin.R's final direct-execution guard inactive.
source(file.path(repair_root, "scripts", "run_mexico_origin.R"), local = FALSE)

dictionary_window <- c("2016-01-01", "2022-12-31")
model_diagnoses <- c("K80", "N20", "N21", "N22", "N23")

parse_dictionary_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    root = repair_root,
    grouped = file.path(repair_root, "data", "mx_origin_grouped_2015_2023.csv.gz"),
    calendar = file.path(repair_root, "data", "mx_updated_calendar.csv.gz"),
    outdir = NULL,
    workers = 8L,
    max_units = Inf
  )
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--") || i == length(args)) stop(sprintf("unexpected argument: %s", token), call. = FALSE)
    key <- sub("^--", "", token)
    value <- args[[i + 1L]]
    if (key == "root") out$root <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "grouped") out$grouped <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "calendar") out$calendar <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "outdir") out$outdir <- normalizePath(value, winslash = "/", mustWork = FALSE)
    else if (key == "workers") out$workers <- as.integer(value)
    else if (key == "max-units") out$max_units <- as.integer(value)
    else stop(sprintf("unknown option: %s", token), call. = FALSE)
    i <- i + 2L
  }
  if (is.null(out$outdir)) stop("--outdir is required", call. = FALSE)
  if (out$workers < 1L || out$workers > 8L) stop("--workers must be 1-8", call. = FALSE)
  if (!is.infinite(out$max_units) && out$max_units < 1L) stop("--max-units must be positive", call. = FALSE)
  out
}

dictionary_source_year_coverage <- function(root, catalog_status_path) {
  if (!file.exists(catalog_status_path)) stop(sprintf("annual dictionary evidence not found: %s", catalog_status_path), call. = FALSE)
  catalog <- jsonlite::fromJSON(catalog_status_path, simplifyVector = FALSE)
  coverage <- catalog$year_coverage
  catalog_years <- as.integer(vapply(coverage, function(x) as.character(x$year), character(1)))
  catalog_status <- vapply(coverage, function(x) as.character(x$status), character(1))
  verified_years <- sort(catalog_years[grepl("^verified", catalog_status)])
  if (!identical(verified_years, 2016:2023)) {
    stop(sprintf("annual catalog evidence years are not exactly 2016-2023: %s", paste(verified_years, collapse = ",")), call. = FALSE)
  }
  verified_captures <- catalog$verified_captures
  capture_coverage <- vapply(verified_captures, function(x) as.character(x$coverage), character(1))
  eligible_capture <- capture_coverage %in% as.character(2016:2022) | grepl("^2023-", capture_coverage)
  capture_mappings <- verified_captures[eligible_capture]
  mapping_ok <- vapply(capture_mappings, function(x) {
    identical(as.character(x$code_1_outpatient_consultation), "CONSULTA EXTERNA") &&
      identical(as.character(x$code_2_urgencias), "URGENCIAS")
  }, logical(1))
  if (!all(mapping_ok) || !all(as.character(2016:2022) %in% capture_coverage) ||
      !any(grepl("^2023-", capture_coverage))) {
    stop("verified annual catalog captures do not consistently label PROCED=1/2 for source years 2016-2023", call. = FALSE)
  }

  coverage_rows <- list()
  for (source_year in 2015:2023) {
    path <- file.path(root, "data", sprintf("mx_origin_source_year_%d.csv.gz", source_year))
    if (!file.exists(path)) stop(sprintf("source-year grouped cache missing: %s", path), call. = FALSE)
    d <- fread_gzip(path, select = c("date", "dx3", "n"),
                    colClasses = list(character = c("date", "dx3")), showProgress = FALSE)
    d[, date := as.Date(date)]
    d[, dx3 := toupper(trimws(dx3))]
    d[, n := suppressWarnings(as.numeric(n))]
    if (anyNA(d$date) || any(!is.finite(d$n) | d$n < 0)) {
      stop(sprintf("source-year grouped cache has invalid dates or counts: %s", path), call. = FALSE)
    }
    selected <- d[date >= as.Date(dictionary_window[[1L]]) & date <= as.Date(dictionary_window[[2L]]) & dx3 %in% model_diagnoses]
    if (nrow(selected)) {
      counts <- selected[, .(events = sum(n)), by = .(event_year = as.integer(format(date, "%Y")))]
      counts[, source_year := source_year]
      coverage_rows[[length(coverage_rows) + 1L]] <- counts[, .(source_year, event_year, events)]
    } else {
      coverage_rows[[length(coverage_rows) + 1L]] <- data.table(source_year = source_year, event_year = NA_integer_, events = 0)
    }
  }
  source_counts <- rbindlist(coverage_rows, use.names = TRUE, fill = TRUE)
  included_source_years <- sort(unique(source_counts[events > 0, source_year]))
  if (!identical(included_source_years, 2016:2023)) {
    stop(sprintf("2016-2022 model events are not sourced exclusively from covered years 2016-2023: %s",
                 paste(included_source_years, collapse = ",")), call. = FALSE)
  }
  if (any(source_counts[source_year == 2015, events] > 0)) {
    stop("unresolved 2015 source-year file contributes events inside the dictionary-supported window", call. = FALSE)
  }
  list(
    verified_source_years = verified_years,
    unresolved_source_years = 2015L,
    catalog_status_sha256 = sha256_file(catalog_status_path),
    source_year_counts = source_counts,
    source_year_cache_hashes = stats::setNames(
      vapply(2015:2023, function(y) sha256_file(file.path(root, "data", sprintf("mx_origin_source_year_%d.csv.gz", y))), character(1)),
      as.character(2015:2023)
    )
  )
}

unit_model_worker_dictionary <- function(unit) {
  d <- panels_by_state[[unit]]
  fail <- function(message, n_days = sum(is.finite(d$E3))) {
    list(status = "failed", unit = unit, n_days = n_days, failure = message)
  }
  target <- tryCatch(unit_target(d), error = function(e) e)
  if (inherits(target, "error")) {
    result <- fail(conditionMessage(target))
    return(list(unit = unit, primary = result, strict = result, shortstay = result))
  }
  primary <- tryCatch(summarize_primary_fit(d, unit, target), error = function(e) fail(conditionMessage(e), target$exposure_days))
  strict <- tryCatch(summarize_strict_ed_fit(d, unit, target), error = function(e) fail(conditionMessage(e), target$exposure_days))
  short_panel <- shortstay_panels[[unit]]
  shortstay <- tryCatch(
    summarize_primary_fit(short_panel, unit, unit_target(short_panel)),
    error = function(e) fail(conditionMessage(e), sum(is.finite(short_panel$E3)))
  )
  list(unit = unit, primary = primary, strict = strict, shortstay = shortstay)
}

main_dictionary_window <- function() {
  args <- parse_dictionary_cli()
  if (!dir.exists(args$root)) stop(sprintf("repair root not found: %s", args$root), call. = FALSE)
  if (dir.exists(args$outdir) && length(list.files(args$outdir, all.files = TRUE, no.. = TRUE))) {
    stop(sprintf("refusing to overwrite non-empty output directory: %s", args$outdir), call. = FALSE)
  }
  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  amendment_path <- file.path(args$root, "admin", "mx_origin_dictionary_window_amendment_20261113.json")
  if (!file.exists(amendment_path)) stop(sprintf("pre-result amendment plan missing: %s", amendment_path), call. = FALSE)
  dictionary_status_path <- file.path(args$root, "admin", "origin_dictionary_crossyear", "crossyear_retrieval_status.json")
  source_year_coverage <- dictionary_source_year_coverage(args$root, dictionary_status_path)

  inputs <- read_model_inputs(args$root, args$grouped, args$calendar)
  source_codes <- inputs$clinical[
    date >= as.Date(dictionary_window[[1L]]) & date <= as.Date(dictionary_window[[2L]]) & dx3 %in% model_diagnoses
  ]
  bad_code1 <- source_codes[origin_raw == "1" & origin_group != "outpatient_consultation_origin"]
  bad_code2 <- source_codes[origin_raw == "2" & origin_group != "ED_origin"]
  ambiguous_codes_misclassified <- source_codes[origin_raw %in% c("4", "5") &
    origin_group %in% c("ED_origin", "outpatient_consultation_origin")]
  if (nrow(bad_code1) || nrow(bad_code2) || nrow(ambiguous_codes_misclassified)) {
    stop("raw PROCED code grouping violates the verified 1/2 map or includes changed code 4/5", call. = FALSE)
  }
  raw_code_counts <- source_codes[, .(records = sum(n)), by = .(origin_raw, origin_group, type_service)]
  setorder(raw_code_counts, origin_raw, origin_group, type_service)
  fwrite(raw_code_counts, file.path(args$outdir, "mx_origin_raw_code_counts_2016_2022.csv"), na = "")
  fwrite(source_year_coverage$source_year_counts,
         file.path(args$outdir, "mx_origin_source_year_coverage_2016_2022.csv"), na = "")
  file.copy(dictionary_status_path, file.path(args$outdir, "crossyear_retrieval_status.json"), overwrite = TRUE)

  primary_panel <- daily_outcomes(inputs$clinical, inputs$calendar, event_window = dictionary_window)
  short_panel <- daily_outcomes(inputs$clinical, inputs$calendar, event_window = dictionary_window, los_filter = "1to30")
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
  cat(sprintf("Running %d Mexico resident states for admission dates %s to %s with %d worker(s)\n",
              length(state_codes), dictionary_window[[1L]], dictionary_window[[2L]], min(args$workers, length(state_codes))))

  results <- if (min(args$workers, length(state_codes)) == 1L) {
    lapply(state_codes, unit_model_worker_dictionary)
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
      "panels_by_state", "shortstay_panels", "unit_model_worker_dictionary", "unit_target",
      "summarize_primary_fit", "summarize_strict_ed_fit", "formal_basis_test",
      "basis_endpoint_vector", "linear_effect"
    ), envir = .GlobalEnv)
    parallel::parLapply(cl, state_codes, unit_model_worker_dictionary)
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
      RRR = RRR_N20_only_vs_K80_ED, lo = lo_N20_only_vs_K80_ED, hi = hi_N20_only_vs_K80_ED, failure)]
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

  main_meta <- primary_rows[, .(model = "primary_ED_minus_outpatient", state, status,
    estimate = logRRR_ED_minus_outpatient, variance = var_ED_minus_outpatient_hac14)]
  sensitivity_meta <- short_rows[, .(model = "sensitivity_recorded_los_1to30", state, status,
    estimate = logRRR_ED_minus_outpatient, variance = var_ED_minus_outpatient_hac14)]
  supportive_meta <- support_rows[, .(model = analysis, state, status, estimate = logRRR, variance = variance_hac14)]
  meta_input <- rbindlist(list(main_meta, sensitivity_meta, supportive_meta), fill = TRUE)
  meta_duplicates <- meta_input[, .N, by = .(model, state)][N != 1L]
  if (nrow(meta_duplicates)) stop("duplicate state contrast rows detected before meta-analysis", call. = FALSE)
  meta_names <- unique(meta_input$model)
  pooled <- rbindlist(lapply(meta_names, function(nm) pool_effect(meta_input, "estimate", "variance", nm)), fill = TRUE)
  fwrite(pooled, file.path(args$outdir, "mx_origin_pooled.csv"), na = "")

  support_counts <- inputs$clinical[
    date >= as.Date(dictionary_window[[1L]]) & date <= as.Date(dictionary_window[[2L]]) & dx3 %in% model_diagnoses,
    .(n = sum(n)), by = .(origin_group, origin_raw, type_service)
  ]
  support_counts[, origin_raw := as.character(origin_raw)]
  audit <- list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    amendment_plan_path = amendment_path,
    amendment_plan_sha256 = sha256_file(amendment_path),
    model_script_path = script_path,
    model_script_sha256 = sha256_file(script_path),
    reused_model_script_path = file.path(args$root, "scripts", "run_mexico_origin.R"),
    reused_model_script_sha256 = sha256_file(file.path(args$root, "scripts", "run_mexico_origin.R")),
    specification_path = file.path(args$root, "admin", "mx_origin_specification.json"),
    specification_sha256 = sha256_file(file.path(args$root, "admin", "mx_origin_specification.json")),
    dictionary_status_path = dictionary_status_path,
    dictionary_status_sha256 = source_year_coverage$catalog_status_sha256,
    dictionary_verified_source_years = source_year_coverage$verified_source_years,
    unresolved_dictionary_source_years_excluded = source_year_coverage$unresolved_source_years,
    source_year_cache_sha256 = source_year_coverage$source_year_cache_hashes,
    source_year_event_counts = source_year_coverage$source_year_counts,
    input_grouped_path = inputs$grouped_path,
    input_grouped_sha256 = sha256_file(inputs$grouped_path),
    input_calendar_path = inputs$calendar_path,
    input_calendar_sha256 = sha256_file(inputs$calendar_path),
    input_grouped_rows = nrow(inputs$clinical),
    event_date_field = "INGRE",
    event_date_window = dictionary_window,
    source_event_year_dx3_totals_2016_2022 = inputs$clinical[
      date >= as.Date(dictionary_window[[1L]]) & date <= as.Date(dictionary_window[[2L]]) & dx3 %in% model_diagnoses,
      .(n = sum(n)), by = .(date_year = format(date, "%Y"), dx3)
    ],
    origin_and_service_code_coverage = support_counts,
    model_years = 2016:2022,
    models = c("primary broad four-endpoint ED-origin minus outpatient-consultation-origin contrast",
               "strict ED-origin N20-only versus K80 contrast",
               "broad four-endpoint sensitivity restricted to recorded DIAS_ESTA 1-30 days"),
    unit_thresholds = list(primary = ">=30 events in each of four cells and >=400 E3-eligible days",
                           strict_N20_only_ED = ">=30 events in each of two cells and >=400 E3-eligible days",
                           los_sensitivity = "same four-cell thresholds"),
    exposure = list(field = "precomputed E3 from mx_updated_calendar.csv.gz", basis_df = 3L,
                    reference = "state-specific P95 versus median from the 2016-2022 all-origin eligible calendar days; recomputed for this window and shared across outcomes"),
    model = list(family = "quasi-Poisson", nuisance = "separate endpoint-specific day-of-week, month, and 12-df natural spline of time",
                 offset = "none", covariance = "calendar-date Bartlett HAC lag 14",
                 pooling = c("DerSimonian-Laird", "REML with Hartung-Knapp")),
    comparison = list(main_cohort_unchanged = TRUE,
                      original_2015_2022_comparator = "results/mx_origin/full_20261112T171633Z",
                      selection_policy = "report estimates in either direction, including null/adverse findings; no selection against known prior estimate"),
    interpretation_limits = c("PROCED=2 is an ED-origin proxy, not proof of disease onset, acuity, triage severity, or mechanism.",
                              "PROCED=1 denotes outpatient-consultation origin; it is not elective.",
                              "Codes 4 and 5 changed labels across annual dictionaries and are excluded from the two modeled origin groups.",
                              "TIPSERV is retained as a raw service code and is not an urgency/elective classifier.",
                              "LOS sensitivity denotes recorded DIAS_ESTA 1-30 days, not hospital certification."),
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
  cat("Mexico dictionary-supported origin analysis finished; results:", args$outdir, "\n")
  cat("Primary states fit:", sum(primary_rows$status == "ok"), "/", nrow(primary_rows), "\n")
  cat("Primary model failures/ineligible:", sum(primary_rows$status != "ok"), "\n")
  cat("Primary event-date window: 2016-2022; within-window state-specific references; no offset; HAC14\n")
}

if (sys.nframe() == 0L) main_dictionary_window()
