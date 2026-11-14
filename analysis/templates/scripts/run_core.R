#!/usr/bin/env Rscript
# Core unit-level repaired analyses. This script consumes canonical calendar
# files whose weather exposures were built before any outcome/period filtering.

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_arg)) stop("cannot resolve run_core.R location from --file", call. = FALSE)
.script_path <- normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/", mustWork = TRUE)
.script_dir <- dirname(.script_path)
.repair_dir <- normalizePath(file.path(.script_dir, ".."), winslash = "/", mustWork = TRUE)
.repo_dir <- normalizePath(file.path(.repair_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "repair_engine.R"), local = FALSE)

parse_args <- function(args) {
  opts <- list(stage = "primary", workers = 8L, save_fits = TRUE,
               brazil = file.path(.repair_dir, "data", "br_calendar.csv.gz"),
               mexico = file.path(.repair_dir, "data", "mx_calendar.csv.gz"),
               outdir = file.path(.repair_dir, "results"))
  i <- 1L
  while (i <= length(args)) {
    token <- args[i]
    if (identical(token, "--save-fits")) {
      opts$save_fits <- TRUE; i <- i + 1L; next
    }
    if (identical(token, "--no-save-fits")) {
      opts$save_fits <- FALSE; i <- i + 1L; next
    }
    if (!startsWith(token, "--")) stop(sprintf("unexpected argument '%s'", token), call. = FALSE)
    keyval <- sub("^--", "", token)
    if (grepl("=", keyval, fixed = TRUE)) {
      bits <- strsplit(keyval, "=", fixed = TRUE)[[1]]
      key <- bits[1]; value <- paste(bits[-1], collapse = "=")
      i <- i + 1L
    } else {
      key <- keyval
      if (i == length(args)) stop(sprintf("missing value for --%s", key), call. = FALSE)
      value <- args[i + 1L]; i <- i + 2L
    }
    if (!key %in% c("stage", "workers", "brazil", "mexico", "outdir")) {
      stop(sprintf("unknown option --%s", key), call. = FALSE)
    }
    opts[[key]] <- value
  }
  opts$stage <- match.arg(opts$stage, c("primary", "sensitivity"))
  opts$workers <- suppressWarnings(as.integer(opts$workers))
  if (length(opts$workers) != 1L || is.na(opts$workers) || opts$workers < 1L) {
    stop("--workers must be a positive integer", call. = FALSE)
  }
  opts$brazil <- normalizePath(opts$brazil, winslash = "/", mustWork = TRUE)
  opts$mexico <- normalizePath(opts$mexico, winslash = "/", mustWork = TRUE)
  opts$outdir <- path.expand(opts$outdir)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", opts$outdir)) opts$outdir <- file.path(.repo_dir, opts$outdir)
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
  opts$outdir <- normalizePath(opts$outdir, winslash = "/", mustWork = FALSE)
  opts
}

.opts <- parse_args(commandArgs(trailingOnly = TRUE))
.result_dir <- file.path(.opts$outdir, "core")
.fit_dir <- file.path(.opts$outdir, "fits")
dir.create(.result_dir, recursive = TRUE, showWarnings = FALSE)
if (.opts$stage == "primary") dir.create(.fit_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c("unit", "date", "E3", "E7", "E14", "E21", "Emax3", "tmean", "tmax",
                      "population", "y_K80", "y_K81", "y_urinary", "y_sum3")

read_calendar <- function(path, country) {
  message(sprintf("Reading canonical %s calendar", country))
  numeric_columns <- setdiff(required_columns, c("unit", "date"))
  classes <- c(unit = "character", date = "character",
               stats::setNames(rep("numeric", length(numeric_columns)), numeric_columns))
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  d <- utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("", "NA"), colClasses = classes)
  missing <- setdiff(required_columns, names(d))
  if (length(missing)) stop(sprintf("%s calendar is missing columns: %s", country,
                                    paste(missing, collapse = ", ")), call. = FALSE)
  d$unit <- as.character(d$unit)
  d$date <- as.Date(d$date)
  if (anyNA(d$unit) || any(!nzchar(d$unit))) stop(sprintf("%s calendar has empty unit IDs", country), call. = FALSE)
  if (anyNA(d$date)) stop(sprintf("%s calendar has invalid or missing dates", country), call. = FALSE)
  dup <- duplicated(data.frame(unit = d$unit, date = d$date))
  if (any(dup)) stop(sprintf("%s calendar has duplicate unit/date rows", country), call. = FALSE)
  d
}

.br <- read_calendar(.opts$brazil, "BR")
.mx <- read_calendar(.opts$mexico, "MX")
.run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
.input_hashes <- c(BR = unname(tools::md5sum(.opts$brazil)), MX = unname(tools::md5sum(.opts$mexico)))
.protocol_files <- c(file.path(.script_dir, "repair_engine.R"), .script_path,
                     file.path(.repair_dir, "admin", "analysis_plan.json"))
.protocol_files <- .protocol_files[file.exists(.protocol_files)]
.protocol_hashes <- tools::md5sum(.protocol_files)
.protocol_signature <- paste(c(R.version.string,
  paste(paste(basename(names(.protocol_hashes)), unname(.protocol_hashes), sep = "="), collapse = ";")),
  collapse = "|")

empty_estimates <- function() data.frame(country = character(), unit = character(), scenario = character(),
  exposure = character(), nuisance = character(), reference = character(), contrast = character(),
  covariance = character(), estimate = numeric(), variance = numeric(), se = numeric(), RR = numeric(),
  lo = numeric(), hi = numeric(), n_complete = integer(), events_first = numeric(),
  events_second = numeric(), events_third = numeric(), p95 = numeric(), median = numeric(),
  stringsAsFactors = FALSE)
empty_failures <- function() data.frame(country = character(), unit = character(), scenario = character(),
  phase = character(), reason = character(), n_complete = integer(), events = character(),
  stringsAsFactors = FALSE)
empty_audit <- function() data.frame(country = character(), unit = character(), scenario = character(),
  converged = logical(), boundary = logical(), rank = integer(), columns = integer(), n = integer(),
  events = numeric(), warnings = character(), formula = character(), stringsAsFactors = FALSE)
empty_cache_manifest <- function() data.frame(run_id = character(), country = character(), unit = character(),
  scenario = character(), cache_path = character(), input_md5 = character(), protocol_signature = character(),
  n_complete = integer(), p95 = numeric(), median = numeric(), stringsAsFactors = FALSE)
bind_records <- function(xs, empty) {
  xs <- Filter(function(x) !is.null(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(empty())
  do.call(rbind, xs)
}
atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop(sprintf("failed to replace output '%s'", path), call. = FALSE)
}

exposure_targets <- function(d, exposure, outcomes = NULL, offset = NULL) {
  if (!exposure %in% names(d)) stop(sprintf("missing exposure '%s'", exposure), call. = FALSE)
  if (is.null(outcomes)) {
    x <- d[[exposure]]
    x <- x[is.finite(x)]
  } else {
    eligible <- eligible_rows(d, outcomes, exposure, offset)
    x <- d[[exposure]][eligible$keep]
  }
  if (length(x) < 4L) stop(sprintf("fewer than four finite full-period values for '%s'", exposure), call. = FALSE)
  list(basis = new_temperature_basis(x, df = 3L), p95 = unname(stats::quantile(x, 0.95, type = 7)),
       median = stats::median(x))
}

period_filter <- function(d, range = NULL, exclude_years = integer()) {
  keep <- rep(TRUE, nrow(d))
  if (!is.null(range)) keep <- keep & d$date >= as.Date(range[1]) & d$date <= as.Date(range[2])
  if (length(exclude_years)) keep <- keep & !as.integer(format(d$date, "%Y")) %in% exclude_years
  d[keep, , drop = FALSE]
}

eligible_rows <- function(d, outcomes, exposure, offset) {
  good <- is.finite(d[[exposure]]) & !is.na(d$date)
  for (nm in outcomes) good <- good & is.finite(d[[nm]]) & d[[nm]] >= 0
  if (!is.null(offset)) good <- good & is.finite(d[[offset]]) & d[[offset]] > 0
  events <- if (any(good)) vapply(outcomes, function(nm) sum(d[[nm]][good]), numeric(1)) else
    setNames(rep(0, length(outcomes)), outcomes)
  list(keep = good, n = sum(good), events = events)
}

safe_unit <- function(unit) paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(unit)))), collapse = "")

run_one_unit <- function(unit_data, country, unit, scenario, exposure, outcomes, nuisance, offset,
                         range = NULL, exclude_years = integer(), reference = "full",
                         full_covariance = TRUE, save_fit = FALSE) {
  failure <- NULL; audit <- NULL; estimates <- empty_estimates()
  cache_path <- NA_character_
  cache_manifest <- empty_cache_manifest()
  n_complete <- NA_integer_; event_text <- ""
  tryCatch({
    all_targets <- exposure_targets(unit_data, exposure, outcomes = outcomes, offset = offset)
    dfit <- period_filter(unit_data, range = range, exclude_years = exclude_years)
    targets <- if (reference == "full") all_targets else {
      period_targets <- exposure_targets(dfit, exposure, outcomes = outcomes, offset = offset)
      list(basis = all_targets$basis, p95 = period_targets$p95, median = period_targets$median)
    }
    eligible <- eligible_rows(dfit, outcomes, exposure, offset)
    event_text <- paste(paste(names(eligible$events), eligible$events, sep = "="), collapse = ";")
    if (eligible$n < 250L || any(eligible$events < 40)) {
      n_complete <- as.integer(eligible$n)
      failure <- data.frame(country = country, unit = unit, scenario = scenario,
        phase = "eligibility", reason = sprintf("eligible complete dates=%d; required >=250 dates and >=40 events per endpoint", eligible$n),
        n_complete = n_complete, events = event_text, stringsAsFactors = FALSE)
    } else {
      n_complete <- as.integer(eligible$n)
      fit <- fit_stacked_contrast(dfit, outcomes = outcomes, exposure = exposure, date = "date",
        offset = offset, basis = targets$basis, p95 = targets$p95, median = targets$median,
        nuisance = nuisance, full_covariance = full_covariance)
      counts <- eligible$events
      if (length(counts) < 3L) counts <- c(counts, rep(NA_real_, 3L - length(counts)))
      rows <- list()
      for (contrast_name in names(fit$contrasts)) {
        item <- fit$contrasts[[contrast_name]]
        for (cov_name in names(item$coefficient_contrast$by_covariance)) {
          q <- item$coefficient_contrast$by_covariance[[cov_name]]
          rows[[length(rows) + 1L]] <- data.frame(country = country, unit = unit, scenario = scenario,
            exposure = exposure, nuisance = nuisance, reference = reference, contrast = contrast_name,
            covariance = cov_name, estimate = unname(q["estimate"]), variance = unname(q["variance"]),
            se = unname(q["se"]), RR = unname(q["RR"]), lo = unname(q["lo"]), hi = unname(q["hi"]),
            n_complete = n_complete, events_first = unname(counts[1]),
            events_second = unname(counts[2]), events_third = unname(counts[3]),
            p95 = fit$p95, median = fit$median, stringsAsFactors = FALSE)
        }
      }
      estimates <- do.call(rbind, rows)
      a <- fit$audit
      audit <- data.frame(country = country, unit = unit, scenario = scenario,
        converged = a$converged, boundary = a$boundary, rank = a$rank, columns = a$columns,
        n = a$n, events = a$events, warnings = paste(a$warnings, collapse = " | "),
        formula = a$formula, stringsAsFactors = FALSE)
      if (save_fit) {
        cache_path <- file.path(.fit_dir, paste0("primary_", country, "_unit_", safe_unit(unit), ".rds"))
        cache <- list(engine_fit = fit, unit = unit, country = country, scenario = scenario,
          exposure = exposure, outcomes = outcomes, nuisance = nuisance, offset = offset,
          p95 = fit$p95, median = fit$median, reference = reference, run_id = .run_id,
          input_md5 = unname(.input_hashes[[country]]), protocol_signature = .protocol_signature)
        tmp <- tempfile(pattern = paste0(basename(cache_path), "."), tmpdir = .fit_dir, fileext = ".rds")
        saveRDS(cache, tmp, compress = "gzip")
        if (file.exists(cache_path)) unlink(cache_path)
        if (!file.rename(tmp, cache_path)) stop(sprintf("failed to write fit cache for unit %s", unit), call. = FALSE)
        cache_manifest <- data.frame(run_id = .run_id, country = country, unit = unit,
          scenario = scenario, cache_path = normalizePath(cache_path, winslash = "/", mustWork = TRUE),
          input_md5 = unname(.input_hashes[[country]]), protocol_signature = .protocol_signature,
          n_complete = n_complete, p95 = fit$p95, median = fit$median, stringsAsFactors = FALSE)
      }
    }
  }, error = function(e) {
    estimates <<- empty_estimates()
    audit <<- NULL
    cache_manifest <<- empty_cache_manifest()
    if (is.null(failure)) failure <<- data.frame(country = country, unit = unit, scenario = scenario,
      phase = "fit", reason = conditionMessage(e), n_complete = n_complete, events = event_text, stringsAsFactors = FALSE)
  })
  list(estimates = estimates, failure = failure, audit = audit, cache_path = cache_path,
       cache_manifest = cache_manifest)
}

run_scenario <- function(calendar, country, scenario, exposure, outcomes, nuisance, offset,
                         range = NULL, exclude_years = integer(), reference = "full",
                         full_covariance = TRUE, save_fit = FALSE, workers = 1L,
                         progress = TRUE) {
  units <- sort(unique(calendar$unit))
  if (progress) message(sprintf("%s: %s (%d units)", country, scenario, length(units)))
  worker <- function(unit) {
    z <- calendar[calendar$unit == unit, , drop = FALSE]
    run_one_unit(z, country, unit, scenario, exposure, outcomes, nuisance, offset,
      range = range, exclude_years = exclude_years, reference = reference,
      full_covariance = full_covariance, save_fit = save_fit)
  }
  if (.Platform$OS.type == "windows" || workers <= 1L || length(units) <= 1L) {
    results <- lapply(units, worker)
  } else {
    results <- parallel::mclapply(units, worker, mc.cores = min(as.integer(workers), length(units)),
                                  mc.preschedule = FALSE, mc.set.seed = FALSE)
  }
  estimates <- bind_records(lapply(results, `[[`, "estimates"), empty_estimates)
  failures <- bind_records(lapply(results, `[[`, "failure"), empty_failures)
  audits <- bind_records(lapply(results, `[[`, "audit"), empty_audit)
  caches <- bind_records(lapply(results, `[[`, "cache_manifest"), empty_cache_manifest)
  list(estimates = estimates, failures = failures, audits = audits, caches = caches, units = length(units))
}

pool_estimates <- function(estimates, country, scenario) {
  if (!nrow(estimates)) {
    return(data.frame(country = character(), scenario = character(), contrast = character(),
      covariance = character(), method = character(), k = integer(), estimate = numeric(), se = numeric(),
      RR = numeric(), lo = numeric(), hi = numeric(), I2 = numeric(), tau2 = numeric(),
      concordance = numeric(), failure = character(), stringsAsFactors = FALSE))
  }
  groups <- split(estimates, interaction(estimates$contrast, estimates$covariance,
                                         drop = TRUE, lex.order = TRUE))
  out <- list()
  for (g in groups) for (method in c("DL", "REML_HK")) {
    p <- pool_meta(g$estimate, g$variance, method = method)
    out[[length(out) + 1L]] <- data.frame(country = country, scenario = scenario,
      contrast = g$contrast[1], covariance = g$covariance[1], method = method,
      k = p$k, estimate = p$estimate, se = p$se, RR = p$RR,
      lo = p$RR95CI[1], hi = p$RR95CI[2], I2 = p$I2, tau2 = p$tau2,
      concordance = p$concordance, failure = if (is.null(p$failure)) "" else p$failure,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

write_scenario <- function(result, country, scenario) {
  stem <- gsub("[^A-Za-z0-9_.-]", "_", scenario)
  atomic_write_csv(result$estimates, file.path(.result_dir, paste0(stem, "_", country, "_per_unit.csv")))
  atomic_write_csv(pool_estimates(result$estimates, country, scenario),
                   file.path(.result_dir, paste0(stem, "_", country, "_pooled.csv")))
  atomic_write_csv(result$failures, file.path(.result_dir, paste0(stem, "_", country, "_failure_audit.csv")))
  atomic_write_csv(result$audits, file.path(.result_dir, paste0(stem, "_", country, "_fit_audit.csv")))
}

append_cache_manifest <- function(rows) {
  path <- file.path(.fit_dir, "primary_cache_manifest.csv")
  current <- if (file.exists(path)) utils::read.csv(path, stringsAsFactors = FALSE) else empty_cache_manifest()
  if (!"run_id" %in% names(current) || !nrow(current)) current <- empty_cache_manifest()
  if (nrow(rows)) current <- rbind(current, rows)
  atomic_write_csv(current, path)
}

run_primary <- function(calendar, country) {
  scenario <- "primary_threeoutcome_E3_separate"
  result <- run_scenario(calendar, country, scenario, exposure = "E3",
    outcomes = c("y_K80", "y_K81", "y_urinary"), nuisance = "separate",
    offset = if (country == "BR") "population" else NULL,
    full_covariance = TRUE, save_fit = isTRUE(.opts$save_fits), workers = .opts$workers)
  write_scenario(result, country, scenario)
  if (isTRUE(.opts$save_fits)) append_cache_manifest(result$caches)
  result
}

run_two_outcome_sensitivity <- function(calendar, country, scenario, exposure = "E3",
    nuisance = "separate", range = NULL, exclude_years = integer(), reference = "full",
    full_covariance = FALSE) {
  result <- run_scenario(calendar, country, scenario, exposure = exposure,
    outcomes = c("y_K80", "y_urinary"), nuisance = nuisance,
    offset = if (country == "BR") "population" else NULL,
    range = range, exclude_years = exclude_years, reference = reference,
    full_covariance = full_covariance, workers = .opts$workers)
  write_scenario(result, country, scenario)
  result
}

run_primary_stage <- function() {
  atomic_write_csv(empty_cache_manifest(), file.path(.fit_dir, "primary_cache_manifest.csv"))
  br <- run_primary(.br, "BR")
  mx <- run_primary(.mx, "MX")
  invisible(list(BR = br, MX = mx))
}

run_mexico_leave_one_state_out <- function() {
  path <- file.path(.result_dir, "primary_threeoutcome_E3_separate_MX_per_unit.csv")
  out_path <- file.path(.result_dir, "sensitivity_MX_leave_one_state_out.csv")
  failure_path <- file.path(.result_dir, "sensitivity_MX_leave_one_state_out_failure_audit.csv")
  if (!file.exists(path)) {
    atomic_write_csv(data.frame(reason = "primary MX per-unit results are unavailable; run --stage primary first"), failure_path)
    return(invisible(NULL))
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$covariance == "hac14" & d$contrast == "y_urinary_vs_y_K80", , drop = FALSE]
  states <- sort(unique(d$unit))
  if (length(states) < 2L) {
    atomic_write_csv(data.frame(reason = "fewer than two successful MX states in the primary result"), failure_path)
    return(invisible(NULL))
  }
  rows <- list()
  for (state in states) {
    keep <- d[d$unit != state, , drop = FALSE]
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(keep$estimate, keep$variance, method = method)
      rows[[length(rows) + 1L]] <- data.frame(country = "MX", scenario = "sensitivity_MX_leave_one_state_out",
        omitted_unit = state, contrast = "y_urinary_vs_y_K80", covariance = "hac14", method = method,
        k = p$k, estimate = p$estimate, se = p$se, RR = p$RR, lo = p$RR95CI[1], hi = p$RR95CI[2],
        I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
        failure = if (is.null(p$failure)) "" else p$failure, stringsAsFactors = FALSE)
    }
  }
  atomic_write_csv(do.call(rbind, rows), out_path)
  if (file.exists(failure_path)) unlink(failure_path)
}

run_sensitivity_stage <- function() {
  # Separate disease-specific year-month adjustment for the primary contrast.
  run_two_outcome_sensitivity(.br, "BR", "sensitivity_primary_yearmonth", nuisance = "yearmonth")
  run_two_outcome_sensitivity(.mx, "MX", "sensitivity_primary_yearmonth", nuisance = "yearmonth")

  # Lag 0-2, 0-6, 0-13 and 0-20-day precomputed weather averages.
  for (exposure in c("E3", "E7", "E14", "E21")) {
    run_two_outcome_sensitivity(.br, "BR", paste0("sensitivity_lag_", exposure), exposure = exposure)
    run_two_outcome_sensitivity(.mx, "MX", paste0("sensitivity_lag_", exposure), exposure = exposure)
  }

  # Fixed full-period temperature target/basis with period-specific target sensitivity.
  periods <- list(
    BR = list(c("2006-01-01", "2016-12-31"), c("2017-01-01", "2020-12-31")),
    MX = list(c("2015-01-01", "2018-12-31"), c("2020-01-01", "2022-12-31"), c("2022-01-01", "2022-12-31")))
  for (country in c("BR", "MX")) {
    calendar <- if (country == "BR") .br else .mx
    for (r in periods[[country]]) {
      label <- paste0("sensitivity_temporal_", gsub("-", "", r[1]), "_", gsub("-", "", r[2]))
      run_two_outcome_sensitivity(calendar, country, paste0(label, "_fixed_full_ref"),
        range = r, reference = "full")
      run_two_outcome_sensitivity(calendar, country, paste0(label, "_period_ref"),
        range = r, reference = "period")
    }
  }

  # Original Brazil eight-cell grid using the corrected lagged exposures.
  for (exposure in c("E3", "Emax3")) for (nuisance in c("shared", "yearmonth")) {
    for (period in list(full = c("2006-01-01", "2020-12-31"),
                        sub = c("2008-01-01", "2015-12-31"))) {
      period_name <- if (identical(period, c("2006-01-01", "2020-12-31"))) "full" else "sub"
      scenario <- paste("sensitivity_original8", exposure, nuisance, period_name, sep = "_")
      run_two_outcome_sensitivity(.br, "BR", scenario, exposure = exposure,
        nuisance = nuisance, range = period, reference = "full")
    }
  }

  # Exclude Mexico's 2020-2021 observations and re-estimate each state's model.
  result <- run_scenario(.mx, "MX", "sensitivity_MX_exclude_2020_2021", exposure = "E3",
    outcomes = c("y_K80", "y_urinary"), nuisance = "separate", offset = NULL,
    exclude_years = c(2020L, 2021L), full_covariance = FALSE, workers = .opts$workers)
  write_scenario(result, "MX", "sensitivity_MX_exclude_2020_2021")
  run_mexico_leave_one_state_out()
}

if (.opts$stage == "primary") run_primary_stage() else run_sensitivity_stage()
message(sprintf("Core %s stage completed; outputs: %s", .opts$stage, .result_dir))
