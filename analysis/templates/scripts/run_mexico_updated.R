#!/usr/bin/env Rscript

# Mexico complete-source cohort repair and prespecified model runs.
# This script writes only data/mx_updated_* and results/mx_updated/.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(root = NULL, stage = "all", workers = 8L)
for (a in args) {
  if (grepl("^--root=", a)) opt$root <- sub("^--root=", "", a)
  else if (grepl("^--stage=", a)) opt$stage <- sub("^--stage=", "", a)
  else if (grepl("^--workers=", a)) opt$workers <- as.integer(sub("^--workers=", "", a))
  else stop(sprintf("unknown argument: %s", a), call. = FALSE)
}
if (!opt$stage %in% c("build", "smoke", "models", "all", "repair-exports")) stop("--stage must be build, smoke, models, all, or repair-exports")
if (is.na(opt$workers) || opt$workers < 1L) stop("--workers must be a positive integer")
if (is.null(opt$root)) {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(file_arg)) stop("cannot infer project root; supply --root=...")
  opt$root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                            winslash = "/", mustWork = TRUE)
}
root <- normalizePath(opt$root, winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data")
script_dir <- file.path(root, "scripts")
result_dir <- file.path(root, "results", "mx_updated")
fit_dir <- file.path(result_dir, "fits")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

# Require splines explicitly before loading or predicting any cached ns object.
if (!requireNamespace("splines", quietly = TRUE)) stop("required namespace 'splines' is unavailable")
if (!requireNamespace("data.table", quietly = TRUE)) stop("required package 'data.table' is unavailable")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("required package 'jsonlite' is unavailable")
if (!requireNamespace("metafor", quietly = TRUE)) stop("required package 'metafor' is unavailable")
data.table <- asNamespace("data.table")
fread <- data.table$fread
fwrite <- data.table$fwrite
rbindlist <- data.table$rbindlist
copy <- data.table$copy
setDT <- data.table$setDT
setorder <- data.table$setorder
setnames <- data.table$setnames
setkeyv <- data.table$setkeyv

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
engine_path <- file.path(script_dir, "repair_engine.R")
if (!file.exists(engine_path)) stop(sprintf("repair engine missing: %s", engine_path))
source(engine_path, local = FALSE)

in_path <- function(name) file.path(data_dir, name)
sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- suppressWarnings(system2("sha256sum", shQuote(normalizePath(path, winslash = "/")),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("sha256sum failed for %s: %s", path, paste(out, collapse = " ")))
  sub("\\s+.*$", "", out[1L])
}
read_table <- function(path) {
  if (!file.exists(path)) stop(sprintf("required input file is missing: %s", path), call. = FALSE)
  cmd <- if (grepl("\\.gz$", path, ignore.case = TRUE)) sprintf("gzip -dc -- %s", shQuote(path)) else NULL
  h <- if (is.null(cmd)) fread(path, keepLeadingZeros = TRUE, showProgress = FALSE,
                               check.names = FALSE, na.strings = c("", "NA")) else
    fread(cmd = cmd, keepLeadingZeros = TRUE, showProgress = FALSE,
          check.names = FALSE, na.strings = c("", "NA"))
  setDT(h)
  for (nm in intersect(c("unit", "state", "source_role", "date"), names(h))) h[[nm]] <- as.character(h[[nm]])
  h
}
atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  fwrite(x, tmp, na = "", quote = "auto")
  if (!file.rename(tmp, path)) stop(sprintf("atomic CSV replacement failed: %s", path), call. = FALSE)
  invisible(path)
}
atomic_gz_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_gz <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".gz")
  on.exit(unlink(tmp_gz), add = TRUE)
  con <- gzfile(tmp_gz, open = "wt")
  tryCatch(write.csv(as.data.frame(x), con, row.names = FALSE, na = ""),
           finally = close(con))
  if (!file.exists(tmp_gz) || file.info(tmp_gz)$size <= 0) stop(sprintf("gzip failed writing %s", path), call. = FALSE)
  if (!file.rename(tmp_gz, path)) stop(sprintf("atomic gzip replacement failed: %s", path), call. = FALSE)
  invisible(path)
}
atomic_json <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = NA)
  if (!file.rename(tmp, path)) stop(sprintf("atomic JSON replacement failed: %s", path), call. = FALSE)
  invisible(path)
}
key_duplicates <- function(d, cols) anyDuplicated(d[, ..cols]) > 0L
check_columns <- function(d, cols, label) {
  missing <- setdiff(cols, names(d))
  if (length(missing)) stop(sprintf("%s lacks required columns: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
}
urinary_parts <- c("n20_only", "n21_only", "n22_other_disease", "n23_unspecified_renal_colic")
health_cols <- c("k80", "k81", urinary_parts)
outcome_cols <- c("y_K80", "y_K81", "y_urinary")
all_outcomes <- outcome_cols
outcomes_from_health <- function(d) {
  for (nm in c("k80", "k81", urinary_parts)) if (!is.numeric(d[[nm]])) d[[nm]] <- as.numeric(d[[nm]])
  d[, `:=`(y_K80 = as.numeric(k80), y_K81 = as.numeric(k81),
           y_urinary = as.numeric(n20_only) + as.numeric(n21_only) +
             as.numeric(n22_other_disease) + as.numeric(n23_unspecified_renal_colic))]
  d[, y_sum3 := y_K80 + y_K81 + y_urinary]
  d
}
add_calendar_fields <- function(d) {
  d[, year := as.integer(format(date, "%Y"))]
  d[, month := as.integer(format(date, "%m"))]
  d[, dow := as.integer(format(date, "%u"))]
  d[, country := "MX"]
  d
}
recompute_E3 <- function(d) {
  setorder(d, unit, date)
  d[, E3_check := calendar_ma(tmean, date, k = 3L), by = unit]
  if ("E3" %in% names(d)) {
    ok <- is.finite(d$E3) & is.finite(d$E3_check)
    bad_na <- is.na(d$E3) != is.na(d$E3_check)
    if (any(bad_na) || (any(ok) && max(abs(d$E3[ok] - d$E3_check[ok])) > 1e-8)) {
      stop("legacy E3 does not match complete-calendar lag0-2 exposure", call. = FALSE)
    }
    d[, E3 := E3_check]
  } else d[, E3 := E3_check]
  d[, E3_check := NULL]
  d
}
read_inputs <- function() {
  paths <- list(
    old_calendar = in_path("mx_calendar.csv.gz"),
    complete_source = in_path("mx_complete_source_state_day_2015_2022.csv"),
    source_2019 = in_path("mx_complete_source_2019_state_day.csv"),
    clinical_roles = in_path("mx_clinical_state_day.csv"),
    weather_2023 = in_path("mx_weather_2023.csv"),
    plan = file.path(root, "admin", "analysis_plan.json"),
    amendment = file.path(root, "admin", "cohort_amendment_2019.json"),
    engine = engine_path,
    script = file.path(script_dir, "run_mexico_updated.R")
  )
  for (p in paths) if (!file.exists(p)) stop(sprintf("required input is missing: %s", p), call. = FALSE)
  hashes <- vapply(paths, sha256, character(1))
  old <- read_table(paths$old_calendar)
  check_columns(old, c("date", "unit", "tmean", "E3", "y_K80", "y_K81", "y_urinary", "y_sum3"), "old Mexico calendar")
  old[, date := as.Date(date)]
  if (anyNA(old$date) || anyNA(old$unit) || key_duplicates(old, c("unit", "date"))) stop("old Mexico calendar has invalid/duplicate unit dates")
  if (nrow(old) != 93504L || length(unique(old$unit)) != 32L || min(old$date) != as.Date("2015-01-01") || max(old$date) != as.Date("2022-12-31")) {
    stop("old weather calendar does not match the expected complete 32-state 2015-2022 grid")
  }
  old[, unit := as.character(unit)]
  old <- recompute_E3(old)

  health <- read_table(paths$complete_source)
  check_columns(health, c("source_role", "state", "date", health_cols), "complete-source 2015-2022 grid")
  health[, `:=`(state = as.character(state), date = as.Date(date))]
  if (nrow(health) != 93504L || anyNA(health$date) || key_duplicates(health, c("state", "date"))) {
    stop("complete-source 2015-2022 health grid is incomplete or duplicated")
  }
  if (length(unique(health$source_role)) != 1L || unique(health$source_role) != "preferred_complete_source_2015_2022") {
    stop("complete-source health role does not identify the preferred all-source cohort")
  }
  if (min(health$date) != as.Date("2015-01-01") || max(health$date) != as.Date("2022-12-31") || length(unique(health$state)) != 32L) {
    stop("complete-source health dates or states do not span the expected grid")
  }
  health <- outcomes_from_health(health)

  # Check that the 2019 extract and recovered role are exactly the same cases as
  # the 2019 slice in the preferred complete-source series, without addition.
  h19 <- health[format(date, "%Y") == "2019"]
  src19 <- read_table(paths$source_2019)
  check_columns(src19, c("state", "date", health_cols), "controlled all-source 2019 extract")
  src19[, `:=`(state = as.character(state), date = as.Date(date))]
  role <- read_table(paths$clinical_roles)
  check_columns(role, c("source_role", "state", "date", health_cols), "clinical source-role state-day table")
  role[, `:=`(state = as.character(state), date = as.Date(date))]
  role19 <- role[source_role == "recovered_additional_2019" & format(date, "%Y") == "2019"]
  check_columns(role, c("state", "date", health_cols), "clinical source-role state-day table")
  if (nrow(h19) != 11680L || nrow(src19) != 11680L || nrow(role19) != 11680L ||
      key_duplicates(src19, c("state", "date")) || key_duplicates(role19, c("state", "date"))) {
    stop("2019 complete-source coverage is not exactly one 32-state daily year")
  }
  comparison <- function(x) {
    z <- x[, c("state", "date", health_cols), with = FALSE]
    setorder(z, state, date)
    z
  }
  if (!identical(comparison(h19), comparison(src19))) {
    stop("2019 complete-source grid disagrees with the all-source 2019 extract")
  }
  # The recovered 2019 file is one source year; the preferred complete-source
  # table also includes later source-year files whose INGRE dates fall in 2019.
  # Require exact state-day coverage and non-negative additions, not equality.
  role_keys <- role19[, .(state, date)]
  full_keys <- h19[, .(state, date)]
  setorder(role_keys, state, date); setorder(full_keys, state, date)
  if (!identical(role_keys, full_keys)) stop("recovered 2019 role does not cover the same 32-state event-date grid")
  role_compare <- merge(comparison(h19), comparison(role19), by = c("state", "date"),
                        all = TRUE, sort = TRUE, suffixes = c("_preferred", "_recovered"))
  role_deltas <- lapply(health_cols, function(nm) {
    delta <- role_compare[[paste0(nm, "_preferred")]] - role_compare[[paste0(nm, "_recovered")]]
    if (any(!is.finite(delta)) || any(delta < 0)) {
      stop(sprintf("preferred all-source 2019 %s has missing/negative cross-source additions", nm))
    }
    data.frame(field = nm, preferred_total = sum(role_compare[[paste0(nm, "_preferred")]]),
      recovered_2019_file_total = sum(role_compare[[paste0(nm, "_recovered")]]),
      later_source_year_additions = sum(delta), dates_with_additions = sum(delta > 0), stringsAsFactors = FALSE)
  })
  role_delta_audit <- rbindlist(role_deltas)
  old19 <- old[format(date, "%Y") == "2019"]
  if (any(vapply(outcome_cols, function(nm) any(is.finite(old19[[nm]])), logical(1)))) {
    stop("old baseline already has finite 2019 outcomes; refusing possible recovered-year double count")
  }

  # Strong monotonic source-completeness check: every original baseline year's
  # annual outcome totals must be retained or increased by the complete-source
  # aggregation. 2019 is deliberately excluded because it was previously absent.
  original_years <- c(2015:2018, 2020:2022)
  checks <- list()
  for (yy in original_years) for (nm in outcome_cols) {
    old_total <- sum(old[year == yy][[nm]], na.rm = TRUE)
    new_total <- sum(health[as.integer(format(date, "%Y")) == yy][[nm]], na.rm = TRUE)
    checks[[length(checks) + 1L]] <- data.frame(year = yy, outcome = nm, old_total = old_total,
      new_total = new_total, difference = new_total - old_total, nondecreasing = new_total >= old_total)
  }
  count_check <- rbindlist(checks)
  if (any(!count_check$nondecreasing)) stop("complete-source annual health totals fell below original counts")

  new_weather <- copy(old)
  new_weather[, c(outcome_cols, "y_sum3") := NULL]
  new <- merge(new_weather, health, by.x = c("unit", "date"), by.y = c("state", "date"),
               all.x = TRUE, sort = FALSE, suffixes = c("", "_health"))
  if (nrow(new) != nrow(old) || anyNA(new$y_K80) || anyNA(new$y_K81) || anyNA(new$y_urinary) ||
      key_duplicates(new, c("unit", "date"))) stop("complete-source health merge failed to preserve the old weather grid")
  # Health rows are all-source counts. The 2019-year role marker exposes its origin.
  if ("source_role" %in% names(new)) setnames(new, "source_role", "health_source_role")
  new <- recompute_E3(new)
  setorder(new, unit, date)
  add_calendar_fields(new)

  # Controlled-fill comparator: preserve the frozen baseline on non-2019 dates;
  # replace the entire absent 2019 year exactly once with the all-source extract.
  controlled <- copy(old)
  h19_out <- h19[, .(unit = state, date, y_K80, y_K81, y_urinary, y_sum3,
                     k80, k81, n20_only, n21_only, n22_other_disease,
                     n23_unspecified_renal_colic, source_role)]
  controlled <- merge(controlled, h19_out, by = c("unit", "date"), all.x = TRUE,
                      sort = FALSE, suffixes = c(".old", ""))
  is19 <- format(controlled$date, "%Y") == "2019"
  for (nm in c(outcome_cols, "y_sum3")) controlled[, (nm) := ifelse(is19, get(nm), get(paste0(nm, ".old")))]
  for (nm in c(outcome_cols, "y_sum3")) controlled[, (paste0(nm, ".old")) := NULL]
  controlled <- recompute_E3(controlled)
  controlled <- add_calendar_fields(controlled)

  updated_path <- in_path("mx_updated_calendar.csv.gz")
  controlled_path <- in_path("mx_updated_controlled_2019_calendar.csv.gz")
  atomic_gz_csv(new, updated_path)
  atomic_gz_csv(controlled, controlled_path)
  if (!identical(unname(hashes["old_calendar"]), sha256(paths$old_calendar))) stop("old mx_calendar changed during canonical rebuild")

  audit <- list(
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    cohort_amendment = jsonlite::fromJSON(paths$amendment, simplifyVector = FALSE),
    cohort = "Mexico preferred complete-source 2015-2022; event dates include full recovered 2019 and cross-year admission dates by INGRE",
    baseline_comparator = "original mx_calendar.csv.gz read-only; original 2015-2018 and 2020-2022 health dates retained",
    rows = nrow(new), units = length(unique(new$unit)), first_date = as.character(min(new$date)),
    last_date = as.character(max(new$date)), dates_per_unit = nrow(new) / length(unique(new$unit)),
    E3_rule = "complete actual-calendar lag0-2; calculated before outcome/date filters; first two dates per state are NA",
    source_role = unique(new$health_source_role), year_totals_nondecreasing = count_check,
    old_2019_outcomes_missing = TRUE, complete_source_2019_rows = nrow(h19),
    recovered_2019_is_superseded_not_added = TRUE, later_source_year_2019_additions = role_delta_audit,
    old_2019_outcomes_double_counted = FALSE,
    health_totals_2015_2022 = lapply(c(health_cols, "y_urinary", "y_sum3"), function(nm) {
      if (nm %in% names(health)) sum(health[[nm]], na.rm = TRUE) else NA_real_
    }) |> setNames(c(health_cols, "y_urinary", "y_sum3")),
    source_hash_sha256 = as.list(hashes),
    output_hash_sha256 = list(mx_updated_calendar = sha256(updated_path),
                              mx_updated_controlled_2019_calendar = sha256(controlled_path)),
    output_rows = list(mx_updated_calendar = nrow(new), mx_updated_controlled_2019_calendar = nrow(controlled))
  )
  atomic_json(audit, file.path(data_dir, "mx_updated_calendar_manifest.json"))
  list(old = old, health = health, new = new, controlled = controlled,
       clinical_roles = role, hashes = hashes, count_check = count_check,
       audit = audit, old_calendar_path = paths$old_calendar)
}

make_targets <- function(d, outcomes) {
  check_columns(d, c("E3", outcomes, "date", "unit"), "target calendar")
  good <- is.finite(d$E3) & !is.na(d$date)
  for (nm in outcomes) good <- good & is.finite(d[[nm]]) & d[[nm]] >= 0
  x <- d$E3[good]
  if (sum(good) < 4L || length(unique(x)) < 2L) stop("insufficient complete full-period paired dates for temperature targets")
  list(basis = new_temperature_basis(x, df = 3L),
       p95 = unname(stats::quantile(x, 0.95, type = 7)), median = stats::median(x),
       n_complete = sum(good), events = vapply(outcomes, function(nm) sum(d[[nm]][good]), numeric(1)),
       date_start = as.character(min(d$date[good])), date_end = as.character(max(d$date[good])))
}
safe_unit <- function(unit) paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(unit)))), collapse = "")
scenario_defs <- list()
scenario_summaries <- list()
input_hashes <- NULL
run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")

empty_per_unit <- function() data.frame(country = character(), unit = character(), scenario = character(),
  exposure = character(), nuisance = character(), reference = character(), contrast = character(),
  covariance = character(), estimate = numeric(), variance = numeric(), se = numeric(), RR = numeric(),
  lo = numeric(), hi = numeric(), n_complete = integer(), events_first = numeric(),
  events_second = numeric(), events_third = numeric(), p95 = numeric(), median = numeric(),
  stringsAsFactors = FALSE)
empty_failure <- function() data.frame(country = character(), unit = character(), scenario = character(),
  phase = character(), reason = character(), n_complete = integer(), events = character(), stringsAsFactors = FALSE)
empty_fit_audit <- function() data.frame(country = character(), unit = character(), scenario = character(),
  converged = logical(), boundary = logical(), rank = integer(), columns = integer(), n = integer(),
  events = numeric(), warnings = character(), formula = character(), stringsAsFactors = FALSE)
empty_cache_manifest <- function() data.frame(run_id = character(), country = character(), unit = character(),
  scenario = character(), cache_path = character(), cache_sha256 = character(), input_sha256 = character(),
  engine_sha256 = character(), n_complete = integer(), p95 = numeric(), median = numeric(), stringsAsFactors = FALSE)
bind_rows <- function(xs, empty) {
  xs <- Filter(function(x) !is.null(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(empty())
  rbindlist(xs, fill = TRUE, use.names = TRUE)
}
calc_targets_by_unit <- function(calendar, outcomes, period_ref = NULL) {
  units <- sort(unique(calendar$unit))
  out <- vector("list", length(units)); names(out) <- units
  for (unit in units) {
    idx <- which(calendar[["unit"]] == unit)
    refdata <- if (is.null(period_ref)) calendar[idx] else {
      ref_idx <- which(period_ref[["unit"]] == unit)
      period_ref[ref_idx]
    }
    out[[unit]] <- tryCatch(make_targets(refdata, outcomes), error = function(e) e)
  }
  out
}
target_row <- function(d, outcomes) {
  good <- is.finite(d$E3) & !is.na(d$date)
  for (nm in outcomes) good <- good & is.finite(d[[nm]]) & d[[nm]] >= 0
  ev <- vapply(outcomes, function(nm) sum(d[[nm]][good]), numeric(1))
  list(keep = good, n = sum(good), events = ev)
}
pool_scenario <- function(estimates, scenario) {
  if (!nrow(estimates)) return(data.frame(country = character(), scenario = character(),
    contrast = character(), covariance = character(), method = character(), k = integer(),
    estimate = numeric(), se = numeric(), RR = numeric(), lo = numeric(), hi = numeric(),
    I2 = numeric(), tau2 = numeric(), concordance = numeric(), failure = character(),
    stringsAsFactors = FALSE))
  groups <- split(seq_len(nrow(estimates)), interaction(estimates$contrast, estimates$covariance,
                                         drop = TRUE, lex.order = TRUE))
  out <- list()
  for (idx in groups) {
    g <- estimates[idx, , drop = FALSE]
    if (anyDuplicated(g$unit) || nrow(g) != length(unique(g$unit))) {
      stop(sprintf("%s pooling group %s/%s has duplicate unit rows", scenario, g$contrast[1], g$covariance[1]))
    }
    for (method in c("DL", "REML_HK")) {
    p <- pool_meta(g$estimate, g$variance, method = method)
    if (p$k != nrow(g) || p$k != length(unique(g$unit))) {
      stop(sprintf("%s pooling k does not equal unique units for %s/%s/%s", scenario,
                   g$contrast[1], g$covariance[1], method))
    }
    out[[length(out) + 1L]] <- data.frame(country = "MX", scenario = scenario,
      contrast = g$contrast[1], covariance = g$covariance[1], method = method,
      k = p$k, estimate = p$estimate, se = p$se, RR = p$RR,
      lo = p$RR95CI[1], hi = p$RR95CI[2], I2 = p$I2, tau2 = p$tau2,
      concordance = p$concordance, failure = if (is.null(p$failure)) "" else p$failure,
      stringsAsFactors = FALSE)
    }
  }
  rbindlist(out, fill = TRUE)
}
save_cache <- function(fit, unit, scenario, target, n_complete, inputs) {
  cache_path <- file.path(fit_dir, paste0(scenario, "_MX_unit_", safe_unit(unit), ".rds"))
  cache <- list(format_version = 1L, fit = fit, unit = unit, country = "MX", scenario = scenario,
    exposure = "E3", outcomes = fit$outcomes, nuisance = fit$nuisance, offset = NULL,
    target_basis = target$basis, p95 = target$p95, median = target$median,
    n_complete = n_complete, run_id = run_id, input_hashes_sha256 = inputs,
    engine_hash_sha256 = sha256(engine_path), R = R.version.string,
    cache_scope = "Independent Mexico updated-cohort cache; does not share or reset the frozen-analysis fit manifest")
  tmp <- tempfile(pattern = paste0(".", basename(cache_path), "."), tmpdir = fit_dir, fileext = ".rds")
  saveRDS(cache, tmp, compress = "gzip")
  if (!file.rename(tmp, cache_path)) stop(sprintf("failed to write full fit cache %s", cache_path))
  data.frame(run_id = run_id, country = "MX", unit = unit, scenario = scenario,
    cache_path = normalizePath(cache_path, winslash = "/", mustWork = TRUE), cache_sha256 = sha256(cache_path),
    input_sha256 = inputs["canonical"], engine_sha256 = sha256(engine_path), n_complete = n_complete,
    p95 = target$p95, median = target$median, stringsAsFactors = FALSE)
}
scenario_cache_manifest <- if (file.exists(file.path(fit_dir, "fit_cache_manifest.csv")))
  read.csv(file.path(fit_dir, "fit_cache_manifest.csv"), stringsAsFactors = FALSE) else empty_cache_manifest()
append_cache_manifest <- function(rows) {
  if (!nrow(rows)) return(invisible(NULL))
  scenario_cache_manifest <<- rbind(scenario_cache_manifest, as.data.frame(rows))
  atomic_csv(scenario_cache_manifest, file.path(fit_dir, "fit_cache_manifest.csv"))
}
write_scenario <- function(result, scenario, metadata) {
  stem <- gsub("[^A-Za-z0-9_.-]", "_", scenario)
  atomic_csv(result$estimates, file.path(result_dir, paste0(stem, "_MX_per_unit.csv")))
  atomic_csv(pool_scenario(result$estimates, scenario), file.path(result_dir, paste0(stem, "_MX_pooled.csv")))
  atomic_csv(result$failures, file.path(result_dir, paste0(stem, "_MX_failure_audit.csv")))
  atomic_csv(result$audits, file.path(result_dir, paste0(stem, "_MX_fit_audit.csv")))
  summary <- list(scenario = scenario, run_id = run_id, units_expected = result$units,
    units_fit = length(unique(result$estimates$unit)), units_failed = nrow(result$failures),
    eligible_unit_k = length(unique(result$estimates$unit)), contrasts = unique(result$estimates$contrast),
    covariance_estimators = unique(result$estimates$covariance), metadata = metadata,
    output_sha256 = list(per_unit = sha256(file.path(result_dir, paste0(stem, "_MX_per_unit.csv"))),
      pooled = sha256(file.path(result_dir, paste0(stem, "_MX_pooled.csv"))),
      failures = sha256(file.path(result_dir, paste0(stem, "_MX_failure_audit.csv"))),
      fit_audit = sha256(file.path(result_dir, paste0(stem, "_MX_fit_audit.csv")))))
  scenario_summaries[[scenario]] <<- summary
  update_run_manifest()
  invisible(summary)
}
update_run_manifest <- function() {
  atomic_json(list(run_id = run_id, updated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    R = R.version.string, workers = opt$workers, OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
    source_input_hashes_sha256 = input_hashes, scenarios = scenario_summaries,
    fit_cache_manifest = "fits/fit_cache_manifest.csv"), file.path(result_dir, "run_manifest.json"))
}
run_scenario <- function(calendar, scenario, outcomes = all_outcomes, target_calendar = calendar,
                         range = NULL, exclude_years = integer(), period_specific_ref = FALSE,
                         save_fits = TRUE) {
  units <- sort(unique(calendar$unit))
  targets <- calc_targets_by_unit(calendar, outcomes, period_ref = target_calendar)
  if (!identical(sort(names(targets)), units)) stop(sprintf("%s target states differ from model states", scenario))
  worker <- function(unit) {
    tryCatch({
      idx <- which(calendar[["unit"]] == unit)
      d <- copy(calendar[idx])
      ref <- targets[[unit]]
      if (inherits(ref, "error")) stop(sprintf("full-period paired target failure: %s", conditionMessage(ref)))
      if (period_specific_ref) {
        dt <- d
        if (!is.null(range)) dt <- dt[date >= as.Date(range[1]) & date <= as.Date(range[2])]
        if (length(exclude_years)) dt <- dt[!year %in% exclude_years]
        elig_ref <- target_row(dt, outcomes)
        if (sum(elig_ref$keep) < 4L) stop("period-specific reference has fewer than four common paired dates")
        x <- dt$E3[elig_ref$keep]
        ref <- list(basis = ref$basis, p95 = unname(quantile(x, 0.95, type = 7)), median = median(x),
                    n_complete = sum(elig_ref$keep), events = elig_ref$events)
      }
      dfit <- d
      if (!is.null(range)) dfit <- dfit[date >= as.Date(range[1]) & date <= as.Date(range[2])]
      if (length(exclude_years)) dfit <- dfit[!year %in% exclude_years]
      eligible <- target_row(dfit, outcomes)
      ev_text <- paste(paste(names(eligible$events), eligible$events, sep = "="), collapse = ";")
      if (eligible$n < 250L || any(eligible$events < 40)) {
        return(list(estimates = empty_per_unit(), failures = data.frame(country = "MX", unit = unit,
          scenario = scenario, phase = "eligibility", reason = sprintf("eligible complete dates=%d; required >=250 dates and >=40 events per endpoint", eligible$n),
          n_complete = eligible$n, events = ev_text, stringsAsFactors = FALSE), audit = NULL, cache = NULL))
      }
      fit <- fit_stacked_contrast(dfit, outcomes = outcomes, exposure = "E3", date = "date",
        offset = NULL, basis = ref$basis, p95 = ref$p95, median = ref$median,
        nuisance = "separate", full_covariance = TRUE)
      counts <- eligible$events
      if (length(counts) < 3L) counts <- c(counts, rep(NA_real_, 3L - length(counts)))
      rows <- list()
      for (contrast_name in names(fit$contrasts)) {
        item <- fit$contrasts[[contrast_name]]
        for (cov_name in names(item$coefficient_contrast$by_covariance)) {
          q <- item$coefficient_contrast$by_covariance[[cov_name]]
          rows[[length(rows) + 1L]] <- data.frame(country = "MX", unit = unit, scenario = scenario,
            exposure = "E3", nuisance = "separate", reference = if (period_specific_ref) "period_specific" else "full_2015_2022",
            contrast = contrast_name, covariance = cov_name, estimate = unname(q["estimate"]),
            variance = unname(q["variance"]), se = unname(q["se"]), RR = unname(q["RR"]),
            lo = unname(q["lo"]), hi = unname(q["hi"]), n_complete = eligible$n,
            events_first = unname(counts[1]), events_second = unname(counts[2]), events_third = unname(counts[3]),
            p95 = ref$p95, median = ref$median,
            dB1 = unname(fit$dB[1L]), dB2 = unname(fit$dB[2L]), dB3 = unname(fit$dB[3L]),
            target_n_complete = ref$n_complete, stringsAsFactors = FALSE)
        }
      }
      a <- fit$audit
      audit <- data.frame(country = "MX", unit = unit, scenario = scenario,
        converged = a$converged, boundary = a$boundary, rank = a$rank, columns = a$columns,
        n = a$n, events = a$events, warnings = paste(a$warnings, collapse = " | "), formula = a$formula,
        stringsAsFactors = FALSE)
      cache <- if (save_fits) save_cache(fit, unit, scenario, ref, eligible$n, input_hashes) else NULL
      list(estimates = rbindlist(rows, fill = TRUE), failures = empty_failure(), audit = audit, cache = cache)
    }, error = function(e) list(estimates = empty_per_unit(), failures = data.frame(country = "MX", unit = unit,
      scenario = scenario, phase = "fit", reason = conditionMessage(e), n_complete = NA_integer_,
      events = "", stringsAsFactors = FALSE), audit = NULL, cache = NULL))
  }
  message(sprintf("MX %s: %d states, %d workers", scenario, length(units), min(opt$workers, length(units))))
  results <- if (.Platform$OS.type == "windows" || opt$workers == 1L || length(units) == 1L)
    lapply(units, worker) else parallel::mclapply(units, worker,
      mc.cores = min(as.integer(opt$workers), length(units)), mc.preschedule = FALSE, mc.set.seed = FALSE)
  result <- list(estimates = bind_rows(lapply(results, `[[`, "estimates"), empty_per_unit),
    failures = bind_rows(lapply(results, `[[`, "failures"), empty_failure),
    audits = bind_rows(lapply(results, `[[`, "audit"), empty_fit_audit),
    caches = bind_rows(lapply(results, `[[`, "cache"), empty_cache_manifest), units = length(units))
  write_scenario(result, scenario, list(outcomes = outcomes, reference = if (period_specific_ref) "period_specific_p95_median_on_common_eligible_dates; 2015-2022 spline basis" else "unit-specific full paired cohort P95/median/basis; fixed before period filters",
    range = range, exclude_years = exclude_years, minimum_complete_dates = 250L,
    minimum_events_each_endpoint = 40L, offset = NULL, nuisance = "separate",
    target_n_by_unit = lapply(targets, function(z) if (inherits(z, "error")) NULL else z$n_complete)))
  append_cache_manifest(result$caches)
  result
}

run_loso <- function(primary) {
  d <- primary[primary$scenario == "primary_complete_source_2015_2022" &
               primary$contrast == "y_urinary_vs_y_K80" & primary$covariance == "hac14", , drop = FALSE]
  if (anyDuplicated(d$unit) || nrow(d) != length(unique(d$unit))) {
    stop("leave-one-state-out source must contain exactly one urinary-vs-K80 HAC14 row per state")
  }
  rows <- list()
  all_states <- sort(unique(d$unit))
  for (omitted in all_states) {
    keep <- d[d$unit != omitted, , drop = FALSE]
    if (anyDuplicated(keep$unit) || nrow(keep) != length(unique(keep$unit))) stop("LOSO exclusion left duplicate state rows")
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(keep$estimate, keep$variance, method = method)
      if (p$k != length(unique(keep$unit))) stop("LOSO pooling k does not equal unique included states")
      rows[[length(rows) + 1L]] <- data.frame(country = "MX", scenario = "sensitivity_leave_one_state_out",
        omitted_unit = omitted, contrast = "y_urinary_vs_y_K80", covariance = "hac14", method = method,
        k = p$k, estimate = p$estimate, se = p$se, RR = p$RR, lo = p$RR95CI[1], hi = p$RR95CI[2],
        I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
        failure = if (is.null(p$failure)) "" else p$failure, stringsAsFactors = FALSE)
    }
  }
  atomic_csv(rbindlist(rows, fill = TRUE), file.path(result_dir, "sensitivity_leave_one_state_out_MX.csv"))
  scenario_summaries[["sensitivity_leave_one_state_out"]] <<- list(scenario = "sensitivity_leave_one_state_out",
    run_id = run_id, states_omitted = length(all_states), primary_input_sha256 = sha256(file.path(result_dir,
      "primary_complete_source_2015_2022_MX_per_unit.csv")), output_sha256 = sha256(file.path(result_dir,
      "sensitivity_leave_one_state_out_MX.csv")))
  update_run_manifest()
}

repair_exports <- function() {
  archive_dir <- file.path(result_dir, "invalid_exports_pre_dBrepair")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  per_unit_paths <- sort(list.files(result_dir, pattern = "_MX_per_unit\\.csv$", full.names = TRUE))
  if (length(per_unit_paths) != 9L) stop(sprintf("expected 9 per-unit scenario tables; found %d", length(per_unit_paths)))
  archival_paths <- c(per_unit_paths,
    sub("_MX_per_unit\\.csv$", "_MX_pooled.csv", per_unit_paths),
    file.path(result_dir, "sensitivity_leave_one_state_out_MX.csv"),
    file.path(result_dir, "run_manifest.json"))
  if (any(!file.exists(archival_paths))) stop("cannot archive incomplete pre-repair exports")
  before <- lapply(archival_paths, function(path) {
    dest <- file.path(archive_dir, basename(path))
    current_hash <- sha256(path)
    if (!file.exists(dest) && !file.copy(path, dest, overwrite = FALSE)) stop(sprintf("failed to archive invalid export %s", path))
    archive_hash <- sha256(dest)
    list(path = basename(path), sha256_before_repair = current_hash, archived_sha256 = archive_hash,
         bytes_before_repair = file.info(path)$size, archive_matches_current = identical(current_hash, archive_hash))
  })
  names(before) <- basename(archival_paths)
  all_repaired <- list()
  for (path in per_unit_paths) {
    d <- read_table(path)
    keys <- c("country", "unit", "scenario", "exposure", "nuisance", "reference", "contrast", "covariance")
    check_columns(d, c(keys, "dB", "estimate", "variance", "se", "RR", "lo", "hi",
                       "n_complete", "events_first", "events_second", "events_third", "p95", "median"), basename(path))
    group_id <- interaction(d$country, d$unit, d$scenario, d$exposure, d$nuisance,
      d$reference, d$contrast, d$covariance, drop = TRUE, lex.order = TRUE)
    groups <- split(seq_len(nrow(d)), group_id)
    collapsed <- vector("list", length(groups))
    gi <- 0L
    original_rows <- nrow(d)
    invariant_columns <- setdiff(names(d), "dB")
    for (idx in groups) {
      g <- d[idx]
      if (nrow(g) != 3L) stop(sprintf("expected exactly three spline-component export rows for %s; found %d", basename(path), nrow(g)))
      for (nm in invariant_columns) {
        if (data.table::uniqueN(g[[nm]]) != 1L) stop(sprintf("duplicate export rows disagree in %s for %s/%s/%s/%s",
          nm, g$scenario[1L], g$unit[1L], g$contrast[1L], g$covariance[1L]))
      }
      db <- as.numeric(g$dB)
      if (length(db) != 3L || any(!is.finite(db))) stop("spline dB component vector is invalid")
      z <- g[1L]
      z[, `:=`(dB1 = db[1L], dB2 = db[2L], dB3 = db[3L], dB = NULL)]
      collapsed[[gi <- gi + 1L]] <- z
    }
    fixed <- rbindlist(collapsed, use.names = TRUE, fill = TRUE)
    dBkeys <- c("country", "unit", "scenario", "contrast", "covariance")
    if (anyDuplicated(fixed[, dBkeys, with = FALSE])) stop("repaired per-unit key is still duplicated")
    db_audit <- fixed[, .(n_dB1 = data.table::uniqueN(dB1), n_dB2 = data.table::uniqueN(dB2),
      n_dB3 = data.table::uniqueN(dB3)), by = .(scenario, unit)]
    if (any(db_audit$n_dB1 != 1L | db_audit$n_dB2 != 1L | db_audit$n_dB3 != 1L)) {
      stop(sprintf("spline dB component tuple differs across contrasts/covariance in %s", basename(path)))
    }
    if (nrow(fixed) != 12L * length(unique(fixed$unit)) && !grepl("controlled", basename(path))) {
      stop(sprintf("per-unit repaired row count is inconsistent with 3 contrasts x 4 covariance for %s", basename(path)))
    }
    if (grepl("controlled", basename(path)) && nrow(fixed) != 4L * length(unique(fixed$unit))) {
      stop("controlled-fill repaired row count is inconsistent with 1 contrast x 4 covariance")
    }
    cols <- names(fixed)
    cols <- cols[cols != "dB1" & cols != "dB2" & cols != "dB3"]
    after_median <- match("median", cols)
    cols <- append(cols, c("dB1", "dB2", "dB3"), after = after_median)
    data.table::setcolorder(fixed, cols)
    atomic_csv(fixed, path)
    stem <- sub("_MX_per_unit\\.csv$", "", basename(path))
    pooled <- pool_scenario(fixed, stem)
    if (nrow(pooled)) for (gname in unique(paste(pooled$contrast, pooled$covariance, sep = "\r"))) {
      g <- pooled[paste(pooled$contrast, pooled$covariance, sep = "\r") == gname, , drop = FALSE]
      original_g <- fixed[fixed$contrast == g$contrast[1L] & fixed$covariance == g$covariance[1L], , drop = FALSE]
      expected_k <- length(unique(original_g$unit))
      if (any(as.integer(g$k) != expected_k)) stop(sprintf("repaired pooling k mismatch in %s", stem))
    }
    pooled_path <- file.path(result_dir, paste0(stem, "_MX_pooled.csv"))
    atomic_csv(pooled, pooled_path)
    all_repaired[[stem]] <- list(path = basename(path), rows_before = original_rows, rows_after = nrow(fixed),
      unit_count = length(unique(fixed$unit)), duplicate_rows_collapsed = original_rows - nrow(fixed),
      sha256_per_unit = sha256(path), sha256_pooled = sha256(pooled_path), pooled_groups = nrow(pooled))
  }
  primary_path <- file.path(result_dir, "primary_complete_source_2015_2022_MX_per_unit.csv")
  primary <- read_table(primary_path)
  loso <- primary[contrast == "y_urinary_vs_y_K80" & covariance == "hac14"]
  if (nrow(loso) != length(unique(loso$unit)) || nrow(loso) != 32L) stop("repaired primary HAC14 source is not one row x 32 states")
  loso_rows <- list()
  for (omitted in sort(unique(loso$unit))) {
    keep <- loso[unit != omitted]
    if (nrow(keep) != length(unique(keep$unit)) || nrow(keep) != 31L) stop("repaired LOSO cohort is not 31 unique states")
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(keep$estimate, keep$variance, method = method)
      if (p$k != 31L) stop("repaired LOSO pooled k does not equal 31")
      loso_rows[[length(loso_rows) + 1L]] <- data.frame(country = "MX", scenario = "sensitivity_leave_one_state_out",
        omitted_unit = omitted, contrast = "y_urinary_vs_y_K80", covariance = "hac14", method = method,
        k = p$k, estimate = p$estimate, se = p$se, RR = p$RR, lo = p$RR95CI[1], hi = p$RR95CI[2],
        I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
        failure = if (is.null(p$failure)) "" else p$failure, stringsAsFactors = FALSE)
    }
  }
  loso_path <- file.path(result_dir, "sensitivity_leave_one_state_out_MX.csv")
  atomic_csv(rbindlist(loso_rows, fill = TRUE), loso_path)

  manifest_path <- file.path(result_dir, "run_manifest.json")
  run_manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  for (stem in names(all_repaired)) {
    entry <- all_repaired[[stem]]
    s <- sub("_MX_per_unit\\.csv$", "", entry$path)
    run_manifest$scenarios[[s]]$units_fit <- entry$unit_count
    run_manifest$scenarios[[s]]$eligible_unit_k <- entry$unit_count
    run_manifest$scenarios[[s]]$output_sha256$per_unit <- entry$sha256_per_unit
    run_manifest$scenarios[[s]]$output_sha256$pooled <- entry$sha256_pooled
  }
  run_manifest$scenarios[["sensitivity_leave_one_state_out"]]$states_omitted <- 32L
  run_manifest$scenarios[["sensitivity_leave_one_state_out"]]$eligible_unit_k <- 31L
  run_manifest$scenarios[["sensitivity_leave_one_state_out"]]$output_sha256 <- sha256(loso_path)
  run_manifest$updated_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  atomic_json(run_manifest, manifest_path)
  audit <- list(status = "resolved", method = "aggregate-only export repair; GLM fits and full fit caches were not rerun",
    repaired_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    defect = "fit$dB was a length-three spline component vector; data.frame recycled scalar estimates into three exported rows per state/contrast/covariance",
    repair = "verified the three duplicated rows had identical model statistics, retained dB1/dB2/dB3 in order, reduced to one unique row per state key, recomputed all pools and leave-one-state-out estimates",
    archived_invalid_export_sha256 = before, scenarios = all_repaired,
    loso_sha256 = sha256(loso_path), run_manifest_sha256 = sha256(manifest_path),
    assertions = list(three_identical_statistic_rows_per_original_group = TRUE,
      unique_unit_per_scenario_contrast_covariance_after_repair = TRUE,
      pooled_k_equals_unique_unit_count = TRUE, loso_k_equals_31 = TRUE))
  atomic_json(audit, file.path(result_dir, "export_repair_audit.json"))
  message(sprintf("Updated-MX export repair complete: %d per-unit tables; pooled k now equals unique state count; full model fits/caches unchanged", length(all_repaired)))
  invisible(audit)
}

build_2023 <- function(base, health_roles, hashes) {
  weather_path <- in_path("mx_weather_2023.csv")
  if (!file.exists(weather_path)) stop("mx_weather_2023.csv has not been uploaded; 2023 refits cannot run")
  weather23 <- read_table(weather_path)
  check_columns(weather23, c("state", "year", "month", "day", "tmean"), "2023 Mexico weather")
  weather23[, `:=`(state = as.character(state), year = as.integer(year), month = as.integer(month), day = as.integer(day))]
  weather23[, date := as.Date(sprintf("%04d-%02d-%02d", year, month, day))]
  if (nrow(weather23) != 11680L || anyNA(weather23$date) || any(!is.finite(weather23$tmean)) ||
      any(weather23$year != 2023L) || key_duplicates(weather23, c("state", "date"))) stop("2023 weather grid is not 32 complete states by 365 days")
  clin23 <- copy(health_roles[source_role == "extension_2023" & format(date, "%Y") == "2023"])
  if (nrow(clin23) != 11680L || key_duplicates(clin23, c("state", "date"))) stop("2023 clinical role lacks one actual-event-date row per state/day")
  # Complete weather history is joined first so E3 for 1 January uses 30-31 Dec 2022.
  wxhist <- rbindlist(list(base[, .(unit, date, tmean)], weather23[, .(unit = state, date, tmean)]), use.names = TRUE)
  setorder(wxhist, unit, date)
  if (key_duplicates(wxhist, c("unit", "date"))) stop("duplicate weather state-date in 2015-2023 concatenation")
  wxhist[, E3 := calendar_ma(tmean, date, k = 3L), by = unit]
  wx23 <- wxhist[format(date, "%Y") == "2023", .(unit, date, tmean, E3)]
  clin23 <- outcomes_from_health(clin23)
  clin23 <- clin23[, .(unit = state, date, source_role, k80, k81, n20_only, n21_only,
    n22_other_disease, n23_unspecified_renal_colic, y_K80, y_K81, y_urinary, y_sum3)]
  out <- merge(wx23, clin23, by = c("unit", "date"), all = TRUE, sort = FALSE)
  if (nrow(out) != 11680L || anyNA(out$y_K80) || anyNA(out$y_K81) || anyNA(out$y_urinary) || anyNA(out$E3)) {
    stop("2023 weather/clinical join is incomplete")
  }
  out <- add_calendar_fields(out)
  path <- in_path("mx_updated_2023_calendar.csv.gz")
  atomic_gz_csv(out, path)
  manifest <- list(created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    source_role = "extension_2023 filtered to actual INGRE event dates in calendar year 2023",
    rows = nrow(out), units = length(unique(out$unit)), first_date = as.character(min(out$date)), last_date = as.character(max(out$date)),
    E3_rule = "calculated over complete 2015-2023 actual weather calendar before filtering; Jan 2023 includes Dec 2022 lags",
    input_hashes_sha256 = as.list(c(hashes, weather23 = sha256(weather_path))), output_sha256 = sha256(path))
  atomic_json(manifest, file.path(data_dir, "mx_updated_2023_calendar_manifest.json"))
  out
}

# Aggregate-only correction for the already fitted runs; deliberately bypasses
# model input loading and never opens or rewrites the full RDS fit caches.
if (opt$stage == "repair-exports") {
  repair_exports()
  quit(save = "no", status = 0L)
}

# Build and audit canonical data before any 2019/new-cohort effects are fitted.
if (opt$stage %in% c("build", "smoke", "all", "models")) {
  inputs <- read_inputs()
  old <- inputs$old
  health <- inputs$health
  updated <- inputs$new
  controlled <- inputs$controlled
  input_hashes <- c(inputs$hashes, canonical = sha256(in_path("mx_updated_calendar.csv.gz")),
                    controlled2019 = sha256(in_path("mx_updated_controlled_2019_calendar.csv.gz")))
  message(sprintf("Canonical source verified: %d rows, %d states, %s to %s; full 2019 rows=%d",
    nrow(updated), length(unique(updated$unit)), min(updated$date), max(updated$date), nrow(health[format(date, "%Y") == "2019"])))
  if (opt$stage == "build") quit(save = "no", status = 0L)

  if (opt$stage %in% c("smoke", "all", "models")) {
    # One real-state full-covariance fit smoke tests the generic engine call,
    # target transfer and cache serialization without suppressing per-state errors.
    u <- sort(unique(updated$unit))[1L]
    d <- copy(updated[unit == u])
    target <- make_targets(d, all_outcomes)
    smoke <- fit_stacked_contrast(d, outcomes = all_outcomes, exposure = "E3", date = "date",
      offset = NULL, basis = target$basis, p95 = target$p95, median = target$median,
      nuisance = "separate", full_covariance = TRUE)
    if (!isTRUE(smoke$audit$converged) || length(smoke$contrasts) != 3L ||
        !all(c("quasi", "day_hc0", "hac14", "hac28") %in% names(smoke$covariance))) {
      stop("actual one-state model smoke did not return three full covariance contrasts")
    }
    message(sprintf("One-state model smoke passed: %s; complete dates=%d; coefficients=%d",
      u, target$n_complete, smoke$audit$columns))
    if (opt$stage == "smoke") quit(save = "no", status = 0L)
  }

  # Install the independent cache manifest only inside mx_updated/fits.
  if (!file.exists(file.path(fit_dir, "fit_cache_manifest.csv"))) atomic_csv(empty_cache_manifest(), file.path(fit_dir, "fit_cache_manifest.csv"))
  scenario_summaries <- list()
  update_run_manifest()

  # Preferred repaired cohort; full paired targets and spline basis are computed
  # on all common-eligible 2015-2022 dates before any scenario restrictions.
  primary <- run_scenario(updated, "primary_complete_source_2015_2022", outcomes = all_outcomes,
                          target_calendar = updated, save_fits = TRUE)
  # Frozen original gap cohort is read-only and retained as a direct comparator.
  original <- run_scenario(old, "comparator_original_gap_2015_2018_2020_2022", outcomes = all_outcomes,
                           target_calendar = old, save_fits = FALSE)
  # Controlled fill: old baseline health on every other year plus the complete
  # all-source 2019 block once. The already-computed complete-weather E3 is kept.
  run_scenario(controlled, "sensitivity_controlled_2019_fill", outcomes = c("y_K80", "y_urinary"),
               target_calendar = controlled, save_fits = TRUE)
  run_scenario(updated, "sensitivity_exclude_2020_2021", outcomes = all_outcomes,
               target_calendar = updated, exclude_years = c(2020L, 2021L), save_fits = FALSE)
  run_scenario(updated, "temporal_2015_2018_fixed_full_reference", outcomes = all_outcomes,
               target_calendar = updated, range = c("2015-01-01", "2018-12-31"), save_fits = FALSE)
  run_scenario(updated, "temporal_2020_2022_fixed_full_reference", outcomes = all_outcomes,
               target_calendar = updated, range = c("2020-01-01", "2022-12-31"), save_fits = FALSE)
  run_loso(primary$estimates)

  # 2023 is an independent supplementary temporal refit, not an extrapolation
  # or forecast. It is run only after all-source actual-date health and weather
  # inputs pass exact state-day coverage checks.
  updated23 <- build_2023(old, inputs$clinical_roles, inputs$hashes)
  # Fit reference basis/targets are held at 2015-2022 values; only model dates
  # are 2023. A second fit uses period-specific 2023 P95/median with that basis.
  run_scenario(updated23, "temporal_replication_2023_fixed_2015_2022_reference", outcomes = all_outcomes,
               target_calendar = updated, save_fits = TRUE)
  run_scenario(updated23, "temporal_replication_2023_jan_nov_boundary_fixed_reference", outcomes = all_outcomes,
               target_calendar = updated, range = c("2023-01-01", "2023-11-30"), save_fits = FALSE)
  run_scenario(updated23, "temporal_replication_2023_period_specific_reference", outcomes = all_outcomes,
               target_calendar = updated, period_specific_ref = TRUE, save_fits = FALSE)
  message(sprintf("Updated Mexico scenario run finished: %d scenario summaries, results in %s", length(scenario_summaries), result_dir))
}
