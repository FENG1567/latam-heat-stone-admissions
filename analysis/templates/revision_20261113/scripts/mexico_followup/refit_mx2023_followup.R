#!/usr/bin/env Rscript

# Bounded Mexico 2023 next-discharge-year follow-up and source-row sensitivities.
# Writes only to revision_20261113/results/mexico_followup. The registered
# calendars and fit caches are read-only inputs.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(root = NULL, mode = "A")
for (a in args) {
  if (grepl("^--root=", a)) opt$root <- sub("^--root=", "", a)
  else if (grepl("^--mode=", a)) opt$mode <- toupper(sub("^--mode=", "", a))
  else stop(sprintf("unknown argument: %s", a), call. = FALSE)
}
if (!opt$mode %in% c("A", "B")) stop("--mode must be A or B", call. = FALSE)
if (is.null(opt$root)) {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(file_arg)) stop("supply --root=/path/to/revision_20261113")
  opt$root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                            winslash = "/", mustWork = TRUE)
}
revision_root <- normalizePath(opt$root, winslash = "/", mustWork = TRUE)
base_root <- normalizePath(file.path(revision_root, ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(revision_root, "results", "mexico_followup")
input_dir <- file.path(revision_root, "data", "mexico_followup")
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]),
                             winslash = "/", mustWork = TRUE)
engine_path <- file.path(base_root, "scripts", "repair_engine.R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", BLAS_NUM_THREADS = "1")
local_library <- file.path(base_root, "Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
for (pkg in c("data.table", "jsonlite", "metafor", "splines")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("required package '%s' is unavailable", pkg))
}
DT <- asNamespace("data.table")
fread <- DT$fread
fwrite <- DT$fwrite
if (!file.exists(engine_path)) stop(sprintf("repair engine missing: %s", engine_path))
source(engine_path, local = FALSE)

outcomes <- c("y_K80", "y_K81", "y_urinary")
endpoint_set <- c("K80", "K81", "N20", "N21", "N22", "N23")
old_2023_scenario <- "temporal_replication_2023_fixed_2015_2022_reference"
output_hashes <- character()
input_hashes <- character()

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- suppressWarnings(system2("sha256sum", shQuote(normalizePath(path, winslash = "/")),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("sha256sum failed for %s: %s", path, paste(out, collapse = " ")))
  sub("\\s+.*$", "", out[1L])
}
read_csv <- function(path, gzip = grepl("\\.gz$", path, ignore.case = TRUE)) {
  if (!file.exists(path)) stop(sprintf("required input file is missing: %s", path), call. = FALSE)
  x <- if (gzip) fread(cmd = sprintf("gzip -dc -- %s", shQuote(path)), showProgress = FALSE,
                       check.names = FALSE, na.strings = c("", "NA"), keepLeadingZeros = TRUE) else
    fread(path, showProgress = FALSE, check.names = FALSE, na.strings = c("", "NA"), keepLeadingZeros = TRUE)
  as.data.frame(x, stringsAsFactors = FALSE)
}
normalize_state <- function(x) {
  token <- trimws(as.character(x))
  if (anyNA(token) || any(!grepl("^[0-9]{1,2}$", token))) stop("state codes must be integer tokens 1-32")
  number <- suppressWarnings(as.integer(token))
  if (anyNA(number) || any(number < 1L | number > 32L)) stop("state codes must be within 01-32")
  sprintf("%02d", number)
}
write_csv <- function(x, name) {
  path <- file.path(out_dir, name)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = out_dir, fileext = ".tmp")
  fwrite(as.data.frame(x), tmp, na = "", quote = "auto")
  if (!file.rename(tmp, path)) stop(sprintf("atomic CSV replacement failed: %s", path), call. = FALSE)
  output_hashes[[name]] <<- sha256(path)
  path
}
write_json <- function(x, name) {
  path <- file.path(out_dir, name)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = out_dir, fileext = ".tmp")
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = NA)
  if (!file.rename(tmp, path)) stop(sprintf("atomic JSON replacement failed: %s", path), call. = FALSE)
  output_hashes[[name]] <<- sha256(path)
  path
}
read_calendar <- function(path) {
  x <- read_csv(path)
  required <- c("unit", "date", "E3", outcomes)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(sprintf("calendar %s lacks: %s", basename(path), paste(missing, collapse = ", ")))
  x$unit <- normalize_state(x$unit)
  x$date <- as.Date(x$date)
  if (anyNA(x$date) || anyDuplicated(paste(x$unit, x$date))) stop(sprintf("invalid or duplicate dates in %s", basename(path)))
  for (nm in c("E3", outcomes)) x[[nm]] <- as.numeric(x[[nm]])
  x
}
validate_calendar <- function(x, year, expected_days, allow_two_leading_e3_na = FALSE) {
  units <- sort(unique(x$unit))
  if (length(units) != 32L) stop(sprintf("expected 32 Mexico states; found %d", length(units)))
  per_unit <- table(x$unit)
  if (any(per_unit != expected_days)) stop(sprintf("%d calendar rows per state, expected %d", expected_days, expected_days))
  if (any(!format(x$date, "%Y") %in% as.character(year))) stop(sprintf("calendar contains dates outside %s", paste(year, collapse = ",")))
  if (any(is.infinite(x$E3))) stop(sprintf("calendar has infinite E3 for %s", paste(year, collapse = ",")))
  if (anyNA(x$E3)) {
    if (!allow_two_leading_e3_na) stop(sprintf("calendar has missing E3 for %s", paste(year, collapse = ",")))
    for (u in units) {
      z <- x[x$unit == u, , drop = FALSE]
      z <- z[order(z$date), , drop = FALSE]
      expected_na <- seq_len(nrow(z)) <= 2L
      if (!identical(is.na(z$E3), expected_na)) stop(sprintf("unexpected missing E3 pattern for state %s", u))
    }
  }
  invisible(units)
}
read_source_agg <- function() {
  p <- file.path(input_dir, "all_daily.csv.gz")
  manifest_path <- file.path(input_dir, "complete.json")
  identity_path <- file.path(input_dir, "count_identity_verified.json")
  if (!file.exists(manifest_path) || !file.exists(identity_path)) stop("parser-repaired source manifests are missing")
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  if (!identical(manifest$status, "complete")) stop("parser-repaired aggregate manifest is incomplete")
  listed <- manifest$outputs[["all_daily.csv.gz"]]
  if (is.null(listed) || !identical(listed, sha256(p))) stop("all_daily hash does not match parser-repaired manifest")
  identity <- jsonlite::fromJSON(identity_path, simplifyVector = FALSE)
  if (!identical(identity$status, "exact_identity_verified") ||
      !identical(identity$source_aggregate_sha256, sha256(p))) stop("registered primary count identity evidence is missing or mismatched")
  a <- read_csv(p, gzip = TRUE)
  req <- c("source_year", "discharge_year", "state", "date", "endpoint", "parser", "n")
  missing <- setdiff(req, names(a))
  if (length(missing)) stop(sprintf("source aggregate lacks: %s", paste(missing, collapse = ", ")))
  a$source_year <- as.integer(a$source_year)
  a$discharge_year <- as.integer(a$discharge_year)
  a$state <- normalize_state(a$state)
  a$date <- as.Date(a$date)
  a$endpoint <- as.character(a$endpoint)
  a$parser <- as.character(a$parser)
  a$n <- as.numeric(a$n)
  if (any(!a$parser %in% c("format_repaired", "original_width_accepted"))) stop("unexpected source parser label")
  # The format_repaired projection is the frozen selected-column reconstruction
  # that was independently proven identical to the registered daily counts.
  a <- a[a$parser == "format_repaired", , drop = FALSE]
  if (anyNA(a$date) || anyNA(a$source_year) || anyNA(a$n) || any(a$n < 0) || any(a$n != floor(a$n))) {
    stop("source aggregate has invalid dates, source years, or counts")
  }
  if (any(!a$endpoint %in% endpoint_set)) stop("source aggregate has an unexpected endpoint code")
  if (any(!a$state %in% sprintf("%02d", 1:32))) stop("source aggregate has a non-Mexico state code")
  if (anyNA(a$discharge_year) || anyDuplicated(paste(a$source_year, a$discharge_year, a$state, a$date, a$endpoint))) {
    stop("source aggregate must be unique by source_year/discharge_year/state/date/endpoint")
  }
  a
}
daily_outcomes <- function(a, state_col = "state", grid = NULL) {
  if (is.null(grid)) {
    if (!nrow(a)) stop("cannot form a daily grid from an empty source subset")
    key <- paste(a[[state_col]], as.character(a$date), sep = "|")
    keys <- unique(key)
    out <- data.frame(state = sub("\\|.*$", "", keys),
                      date = as.Date(sub("^.*\\|", "", keys)), stringsAsFactors = FALSE)
  } else {
    grid_state <- if ("state" %in% names(grid)) as.character(grid$state) else as.character(grid$unit)
    out <- data.frame(state = grid_state, date = as.Date(grid$date), stringsAsFactors = FALSE)
    if (anyNA(out$date) || anyDuplicated(paste(out$state, out$date))) stop("daily outcome grid has invalid or duplicate keys")
  }
  for (ep in endpoint_set) {
    z <- a[a$endpoint == ep, , drop = FALSE]
    if (nrow(z)) {
      agg <- aggregate(z$n, by = list(state = z[[state_col]], date = z$date), FUN = sum)
      names(agg)[3L] <- ep
      out <- merge(out, agg, by = c("state", "date"), all.x = TRUE, sort = FALSE)
    } else out[[ep]] <- 0
    out[[ep]][is.na(out[[ep]])] <- 0
  }
  out$y_K80 <- out$K80
  out$y_K81 <- out$K81
  out$y_urinary <- out$N20 + out$N21 + out$N22 + out$N23
  out$y_sum3 <- out$y_K80 + out$y_K81 + out$y_urinary
  out
}
assert_daily_equal <- function(calendar, aggregate_daily, label) {
  agg_key <- paste(aggregate_daily$state, aggregate_daily$date, sep = "|")
  cal_key <- paste(calendar$unit, calendar$date, sep = "|")
  idx <- match(cal_key, agg_key)
  if (anyNA(idx) || nrow(calendar) != nrow(aggregate_daily)) {
    stop(sprintf("%s daily grid differs from calendar (calendar rows=%d; aggregate rows=%d)",
                 label, nrow(calendar), nrow(aggregate_daily)))
  }
  rows <- lapply(outcomes, function(nm) {
    d <- as.numeric(calendar[[nm]]) - as.numeric(aggregate_daily[[nm]][idx])
    data.frame(outcome = nm, rows_checked = length(d), differing_rows = sum(d != 0),
      max_abs_difference = if (length(d)) max(abs(d)) else 0,
      calendar_total = sum(calendar[[nm]]), aggregate_total = sum(aggregate_daily[[nm]]),
      stringsAsFactors = FALSE)
  })
  check <- do.call(rbind, rows)
  check$comparison <- label
  check$exact_reproduction <- check$differing_rows == 0L
  check
}
make_targets <- function(d) {
  good <- is.finite(d$E3) & !is.na(d$date)
  for (nm in outcomes) good <- good & is.finite(d[[nm]]) & d[[nm]] >= 0
  x <- d$E3[good]
  if (sum(good) < 4L || length(unique(x)) < 2L) stop("insufficient complete target dates")
  list(basis = new_temperature_basis(x, df = 3L),
       p95 = unname(stats::quantile(x, 0.95, type = 7)), median = stats::median(x),
       n_complete = sum(good), events = vapply(outcomes, function(nm) sum(d[[nm]][good]), numeric(1)))
}
attempt_row <- function(unit, scenario, n, events, eligible, reason = "") {
  data.frame(unit = unit, scenario = scenario, attempted = TRUE, eligible = eligible,
    n_complete = n, events_K80 = unname(events[1]), events_K81 = unname(events[2]),
    events_urinary = unname(events[3]), reason = reason, stringsAsFactors = FALSE)
}
empty_per_unit <- function() data.frame(country = character(), unit = character(), scenario = character(),
  exposure = character(), nuisance = character(), reference = character(), contrast = character(),
  covariance = character(), estimate = numeric(), variance = numeric(), se = numeric(), RR = numeric(),
  lo = numeric(), hi = numeric(), n_complete = integer(), events_first = numeric(),
  events_second = numeric(), events_third = numeric(), p95 = numeric(), median = numeric(),
  dB1 = numeric(), dB2 = numeric(), dB3 = numeric(), target_n_complete = integer(), stringsAsFactors = FALSE)
fit_cohort <- function(calendar, target_calendar, scenario, file_prefix, minimum_dates = 250L,
                       minimum_events = 40L, reference_label = "full_2015_2022") {
  units <- sort(unique(calendar$unit))
  if (!identical(units, sort(unique(target_calendar$unit)))) stop(sprintf("%s target/model unit mismatch", scenario))
  per_unit <- list(); failures <- list(); audits <- list(); attempts <- list()
  for (unit in units) {
    d <- calendar[calendar$unit == unit, , drop = FALSE]
    refdata <- target_calendar[target_calendar$unit == unit, , drop = FALSE]
    ref <- tryCatch(make_targets(refdata), error = function(e) e)
    if (inherits(ref, "error")) {
      attempts[[length(attempts) + 1L]] <- attempt_row(unit, scenario, NA_integer_, c(NA, NA, NA), FALSE,
        paste("target failure:", conditionMessage(ref)))
      failures[[length(failures) + 1L]] <- data.frame(country = "MX", unit = unit, scenario = scenario,
        phase = "target", reason = conditionMessage(ref), n_complete = NA_integer_, events = "", stringsAsFactors = FALSE)
      next
    }
    good <- is.finite(d$E3) & !is.na(d$date)
    for (nm in outcomes) good <- good & is.finite(d[[nm]]) & d[[nm]] >= 0
    ev <- vapply(outcomes, function(nm) sum(d[[nm]][good]), numeric(1))
    n <- sum(good)
    if (n < minimum_dates || any(ev < minimum_events)) {
      reason <- sprintf("eligible complete dates=%d; required >=%d dates and >=%d events per endpoint", n, minimum_dates, minimum_events)
      attempts[[length(attempts) + 1L]] <- attempt_row(unit, scenario, n, ev, FALSE, reason)
      failures[[length(failures) + 1L]] <- data.frame(country = "MX", unit = unit, scenario = scenario,
        phase = "eligibility", reason = reason, n_complete = n,
        events = paste(paste(names(ev), ev, sep = "="), collapse = ";"), stringsAsFactors = FALSE)
      next
    }
    attempts[[length(attempts) + 1L]] <- attempt_row(unit, scenario, n, ev, TRUE)
    fit <- tryCatch(fit_stacked_contrast(d, outcomes = outcomes, exposure = "E3", date = "date", offset = NULL,
      basis = ref$basis, p95 = ref$p95, median = ref$median, nuisance = "separate", full_covariance = TRUE),
      error = function(e) e)
    if (inherits(fit, "error")) {
      failures[[length(failures) + 1L]] <- data.frame(country = "MX", unit = unit, scenario = scenario,
        phase = "fit", reason = conditionMessage(fit), n_complete = n, events = "", stringsAsFactors = FALSE)
      attempts[[length(attempts)]]$reason <- paste("fit failure:", conditionMessage(fit))
      next
    }
    rows <- list()
    for (contrast_name in names(fit$contrasts)) {
      item <- fit$contrasts[[contrast_name]]
      for (cov_name in names(item$coefficient_contrast$by_covariance)) {
        q <- item$coefficient_contrast$by_covariance[[cov_name]]
        rows[[length(rows) + 1L]] <- data.frame(country = "MX", unit = unit, scenario = scenario,
          exposure = "E3", nuisance = "separate", reference = reference_label,
          contrast = contrast_name, covariance = cov_name, estimate = unname(q["estimate"]),
          variance = unname(q["variance"]), se = unname(q["se"]), RR = unname(q["RR"]),
          lo = unname(q["lo"]), hi = unname(q["hi"]), n_complete = n,
          events_first = unname(ev[1]), events_second = unname(ev[2]), events_third = unname(ev[3]),
          p95 = ref$p95, median = ref$median, dB1 = unname(fit$dB[1L]), dB2 = unname(fit$dB[2L]),
          dB3 = unname(fit$dB[3L]), target_n_complete = ref$n_complete, stringsAsFactors = FALSE)
      }
    }
    a <- fit$audit
    audits[[length(audits) + 1L]] <- data.frame(country = "MX", unit = unit, scenario = scenario,
      converged = a$converged, boundary = a$boundary, rank = a$rank, columns = a$columns,
      n = a$n, events = a$events, warnings = paste(a$warnings, collapse = " | "),
      formula = a$formula, stringsAsFactors = FALSE)
    per_unit[[length(per_unit) + 1L]] <- do.call(rbind, rows)
  }
  estimates <- if (length(per_unit)) do.call(rbind, per_unit) else empty_per_unit()
  fail <- if (length(failures)) do.call(rbind, failures) else data.frame(country = character(), unit = character(),
    scenario = character(), phase = character(), reason = character(), n_complete = integer(), events = character())
  audit <- if (length(audits)) do.call(rbind, audits) else data.frame(country = character(), unit = character(),
    scenario = character(), converged = logical(), boundary = logical(), rank = integer(), columns = integer(),
    n = integer(), events = numeric(), warnings = character(), formula = character())
  attempt <- if (length(attempts)) do.call(rbind, attempts) else data.frame()
  pooled <- pool_estimates(estimates, scenario)
  write_csv(estimates, paste0(file_prefix, "_per_unit.csv"))
  write_csv(pooled, paste0(file_prefix, "_pooled.csv"))
  write_csv(fail, paste0(file_prefix, "_failure_audit.csv"))
  write_csv(audit, paste0(file_prefix, "_fit_audit.csv"))
  write_csv(attempt, paste0(file_prefix, "_attempts.csv"))
  data.frame(scenario = scenario, units_attempted = length(units), units_eligible = sum(attempt$eligible),
    units_fit = length(unique(estimates$unit)), units_failed = nrow(fail), criteria_days = minimum_dates,
    criteria_events_each_endpoint = minimum_events, stringsAsFactors = FALSE)
}
pool_estimates <- function(estimates, scenario) {
  if (!nrow(estimates)) return(data.frame(country = character(), scenario = character(), contrast = character(),
    covariance = character(), method = character(), k = integer(), estimate = numeric(), se = numeric(), RR = numeric(),
    lo = numeric(), hi = numeric(), I2 = numeric(), tau2 = numeric(), concordance = numeric(), failure = character()))
  groups <- unique(estimates[c("contrast", "covariance")])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- estimates[estimates$contrast == groups$contrast[i] & estimates$covariance == groups$covariance[i], , drop = FALSE]
    if (anyDuplicated(g$unit)) stop(sprintf("%s pooling has duplicate unit rows", scenario))
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(g$estimate, g$variance, method = method)
      rows[[length(rows) + 1L]] <- data.frame(country = "MX", scenario = scenario, contrast = g$contrast[1],
        covariance = g$covariance[1], method = method, k = p$k, estimate = p$estimate, se = p$se, RR = p$RR,
        lo = p$RR95CI[1], hi = p$RR95CI[2], I2 = p$I2, tau2 = p$tau2,
        concordance = p$concordance, failure = if (is.null(p$failure)) "" else p$failure, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
read_existing <- function(path, scenario = NULL) {
  x <- read_csv(path, gzip = FALSE)
  if (!is.null(scenario) && "scenario" %in% names(x)) x <- x[x$scenario == scenario, , drop = FALSE]
  x
}
compare_tables <- function(new, old, keys, metrics, label) {
  if (!all(keys %in% names(new)) || !all(keys %in% names(old))) stop(sprintf("%s comparison lacks keys", label))
  nk <- do.call(paste, c(new[keys], sep = "|")); ok <- do.call(paste, c(old[keys], sep = "|"))
  if (anyDuplicated(nk) || anyDuplicated(ok) || !setequal(nk, ok)) stop(sprintf("%s keys do not match", label))
  oi <- match(nk, ok)
  rows <- lapply(metrics, function(nm) {
    if (!nm %in% names(new) || !nm %in% names(old)) stop(sprintf("%s comparison lacks metric %s", label, nm))
    a <- as.numeric(new[[nm]]); b <- as.numeric(old[[nm]][oi]); delta <- a - b
    data.frame(metric = nm, rows_compared = length(delta), differing_rows = sum(delta != 0),
      max_abs_difference = max(abs(delta), na.rm = TRUE), stringsAsFactors = FALSE)
  })
  check <- do.call(rbind, rows)
  check$comparison <- label
  check
}
candidate_adjustment <- function(calendar, candidates, sources, years, label) {
  chosen <- candidates[candidates$source_year %in% sources &
    format(candidates$date, "%Y") %in% as.character(years), , drop = FALSE]
  agg <- data.frame(state = character(), date = as.Date(character()), endpoint = character(), removed = numeric())
  if (nrow(chosen)) {
    agg <- aggregate(chosen$candidate_repeat_extra,
      by = list(state = chosen$state, date = chosen$date, endpoint = chosen$endpoint), FUN = sum)
    names(agg)[4L] <- "removed"
  }
  result <- apply_endpoint_adjustment(calendar, agg, label)
  changed <- rowSums(result$daily_adjustment[c("subtract_K80", "subtract_K81", "subtract_urinary")]) > 0
  result$daily_adjustment <- result$daily_adjustment[changed, , drop = FALSE]
  result
}
source_filter_calendar <- function(calendar, source_daily, label) {
  map <- paste(source_daily$state, source_daily$date, sep = "|")
  idx <- match(paste(calendar$unit, calendar$date, sep = "|"), map)
  if (anyNA(idx) || length(idx) != nrow(calendar)) stop(sprintf("%s source aggregate lacks calendar dates", label))
  out <- calendar
  for (nm in outcomes) out[[nm]] <- source_daily[[nm]][idx]
  out$y_sum3 <- out$y_K80 + out$y_K81 + out$y_urinary
  out
}
apply_endpoint_adjustment <- function(calendar, removals, label) {
  k <- paste(removals$state, removals$date, sep = "|")
  rkey <- paste(calendar$unit, calendar$date, sep = "|")
  add <- matrix(0, nrow = nrow(calendar), ncol = 3L, dimnames = list(NULL, outcomes))
  for (ep in endpoint_set) {
    hit <- which(removals$endpoint == ep)
    if (!length(hit)) next
    jj <- match(k[hit], rkey)
    if (anyNA(jj)) stop(sprintf("%s candidate rows do not map to analysis calendar", label))
    col <- if (ep == "K80") "y_K80" else if (ep == "K81") "y_K81" else "y_urinary"
    for (j in seq_along(hit)) add[jj[j], col] <- add[jj[j], col] + removals$removed[hit[j]]
  }
  out <- calendar
  count_rows <- list()
  for (nm in outcomes) {
    out[[nm]] <- as.numeric(out[[nm]]) - add[, nm]
    if (any(out[[nm]] < 0)) stop(sprintf("%s subtraction would make %s counts negative", label, nm))
    count_rows[[length(count_rows) + 1L]] <- data.frame(sensitivity = label, outcome = nm,
      candidate_events_subtracted = sum(add[, nm]), original_events = sum(calendar[[nm]]),
      bounded_events = sum(out[[nm]]), dates_affected = sum(add[, nm] > 0),
      stringsAsFactors = FALSE)
  }
  out$y_sum3 <- out$y_K80 + out$y_K81 + out$y_urinary
  list(calendar = out, counts = do.call(rbind, count_rows), daily_adjustment = data.frame(
    state = calendar$unit, date = calendar$date, subtract_K80 = add[, "y_K80"],
    subtract_K81 = add[, "y_K81"], subtract_urinary = add[, "y_urinary"], stringsAsFactors = FALSE))
}
write_reproduction_check <- function(new, old, name) {
  keys <- c("unit", "contrast", "covariance")
  metrics <- c("estimate", "variance", "RR", "lo", "hi", "n_complete", "events_first", "events_second",
               "events_third", "p95", "median", "dB1", "dB2", "dB3", "target_n_complete")
  check <- compare_tables(new, old, keys, metrics, name)
  write_csv(check, paste0(name, "_check.csv"))
  if (any(check$max_abs_difference > 1e-10)) stop(sprintf("%s failed exact numeric reproduction (tolerance 1e-10)", name))
  check
}
write_delta <- function(new, reference, keys, metrics, label, file) {
  take <- unique(c(keys, metrics))
  if (!all(take %in% names(new)) || !all(take %in% names(reference))) stop(sprintf("%s delta inputs lack required columns", label))
  n <- new[take]; r <- reference[take]
  names(n)[match(metrics, names(n))] <- paste0(metrics, "_sensitivity")
  names(r)[match(metrics, names(r))] <- paste0(metrics, "_reference")
  z <- merge(r, n, by = keys, all = TRUE, sort = FALSE)
  has_r <- !is.na(z[[paste0(metrics[1L], "_reference")]])
  has_n <- !is.na(z[[paste0(metrics[1L], "_sensitivity")]])
  z$membership <- ifelse(has_r & has_n, "both_fitted", ifelse(has_r, "reference_only", "sensitivity_only"))
  for (nm in metrics) z[[paste0(nm, "_difference")]] <- z[[paste0(nm, "_sensitivity")]] - z[[paste0(nm, "_reference")]]
  z$sensitivity <- label
  write_csv(z, file)
  z
}
count_table <- function(a, label, source_years, years) {
  z <- a[a$source_year %in% source_years & format(a$date, "%Y") %in% as.character(years), , drop = FALSE]
  if (!nrow(z)) return(data.frame())
  z$admission_year <- as.integer(format(z$date, "%Y"))
  agg <- aggregate(z$n, by = list(source_year = z$source_year, discharge_year = z$discharge_year,
    admission_year = z$admission_year, endpoint = z$endpoint), FUN = sum)
  names(agg)[5L] <- "events"
  states <- aggregate(z$state, by = list(source_year = z$source_year, discharge_year = z$discharge_year,
    admission_year = z$admission_year, endpoint = z$endpoint), FUN = function(v) length(unique(v)))
  days <- aggregate(paste(z$state, z$date), by = list(source_year = z$source_year, discharge_year = z$discharge_year,
    admission_year = z$admission_year, endpoint = z$endpoint), FUN = function(v) length(unique(v)))
  agg$states <- states[[5L]]; agg$state_days <- days[[5L]]; agg$cohort_label <- label
  agg
}
write_manifest <- function(name, status, details = list()) {
  payload <- c(list(status = status, mode = opt$mode, created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    R = R.version.string, revision_root = revision_root, base_root = base_root,
    workers = 1L, BLAS_threads = 1L, source_engine = engine_path,
    input_sha256 = as.list(input_hashes), output_sha256 = as.list(output_hashes)), details)
  write_json(payload, name)
}

main_A <- function() {
  primary_path <- file.path(base_root, "data", "mx_updated_calendar.csv.gz")
  calendar23_path <- file.path(base_root, "data", "mx_updated_2023_calendar.csv.gz")
  source_path <- file.path(input_dir, "all_daily.csv.gz")
  parser_manifest_path <- file.path(input_dir, "complete.json")
  identity_path <- file.path(input_dir, "count_identity_verified.json")
  old_per_path <- file.path(base_root, "results", "mx_updated", paste0(old_2023_scenario, "_MX_per_unit.csv"))
  old_pool_path <- file.path(base_root, "results", "mx_updated", paste0(old_2023_scenario, "_MX_pooled.csv"))
  for (p in c(primary_path, calendar23_path, source_path, parser_manifest_path, identity_path, old_per_path, old_pool_path)) {
    if (!file.exists(p)) stop(sprintf("required A input is missing: %s", p))
    input_hashes[basename(p)] <<- sha256(p)
  }
  input_hashes["repair_engine.R"] <<- sha256(engine_path)
  input_hashes["refit_mx2023_followup.R"] <<- sha256(script_path)
  primary <- read_calendar(primary_path)
  cal23_existing <- read_calendar(calendar23_path)
  validate_calendar(primary, 2015:2022, 2922L, allow_two_leading_e3_na = TRUE)
  validate_calendar(cal23_existing, 2023, 365L)
  source <- read_source_agg()
  annual_target <- source[source$source_year %in% 2023:2024, , drop = FALSE]
  annual_discharge_check <- as.data.frame(table(annual_target$source_year, annual_target$discharge_year), stringsAsFactors = FALSE)
  names(annual_discharge_check) <- c("source_year", "discharge_year", "aggregate_rows")
  annual_discharge_check$source_year <- as.integer(as.character(annual_discharge_check$source_year))
  annual_discharge_check$discharge_year <- as.integer(as.character(annual_discharge_check$discharge_year))
  annual_discharge_check <- annual_discharge_check[annual_discharge_check$aggregate_rows > 0, , drop = FALSE]
  if (nrow(annual_target) && any(annual_target$source_year != annual_target$discharge_year)) {
    stop("2023/2024 target rows include off-year discharge records; see source audit before follow-up")
  }
  write_csv(annual_discharge_check, "a_source_2023_2024_discharge_year_check.csv")
  source23 <- source[source$source_year == 2023L & source$discharge_year == 2023L &
    format(source$date, "%Y") == "2023", , drop = FALSE]
  source24_23 <- source[source$source_year == 2024L & source$discharge_year == 2024L &
    format(source$date, "%Y") == "2023", , drop = FALSE]
  daily23 <- daily_outcomes(source23, grid = cal23_existing[c("unit", "date")])
  baseline23 <- cal23_existing
  rebuilt_idx <- match(paste(baseline23$unit, baseline23$date, sep = "|"), paste(daily23$state, daily23$date, sep = "|"))
  if (anyNA(rebuilt_idx)) stop("source-year 2023 aggregate lacks an existing 2023 state-day")
  baseline23$y_K80 <- daily23$y_K80[rebuilt_idx]
  baseline23$y_K81 <- daily23$y_K81[rebuilt_idx]
  baseline23$y_urinary <- daily23$y_urinary[rebuilt_idx]
  baseline23$y_sum3 <- baseline23$y_K80 + baseline23$y_K81 + baseline23$y_urinary
  source_check <- assert_daily_equal(cal23_existing, daily23, "source-year 2023 baseline")
  write_csv(source_check, "a_source_2023_calendar_count_check.csv")
  if (any(!source_check$exact_reproduction)) stop("source-year 2023 aggregate differs from existing 2023 calendar; follow-up fit withheld")
  add24 <- daily_outcomes(source24_23, grid = cal23_existing[c("unit", "date")])
  extended23 <- baseline23
  k <- match(paste(extended23$unit, extended23$date, sep = "|"), paste(add24$state, add24$date, sep = "|"))
  if (anyNA(k)) stop("source-year 2024 admission-date-2023 aggregate lacks extension state-days")
  for (nm in outcomes) extended23[[nm]] <- baseline23[[nm]] + add24[[nm]][k]
  extended23$y_sum3 <- extended23$y_K80 + extended23$y_K81 + extended23$y_urinary
  if (!identical(as.numeric(extended23$E3), as.numeric(cal23_existing$E3))) stop("extension unexpectedly changed actual-calendar E3")
  if (!identical(as.character(extended23$date), as.character(cal23_existing$date)) ||
      !identical(as.character(extended23$unit), as.character(cal23_existing$unit))) stop("extension changed calendar keys")
  source24 <- source[source$source_year == 2024L & source$discharge_year == 2024L, , drop = FALSE]
  prior24 <- count_table(source24, "2024 discharge source; admission year 2015-2022 (reported separately)",
    source_years = 2024L, years = 2015:2022)
  follow24 <- count_table(source24, "2024 discharge source; admission year 2023 (next-year follow-up)",
    source_years = 2024L, years = 2023L)
  if (nrow(prior24)) write_csv(prior24, "a_source_2024_prior_admission_counts_2015_2022.csv")
  else write_csv(data.frame(source_year = integer(), admission_year = integer(), endpoint = character(), events = numeric(), states = integer(), state_days = integer(), cohort_label = character()), "a_source_2024_prior_admission_counts_2015_2022.csv")
  write_csv(follow24, "a_source_2024_admission_year_2023_followup_counts.csv")
  count_rows <- list()
  for (nm in outcomes) {
    count_rows[[length(count_rows) + 1L]] <- data.frame(cohort = "source2023_baseline", outcome = nm,
      events = sum(baseline23[[nm]]), stringsAsFactors = FALSE)
    count_rows[[length(count_rows) + 1L]] <- data.frame(cohort = "source2024_followup_added_to_2023", outcome = nm,
      events = sum(extended23[[nm]]) - sum(baseline23[[nm]]), stringsAsFactors = FALSE)
    count_rows[[length(count_rows) + 1L]] <- data.frame(cohort = "source2023_plus_2024_extended_2023", outcome = nm,
      events = sum(extended23[[nm]]), stringsAsFactors = FALSE)
  }
  write_csv(do.call(rbind, count_rows), "a_cohort_event_counts.csv")
  base_summary <- fit_cohort(baseline23, primary, old_2023_scenario, "a_baseline_2023_rebuilt",
    minimum_dates = 250L, minimum_events = 40L)
  new_per <- read_csv(file.path(out_dir, "a_baseline_2023_rebuilt_per_unit.csv"), gzip = FALSE)
  old_per <- read_existing(old_per_path, old_2023_scenario)
  reproduction <- write_reproduction_check(new_per, old_per, "a_baseline_reproduction_per_unit")
  old_fail <- read_existing(file.path(base_root, "results", "mx_updated", paste0(old_2023_scenario, "_MX_failure_audit.csv")), old_2023_scenario)
  new_fail <- read_csv(file.path(out_dir, "a_baseline_2023_rebuilt_failure_audit.csv"), gzip = FALSE)
  old_failed <- sort(as.character(old_fail$unit)); new_failed <- sort(as.character(new_fail$unit))
  if (!identical(old_failed, new_failed)) stop("baseline 2023 eligibility membership differs from the registered 30-state fit")
  if (nrow(old_per) == 0L || length(unique(old_per$unit)) != 30L || length(old_failed) != 2L) {
    stop(sprintf("baseline membership expected 30 fitted and 2 ineligible states; got %d and %d",
      length(unique(old_per$unit)), length(old_failed)))
  }
  new_pool <- read_csv(file.path(out_dir, "a_baseline_2023_rebuilt_pooled.csv"), gzip = FALSE)
  old_pool <- read_existing(old_pool_path, old_2023_scenario)
  pooled_check <- compare_tables(new_pool, old_pool, c("contrast", "covariance", "method"),
    c("estimate", "se", "RR", "lo", "hi", "k"), "a_baseline_reproduction_pooled")
  write_csv(pooled_check, "a_baseline_reproduction_pooled_check.csv")
  if (any(pooled_check$max_abs_difference > 1e-10)) stop("baseline pooled results do not reproduce the registered 2023 fit")
  baseline_rrr <- new_pool[new_pool$contrast == "y_urinary_vs_y_K80" & new_pool$covariance == "hac14" & new_pool$method == "DL", , drop = FALSE]
  if (nrow(baseline_rrr) != 1L) stop("baseline direct urinary-vs-K80 HAC14/DL row was not unique")
  if (abs(baseline_rrr$RR - 1.100) > 0.002) stop(sprintf("baseline direct RRR %.6f does not reproduce the expected approximately 1.100", baseline_rrr$RR))
  write_csv(base_summary, "a_baseline_2023_cohort_summary.csv")
  extended_summary <- fit_cohort(extended23, primary,
    "temporal_replication_2023_followup_2024_fixed_2015_2022_reference",
    "a_extended_2023_2024", minimum_dates = 250L, minimum_events = 40L)
  write_csv(extended_summary, "a_extended_2023_2024_cohort_summary.csv")
  write_manifest("a_manifest.json", "complete", list(
    input_sources = list(source_aggregate = "parser_repaired/all_daily.csv.gz, parser=format_repaired; source year/discharge year 2023 baseline and 2024 rows admitted in 2023",
      reference = "registered complete-source 2015-2022 mx_updated_calendar.csv.gz; per-state basis/P95/shared median fixed as in prior 2023 refit",
      weather_exposure = "existing 2023 actual-date E3 preserved exactly; it was built on full 2015-2023 weather history including 2022 lags",
      source_year_2024_prior_admissions = "counts for admission years 2015-2022 reported separately and excluded from the 2023 follow-up fit"),
    eligibility = list(all_32_states_attempted = TRUE, minimum_complete_dates = 250L,
      minimum_events_each_of_three_endpoints = 40L, model = "three-outcome joint quasi-Poisson; separate calendar nuisance terms; HAC14/HAC28; DL and REML_HK pooling"),
    baseline_reproduction = list(states_fit = 30L, states_ineligible = 2L,
      max_per_unit_difference = max(reproduction$max_abs_difference),
      max_pooled_difference = max(pooled_check$max_abs_difference),
      RRR_urinary_vs_K80_HAC14_DL = baseline_rrr$RR,
      registered_RRR = 1.100, registered_RRR_numeric_check = abs(baseline_rrr$RR - 1.100) < 0.002)))
}

main_B <- function() {
  primary_path <- file.path(base_root, "data", "mx_updated_calendar.csv.gz")
  calendar23_path <- file.path(base_root, "data", "mx_updated_2023_calendar.csv.gz")
  agg_path <- file.path(input_dir, "all_daily.csv.gz")
  candidate_path <- file.path(input_dir, "all_duplicates.csv")
  candidate_manifest_path <- file.path(input_dir, "complete.json")
  identity_path <- file.path(input_dir, "count_identity_verified.json")
  principal_path <- agg_path
  registered_primary_path <- file.path(base_root, "results", "mx_updated", "primary_complete_source_2015_2022_MX_per_unit.csv")
  registered_23_path <- file.path(base_root, "results", "mx_updated", paste0(old_2023_scenario, "_MX_per_unit.csv"))
  for (p in c(primary_path, calendar23_path, agg_path, candidate_path, candidate_manifest_path,
              identity_path, registered_primary_path, registered_23_path)) {
    if (!file.exists(p)) stop(sprintf("required B input is missing: %s", p))
    input_hashes[basename(p)] <<- sha256(p)
  }
  input_hashes["repair_engine.R"] <<- sha256(engine_path)
  input_hashes["refit_mx2023_followup.R"] <<- sha256(script_path)
  extended_base_path <- file.path(out_dir, "a_extended_2023_2024_per_unit.csv")
  extended_pool_path <- file.path(out_dir, "a_extended_2023_2024_pooled.csv")
  if (!file.exists(extended_base_path) || !file.exists(extended_pool_path)) stop("A extension outputs are missing; complete mode A before mode B")
  input_hashes["a_extended_2023_2024_per_unit.csv"] <<- sha256(extended_base_path)
  input_hashes["a_extended_2023_2024_pooled.csv"] <<- sha256(extended_pool_path)
  extended_base <- read_csv(extended_base_path, gzip = FALSE)
  extended_base_pool <- read_csv(extended_pool_path, gzip = FALSE)
  candidate_manifest <- jsonlite::fromJSON(candidate_manifest_path, simplifyVector = FALSE)
  if (!identical(candidate_manifest$status, "complete")) stop("parser-repaired source manifest is not complete")
  expected_candidate_hash <- candidate_manifest$outputs[["all_duplicates.csv"]]
  if (is.null(expected_candidate_hash) || !identical(expected_candidate_hash, sha256(candidate_path))) stop("all_duplicates hash does not match parser-repaired manifest")
  identity <- jsonlite::fromJSON(identity_path, simplifyVector = FALSE)
  if (!identical(identity$status, "exact_identity_verified") || !identical(identity$source_aggregate_sha256, sha256(agg_path))) {
    stop("registered primary count identity evidence is missing or mismatched")
  }
  primary <- read_calendar(primary_path)
  calendar23 <- read_calendar(calendar23_path)
  validate_calendar(calendar23, 2023, 365L)
  source <- read_source_agg()
  cand <- read_csv(candidate_path, gzip = FALSE)
  if (!all(c("source_year", "discharge_year", "state", "date", "endpoint", "n") %in% names(cand))) {
    stop("all_duplicates candidate table has unexpected schema")
  }
  cand$source_year <- as.integer(cand$source_year); cand$state <- normalize_state(cand$state)
  cand$discharge_year <- as.integer(cand$discharge_year)
  cand$date <- as.Date(cand$date); cand$endpoint <- as.character(cand$endpoint)
  cand$candidate_repeat_extra <- as.numeric(cand$n)
  if (anyNA(cand$date) || anyNA(cand$candidate_repeat_extra) || any(cand$candidate_repeat_extra < 0) ||
      any(cand$candidate_repeat_extra != floor(cand$candidate_repeat_extra)) ||
      any(!cand$endpoint %in% endpoint_set)) stop("candidate repeat table contains invalid dates/codes/counts")
  all_source <- source
  all_source <- all_source[all_source$parser == "format_repaired", , drop = FALSE]
  # Reconstruct the registered primary daily outcomes before either sensitivity.
  prim_src <- all_source[all_source$source_year %in% 2015:2023 &
    format(all_source$date, "%Y") %in% as.character(2015:2022), , drop = FALSE]
  prim_daily <- daily_outcomes(prim_src, grid = primary[c("unit", "date")])
  prim_reconstruction <- assert_daily_equal(primary, prim_daily, "all-source aggregate vs registered primary")
  write_csv(prim_reconstruction, "b_primary_aggregate_reproduction_check.csv")
  if (any(!prim_reconstruction$exact_reproduction)) stop("all-source aggregate does not reproduce registered primary; B sensitivities withheld")
  # Verify that 2023/2024 source aggregate rows reconcile with A's event base.
  ext_src <- all_source[all_source$source_year %in% 2023:2024 &
    all_source$discharge_year == all_source$source_year & format(all_source$date, "%Y") == "2023", , drop = FALSE]
  ext_daily <- daily_outcomes(ext_src, grid = calendar23[c("unit", "date")])
  ext_base <- calendar23
  ext23_idx <- match(paste(ext_base$unit, ext_base$date, sep = "|"), paste(ext_daily$state, ext_daily$date, sep = "|"))
  if (anyNA(ext23_idx)) stop("all-source aggregate lacks 2023 extension state-days")
  for (nm in outcomes) ext_base[[nm]] <- ext_daily[[nm]][ext23_idx]
  ext_base$y_sum3 <- ext_base$y_K80 + ext_base$y_K81 + ext_base$y_urinary
  a_extended <- calendar23
  source2024_23 <- source[source$source_year == 2024L & source$discharge_year == 2024L & format(source$date, "%Y") == "2023", , drop = FALSE]
  daily2024_23 <- daily_outcomes(source2024_23, grid = calendar23[c("unit", "date")])
  j24 <- match(paste(a_extended$unit, a_extended$date, sep = "|"), paste(daily2024_23$state, daily2024_23$date, sep = "|"))
  source2023_23 <- source[source$source_year == 2023L & source$discharge_year == 2023L & format(source$date, "%Y") == "2023", , drop = FALSE]
  daily2023_23 <- daily_outcomes(source2023_23, grid = calendar23[c("unit", "date")])
  j23 <- match(paste(a_extended$unit, a_extended$date, sep = "|"), paste(daily2023_23$state, daily2023_23$date, sep = "|"))
  if (anyNA(j24) || anyNA(j23)) stop("source 2023/2024 aggregate lacks 2023 extension state-days")
  for (nm in outcomes) a_extended[[nm]] <- daily2023_23[[nm]][j23] + daily2024_23[[nm]][j24]
  a_extended$y_sum3 <- a_extended$y_K80 + a_extended$y_K81 + a_extended$y_urinary
  ext_check <- assert_daily_equal(a_extended, ext_daily, "all-source aggregate vs A source 2023+2024 extension")
  write_csv(ext_check, "b_extension_aggregate_reproduction_check.csv")
  if (any(!ext_check$exact_reproduction)) stop("all-source aggregate does not reproduce A extension; B sensitivities withheld")
  # Add discharge-year strata without exposing state/date rows.
  z <- all_source[all_source$source_year %in% 2015:2024 & format(all_source$date, "%Y") %in% as.character(2015:2024), , drop = FALSE]
  z$admission_year <- as.integer(format(z$date, "%Y"))
  sourceyear_discharge_trace <- aggregate(z$n, by = list(source_year = z$source_year,
    discharge_year = z$discharge_year, admission_year = z$admission_year, endpoint = z$endpoint), FUN = sum)
  names(sourceyear_discharge_trace)[5L] <- "events"
  sourceyear_discharge_trace$source_year_matches_discharge_year <-
    sourceyear_discharge_trace$source_year == sourceyear_discharge_trace$discharge_year
  write_csv(sourceyear_discharge_trace, "b_source_discharge_year_trace.csv")
  registered_primary <- read_existing(registered_primary_path, "primary_complete_source_2015_2022")
  reference_rows <- registered_primary[registered_primary$contrast == "y_urinary_vs_y_K80" &
    registered_primary$covariance == "hac14", , drop = FALSE]
  targets <- do.call(rbind, lapply(sort(unique(primary$unit)), function(unit) {
    t <- make_targets(primary[primary$unit == unit, , drop = FALSE])
    data.frame(unit = unit, p95 = t$p95, median = t$median, target_n_complete = t$n_complete)
  }))
  ref_check <- compare_tables(targets, reference_rows, c("unit"), c("p95", "median", "target_n_complete"),
    "b_registered_primary_fixed_targets")
  write_csv(ref_check, "b_registered_primary_integrity_check.csv")
  if (any(ref_check$max_abs_difference > 1e-10)) stop("primary fixed targets differ from registered complete-source targets")
  # Candidate-repeat upper bound: only exact non-ID matches within annual files.
  cand_primary <- candidate_adjustment(primary, cand, sources = 2015:2023, years = 2015:2022,
    label = "primary_2015_2022_nonID_repeat_upper_bound")
  cand_extended <- candidate_adjustment(a_extended, cand, sources = 2023:2024, years = 2023,
    label = "extended_2023_nonID_repeat_upper_bound")
  # Source-year/discharge-year concordance is an independent source-selection
  # sensitivity, not an episode or patient deduplication.
  concord_primary_src <- all_source[all_source$source_year %in% 2015:2023 &
    all_source$discharge_year == all_source$source_year &
    format(all_source$date, "%Y") %in% as.character(2015:2022), , drop = FALSE]
  concord_primary_daily <- daily_outcomes(concord_primary_src, grid = primary[c("unit", "date")])
  concord_primary <- source_filter_calendar(primary, concord_primary_daily, "primary source_year=discharge_year")
  concord_extended_src <- all_source[all_source$source_year %in% 2023:2024 &
    all_source$discharge_year == all_source$source_year & format(all_source$date, "%Y") == "2023", , drop = FALSE]
  concord_extended_daily <- daily_outcomes(concord_extended_src, grid = calendar23[c("unit", "date")])
  concord_extended <- source_filter_calendar(a_extended, concord_extended_daily, "extension source_year=discharge_year")
  candidate_totals <- rbind(cand_primary$counts, cand_extended$counts)
  write_csv(candidate_totals, "b_candidate_repeat_subtraction_counts.csv")
  write_csv(rbind(cand_primary$daily_adjustment, cand_extended$daily_adjustment), "b_candidate_repeat_daily_adjustments.csv")
  cand$admission_year <- as.integer(format(cand$date, "%Y"))
  cand_source_counts <- aggregate(cand$candidate_repeat_extra,
    by = list(source_year = cand$source_year, discharge_year = cand$discharge_year,
      admission_year = cand$admission_year, endpoint = cand$endpoint), FUN = sum)
  names(cand_source_counts)[5L] <- "nonID_repeat_candidate_extras"
  write_csv(cand_source_counts, "b_candidate_repeat_source_counts.csv")
  concord_totals <- rbind(data.frame(sensitivity = "primary_2015_2022_source_year_equals_discharge_year",
    outcome = outcomes, registered_events = colSums(primary[outcomes]), concordant_events = colSums(concord_primary[outcomes]),
    events_excluded = colSums(primary[outcomes]) - colSums(concord_primary[outcomes]), stringsAsFactors = FALSE),
    data.frame(sensitivity = "extended_2023_source_year_equals_discharge_year",
    outcome = outcomes, registered_events = colSums(a_extended[outcomes]), concordant_events = colSums(concord_extended[outcomes]),
    events_excluded = colSums(a_extended[outcomes]) - colSums(concord_extended[outcomes]), stringsAsFactors = FALSE))
  write_csv(concord_totals, "b_sourceyear_concordance_event_counts.csv")
  # No B model is started until the registered primary is proven exactly reproducible.
  scenarios <- list(
    list(cal = cand_primary$calendar, ref = primary, scenario = "primary_repeat_bound_2015_2022",
      prefix = "b_primary_repeat_bound", label = "full_2015_2022_fixed_reference"),
    list(cal = cand_extended$calendar, ref = primary, scenario = "extended_2023_repeat_bound",
      prefix = "b_extended_2023_repeat_bound", label = "full_2015_2022_fixed_reference"),
    list(cal = concord_primary, ref = primary, scenario = "primary_source_year_equals_discharge_year",
      prefix = "b_primary_sourceyear_concordant", label = "full_2015_2022_fixed_reference"),
    list(cal = concord_extended, ref = primary, scenario = "extended_2023_source_year_equals_discharge_year",
      prefix = "b_extended_2023_sourceyear_concordant", label = "full_2015_2022_fixed_reference")
  )
  summaries <- list()
  for (s in scenarios) {
    summaries[[length(summaries) + 1L]] <- fit_cohort(s$cal, s$ref, s$scenario, s$prefix,
      minimum_dates = 250L, minimum_events = 40L, reference_label = s$label)
  }
  write_csv(do.call(rbind, summaries), "b_sensitivity_cohort_summaries.csv")
  registered_primary_pool <- read_existing(file.path(base_root, "results", "mx_updated",
    "primary_complete_source_2015_2022_MX_pooled.csv"), "primary_complete_source_2015_2022")
  sensitivity_specs <- list(
    list(prefix = "b_primary_repeat_bound", baseline = registered_primary, baseline_pool = registered_primary_pool, label = "primary_nonID_repeat_upper_bound"),
    list(prefix = "b_primary_sourceyear_concordant", baseline = registered_primary, baseline_pool = registered_primary_pool, label = "primary_source_year_equals_discharge_year"),
    list(prefix = "b_extended_2023_repeat_bound", baseline = extended_base, baseline_pool = extended_base_pool, label = "extended_2023_nonID_repeat_upper_bound"),
    list(prefix = "b_extended_2023_sourceyear_concordant", baseline = extended_base, baseline_pool = extended_base_pool, label = "extended_2023_source_year_equals_discharge_year")
  )
  per_unit_deltas <- list(); pooled_deltas <- list()
  for (sp in sensitivity_specs) {
    sens <- read_csv(file.path(out_dir, paste0(sp$prefix, "_per_unit.csv")), gzip = FALSE)
    pool <- read_csv(file.path(out_dir, paste0(sp$prefix, "_pooled.csv")), gzip = FALSE)
    per_unit_deltas[[length(per_unit_deltas) + 1L]] <- write_delta(sens, sp$baseline,
      c("unit", "contrast", "covariance"), c("estimate", "RR", "lo", "hi"), sp$label,
      paste0(sp$prefix, "_delta_per_unit.csv"))
    pooled_deltas[[length(pooled_deltas) + 1L]] <- write_delta(pool, sp$baseline_pool,
      c("contrast", "covariance", "method"), c("estimate", "RR", "lo", "hi", "k"), sp$label,
      paste0(sp$prefix, "_delta_pooled.csv"))
  }
  write_csv(do.call(rbind, per_unit_deltas), "b_sensitivity_deltas_per_unit.csv")
  write_csv(do.call(rbind, pooled_deltas), "b_sensitivity_deltas_pooled.csv")
  write_manifest("b_manifest.json", "complete", list(
  candidate_repeat_definition = "all_duplicates.csv groups exact non-ID row repetitions within each annual source; candidate extra rows are a conservative upper bound, not confirmed duplicate episodes/people",
    candidate_sources_primary = "source_year 2015-2023; admission dates 2015-2022",
    candidate_sources_extended = "source_year 2023-2024; admission dates in 2023",
    source_year_concordance_definition = "retained only aggregate rows where source_year equals discharge_year; independent source-selection sensitivity, not deduplication",
    primary_reconstruction_gate = list(max_count_difference = max(prim_reconstruction$max_abs_difference),
      exact_reproduction = TRUE, outcome_rows_checked = sum(prim_reconstruction$rows_checked)),
    fixed_reference = "same per-state 2015-2022 E3 natural spline basis/P95/shared median as registered complete-source Mexico fit",
    eligibility = list(all_32_states_attempted = TRUE, minimum_complete_dates = 250L,
      minimum_events_each_of_three_endpoints = 40L),
    source_trace = "b_source_discharge_year_trace.csv reports source-year/discharge-year/admission-year/endpoint counts; source_2024 admission years 2015-2022 remain separately tabulated in A output"))
}

tryCatch({
  if (opt$mode == "A") main_A() else main_B()
}, error = function(e) {
  try(write_manifest(if (opt$mode == "A") "a_manifest.json" else "b_manifest.json", "failed",
                     list(error = conditionMessage(e))), silent = TRUE)
  stop(e)
})

