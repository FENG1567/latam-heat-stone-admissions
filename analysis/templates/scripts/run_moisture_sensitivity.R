#!/usr/bin/env Rscript
# Locked, bounded moisture sensitivity run for the repaired BR/MX cohorts.
# The frozen engine is sourced read-only; this runner owns only its separate
# results/moisture_sensitivity directory.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(root = NULL, stage = "validate", workers = 4L, source_qa = NULL)
for (a in args) {
  if (grepl("^--root=", a)) opt$root <- sub("^--root=", "", a)
  else if (grepl("^--stage=", a)) opt$stage <- sub("^--stage=", "", a)
  else if (grepl("^--workers=", a)) opt$workers <- as.integer(sub("^--workers=", "", a))
  else if (grepl("^--source-qa=", a)) opt$source_qa <- sub("^--source-qa=", "", a)
  else stop(sprintf("unknown argument: %s", a), call. = FALSE)
}
if (!opt$stage %in% c("validate", "fit")) stop("--stage must be validate or fit", call. = FALSE)
if (length(opt$workers) != 1L || is.na(opt$workers) || opt$workers < 1L) stop("--workers must be a positive integer")
if (is.null(opt$root)) {
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(file_arg)) stop("cannot infer repair root; supply --root=...")
  opt$root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                            winslash = "/", mustWork = TRUE)
}
root <- normalizePath(opt$root, winslash = "/", mustWork = TRUE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]),
                             winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_path)
engine_path <- file.path(script_dir, "repair_engine.R")
plan_path <- file.path(root, "admin", "moisture_sensitivity_plan_20261113.json")
out_dir <- file.path(root, "results", "moisture_sensitivity")
qa_dir <- file.path(out_dir, "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("required package 'data.table' is unavailable")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("required package 'jsonlite' is unavailable")
if (!requireNamespace("metafor", quietly = TRUE)) stop("required package 'metafor' is unavailable")
if (!requireNamespace("splines", quietly = TRUE)) stop("required namespace 'splines' is unavailable")
data.table::setDTthreads(1L)
source(engine_path, local = FALSE)

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  cmd <- if (nzchar(Sys.which("sha256sum"))) "sha256sum" else "shasum"
  out <- if (identical(cmd, "sha256sum"))
    suppressWarnings(system2(cmd, shQuote(normalizePath(path, winslash = "/")), stdout = TRUE, stderr = TRUE)) else
      suppressWarnings(system2(cmd, c("-a", "256", shQuote(normalizePath(path, winslash = "/"))), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("SHA-256 command failed for %s: %s", path, paste(out, collapse = " ")))
  sub("\\s+.*$", "", out[1L])
}
atomic_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop(sprintf("failed to atomically write %s", path))
}
atomic_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop(sprintf("failed to atomically write %s", path))
}
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(x, tmp, compress = "gzip")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop(sprintf("failed to atomically write scalar cache %s", path))
}
hash_text <- function(x) {
  path <- tempfile(pattern = "moisture_hash_")
  on.exit(unlink(path), add = TRUE)
  writeLines(enc2utf8(paste(x, collapse = "\n")), path, useBytes = TRUE)
  sha256(path)
}
near <- function(x, y, tol = 1e-9) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol, check.attributes = FALSE))
stopf <- function(...) stop(sprintf(...), call. = FALSE)

model_spec <- list(
  status = "locked_before_effect_fits",
  settings = c("matched_unadjusted", "matched_dew3_ns3", "matched_rh3_ns3"),
  common_sample = "For each unit, the intersection of original primary eligible paired dates and dates with finite dew3 and rh3; use this same sample in all three settings.",
  outcome_order = c("y_K80", "y_K81", "y_urinary"),
  temperature_reference = "Full original eligible paired E3 calendar, unit-specific type-7 P95, median and ns(df=3) basis; unchanged by moisture date loss.",
  nuisance = "Outcome-specific disease*(B1+B2+B3+dow+month+ns(time,12)); moisture settings add disease*ns(moisture,3).",
  covariance = "Joint daily-score Bartlett HAC14; no quasi-dispersion multiplier; no other covariance output.",
  pooling = c("DL", "REML_HK"),
  contrasts = c("y_urinary_vs_y_K80", "y_urinary_vs_y_K81", "y_K81_vs_y_K80"),
  eligibility = "At least 250 matched dates and 40 events in each outcome.",
  coverage_exclusion_BR = c("29020", "16001", "41001"),
  source_qa_gate = "passed status and exact SHA-256 match for both transferred moisture calendars.",
  estimate_interpretation = "Moisture-adjusted relative associations under the additive exposure model and conditioning assumptions; no temperature-by-moisture interaction."
)
source_hashes <- c(runner = sha256(script_path), engine = sha256(engine_path), plan = sha256(plan_path))
model_lock_path <- file.path(qa_dir, "model_lock_manifest.json")

# The validation stage uses no study inputs and writes no study effect estimates.
fit_stacked_moisture <- function(d, outcome_names, basis, p95, median,
                                offset_col = NULL, moisture_col = NULL) {
  if (!is.data.frame(d) || !all(c("date", "E3", outcome_names) %in% names(d))) stop("invalid unit data for stacked fit")
  dates <- as.Date(d$date)
  if (anyNA(dates) || anyDuplicated(dates)) stop("fit data dates must be valid and unique")
  keep <- is.finite(d$E3) & is.finite(d$time_years)
  for (nm in outcome_names) keep <- keep & is.finite(d[[nm]]) & d[[nm]] >= 0
  if (!is.null(offset_col)) keep <- keep & is.finite(d[[offset_col]]) & d[[offset_col]] > 0
  if (!is.null(moisture_col)) keep <- keep & is.finite(d[[moisture_col]])
  if (!any(keep)) stop("no complete rows remain for stacked moisture fit")
  d <- d[keep, , drop = FALSE]
  dates <- as.Date(d$date)
  e <- as.numeric(d$E3)
  bb <- .basis_at(basis, e)
  dow <- factor(as.integer(format(dates, "%u")), levels = 1:7)
  month <- factor(as.integer(format(dates, "%m")), levels = 1:12)
  time <- as.numeric(d$time_years)
  off <- if (is.null(offset_col)) rep(0, nrow(d)) else log(as.numeric(d[[offset_col]]))
  common <- data.frame(B1 = bb[, 1L], B2 = bb[, 2L], B3 = bb[, 3L], dow = dow, month = month,
    time = time, off = off, date_actual = dates)
  if (!is.null(moisture_col)) {
    mm <- splines::ns(as.numeric(d[[moisture_col]]), df = 3L)
    if (ncol(mm) != 3L || any(!is.finite(mm))) stop("moisture spline did not produce three finite columns")
    common$M1 <- mm[, 1L]; common$M2 <- mm[, 2L]; common$M3 <- mm[, 3L]
  }
  st <- do.call(rbind, lapply(seq_along(outcome_names), function(j) {
    z <- common
    z$y <- as.numeric(d[[outcome_names[j]]])
    z$disease <- factor(outcome_names[j], levels = outcome_names)
    z
  }))
  rownames(st) <- NULL
  rhs <- if (is.null(moisture_col))
    "disease * (B1 + B2 + B3 + dow + month + splines::ns(time, 12))" else
      "disease * (B1 + B2 + B3 + dow + month + splines::ns(time, 12) + M1 + M2 + M3)"
  formula <- if (is.null(offset_col)) stats::as.formula(paste("y ~", rhs)) else
    stats::as.formula(paste("y ~", rhs, "+ offset(off)"))
  audited <- .fit_glm_audited(formula, st)
  fit <- audited$fit; X <- audited$X
  # Set moisture spline coordinates to zero for both temperature targets.
  # This explicitly holds the additive moisture covariates fixed in every
  # outcome-specific temperature contrast.
  target_template <- st
  if (!is.null(moisture_col)) target_template[, c("M1", "M2", "M3")] <- 0
  effects <- lapply(outcome_names, function(endpoint) {
    xhi <- .design_at(fit, target_template, endpoint = endpoint, basis = basis, temp_value = p95)
    xmed <- .design_at(fit, target_template, endpoint = endpoint, basis = basis, temp_value = median)
    as.numeric(xhi - xmed)
  })
  names(effects) <- outcome_names
  pairs <- utils::combn(seq_along(outcome_names), 2L, simplify = FALSE)
  labels <- vapply(pairs, function(p) paste0(outcome_names[p[2L]], "_vs_", outcome_names[p[1L]]), character(1))
  Ls <- lapply(pairs, function(p) setNames(effects[[p[2L]]] - effects[[p[1L]]], colnames(X)))
  names(Ls) <- labels
  V <- .score_covariances(fit, X, y = fit$y, mu = fit$fitted.values,
                          dates = st$date_actual, hac_lags = 14L, full_covariance = TRUE)
  if (!"hac14" %in% names(V$hac)) stop("score engine did not return the required HAC14 covariance")
  out <- lapply(names(Ls), function(label) {
    L <- Ls[[label]]
    cc <- .contrast_summary(L, fit, list(hac14 = V$hac$hac14))
    q <- cc$by_covariance$hac14
    list(label = label, L = L, estimate = unname(q["estimate"]), variance = unname(q["variance"]),
      se = unname(q["se"]), RR = unname(q["RR"]), lo = unname(q["lo"]), hi = unname(q["hi"]))
  })
  names(out) <- labels
  dB <- as.numeric(.basis_at(basis, p95) - .basis_at(basis, median))
  list(fit = fit, design = X, score = V, contrasts = out, dB = dB,
    p95 = p95, median = median, audit = audited$audit, retained_dates = dates)
}

manual_bartlett_hac <- function(U, dates, lag = 14L) {
  dates <- as.Date(dates)
  dnum <- as.integer(dates)
  out <- crossprod(U)
  if (lag > 0L) for (ell in seq_len(lag)) {
    right <- match(as.character(as.Date(dates - ell, origin = "1970-01-01")), as.character(dates))
    use <- which(!is.na(right))
    if (length(use)) {
      cross <- crossprod(U[use, , drop = FALSE], U[right[use], , drop = FALSE])
      w <- (lag + 1L - ell) / (lag + 1L)
      out <- out + w * (cross + t(cross))
    }
  }
  (out + t(out)) / 2
}
scalar_cache_path <- function(cache_dir, country, unit) {
  unit_hex <- paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(unit)))), collapse = "")
  file.path(cache_dir, paste0(country, "_unit_", unit_hex, ".rds"))
}
read_scalar_cache <- function(path, cache_key, country, unit) {
  if (!file.exists(path)) return(NULL)
  z <- tryCatch(readRDS(path), error = function(e) stopf("scalar cache %s cannot be read: %s", basename(path), conditionMessage(e)))
  if (!identical(z$cache_key, cache_key) || !identical(z$country, country) || !identical(z$unit, unit)) {
    stopf("scalar cache key/unit mismatch at %s", path)
  }
  expected <- c("matched_unadjusted", "matched_dew3_ns3", "matched_rh3_ns3")
  if (!is.list(z$result) || !identical(sort(names(z$result$results)), sort(expected))) stopf("scalar cache payload is malformed at %s", path)
  for (scenario in expected) {
    p <- z$result$results[[scenario]]
    if (!is.list(p) || !all(c("estimates", "failure", "audit") %in% names(p))) stopf("scalar cache scenario %s is malformed", scenario)
    e <- p$estimates
    if (nrow(e) && (nrow(e) != 3L || any(e$unit != unit) || any(e$scenario != scenario) ||
        any(e$covariance != "hac14") || anyDuplicated(e[c("unit", "contrast", "covariance")]))) {
      stopf("scalar cache unit/contrast keys are invalid for %s/%s/%s", country, unit, scenario)
    }
    if (nrow(e) && !setequal(e$contrast, c("y_K81_vs_y_K80", "y_urinary_vs_y_K80", "y_urinary_vs_y_K81"))) {
      stopf("scalar cache has incomplete contrasts for %s/%s/%s", country, unit, scenario)
    }
  }
  z
}
write_scalar_cache <- function(path, cache_key, country, unit, result, metadata) {
  atomic_rds(list(cache_key = cache_key, country = country, unit = unit,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    metadata = metadata, result = result), path)
  invisible(path)
}
run_scalar_cached <- function(path, cache_key, country, unit, compute) {
  cached <- read_scalar_cache(path, cache_key, country, unit)
  if (!is.null(cached)) return(list(result = cached$result, cache_path = path,
    cache_reused = TRUE, metadata = cached$metadata))
  fresh <- compute()
  if (!is.list(fresh) || !all(c("result", "metadata") %in% names(fresh))) stop("scalar cache compute callback returned an invalid payload")
  write_scalar_cache(path, cache_key, country, unit, fresh$result, fresh$metadata)
  list(result = fresh$result, cache_path = path, cache_reused = FALSE, metadata = fresh$metadata)
}
synthetic_validation <- function() {
  set.seed(20261113)
  n <- 720L
  date <- seq(as.Date("2017-01-01"), by = "day", length.out = n)
  time <- as.numeric(date - min(date)) / 365.25
  tmean <- 24 + 5 * sin(2 * pi * as.numeric(format(date, "%j")) / 365.25) + 0.6 * sin(seq_len(n) / 13)
  dew3 <- 16 + 4 * sin(2 * pi * as.numeric(format(date, "%j")) / 365.25 + 0.35) + rnorm(n, 0, 0.7)
  rh3 <- pmin(98, pmax(5, 60 + 16 * sin(seq_len(n) / 51) + rnorm(n, 0, 3)))
  season <- 0.15 * sin(2 * pi * as.numeric(format(date, "%j")) / 365.25)
  d <- data.frame(date = date, E3 = calendar_ma(tmean, date, 3L), time_years = time,
    dew3 = dew3, rh3 = rh3,
    y_K80 = stats::rpois(n, exp(1.1 + 0.02 * (tmean - 24) + season)),
    y_K81 = stats::rpois(n, exp(1.4 + 0.03 * (tmean - 24) + season)),
    y_urinary = stats::rpois(n, exp(1.7 + 0.04 * (tmean - 24) + season)))
  good <- is.finite(d$E3)
  basis <- new_temperature_basis(d$E3[good], df = 3L)
  p95 <- unname(stats::quantile(d$E3[good], 0.95, type = 7)); med <- stats::median(d$E3[good])
  f <- fit_stacked_moisture(d, c("y_K80", "y_K81", "y_urinary"), basis, p95, med,
    moisture_col = "dew3")
  if (!f$audit$converged || f$audit$rank != f$audit$columns || length(f$contrasts) != 3L) stop("synthetic moisture model did not pass convergence/rank checks")
  # Independently aggregate scores by calendar date and compute the Bartlett
  # HAC meat over exact daily offsets, then verify H^-1 meat H^-1.
  U0 <- rowsum(f$design * as.numeric(f$fit$y - f$fit$fitted.values),
               group = as.character(f$score$dates), reorder = FALSE)
  dates <- sort(unique(as.Date(rownames(U0))))
  U <- U0[match(as.character(dates), rownames(U0)), , drop = FALSE]
  rownames(U) <- as.character(dates)
  meat <- manual_bartlett_hac(U, dates, 14L)
  v_direct <- f$score$Hinv %*% meat %*% f$score$Hinv
  v_engine <- f$score$hac$hac14
  if (!near(v_direct, v_engine, tol = 1e-8)) stopf("synthetic direct HAC14 does not match immutable engine covariance (max abs delta %.12g; direct %.12g; engine %.12g)",
    max(abs(v_direct - v_engine)), v_direct[1L, 1L], v_engine[1L, 1L])
  beta <- stats::coef(f$fit)
  checks <- list()
  for (label in names(f$contrasts)) {
    q <- f$contrasts[[label]]; L <- q$L
    analytic_est <- as.numeric(L %*% beta)
    analytic_var <- as.numeric(t(L) %*% v_direct %*% L)
    if (!near(analytic_est, q$estimate, 1e-9) || !near(analytic_var, q$variance, 1e-8)) {
      stopf("synthetic analytic L beta / L V L assertion failed for %s", label)
    }
    checks[[label]] <- list(analytic_L_beta = analytic_est, reported_estimate = q$estimate,
      analytic_LVL = analytic_var, reported_variance = q$variance,
      estimate_abs_delta = abs(analytic_est - q$estimate), variance_abs_delta = abs(analytic_var - q$variance))
  }
  # Simulate interruption after the first completed unit, then resume. The
  # first unit must be reused unchanged and only the second unit computed.
  cache_dir <- tempfile(pattern = "moisture_scalar_cache_")
  dir.create(cache_dir, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  cache_key <- hash_text(c("synthetic-validation", source_hashes["runner"], source_hashes["engine"]))
  calls <- 0L
  stub_payload <- function(unit, failing = FALSE) {
    scenarios <- c("matched_unadjusted", "matched_dew3_ns3", "matched_rh3_ns3")
    contrasts <- c("y_K81_vs_y_K80", "y_urinary_vs_y_K80", "y_urinary_vs_y_K81")
    results <- setNames(lapply(scenarios, function(s) {
      e <- if (failing) data.frame() else data.frame(unit = unit, scenario = s, covariance = "hac14",
        contrast = contrasts, estimate = 0.1, variance = 0.01)
      list(estimates = e, failure = if (failing) data.frame(unit = unit, scenario = s,
        phase = "fit", reason = "synthetic rank-deficiency audit") else data.frame(), audit = data.frame())
    }), scenarios)
    list(results = results, matched_n = 720L, common = !failing)
  }
  compute_stub <- function(unit, failing = FALSE) {
    function() {
      calls <<- calls + 1L
      list(result = stub_payload(unit, failing), metadata = list(n_dates = 720L))
    }
  }
  unit_one <- scalar_cache_path(cache_dir, "QA", "u01")
  first <- run_scalar_cached(unit_one, cache_key, "QA", "u01", compute_stub("u01"))
  first_hash <- sha256(unit_one)
  # The second invocation resumes unit u01, computes missing u02, and records
  # a rank-failure result for u03. A third invocation reuses all three.
  resumed <- run_scalar_cached(unit_one, cache_key, "QA", "u01", compute_stub("u01"))
  unit_two <- scalar_cache_path(cache_dir, "QA", "u02")
  new_unit <- run_scalar_cached(unit_two, cache_key, "QA", "u02", compute_stub("u02"))
  unit_three <- scalar_cache_path(cache_dir, "QA", "u03")
  failed_unit <- run_scalar_cached(unit_three, cache_key, "QA", "u03", compute_stub("u03", TRUE))
  test_units <- c("u01", "u02", "u03")
  test_paths <- list(u01 = unit_one, u02 = unit_two, u03 = unit_three)
  again <- lapply(test_units, function(u) run_scalar_cached(test_paths[[u]], cache_key, "QA", u, compute_stub(u)))
  changed_key <- hash_text(c(cache_key, "changed-input"))
  changed <- scalar_cache_path(file.path(cache_dir, substr(changed_key, 1L, 16L)), "QA", "u01")
  changed_run <- run_scalar_cached(changed, changed_key, "QA", "u01", compute_stub("u01"))
  mismatch_rejected <- inherits(try(read_scalar_cache(unit_one, changed_key, "QA", "u01"), silent = TRUE), "try-error")
  if (first$cache_reused || !resumed$cache_reused || new_unit$cache_reused ||
      !failed_unit$result$results$matched_unadjusted$failure$reason[1L] %in% "synthetic rank-deficiency audit" ||
      any(!vapply(again, function(x) x$cache_reused, logical(1))) || changed_run$cache_reused ||
      !mismatch_rejected || !identical(first_hash, sha256(unit_one)) || calls != 4L) {
    stop("scalar-cache interrupted-resume/source-fingerprint fixture failed")
  }
  list(status = "passed", fixture = "720-day three-outcome count fixture with dew3 spline adjustment",
    generated_calendar_dates = n, fitted_dates = length(f$retained_dates), stacked_rows = nrow(f$design),
    rank = f$audit$rank, columns = f$audit$columns,
    score_covariance_direct_vs_engine_max_abs_delta = max(abs(v_direct - v_engine)),
    covariance_formula = "H^-1 [U'U + sum_l w_l(U_t'U_{t-l}+transpose)] H^-1, Bartlett lag 14",
    contrasts = checks,
    scalar_cache_resume = list(status = "passed", completed_unit_reused = "u01", resumed_unit_computed = "u02",
      rank_failure_cached = "u03", changed_fingerprint_refit = TRUE, changed_key_mismatch_rejected = mismatch_rejected,
      completed_unit_sha256_unchanged = TRUE, compute_callbacks = calls, scalar_cache_payload_only = TRUE))
}

write_lock <- function() {
  lock <- list(source_sha256 = as.list(source_hashes), model_specification = model_spec,
    runner = basename(script_path), engine = basename(engine_path), plan = basename(plan_path),
    no_engine_or_global_helper_changes = TRUE, no_full_fit_rds_cache = TRUE)
  atomic_json(lock, model_lock_path)
  lock
}
verify_lock <- function() {
  if (!file.exists(model_lock_path)) stop("model lock missing; run --stage=validate before any effect fit")
  locked <- jsonlite::fromJSON(model_lock_path, simplifyVector = FALSE)
  old <- unlist(locked$source_sha256, use.names = TRUE)
  if (!identical(old[c("runner", "engine", "plan")], source_hashes[c("runner", "engine", "plan")])) {
    stop("model lock source hashes changed after validation")
  }
  invisible(locked)
}

status_path <- file.path(qa_dir, "run_status.json")
write_run_status <- function(status, details = list()) {
  atomic_json(c(list(status = status, updated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    run_id = run_id), details), status_path)
}
run_id <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")

if (opt$stage == "validate") {
  write_run_status("validating_model_and_synthetic_fixture", list(source_sha256 = as.list(source_hashes)))
  synthetic <- synthetic_validation()
  atomic_json(synthetic, file.path(qa_dir, "synthetic_covariance_contrast_validation.json"))
  write_lock()
  write_run_status("validation_passed_no_study_effects_fitted", list(
    model_lock_sha256 = sha256(model_lock_path), synthetic_validation_sha256 = sha256(file.path(qa_dir, "synthetic_covariance_contrast_validation.json"))))
  message(sprintf("Synthetic covariance/contrast validation passed; model lock SHA-256 %s", sha256(model_lock_path)))
  quit(save = "no", status = 0L)
}

verify_lock()
if (is.null(opt$source_qa) || !file.exists(opt$source_qa)) stop("--source-qa must name the completed source QA JSON before fitting")
qa_path <- normalizePath(opt$source_qa, winslash = "/", mustWork = TRUE)
source_qa <- jsonlite::fromJSON(qa_path, simplifyVector = FALSE)
if (!identical(tolower(as.character(source_qa$status)), "passed")) stop("moisture source QA status is not exactly 'passed'; no effect fits started")

data_dir <- file.path(root, "data")
primary_dir <- file.path(root, "results", "core")
mx_dir <- file.path(root, "results", "mx_updated")
source_paths <- c(
  BR_calendar = file.path(data_dir, "br_calendar.csv.gz"),
  MX_calendar = file.path(data_dir, "mx_updated_calendar.csv.gz"),
  BR_moisture = file.path(data_dir, "br_moisture_calendar.csv.gz"),
  MX_moisture = file.path(data_dir, "mx_moisture_calendar_2015_2022.csv.gz"),
  BR_primary = file.path(primary_dir, "primary_threeoutcome_E3_separate_BR_per_unit.csv"),
  BR_primary_failure = file.path(primary_dir, "primary_threeoutcome_E3_separate_BR_failure_audit.csv"),
  BR_primary_audit = file.path(primary_dir, "primary_threeoutcome_E3_separate_BR_fit_audit.csv"),
  BR_primary_manifest = file.path(root, "results", "fits", "primary_cache_manifest.csv"),
  MX_primary = file.path(mx_dir, "primary_complete_source_2015_2022_MX_per_unit.csv"),
  MX_primary_failure = file.path(mx_dir, "primary_complete_source_2015_2022_MX_failure_audit.csv"),
  MX_primary_audit = file.path(mx_dir, "primary_complete_source_2015_2022_MX_fit_audit.csv"),
  MX_primary_manifest = file.path(mx_dir, "fits", "fit_cache_manifest.csv"),
  MX_run_manifest = file.path(mx_dir, "run_manifest.json"),
  primary_core_runner = file.path(script_dir, "run_core.R"),
  original_analysis_plan = file.path(root, "admin", "analysis_plan.json"),
  source_QA = qa_path)
missing_paths <- source_paths[!file.exists(source_paths)]
if (length(missing_paths)) stopf("required input files are missing: %s", paste(names(missing_paths), collapse = ", "))

qa_flat <- unlist(source_qa, recursive = TRUE, use.names = TRUE)
for (nm in c("BR_moisture", "MX_moisture")) {
  base <- basename(source_paths[[nm]])
  idx <- grepl(base, names(qa_flat), fixed = TRUE)
  candidates <- unique(as.character(qa_flat[idx]))
  candidates <- candidates[grepl("^[0-9a-fA-F]{64}$", candidates)]
  if (length(candidates) != 1L || !identical(tolower(candidates), tolower(sha256(source_paths[[nm]])))) {
    stopf("source QA does not contain a unique matching SHA-256 for %s", base)
  }
}
input_hashes <- c(source_sha256 = vapply(source_paths, sha256, character(1)),
  runner_sha256 = source_hashes["runner"], engine_sha256 = source_hashes["engine"],
  plan_sha256 = source_hashes["plan"], model_lock_sha256 = sha256(model_lock_path))
cache_run_key <- hash_text(paste(names(input_hashes), unname(input_hashes), sep = "=\t"))
cache_dir <- file.path(out_dir, "scalar_cache", substr(cache_run_key, 1L, 16L))
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
existing_run_manifest <- file.path(out_dir, "run_manifest.json")
if (file.exists(existing_run_manifest)) {
  previous <- jsonlite::fromJSON(existing_run_manifest, simplifyVector = FALSE)
  if (identical(previous$status, "completed") && !identical(previous$cache_run_key, cache_run_key)) {
    stop("a completed moisture run exists for different input hashes; preserve it and use a new output directory")
  }
  if (identical(previous$status, "completed") && identical(previous$cache_run_key, cache_run_key)) {
    message(sprintf("This source-keyed moisture run is already complete (cache key %s); no files rewritten", substr(cache_run_key, 1L, 16L)))
    quit(save = "no", status = 0L)
  }
}
write_run_status("source_qa_passed_loading_inputs_before_preflight", list(
  input_sha256 = as.list(input_hashes), cache_run_key = cache_run_key, scalar_cache_dir = cache_dir,
  source_qa_path = qa_path, source_qa_sha256 = sha256(qa_path), workers_requested = opt$workers))

read_table <- function(path, select = NULL, character_cols = c("unit")) {
  common <- list(colClasses = list(character = character_cols), na.strings = c("", "NA", "NaN"),
    showProgress = FALSE, nThread = 1L, data.table = FALSE)
  if (!is.null(select)) common$select <- select
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    # `fread(file=*.gz)` delegates decompression to R.utils. Keep the locked
    # server library unchanged and stream gzip through the system gzip binary.
    common$cmd <- paste("gzip -cd", shQuote(path))
  } else common$file <- path
  do.call(data.table::fread, common)
}
outcomes <- c("y_K80", "y_K81", "y_urinary")
expected_contrasts <- c("y_K81_vs_y_K80", "y_urinary_vs_y_K80", "y_urinary_vs_y_K81")
read_country <- function(country) {
  if (country == "BR") {
    cal_path <- source_paths[["BR_calendar"]]; wx_path <- source_paths[["BR_moisture"]]
    cal_cols <- c("unit", "date", "tmean", "E3", "population", outcomes)
    wx_cols <- c("unit", "date", "tmean", "dew3", "rh3", "temperature_valid_mass")
    primary_path <- source_paths[["BR_primary"]]; failure_path <- source_paths[["BR_primary_failure"]]
    audit_path <- source_paths[["BR_primary_audit"]]; manifest_path <- source_paths[["BR_primary_manifest"]]
    offset_col <- "population"
  } else {
    cal_path <- source_paths[["MX_calendar"]]; wx_path <- source_paths[["MX_moisture"]]
    cal_cols <- c("unit", "date", "tmean", "E3", outcomes)
    wx_cols <- c("unit", "date", "tmean", "dew3", "rh3")
    primary_path <- source_paths[["MX_primary"]]; failure_path <- source_paths[["MX_primary_failure"]]
    audit_path <- source_paths[["MX_primary_audit"]]; manifest_path <- source_paths[["MX_primary_manifest"]]
    offset_col <- NULL
  }
  cal <- read_table(cal_path, cal_cols)
  wx <- read_table(wx_path, wx_cols)
  for (d in list(cal, wx)) {
    if (anyNA(d$unit) || any(!nzchar(d$unit))) stopf("%s input has empty unit IDs", country)
    d$date <- as.Date(d$date)
    if (anyNA(d$date)) stopf("%s input has invalid dates", country)
  }
  cal$date <- as.Date(cal$date); wx$date <- as.Date(wx$date)
  if (anyDuplicated(cal[c("unit", "date")]) || anyDuplicated(wx[c("unit", "date")])) stopf("%s input has duplicate unit-date keys", country)
  if (any(is.finite(cal$E3) & cal$E3 < 0)) stopf("%s E3 contains invalid negative values", country)
  if (any(!is.finite(wx$dew3) & !is.na(wx$dew3)) || any(!is.finite(wx$rh3) & !is.na(wx$rh3))) stopf("%s moisture columns contain infinite values", country)
  if (any(is.finite(wx$rh3) & (wx$rh3 < 0 | wx$rh3 > 100))) stopf("%s rh3 is outside [0,100]", country)
  if (country == "BR" && any(is.finite(wx$temperature_valid_mass) &
      (wx$temperature_valid_mass < -1e-6 | wx$temperature_valid_mass > 1 + 1e-6))) {
    stop("BR temperature_valid_mass falls outside [0,1] beyond 1e-6 rounding tolerance")
  }
  key <- function(d) paste(d$unit, as.character(d$date), sep = "\r")
  idx <- match(key(cal), key(wx))
  if (anyNA(idx)) stopf("%s moisture calendar lacks %d canonical unit-date rows", country, sum(is.na(idx)))
  canonical_t <- as.numeric(cal$tmean); moisture_t <- as.numeric(wx$tmean[idx])
  if (!identical(is.finite(canonical_t), is.finite(moisture_t))) stopf("%s canonical and moisture temperature masks differ", country)
  tdiff <- abs(canonical_t[is.finite(canonical_t)] - moisture_t[is.finite(moisture_t)])
  max_tdiff <- if (length(tdiff)) max(tdiff) else NA_real_
  if (is.finite(max_tdiff) && max_tdiff > 1e-6) stopf("%s moisture temperature differs from canonical tmean by %.9g C (tolerance 1e-6)", country, max_tdiff)
  d <- cal
  d$dew3 <- as.numeric(wx$dew3[idx]); d$rh3 <- as.numeric(wx$rh3[idx])
  if (country == "BR") d$temperature_valid_mass <- as.numeric(wx$temperature_valid_mass[idx])
  d$time_years <- as.numeric(d$date - min(d$date)) / 365.25
  primary <- read_table(primary_path)
  primary_failure <- read_table(failure_path)
  primary_audit <- read_table(audit_path)
  primary_manifest <- read_table(manifest_path)
  run_manifest <- if (country == "MX") jsonlite::fromJSON(source_paths[["MX_run_manifest"]], simplifyVector = FALSE) else NULL
  if (nrow(primary_failure) && any(primary_failure$scenario == if (country == "BR") "primary_threeoutcome_E3_separate" else "primary_complete_source_2015_2022")) {
    # Failures are retained by unit below; their presence does not invalidate the other primary units.
  }
  if (any(primary$covariance == "hac14" & !primary$contrast %in% expected_contrasts)) stopf("%s primary HAC14 has unexpected contrast keys", country)
  list(data = d, primary = primary, primary_failure = primary_failure, primary_audit = primary_audit,
    primary_manifest = primary_manifest, run_manifest = run_manifest, offset_col = offset_col,
    moisture_tmean_max_abs_diff_c = max_tdiff)
}

md5 <- function(path) unname(tools::md5sum(path)[1L])
validate_primary_provenance <- function(country, obj) {
  pri <- obj$primary[obj$primary$covariance == "hac14", , drop = FALSE]
  keys <- c("unit", "contrast", "covariance")
  if (anyDuplicated(pri[keys])) stopf("%s primary HAC14 output has duplicate unit/contrast keys", country)
  primary_units <- sort(unique(pri$unit))
  manifest <- obj$primary_manifest
  expected_scenario <- if (country == "BR") "primary_threeoutcome_E3_separate" else "primary_complete_source_2015_2022"
  manifest <- manifest[manifest$scenario == expected_scenario, , drop = FALSE]
  if (anyDuplicated(manifest$unit)) stopf("%s primary cache manifest has duplicate unit keys", country)
  manifest <- manifest[manifest$unit %in% primary_units, , drop = FALSE]
  if (!setequal(sort(unique(manifest$unit)), primary_units)) stopf("%s primary output and cache manifest unit keys differ", country)
  if (country == "BR") {
    protocol_files <- c(engine_path, source_paths[["primary_core_runner"]], source_paths[["original_analysis_plan"]])
    protocol_hashes <- tools::md5sum(protocol_files)
    names(protocol_hashes) <- basename(names(protocol_hashes))
    expected_protocol <- paste(c(R.version.string,
      paste(paste(names(protocol_hashes), unname(protocol_hashes), sep = "="), collapse = ";")), collapse = "|")
    if (any(tolower(manifest$input_md5) != tolower(md5(source_paths[["BR_calendar"]]))) ||
        any(manifest$protocol_signature != expected_protocol)) {
      stop("BR primary cache manifest does not bind the expected calendar, frozen engine, core runner and analysis plan")
    }
  } else {
    if (any(tolower(manifest$input_sha256) != tolower(sha256(source_paths[["MX_calendar"]]))) ||
        any(tolower(manifest$engine_sha256) != tolower(sha256(engine_path)))) {
      stop("MX primary fit manifest does not bind the expected canonical calendar and immutable engine")
    }
    s <- obj$run_manifest$scenarios$primary_complete_source_2015_2022
    if (is.null(s$output_sha256$per_unit) || !identical(tolower(s$output_sha256$per_unit), tolower(sha256(source_paths[["MX_primary"]])))) {
      stop("MX primary run manifest does not bind the current per-unit result table")
    }
  }
  for (unit in primary_units) {
    m <- manifest[manifest$unit == unit, , drop = FALSE]
    p <- pri[pri$unit == unit, , drop = FALSE]
    if (nrow(m) != 1L || nrow(p) != 3L || anyDuplicated(p$contrast) || !setequal(p$contrast, expected_contrasts)) {
      stopf("%s/%s primary provenance keys are incomplete", country, unit)
    }
    if (any(as.integer(p$n_complete) != as.integer(m$n_complete)) ||
        any(!vapply(seq_len(nrow(p)), function(i) near(p$p95[i], m$p95[1L], 1e-9) && near(p$median[i], m$median[1L], 1e-9), logical(1)))) {
      stopf("%s/%s primary provenance dates or temperature targets disagree with its fit manifest", country, unit)
    }
    if (all(c("dB1", "dB2", "dB3") %in% names(m)) && all(c("dB1", "dB2", "dB3") %in% names(p)) &&
        any(!vapply(seq_len(nrow(p)), function(i) near(c(p$dB1[i], p$dB2[i], p$dB3[i]),
          as.numeric(m[1L, c("dB1", "dB2", "dB3")]), 1e-8), logical(1)))) {
      stopf("%s/%s primary fitted temperature basis contrast disagrees with its fit manifest", country, unit)
    }
  }
  invisible(TRUE)
}

baseline_eligible <- function(d, offset_col) {
  keep <- is.finite(d$E3) & !is.na(d$date)
  for (nm in outcomes) keep <- keep & is.finite(d[[nm]]) & d[[nm]] >= 0
  if (!is.null(offset_col)) keep <- keep & is.finite(d[[offset_col]]) & d[[offset_col]] > 0
  list(keep = keep, n = sum(keep), events = vapply(outcomes, function(nm) sum(d[[nm]][keep]), numeric(1)))
}
make_targets <- function(d, offset_col) {
  elig <- baseline_eligible(d, offset_col)
  x <- as.numeric(d$E3[elig$keep])
  if (length(x) < 4L || length(unique(x)) < 2L) stop("full original paired calendar cannot define a three-df temperature basis")
  basis <- new_temperature_basis(x, df = 3L)
  list(basis = basis, p95 = unname(stats::quantile(x, 0.95, type = 7)), median = stats::median(x),
    dB = as.numeric(.basis_at(basis, unname(stats::quantile(x, 0.95, type = 7))) - .basis_at(basis, stats::median(x))),
    n_complete = elig$n, events = elig$events, eligible = elig$keep)
}
empty_estimates <- function() data.frame(country = character(), unit = character(), scenario = character(),
  exposure = character(), nuisance = character(), reference = character(), contrast = character(), covariance = character(),
  estimate = numeric(), variance = numeric(), se = numeric(), RR = numeric(), lo = numeric(), hi = numeric(),
  n_complete = integer(), events_first = numeric(), events_second = numeric(), events_third = numeric(),
  p95 = numeric(), median = numeric(), dB1 = numeric(), dB2 = numeric(), dB3 = numeric(),
  target_n_complete = integer(), sample_type = character(), stringsAsFactors = FALSE)
empty_failures <- function() data.frame(country = character(), unit = character(), scenario = character(),
  phase = character(), reason = character(), n_complete = integer(), events = character(), stringsAsFactors = FALSE)
empty_audit <- function() data.frame(country = character(), unit = character(), scenario = character(),
  converged = logical(), boundary = logical(), rank = integer(), columns = integer(), n = integer(), events = numeric(),
  warnings = character(), formula = character(), moisture_basis = character(), audit_source = character(), stringsAsFactors = FALSE)
failure_row <- function(country, unit, scenario, phase, reason, n = NA_integer_, events = "")
  data.frame(country = country, unit = unit, scenario = scenario, phase = phase, reason = reason,
    n_complete = as.integer(n), events = events, stringsAsFactors = FALSE)
audit_row <- function(country, unit, scenario, a, source = "new_fit", moisture_basis = "none")
  data.frame(country = country, unit = unit, scenario = scenario, converged = isTRUE(a$converged),
    boundary = isTRUE(a$boundary), rank = as.integer(a$rank), columns = as.integer(a$columns),
    n = as.integer(a$n), events = as.numeric(a$events), warnings = paste(a$warnings, collapse = " | "),
    formula = a$formula, moisture_basis = moisture_basis, audit_source = source, stringsAsFactors = FALSE)
event_text <- function(events) paste(paste(names(events), unname(events), sep = "="), collapse = ";")

fit_unit_scenario <- function(country, unit, d, target, offset_col, scenario, moisture_col = NULL, sample_type = "matched") {
  tryCatch({
    f <- fit_stacked_moisture(d, outcomes, target$basis, target$p95, target$median,
      offset_col = offset_col, moisture_col = moisture_col)
    ev <- vapply(outcomes, function(nm) sum(d[[nm]]), numeric(1))
    rows <- lapply(f$contrasts, function(z) data.frame(country = country, unit = unit, scenario = scenario,
      exposure = "E3", nuisance = "separate", reference = if (country == "MX") "full_2015_2022" else "full_2006_2020",
      contrast = z$label, covariance = "hac14", estimate = z$estimate, variance = z$variance, se = z$se,
      RR = z$RR, lo = z$lo, hi = z$hi, n_complete = nrow(d), events_first = ev[1L], events_second = ev[2L],
      events_third = ev[3L], p95 = f$p95, median = f$median, dB1 = f$dB[1L], dB2 = f$dB[2L], dB3 = f$dB[3L],
      target_n_complete = target$n_complete, sample_type = sample_type, stringsAsFactors = FALSE))
    list(estimates = do.call(rbind, rows), failure = empty_failures(),
      audit = audit_row(country, unit, scenario, f$audit, moisture_basis = if (is.null(moisture_col)) "none" else paste0("splines::ns(", moisture_col, ", df=3)",
        "; fitted separately by outcome; values held fixed for temperature contrast")))
  }, error = function(e) list(estimates = empty_estimates(),
    failure = failure_row(country, unit, scenario, "fit", conditionMessage(e), nrow(d), event_text(vapply(outcomes, function(nm) sum(d[[nm]]), numeric(1)))),
    audit = empty_audit(), details = NULL))
}

reuse_primary_unit <- function(country, unit, matched_dates, pair_dates, target, event_counts, primary, primary_failure, primary_audit) {
  # Reuse is allowed only for exact dates, exact counts/events/targets, unique
  # three-row HAC14 keys, and an audited full-rank converged primary fit.
  if (!identical(as.character(matched_dates), as.character(pair_dates))) return(NULL)
  old_fail <- primary_failure[primary_failure$unit == unit, , drop = FALSE]
  old <- primary[primary$unit == unit & primary$covariance == "hac14", , drop = FALSE]
  if (nrow(old_fail) || !nrow(old)) return(NULL)
  if (anyDuplicated(old[c("unit", "contrast", "covariance")]) || nrow(old) != 3L ||
      !setequal(old$contrast, expected_contrasts)) stopf("%s/%s primary reuse keys/count failed", country, unit)
  expected_reference <- if (country == "MX") "full_2015_2022" else "full"
  if (any(old$country != country) || any(old$unit != unit) ||
      any(old$scenario != if (country == "MX") "primary_complete_source_2015_2022" else "primary_threeoutcome_E3_separate") ||
      any(old$exposure != "E3") || any(old$nuisance != "separate") || any(old$reference != expected_reference)) {
    stopf("%s/%s primary reuse model metadata do not match the locked primary setting", country, unit)
  }
  if (any(as.integer(old$n_complete) != length(matched_dates)) || any(!vapply(seq_len(nrow(old)), function(i)
      near(old$p95[i], target$p95, 1e-9) && near(old$median[i], target$median, 1e-9), logical(1)))) {
    stopf("%s/%s primary reuse date-count/temperature-target assertion failed", country, unit)
  }
  if (any(as.numeric(old$events_first) != event_counts[1L] | as.numeric(old$events_second) != event_counts[2L] |
          as.numeric(old$events_third) != event_counts[3L])) stopf("%s/%s primary reuse event assertion failed", country, unit)
  if (all(c("dB1", "dB2", "dB3") %in% names(old)) && any(!vapply(seq_len(nrow(old)), function(i)
      near(c(old$dB1[i], old$dB2[i], old$dB3[i]), target$dB, 1e-8), logical(1)))) {
    stopf("%s/%s primary reuse temperature basis contrast assertion failed", country, unit)
  }
  au <- primary_audit[primary_audit$unit == unit, , drop = FALSE]
  if (nrow(au) != 1L || !isTRUE(as.logical(au$converged[1L])) || as.integer(au$rank[1L]) != as.integer(au$columns[1L])) {
    stopf("%s/%s primary reuse fit audit is missing, nonconverged or rank deficient", country, unit)
  }
  expected_formula <- if (country == "BR")
    "y ~ disease * (B1 + B2 + B3 + dow + month + splines::ns(time, 12)) + offset(off)" else
      "y ~ disease * (B1 + B2 + B3 + dow + month + splines::ns(time, 12))"
  normalize_formula <- function(x) gsub("[[:space:]]+", " ", trimws(as.character(x)))
  if (!identical(as.integer(au$n[1L]), as.integer(length(matched_dates) * length(outcomes))) ||
      !near(as.numeric(au$events[1L]), sum(event_counts), 0) ||
      !identical(normalize_formula(au$formula[1L]), normalize_formula(expected_formula))) {
    stopf("%s/%s primary reuse audit rows/events/formula do not match the locked model", country, unit)
  }
  old$scenario <- "matched_unadjusted"
  old$sample_type <- "original_primary_dates_verified_exact_match"
  old$dB1 <- if ("dB1" %in% names(old)) old$dB1 else target$dB[1L]
  old$dB2 <- if ("dB2" %in% names(old)) old$dB2 else target$dB[2L]
  old$dB3 <- if ("dB3" %in% names(old)) old$dB3 else target$dB[3L]
  old$target_n_complete <- target$n_complete
  au$scenario <- "matched_unadjusted"; au$audit_source <- "reused_primary_hac14"
  au$moisture_basis <- "none"
  list(estimates = old, failure = empty_failures(), audit = au)
}

process_unit <- function(country, unit, country_data, offset_col, primary, primary_failure, primary_audit) {
  d0 <- country_data[country_data$unit == unit, , drop = FALSE]
  target <- tryCatch(make_targets(d0, offset_col), error = function(e) e)
  scen_names <- c("matched_unadjusted", "matched_dew3_ns3", "matched_rh3_ns3")
  if (inherits(target, "error")) {
    return(list(results = setNames(lapply(scen_names, function(s) list(estimates = empty_estimates(),
      failure = failure_row(country, unit, s, "reference", conditionMessage(target)), audit = empty_audit())), scen_names),
      matched_n = 0L, common = FALSE))
  }
  pair <- baseline_eligible(d0, offset_col)
  matched <- pair$keep & is.finite(d0$dew3) & is.finite(d0$rh3)
  d <- d0[matched, , drop = FALSE]
  counts <- vapply(outcomes, function(nm) sum(d[[nm]], na.rm = TRUE), numeric(1))
  if (sum(matched) < 250L || any(counts < 40)) {
    msg <- sprintf("matched complete dates=%d; events=%s; required >=250 dates and >=40 events per outcome",
      sum(matched), event_text(counts))
    out <- setNames(lapply(scen_names, function(s) list(estimates = empty_estimates(),
      failure = failure_row(country, unit, s, "eligibility", msg, sum(matched), event_text(counts)), audit = empty_audit())), scen_names)
    return(list(results = out, matched_n = sum(matched), common = FALSE))
  }
  dates_pair <- as.character(d0$date[pair$keep]); dates_matched <- as.character(d$date)
  unadj <- reuse_primary_unit(country, unit, dates_matched, dates_pair, target, counts,
    primary, primary_failure, primary_audit)
  if (is.null(unadj)) {
    same_dates <- identical(dates_matched, dates_pair)
    old_fail <- primary_failure[primary_failure$unit == unit, , drop = FALSE]
    if (same_dates && nrow(old_fail)) {
      unadj <- list(estimates = empty_estimates(),
        failure = failure_row(country, unit, "matched_unadjusted", "inherited_primary_failure",
          paste(old_fail$phase[1L], old_fail$reason[1L], sep = ": "), sum(matched), event_text(counts)), audit = empty_audit())
    } else {
      unadj <- fit_unit_scenario(country, unit, d, target, offset_col, "matched_unadjusted", NULL,
        if (same_dates) "original_dates_fit_due_primary_row_unavailable" else "moisture_matched_baseline_refit")
    }
  }
  dew <- fit_unit_scenario(country, unit, d, target, offset_col, "matched_dew3_ns3", "dew3")
  rh <- fit_unit_scenario(country, unit, d, target, offset_col, "matched_rh3_ns3", "rh3")
  out <- list(matched_unadjusted = unadj, matched_dew3_ns3 = dew, matched_rh3_ns3 = rh)
  common <- all(vapply(out, function(z) nrow(z$estimates) == 3L && all(z$estimates$covariance == "hac14"), logical(1)))
  list(results = out, matched_n = sum(matched), common = common)
}

run_unit_cached <- function(country, unit, country_data, offset_col, primary, primary_failure,
                            primary_audit, cache_dir, cache_key, input_hashes) {
  path <- scalar_cache_path(cache_dir, country, unit)
  run_scalar_cached(path, cache_key, country, unit, function() {
    result <- process_unit(country, unit, country_data, offset_col, primary, primary_failure, primary_audit)
    d <- country_data[country_data$unit == unit, , drop = FALSE]
    baseline <- baseline_eligible(d, offset_col)
    matched <- baseline$keep & is.finite(d$dew3) & is.finite(d$rh3)
    target <- tryCatch(make_targets(d, offset_col), error = function(e) e)
    metadata <- list(input_sha256 = as.list(input_hashes), baseline_eligible_n = baseline$n,
      matched_n = sum(matched), baseline_eligible_dates_sha256 = hash_text(as.character(d$date[baseline$keep])),
      matched_dates_sha256 = hash_text(as.character(d$date[matched])),
      p95 = if (inherits(target, "error")) NA_real_ else target$p95,
      median = if (inherits(target, "error")) NA_real_ else target$median,
      temperature_dB = if (inherits(target, "error")) NULL else as.list(target$dB),
      common_successful_unit = result$common)
    list(result = result, metadata = metadata)
  })
}

pool_rows <- function(rows, country, scenario, scope) {
  rows <- rows[rows$covariance == "hac14", , drop = FALSE]
  out <- list()
  for (contrast in expected_contrasts) {
    g <- rows[rows$contrast == contrast, , drop = FALSE]
    if (nrow(g) && anyDuplicated(g$unit)) stopf("duplicate pooled unit key for %s/%s/%s", country, scenario, contrast)
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(as.numeric(g$estimate), as.numeric(g$variance), method = method)
      out[[length(out) + 1L]] <- data.frame(country = country, scenario = scenario, sample_scope = scope,
        contrast = contrast, covariance = "hac14", method = method, k = p$k,
        estimate = p$estimate, se = p$se, RR = p$RR, lo = p$RR95CI[1L], hi = p$RR95CI[2L],
        I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
        failure = if (is.null(p$failure)) "" else p$failure, unit_ids = paste(sort(unique(g$unit)), collapse = ";"),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

write_run_status("preflight_pending", list(input_sha256 = as.list(input_hashes), workers_requested = opt$workers))
countries <- list(BR = read_country("BR"), MX = read_country("MX"))
validate_primary_provenance("BR", countries$BR)
validate_primary_provenance("MX", countries$MX)
workers <- opt$workers
run_workers <- min(workers, max(1L, parallel::detectCores(logical = TRUE) - 1L))
preflight_results <- list()
preflight_attempts <- list()
for (country in names(countries)) {
  obj <- countries[[country]]
  units <- sort(unique(obj$data$unit))
  eligible_units <- units[vapply(units, function(u) {
    d <- obj$data[obj$data$unit == u, , drop = FALSE]
    pair <- baseline_eligible(d, obj$offset_col)
    k <- pair$keep & is.finite(d$dew3) & is.finite(d$rh3)
    sum(k) >= 250L && all(vapply(outcomes, function(nm) sum(d[[nm]][k]) >= 40, logical(1)))
  }, logical(1))]
  if (!length(eligible_units)) stopf("%s has no unit eligible for matched-model preflight", country)
  passed <- FALSE
  attempted <- list()
  for (unit in eligible_units) {
    record <- run_unit_cached(country, unit, obj$data, obj$offset_col, obj$primary,
      obj$primary_failure, obj$primary_audit, cache_dir, cache_run_key, input_hashes)
    attempted[[unit]] <- record
    if (isTRUE(record$result$common)) {
      preflight_results[[country]] <- list(unit = unit, record = record)
      passed <- TRUE
      break
    }
  }
  preflight_attempts[[country]] <- attempted
  if (!passed) {
    details <- lapply(attempted, function(record) lapply(record$result$results, function(v) v$failure$reason))
    atomic_json(list(status = "failed", country = country, attempts = names(attempted), failures = details,
      source_qa_sha256 = sha256(qa_path), input_sha256 = as.list(input_hashes)), file.path(qa_dir, "actual_data_preflight.json"))
    write_run_status("preflight_failed_no_full_fit_started", list(country = country, attempted_units = names(attempted),
      cache_run_key = cache_run_key, input_sha256 = as.list(input_hashes), failure_reasons = details))
    stopf("%s had no successful actual-data preflight unit among %d deterministic eligible candidates", country, length(attempted))
  }
}
atomic_json(list(status = "passed", created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  source_qa_sha256 = sha256(qa_path), units = lapply(preflight_results, function(x) x$unit),
  attempted_units = lapply(preflight_attempts, names),
  failed_candidates = lapply(preflight_attempts, function(xs) names(xs)[!vapply(xs, function(z) isTRUE(z$result$common), logical(1))]),
  moisture_tmean_max_abs_diff_c = lapply(countries, `[[`, "moisture_tmean_max_abs_diff_c"),
  settings = model_spec$settings, input_sha256 = as.list(input_hashes)), file.path(qa_dir, "actual_data_preflight.json"))
write_run_status("preflight_passed_fitting_bounded_scenarios", list(
  preflight_units = lapply(preflight_results, function(x) x$unit), attempted_units = lapply(preflight_attempts, names),
  cache_run_key = cache_run_key, scalar_cache_dir = cache_dir,
  workers_requested = opt$workers, workers_used = run_workers))

country_results <- list()
for (country in names(countries)) {
  message(sprintf("%s moisture sensitivity: %d units; %d workers", country,
    length(unique(countries[[country]]$data$unit)), run_workers))
  obj <- countries[[country]]
  units <- sort(unique(obj$data$unit))
  pre_u <- preflight_results[[country]]$unit
  remaining <- setdiff(units, pre_u)
  worker <- function(unit) run_unit_cached(country, unit, obj$data, obj$offset_col,
    obj$primary, obj$primary_failure, obj$primary_audit, cache_dir, cache_run_key, input_hashes)
  fit_records <- list()
  chunks <- split(remaining, ceiling(seq_along(remaining) / 25L))
  completed <- 1L
  for (chunk in seq_along(chunks)) {
    this <- chunks[[chunk]]
    part <- if (.Platform$OS.type == "windows" || run_workers <= 1L || length(this) <= 1L)
      lapply(this, worker) else parallel::mclapply(this, worker,
        mc.cores = min(run_workers, length(this)), mc.preschedule = FALSE, mc.set.seed = FALSE)
    names(part) <- this
    fit_records <- c(fit_records, part)
    completed <- completed + length(this)
    message(sprintf("%s checkpoint: %d/%d unit scalar caches complete", country, completed, length(units)))
    write_run_status(paste0("fitting_", country, "_checkpoint"), list(country = country,
      completed_units = completed, total_units = length(units), cache_run_key = cache_run_key,
      scalar_cache_dir = cache_dir, workers_used = run_workers))
  }
  fit_records[[pre_u]] <- preflight_results[[country]]$record
  fit_records <- fit_records[units]
  cache_manifest_rows <- lapply(names(fit_records), function(u) {
    rec <- fit_records[[u]]
    data.frame(country = country, unit = u, cache_path = normalizePath(rec$cache_path, winslash = "/", mustWork = TRUE),
      cache_sha256 = sha256(rec$cache_path), cache_run_key = cache_run_key, reused_on_this_run = rec$cache_reused,
      cache_bytes = file.info(rec$cache_path)$size, matched_n = rec$metadata$matched_n,
      matched_dates_sha256 = rec$metadata$matched_dates_sha256,
      baseline_eligible_dates_sha256 = rec$metadata$baseline_eligible_dates_sha256,
      p95 = rec$metadata$p95, median = rec$metadata$median, stringsAsFactors = FALSE)
  })
  cache_manifest_rows <- do.call(rbind, cache_manifest_rows)
  all_unit_results <- lapply(fit_records, `[[`, "result")
  by_scenario <- setNames(vector("list", 3L), model_spec$settings)
  for (scenario in model_spec$settings) {
    est <- do.call(rbind, lapply(all_unit_results, function(z) z$results[[scenario]]$estimates))
    fail <- do.call(rbind, lapply(all_unit_results, function(z) z$results[[scenario]]$failure))
    aud <- do.call(rbind, lapply(all_unit_results, function(z) z$results[[scenario]]$audit))
    if (is.null(est)) est <- empty_estimates()
    if (is.null(fail)) fail <- empty_failures()
    if (is.null(aud)) aud <- empty_audit()
    rownames(est) <- rownames(fail) <- rownames(aud) <- NULL
    if (nrow(est) && (anyDuplicated(est[c("unit", "contrast", "covariance")]) ||
        any(est$covariance != "hac14") || any(!est$contrast %in% expected_contrasts))) {
      stopf("%s/%s emitted duplicate or unexpected result keys", country, scenario)
    }
    by_scenario[[scenario]] <- list(estimates = est, failures = fail, audits = aud,
      pooled = pool_rows(est, country, scenario, "all_successful_units"))
  }
  common_units <- Reduce(intersect, lapply(by_scenario, function(z) unique(z$estimates$unit)))
  common_pool <- do.call(rbind, lapply(model_spec$settings, function(s) {
    z <- by_scenario[[s]]$estimates
    pool_rows(z[z$unit %in% common_units, , drop = FALSE], country, s, "common_successful_units")
  }))
  membership <- data.frame(country = country, unit = units, matched_common_success = units %in% common_units,
    matched_dates = vapply(all_unit_results[units], function(z) z$matched_n, integer(1)), stringsAsFactors = FALSE)
  for (s in model_spec$settings) membership[[paste0(s, "_success")]] <- units %in% unique(by_scenario[[s]]$estimates$unit)
  country_results[[country]] <- list(by_scenario = by_scenario, common_pool = common_pool, membership = membership,
    common_units = common_units, all_units = units, scalar_cache_manifest = cache_manifest_rows)
}

# Coverage analysis recalculates country pooling with newly derived weights and
# heterogeneity after excluding the three plan-identified low-coverage units.
coverage_rows <- list()
for (scenario in c("primary_original", model_spec$settings)) {
  for (country in names(country_results)) {
    if (scenario == "primary_original") {
      if (country == "BR") d <- countries$BR$primary[countries$BR$primary$covariance == "hac14", , drop = FALSE] else
        d <- countries$MX$primary[countries$MX$primary$covariance == "hac14", , drop = FALSE]
      if (!nrow(d)) next
      d$scenario <- "primary_original"
      # The original source row order and keys are validated here for the
      # documented geographic exclusion re-pool.
      if (anyDuplicated(d[c("unit", "contrast", "covariance")])) stopf("%s primary HAC14 keys duplicate", country)
      d$sample_type <- "original_primary"
    } else d <- country_results[[country]]$by_scenario[[scenario]]$estimates
    for (scope in c("all_successful_units", "exclude_BR_low_coverage")) {
      if (scope == "exclude_BR_low_coverage" && country != "BR") next
      keep <- if (scope == "exclude_BR_low_coverage") !d$unit %in% model_spec$coverage_exclusion_BR else rep(TRUE, nrow(d))
      coverage_rows[[length(coverage_rows) + 1L]] <- pool_rows(d[keep, , drop = FALSE], country, scenario, scope)
    }
  }
}
coverage_pools <- do.call(rbind, coverage_rows)

# The source companion is retained as a QC table; the exclusion list comes
# from the approved plan and is checked against it rather than refitted.
temp_mass <- countries$BR$data[, c("unit", "temperature_valid_mass"), drop = FALSE]
temp_groups <- split(temp_mass$temperature_valid_mass, temp_mass$unit)
qc_expanded <- do.call(rbind, lapply(names(temp_groups), function(u) {
  x <- as.numeric(temp_groups[[u]])
  data.frame(unit = u, n_finite = sum(is.finite(x)),
    minimum = if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_,
    mean = if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_,
    maximum = if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)
}))
qc_expanded$low_coverage_exclusion <- qc_expanded$unit %in% model_spec$coverage_exclusion_BR

all_pooled <- do.call(rbind, lapply(names(country_results), function(country) {
  z <- country_results[[country]]
  do.call(rbind, c(lapply(model_spec$settings, function(s) z$by_scenario[[s]]$pooled), list(z$common_pool)))
}))
all_membership <- do.call(rbind, lapply(country_results, `[[`, "membership"))
all_cache_manifest <- do.call(rbind, lapply(country_results, `[[`, "scalar_cache_manifest"))
manifest <- list(run_id = run_id, status = "completed", completed_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  R = R.version.string, workers_requested = opt$workers, workers_used = run_workers,
  cache_run_key = cache_run_key, scalar_cache_dir = cache_dir,
  scalar_cache_policy = "Source-hash keyed atomic per-unit scalar result and failure payloads only; no model/design/covariance/fit RDS objects.",
  native_thread_environment = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"))),
  data_table_threads = data.table::getDTthreads(), source_sha256 = as.list(input_hashes),
  source_QA_path = qa_path, source_QA_status = source_qa$status, model_specification = model_spec,
  preflight_units = lapply(preflight_results, function(x) x$unit),
  units = lapply(country_results, function(x) list(expected = length(x$all_units),
    common_successful = length(x$common_units), fitted_by_scenario = lapply(x$by_scenario, function(z) length(unique(z$estimates$unit))),
    failed_by_scenario = lapply(x$by_scenario, function(z) nrow(z$failures)))))

for (country in names(country_results)) for (scenario in model_spec$settings) {
  z <- country_results[[country]]$by_scenario[[scenario]]
  atomic_csv(z$estimates, file.path(out_dir, paste0(scenario, "_", country, "_per_unit.csv")))
  atomic_csv(z$failures, file.path(out_dir, paste0(scenario, "_", country, "_failure_audit.csv")))
  atomic_csv(z$audits, file.path(out_dir, paste0(scenario, "_", country, "_fit_audit.csv")))
}
atomic_csv(all_pooled, file.path(out_dir, "country_pooled_DL_REML_HK.csv"))
atomic_csv(coverage_pools, file.path(out_dir, "coverage_exclusion_pooled_DL_REML_HK.csv"))
atomic_csv(all_membership, file.path(out_dir, "common_successful_unit_comparison.csv"))
atomic_csv(qc_expanded, file.path(out_dir, "brazil_temperature_valid_mass_qc.csv"))
atomic_csv(all_cache_manifest, file.path(out_dir, "scalar_cache_manifest.csv"))
outputs <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
outputs <- outputs[!basename(outputs) %in% c("run_manifest.json")]
manifest$output_sha256 <- as.list(setNames(vapply(outputs, sha256, character(1)), basename(outputs)))
manifest$qa_sha256 <- list(source_QA = sha256(qa_path), synthetic_validation = sha256(file.path(qa_dir, "synthetic_covariance_contrast_validation.json")),
  actual_data_preflight = sha256(file.path(qa_dir, "actual_data_preflight.json")), model_lock = sha256(model_lock_path))
atomic_json(manifest, file.path(out_dir, "run_manifest.json"))
write_run_status("completed", list(run_manifest_sha256 = sha256(file.path(out_dir, "run_manifest.json")),
  cache_run_key = cache_run_key, scalar_cache_manifest_sha256 = sha256(file.path(out_dir, "scalar_cache_manifest.csv")),
  units = manifest$units, output_sha256 = manifest$output_sha256))
message(sprintf("Moisture sensitivity complete: outputs in %s", out_dir))
