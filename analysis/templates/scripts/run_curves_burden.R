# Downstream repaired main-effect, curve and standardized-burden analyses.
# Usage: Rscript run_curves_burden.R main [BR|MX|all]
#        Rscript run_curves_burden.R curves_burden [BR|MX|all]
# This script consumes the primary full-covariance fit cache from run_core.R.
# It writes only aggregate unit/country outputs; no patient-level records are saved.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
options(stringsAsFactors = FALSE, scipen = 999)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
if (is.na(script_file) || !nzchar(script_file)) stop("could not locate this R script", call. = FALSE)
script_dir <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "repair_engine.R"), local = FALSE)
if (!requireNamespace("splines", quietly = TRUE)) stop("splines namespace is required for cached natural-spline predictions", call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1L] %in% c("main", "curves_burden")) {
  stop("usage: Rscript run_curves_burden.R main|curves_burden [BR|MX|all]", call. = FALSE)
}
stage <- args[1L]
country_arg <- if (length(args) >= 2L) toupper(args[2L]) else "ALL"
if (!country_arg %in% c("ALL", "BR", "MX")) stop("country must be BR, MX, or all", call. = FALSE)
countries <- if (country_arg == "ALL") c("BR", "MX") else country_arg
run_id <- format(Sys.time(), "%Y%m%dT%H%M%S")
out_dir <- file.path(project_dir, "results", paste0("downstream_", stage, "_", run_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fit_dir <- file.path(project_dir, "results", "fits")

append_csv <- function(row, path) {
  if (is.null(row) || !nrow(row)) return(invisible(NULL))
  write.table(row, file = path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path),
              na = "", qmethod = "double")
  invisible(NULL)
}

failure_path <- file.path(out_dir, "failure_audit.csv")
audit_failure <- function(country, unit = NA_character_, analysis, endpoint = NA_character_,
                          reason, n_days = NA_integer_, events = NA_real_) {
  append_csv(data.frame(stage = stage, country = country, unit = unit,
                        analysis = analysis, endpoint = endpoint,
                        n_days = n_days, events = events,
                        reason = as.character(reason), run_id = run_id,
                        stringsAsFactors = FALSE), failure_path)
}

calendar_path <- function(country) file.path(project_dir, "data", paste0(tolower(country), "_calendar.csv.gz"))
read_calendar <- function(country) {
  p <- calendar_path(country)
  if (!file.exists(p)) stop(sprintf("prepared calendar is missing: %s", p), call. = FALSE)
  d <- utils::read.csv(gzfile(p), stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = c(unit = "character", date = "character",
                                      E3 = "numeric", E7 = "numeric", E14 = "numeric", E21 = "numeric",
                                      Emax3 = "numeric", tmean = "numeric", tmax = "numeric",
                                      population = "numeric", y_K80 = "numeric", y_K81 = "numeric",
                                      y_urinary = "numeric", y_sum3 = "numeric"))
  required <- c("unit", "date", "E3", "population", "y_K80", "y_K81", "y_urinary", "y_sum3")
  missing <- setdiff(required, names(d))
  if (length(missing)) stop(sprintf("%s calendar lacks columns: %s", country, paste(missing, collapse = ", ")), call. = FALSE)
  d$date <- as.Date(d$date)
  if (anyNA(d$date) || anyNA(d$unit)) stop(sprintf("%s calendar has invalid unit/date values", country), call. = FALSE)
  if (anyDuplicated(paste(d$unit, d$date, sep = "\r"))) stop(sprintf("%s calendar has duplicate unit/date rows", country), call. = FALSE)
  d
}

safe_unit_id <- function(unit) paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(as.character(unit))))), collapse = "")

cache_path <- function(country, unit) {
  file.path(fit_dir, sprintf("primary_%s_unit_%s.rds", country, safe_unit_id(unit)))
}

.input_md5 <- list()
current_cache_manifest <- function() {
  path <- file.path(fit_dir, "primary_cache_manifest.csv")
  if (!file.exists(path)) stop("current successful primary-cache manifest is missing", call. = FALSE)
  m <- read.csv(path, stringsAsFactors = FALSE,
                colClasses = c(unit = "character", run_id = "character"))
  required <- c("run_id", "country", "unit", "scenario", "cache_path", "input_md5", "protocol_signature")
  if (!all(required %in% names(m)) || !nrow(m)) stop("current primary-cache manifest is empty or incomplete", call. = FALSE)
  if (length(unique(m$run_id)) != 1L || anyDuplicated(m[, c("country", "unit", "scenario")])) {
    stop("primary-cache manifest has mixed runs or duplicate units", call. = FALSE)
  }
  m
}

read_primary_cache <- function(path, country, unit = NULL) {
  if (!file.exists(path)) stop(sprintf("primary fit cache is missing: %s", path), call. = FALSE)
  z <- readRDS(path)
  if (!is.list(z) || is.null(z$engine_fit) || is.null(z$unit) || is.null(z$country)) {
    stop("cache does not have the documented primary-fit structure", call. = FALSE)
  }
  if (toupper(as.character(z$country)) != country) stop("cache country metadata does not match its filename", call. = FALSE)
  if (!is.null(unit) && as.character(z$unit) != as.character(unit)) stop("cache unit metadata does not match its filename", call. = FALSE)
  if (!identical(z$scenario, "primary_threeoutcome_E3_separate")) stop("cache scenario is not the primary separate-calendar E3 fit", call. = FALSE)
  m <- current_cache_manifest()
  approved <- m[m$country == country & m$unit == as.character(z$unit) & m$scenario == z$scenario, , drop = FALSE]
  if (nrow(approved) != 1L || is.null(z$run_id) || is.null(z$input_md5) || is.null(z$protocol_signature) ||
      !identical(as.character(z$run_id), approved$run_id) ||
      !identical(as.character(z$input_md5), approved$input_md5) ||
      !identical(as.character(z$protocol_signature), approved$protocol_signature) ||
      normalizePath(path, winslash = "/", mustWork = TRUE) != approved$cache_path) {
    stop("cache is absent from the current successful run or its provenance differs", call. = FALSE)
  }
  if (is.null(.input_md5[[country]])) .input_md5[[country]] <<- unname(tools::md5sum(calendar_path(country)))
  if (!identical(as.character(z$input_md5), .input_md5[[country]])) {
    stop("cache input checksum differs from the current canonical calendar", call. = FALSE)
  }
  ef <- z$engine_fit
  if (!identical(ef$covariance_scope, "full") || is.null(ef$covariance$hac14) ||
      is.null(ef$score$daily_scores) || is.null(ef$score$Hinv)) {
    stop("cache lacks full joint HAC14 covariance or daily score fields", call. = FALSE)
  }
  z
}

eligible_counts <- function(d, endpoint, offset) {
  keep <- is.finite(d[[endpoint]]) & is.finite(d$E3)
  if (!is.null(offset)) keep <- keep & is.finite(d[[offset]]) & d[[offset]] > 0
  c(days = sum(keep), events = sum(d[[endpoint]][keep]))
}

shared_reference <- function(d, country) {
  keep <- is.finite(d$E3) & is.finite(d$y_K80) & is.finite(d$y_K81) & is.finite(d$y_urinary)
  if (country == "BR") keep <- keep & is.finite(d$population) & d$population > 0
  if (!any(keep)) return(c(p95 = NA_real_, median = NA_real_))
  c(p95 = unname(stats::quantile(d$E3[keep], 0.95, type = 7)),
    median = unname(stats::median(d$E3[keep])),
    p10 = unname(stats::quantile(d$E3[keep], 0.10, type = 7)),
    p60 = unname(stats::quantile(d$E3[keep], 0.60, type = 7)))
}

check_cache_reference <- function(cache, d, country) {
  target <- shared_reference(d, country)
  if (any(!is.finite(target)) || !is.finite(cache$p95) || !is.finite(cache$median) ||
      abs(target[["p95"]] - cache$p95) > 1e-8 || abs(target[["median"]] - cache$median) > 1e-8) {
    stop("cached P95/shared median do not match the common complete paired-date E3 reference", call. = FALSE)
  }
  invisible(target)
}

symmetrize <- function(V) (V + t(V)) / 2

contrast_stats <- function(L, beta, V) {
  est <- as.numeric(L %*% beta)
  v <- as.numeric(t(L) %*% V %*% L)
  if (!is.finite(v) || v < -1e-8) stop("contrast covariance is negative or non-finite", call. = FALSE)
  se <- sqrt(max(0, v))
  c(estimate = est, variance = max(0, v), se = se, RR = exp(est),
    lo = exp(est - stats::qnorm(0.975) * se), hi = exp(est + stats::qnorm(0.975) * se))
}

stable_seed <- function(country, unit) {
  x <- utf8ToInt(paste(country, unit, sep = ":"))
  as.integer((sum((x %% 1009L) * seq_along(x)) %% 2147483646) + 1)
}

normal_draws <- function(mu, V, n = 1000L) {
  V <- symmetrize(V)
  ee <- eigen(V, symmetric = TRUE)
  scale <- max(1, max(abs(ee$values)))
  if (min(ee$values) < -1e-8 * scale) stop("HAC exposure-coefficient covariance is not positive semidefinite", call. = FALSE)
  root <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), nrow = length(mu))
  Z <- matrix(stats::rnorm(n * length(mu)), nrow = n, ncol = length(mu))
  sweep(Z %*% t(root), 2L, mu, `+`)
}

meta_row <- function(country, endpoint, contrast, method, yi, vi, region_days = NA_real_,
                     n_regions = NULL, support_regions = NA_integer_, support_days = NA_real_,
                     supporting_region_days = NA_real_) {
  pooled <- pool_meta(yi, vi, method = method)
  if (is.null(n_regions)) n_regions <- pooled$k
  data.frame(country = country, endpoint = endpoint, contrast = contrast, method = method,
             k = pooled$k, I2 = pooled$I2, concordance = pooled$concordance,
             log_estimate = pooled$estimate, log_se = pooled$se, RR = pooled$RR,
             lo = pooled$RR95CI[1L], hi = pooled$RR95CI[2L], tau2 = pooled$tau2,
             contributing_regions = n_regions, contributor_regiondays = region_days,
             supporting_regions = support_regions, support_days_within_1C = support_days,
             supporting_regiondays = supporting_region_days,
             failure = if (is.null(pooled$failure)) "" else pooled$failure,
             stringsAsFactors = FALSE)
}

run_main <- function(country, calendar) {
  unit_ids <- sort(unique(calendar$unit))
  result_path <- file.path(out_dir, "main_unit.csv")
  per_unit <- list()
  for (unit in unit_ids) {
    d <- calendar[calendar$unit == unit, , drop = FALSE]
    p <- cache_path(country, unit)
    cache <- tryCatch(read_primary_cache(p, country, unit), error = function(e) {
      audit_failure(country, unit, "main_reference_cache", reason = conditionMessage(e))
      NULL
    })
    if (is.null(cache)) next
    reference_ok <- tryCatch(check_cache_reference(cache, d, country), error = function(e) e)
    if (inherits(reference_ok, "error")) {
      audit_failure(country, unit, "main_reference", reason = conditionMessage(reference_ok))
      next
    }
    offset <- if (country == "BR") "population" else NULL
    endpoints <- c("y_K80", "y_K81", "y_urinary", "y_sum3")
    set.seed(stable_seed(country, unit))
    for (endpoint in endpoints) {
      counts <- eligible_counts(d, endpoint, offset)
      if (counts[["days"]] < 400 || counts[["events"]] < 100) {
        audit_failure(country, unit, "main_effect", endpoint,
                      sprintf("eligibility threshold not met (<400 eligible days or <100 events)"),
                      counts[["days"]], counts[["events"]])
        next
      }
      fit <- tryCatch(fit_main_effect(data = d, outcome = endpoint, exposure = "E3", date = "date",
                                      offset = offset, p95 = cache$p95, median = cache$median,
                                      basis = cache$engine_fit$basis,
                                      mrt_bounds = unname(reference_ok[c("p10", "p60")]),
                                      hac_lags = 14L), error = function(e) e)
      if (inherits(fit, "error")) {
        audit_failure(country, unit, "main_effect", endpoint, conditionMessage(fit), counts[["days"]], counts[["events"]])
        next
      }
      hac <- fit$hac_vcov$hac14
      if (is.null(hac)) {
        audit_failure(country, unit, "main_effect", endpoint, "engine returned no HAC14 covariance", counts[["days"]], counts[["events"]])
        next
      }
      beta <- stats::coef(fit$fit)
      names(beta) <- names(stats::coef(fit$fit))
      names_hi <- .design_at(fit$fit, fit$data, basis = fit$basis, temp_value = fit$p95)
      names_med <- .design_at(fit$fit, fit$data, basis = fit$basis, temp_value = fit$median)
      L_med <- as.numeric(names_hi - names_med); names(L_med) <- names(beta)
      median_stat <- tryCatch(contrast_stats(L_med, beta, hac), error = function(e) e)
      if (inherits(median_stat, "error")) {
        audit_failure(country, unit, "main_shared_median", endpoint, conditionMessage(median_stat), counts[["days"]], counts[["events"]])
        next
      }
      mrt_stat <- fit$contrast$by_covariance$hac14
      bcols <- c("B1", "B2", "B3")
      Vb <- hac[bcols, bcols, drop = FALSE]
      draw_mrt <- tryCatch({
        draws <- normal_draws(beta[bcols], Vb, n = 1000L)
        grid <- seq(fit$mrt_bounds[1L], fit$mrt_bounds[2L], length.out = 2001L)
        Bgrid <- .basis_at(fit$basis, grid)
        bp95 <- as.numeric(.basis_at(fit$basis, fit$p95))
        risk <- draws %*% t(Bgrid)
        chosen <- max.col(-risk, ties.method = "first")
        mrt_draw <- grid[chosen]
        rr_draw <- exp(rowSums(draws * rep(bp95, each = nrow(draws))) - risk[cbind(seq_len(nrow(risk)), chosen)])
        list(lo = unname(stats::quantile(rr_draw, 0.025, names = FALSE, type = 7)),
             hi = unname(stats::quantile(rr_draw, 0.975, names = FALSE, type = 7)),
             lower_frequency = mean(chosen == 1L), upper_frequency = mean(chosen == length(grid)),
             mrt_lo = unname(stats::quantile(mrt_draw, 0.025, names = FALSE, type = 7)),
             mrt_hi = unname(stats::quantile(mrt_draw, 0.975, names = FALSE, type = 7)))
      }, error = function(e) e)
      if (inherits(draw_mrt, "error")) {
        audit_failure(country, unit, "main_mrt_coefficient_simulation", endpoint,
                      paste("1000-draw asymptotic HAC coefficient simulation failed:", conditionMessage(draw_mrt)),
                      counts[["days"]], counts[["events"]])
        next
      }
      boundary <- if (abs(fit$mrt - fit$mrt_bounds[1L]) < 1e-5) "P10" else
        if (abs(fit$mrt - fit$mrt_bounds[2L]) < 1e-5) "P60" else "interior"
      row <- data.frame(country = country, unit = unit, endpoint = endpoint,
                        eligible_days = counts[["days"]], events = counts[["events"]],
                        p95 = fit$p95, shared_median = fit$median,
                        mrt = fit$mrt, mrt_p10 = fit$mrt_bounds[1L], mrt_p60 = fit$mrt_bounds[2L],
                        mrt_boundary = boundary,
                        mrt_boundary_frequency_simulated = draw_mrt$lower_frequency + draw_mrt$upper_frequency,
                        mrt_p10_boundary_sim_frequency = draw_mrt$lower_frequency,
                        mrt_p60_boundary_sim_frequency = draw_mrt$upper_frequency,
                        mrt_sim_p025 = draw_mrt$mrt_lo, mrt_sim_p975 = draw_mrt$mrt_hi,
                        p95_vs_mrt_logRR = unname(mrt_stat["estimate"]),
                        p95_vs_mrt_RR = unname(mrt_stat["RR"]),
                        p95_vs_mrt_lo = unname(mrt_stat["lo"]), p95_vs_mrt_hi = unname(mrt_stat["hi"]),
                        p95_vs_mrt_sim_lo = draw_mrt$lo, p95_vs_mrt_sim_hi = draw_mrt$hi,
                        p95_vs_shared_median_logRR = unname(median_stat["estimate"]),
                        p95_vs_shared_median_RR = unname(median_stat["RR"]),
                        p95_vs_shared_median_lo = unname(median_stat["lo"]),
                        p95_vs_shared_median_hi = unname(median_stat["hi"]),
                        hac14_variance_p95_vs_mrt = unname(mrt_stat["variance"]),
                        hac14_variance_p95_vs_median = unname(median_stat["variance"]),
                        simulation = "1000 asymptotic joint HAC14 draws of B1:B3; MRT reselected on fixed 2001-point P10-P60 grid; not a refit bootstrap",
                        stringsAsFactors = FALSE)
      append_csv(row, result_path)
      per_unit[[length(per_unit) + 1L]] <- row
    }
  }
  if (length(per_unit)) {
    u <- do.call(rbind, per_unit)
    country_path <- file.path(out_dir, "main_country.csv")
    for (endpoint in unique(u$endpoint)) for (contrast in c("p95_vs_mrt", "p95_vs_shared_median")) {
      est_col <- paste0(contrast, "_logRR")
      var_col <- if (contrast == "p95_vs_mrt") "hac14_variance_p95_vs_mrt" else "hac14_variance_p95_vs_median"
      for (method in c("DL", "REML_HK")) {
        row <- meta_row(country, endpoint, contrast, method, u[[est_col]][u$endpoint == endpoint],
                        u[[var_col]][u$endpoint == endpoint], n_regions = sum(u$endpoint == endpoint))
        append_csv(row, country_path)
      }
    }
  }
}

model_matrix_at <- function(fit, template, basis, temperature) {
  nd <- template
  b <- .basis_at(basis, temperature)
  nd$B1 <- b[, 1L]; nd$B2 <- b[, 2L]; nd$B3 <- b[, 3L]
  tt <- stats::delete.response(stats::terms(fit))
  X <- stats::model.matrix(tt, data = nd, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  X[, names(stats::coef(fit)), drop = FALSE]
}

curve_vector <- function(ef, endpoint, temperature) {
  Xx <- .design_at(ef$fit, ef$data, endpoint = endpoint, basis = ef$basis, temp_value = temperature)
  Xm <- .design_at(ef$fit, ef$data, endpoint = endpoint, basis = ef$basis, temp_value = ef$median)
  L <- as.numeric(Xx - Xm); names(L) <- names(stats::coef(ef$fit))
  L
}

curve_effect <- function(ef, endpoint, temperature) {
  contrast_stats(curve_vector(ef, endpoint, temperature), stats::coef(ef$fit), ef$covariance$hac14)
}

burden_components <- function(ef) {
  fit <- ef$fit; beta <- stats::coef(fit); p <- length(beta)
  endpoints <- c("y_K80", "y_K81", "y_urinary")
  by_endpoint <- list()
  for (endpoint in endpoints) {
    z <- ef$data[ef$data$disease == endpoint, , drop = FALSE]
    if (!nrow(z)) stop(sprintf("no model rows for burden endpoint %s", endpoint), call. = FALSE)
    Xhi <- model_matrix_at(fit, z, ef$basis, ef$p95)
    Xmed <- model_matrix_at(fit, z, ef$basis, ef$median)
    off <- if ("off" %in% names(z)) z$off else rep(0, nrow(z))
    muhi <- exp(as.numeric(Xhi %*% beta) + off)
    mume <- exp(as.numeric(Xmed %*% beta) + off)
    by_endpoint[[endpoint]] <- list(date = as.Date(z$date_actual),
                                    delta = muhi - mume,
                                    derivative = muhi * Xhi - mume * Xmed)
  }
  ref_dates <- by_endpoint[[1L]]$date
  for (endpoint in endpoints) if (!identical(by_endpoint[[endpoint]]$date, ref_dates)) {
    stop("paired endpoint prediction dates are not aligned", call. = FALSE)
  }
  daily_delta <- list(y_K80 = by_endpoint$y_K80$delta,
                      y_K81 = by_endpoint$y_K81$delta,
                      y_urinary = by_endpoint$y_urinary$delta)
  daily_grad <- list(y_K80 = by_endpoint$y_K80$derivative,
                     y_K81 = by_endpoint$y_K81$derivative,
                     y_urinary = by_endpoint$y_urinary$derivative)
  daily_delta$y_sum3 <- daily_delta$y_K80 + daily_delta$y_K81 + daily_delta$y_urinary
  daily_grad$y_sum3 <- daily_grad$y_K80 + daily_grad$y_K81 + daily_grad$y_urinary
  daily_delta$urinary_minus_K80 <- daily_delta$y_urinary - daily_delta$y_K80
  daily_grad$urinary_minus_K80 <- daily_grad$y_urinary - daily_grad$y_K80
  out <- list(date = ref_dates, endpoints = list())
  V <- ef$covariance$hac14
  for (endpoint in names(daily_delta)) {
    g <- colMeans(daily_grad[[endpoint]])
    estimate <- mean(daily_delta[[endpoint]])
    variance <- as.numeric(t(g) %*% V %*% g)
    if (!is.finite(variance) || variance < -1e-8) stop(sprintf("invalid burden gradient variance for %s", endpoint), call. = FALSE)
    se <- sqrt(max(0, variance))
    out$endpoints[[endpoint]] <- list(estimate = estimate, gradient = g,
                                      variance = max(0, variance), se = se,
                                      lo = estimate - stats::qnorm(0.975) * se,
                                      hi = estimate + stats::qnorm(0.975) * se,
                                      daily_delta = daily_delta[[endpoint]])
  }
  U <- ef$score$daily_scores
  score_dates <- as.Date(ef$score$daily_dates)
  Hinv <- ef$score$Hinv
  influence <- list()
  for (endpoint in names(out$endpoints)) {
    g <- out$endpoints[[endpoint]]$gradient
    influence[[endpoint]] <- as.numeric(U %*% as.numeric(crossprod(Hinv, g)))
  }
  out$influence <- influence
  out$influence_dates <- score_dates
  out
}

national_hac14 <- function(psi) {
  n <- length(psi)
  v <- sum(psi^2)
  if (n > 1L) for (lag in seq_len(min(14L, n - 1L))) {
    weight <- 1 - lag / 15
    v <- v + 2 * weight * sum(psi[(lag + 1L):n] * psi[seq_len(n - lag)])
  }
  if (!is.finite(v) || v < -1e-8) stop("synchronized national Bartlett HAC14 variance is invalid", call. = FALSE)
  max(0, v)
}

run_curves_burden <- function(country, calendar, cache_files) {
  curve_path <- file.path(out_dir, "curves_unit.csv")
  audit32_path <- file.path(out_dir, "curves_32C_support_audit.csv")
  burden_unit_path <- file.path(out_dir, "burden_unit.csv")
  daily_path <- file.path(out_dir, "burden_daily_influence_country.csv")
  burden_country_path <- file.path(out_dir, "burden_country.csv")
  curve_unit_rows <- list()
  burden_daily_total <- list()
  included_units <- setNames(integer(5L), c("y_K80", "y_K81", "y_urinary", "y_sum3", "urinary_minus_K80"))
  burden_point_total <- setNames(numeric(5L), names(included_units))
  all_dates <- seq(min(calendar$date), max(calendar$date), by = "day")
  for (endpoint in names(included_units)) burden_daily_total[[endpoint]] <- setNames(rep(0, length(all_dates)), as.character(all_dates))
  curve_temps <- seq(14, 30, by = 0.5)
  cached_units <- sub("^primary_[A-Z]+_unit_", "", basename(cache_files))
  cached_units <- sub("\\.rds$", "", cached_units)
  panel_safe_units <- vapply(unique(calendar$unit), safe_unit_id, character(1))
  for (i in which(!panel_safe_units %in% cached_units)) {
    unit <- unique(calendar$unit)[i]
    audit_failure(country, unit, "curves_burden_cache", reason = "primary full-covariance fit cache is missing")
  }
  for (path in cache_files) {
    cache <- tryCatch(read_primary_cache(path, country), error = function(e) {
      audit_failure(country, NA_character_, "curves_burden_cache", reason = conditionMessage(e))
      NULL
    })
    if (is.null(cache)) next
    unit <- as.character(cache$unit)
    ef <- cache$engine_fit
    panel_unit <- calendar[calendar$unit == unit, , drop = FALSE]
    if (!nrow(panel_unit)) {
      audit_failure(country, unit, "curves_burden", reason = "unit is absent from canonical prepared calendar")
      next
    }
    reference_ok <- tryCatch(check_cache_reference(cache, panel_unit, country), error = function(e) e)
    if (inherits(reference_ok, "error")) {
      audit_failure(country, unit, "curves_burden_reference", reason = conditionMessage(reference_ok))
      next
    }
    modeled <- ef$data[ef$data$disease == ef$outcomes[1L], , drop = FALSE]
    # Recover actual modeled E3 by an exact date join to the prepared paired calendar.
    support_panel <- panel_unit[panel_unit$date %in% as.Date(modeled$date_actual), , drop = FALSE]
    support_panel <- support_panel[match(as.Date(modeled$date_actual), support_panel$date), , drop = FALSE]
    if (nrow(support_panel) != nrow(modeled) || anyNA(support_panel$E3)) {
      audit_failure(country, unit, "curve_support", reason = "modeled paired dates could not be aligned to actual E3 values")
      next
    }
    e <- support_panel$E3
    p2 <- unname(stats::quantile(e, 0.02, type = 7)); p98 <- unname(stats::quantile(e, 0.98, type = 7))
    day_count32 <- sum(abs(e - 32) <= 1)
    append_csv(data.frame(country = country, unit = unit, modeled_E3_days = length(e),
                          local_p2 = p2, local_p98 = p98, days_within_1C_of_32 = day_count32,
                          supported_at_32 = 32 >= p2 && 32 <= p98 && day_count32 >= 30,
                          plotted = FALSE, stringsAsFactors = FALSE), audit32_path)
    unit_curve_rows <- list()
    curve_error <- NULL
    for (temp in curve_temps) {
      near_days <- sum(abs(e - temp) <= 1)
      locally_supported <- temp >= p2 && temp <= p98 && near_days >= 30
      for (endpoint in c("y_K80", "y_urinary")) {
        stat <- tryCatch(curve_effect(ef, endpoint, temp), error = function(err) err)
        if (inherits(stat, "error")) {
          curve_error <- conditionMessage(stat)
          audit_failure(country, unit, "curve_prediction", paste(endpoint, temp, sep = "@"), curve_error,
                        length(e), NA_real_)
          next
        }
        row <- data.frame(country = country, unit = unit, temperature = temp,
                          endpoint = endpoint, contrast = "RR_vs_shared_median",
                          log_estimate = unname(stat["estimate"]), variance = unname(stat["variance"]),
                          RR = unname(stat["RR"]), lo = unname(stat["lo"]), hi = unname(stat["hi"]),
                          locally_supported = locally_supported, local_p2 = p2, local_p98 = p98,
                          eligible_regiondays = length(e), days_within_1C = near_days,
                          stringsAsFactors = FALSE)
        unit_curve_rows[[length(unit_curve_rows) + 1L]] <- row
      }
      ku <- tryCatch(curve_effect(ef, "y_K80", temp), error = function(err) err)
      ur <- tryCatch(curve_effect(ef, "y_urinary", temp), error = function(err) err)
      if (!inherits(ku, "error") && !inherits(ur, "error")) {
        Lku <- curve_vector(ef, "y_K80", temp)
        Lur <- curve_vector(ef, "y_urinary", temp)
        Ldirect <- as.numeric(Lur - Lku); names(Ldirect) <- names(stats::coef(ef$fit))
        direct <- tryCatch(contrast_stats(Ldirect, stats::coef(ef$fit), ef$covariance$hac14), error = function(err) err)
        if (!inherits(direct, "error")) {
          row <- data.frame(country = country, unit = unit, temperature = temp,
                            endpoint = "urinary_vs_K80", contrast = "direct_RRR",
                            log_estimate = unname(direct["estimate"]), variance = unname(direct["variance"]),
                            RR = unname(direct["RR"]), lo = unname(direct["lo"]), hi = unname(direct["hi"]),
                            locally_supported = locally_supported, local_p2 = p2, local_p98 = p98,
                            eligible_regiondays = length(e), days_within_1C = near_days,
                            stringsAsFactors = FALSE)
          unit_curve_rows[[length(unit_curve_rows) + 1L]] <- row
        } else audit_failure(country, unit, "curve_direct_RRR", as.character(temp), conditionMessage(direct), length(e), NA_real_)
      } else {
        audit_failure(country, unit, "curve_direct_RRR", as.character(temp), "K80 or urinary prediction failed", length(e), NA_real_)
      }
    }
    if (length(unit_curve_rows)) {
      ur <- do.call(rbind, unit_curve_rows)
      for (i in seq_len(nrow(ur))) append_csv(ur[i, , drop = FALSE], curve_path)
      curve_unit_rows[[length(curve_unit_rows) + 1L]] <- ur
    }
    # A failed isolated prediction is audited, but successful temperature/endpoint rows remain usable.
    if (!length(unit_curve_rows) || is.null(ef$covariance$hac14)) {
      if (!is.null(curve_error)) audit_failure(country, unit, "curve_unit", reason = curve_error)
    }
    burden <- tryCatch(burden_components(ef), error = function(e) e)
    if (inherits(burden, "error")) {
      audit_failure(country, unit, "standardized_burden", reason = conditionMessage(burden))
      next
    }
    burden_names <- names(burden$endpoints)
    for (endpoint in burden_names) {
      b <- burden$endpoints[[endpoint]]
      row <- data.frame(country = country, unit = unit, endpoint = endpoint,
                        eligible_paired_days = length(burden$date), p95 = cache$p95,
                        shared_median = cache$median, mean_daily_admission_difference = b$estimate,
                        se = b$se, lo = b$lo, hi = b$hi, delta_method_variance = b$variance,
                        covariance = "full joint HAC14; gradient of standardized mean prediction",
                        stringsAsFactors = FALSE)
      append_csv(row, burden_unit_path)
      included_units[endpoint] <- included_units[endpoint] + 1L
      burden_point_total[endpoint] <- burden_point_total[endpoint] + b$estimate
      daily_idx <- match(as.character(burden$influence_dates), names(burden_daily_total[[endpoint]]))
      valid_idx <- !is.na(daily_idx)
      burden_daily_total[[endpoint]][daily_idx[valid_idx]] <-
        burden_daily_total[[endpoint]][daily_idx[valid_idx]] + burden$influence[[endpoint]][valid_idx]
    }
  }
  burden_ok <- sum(included_units) > 0L
  if (!burden_ok) audit_failure(country, analysis = "standardized_burden", reason = "no primary fit cache produced a valid burden estimate")
  if (!length(curve_unit_rows)) {
    audit_failure(country, analysis = "curve_pooling", reason = "no regional curve predictions succeeded")
  }
  if (burden_ok) {
    daily_df <- data.frame(country = country, date = all_dates, stringsAsFactors = FALSE)
    for (endpoint in names(burden_daily_total)) daily_df[[paste0("psi_", endpoint)]] <- unname(burden_daily_total[[endpoint]])
    append_csv(daily_df, daily_path)
    for (endpoint in names(burden_point_total)) {
      psi <- burden_daily_total[[endpoint]]
      var <- tryCatch(national_hac14(psi), error = function(e) e)
      if (inherits(var, "error")) {
        audit_failure(country, NA_character_, "national_burden_HAC14", endpoint, conditionMessage(var))
        next
      }
      se <- sqrt(var)
      row <- data.frame(country = country, endpoint = endpoint,
                        included_region_count = included_units[endpoint],
                        national_point_sum_of_unit_means = burden_point_total[endpoint],
                        synchronized_daily_HAC14_variance = var, se = se,
                        lo = burden_point_total[endpoint] - stats::qnorm(0.975) * se,
                        hi = burden_point_total[endpoint] + stats::qnorm(0.975) * se,
                        uncertainty = "synchronized calendar-day influence sum; Bartlett HAC14 includes cross-region and cross-day covariance",
                        stringsAsFactors = FALSE)
      append_csv(row, burden_country_path)
    }
  }
  if (length(curve_unit_rows)) {
    u <- do.call(rbind, curve_unit_rows)
    pool_path <- file.path(out_dir, "curves_country.csv")
    groups <- unique(u[c("temperature", "endpoint", "contrast")])
    for (j in seq_len(nrow(groups))) {
      gg <- groups[j, , drop = FALSE]
      z_all <- u[u$temperature == gg$temperature & u$endpoint == gg$endpoint & u$contrast == gg$contrast, , drop = FALSE]
      z_supported <- z_all[z_all$locally_supported, , drop = FALSE]
      for (scope in c("all_successful", "support_restricted")) {
        z <- if (scope == "all_successful") z_all else z_supported
        for (method in c("DL", "REML_HK")) {
          row <- meta_row(country, gg$endpoint, paste(gg$contrast, scope, sep = ":"), method,
                          z$log_estimate, z$variance,
                          region_days = sum(z$eligible_regiondays), n_regions = nrow(z),
                          support_regions = nrow(z_supported),
                          support_days = sum(z_supported$days_within_1C),
                          supporting_region_days = sum(z_supported$eligible_regiondays))
          row$temperature <- gg$temperature
          row$curve_scope <- scope
          append_csv(row, pool_path)
        }
      }
    }
  }
}

for (country in countries) {
  calendar <- tryCatch(read_calendar(country), error = function(e) {
    audit_failure(country, analysis = "input_calendar", reason = conditionMessage(e))
    NULL
  })
  if (is.null(calendar)) next
  if (!dir.exists(fit_dir)) {
    audit_failure(country, analysis = "fit_cache", reason = sprintf("fit cache directory is missing: %s", fit_dir))
    next
  }
  if (stage == "main") {
    run_main(country, calendar)
  } else {
    manifest <- current_cache_manifest()
    files <- manifest$cache_path[manifest$country == country & manifest$scenario == "primary_threeoutcome_E3_separate"]
    if (!length(files)) {
      audit_failure(country, analysis = "fit_cache", reason = "no primary full-covariance fit caches found")
      next
    }
    run_curves_burden(country, calendar, files)
  }
}

cat(sprintf("Stage %s finished; outputs: %s\n", stage, normalizePath(out_dir, winslash = "/", mustWork = TRUE)))
