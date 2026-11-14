#!/usr/bin/env Rscript
# Prespecified secondary DLNM timing sensitivity. Weather histories are
# constructed from exact calendar dates before any outcome filtering.

set_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", BLAS_NUM_THREADS = "1")
}
set_threads()

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_path <- if (length(.script_arg)) {
  normalizePath(sub("^--file=", "", .script_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  .source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(.source_file)) stop("cannot resolve run_dlnm.R location", call. = FALSE)
  normalizePath(.source_file, winslash = "/", mustWork = TRUE)
}
.script_dir <- dirname(.script_path)
.repair_dir <- normalizePath(file.path(.script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "repair_engine.R"), local = FALSE)

.dlnm_constants <- list(
  outcomes = c("y_K80", "y_urinary"),
  primary_outcomes = c("y_K80", "y_K81", "y_urinary"),
  temp_df = 3L, lag_df = 3L, max_lag = 21L, hac_lag = 28L,
  min_days = 400L, min_cases = 40L, min_workers = 8L,
  lag_knot = sqrt(21)
)
.cb_names <- as.vector(outer(1:3, 1:3, function(j, k) paste0("CB", j, "_", k)))

parse_args <- function(args) {
  opts <- list(country = "all", workers = 8L,
               brazil = file.path(.repair_dir, "data", "br_calendar.csv.gz"),
               mexico = file.path(.repair_dir, "data", "mx_calendar.csv.gz"),
               outdir = file.path(.repair_dir, "results", "dlnm"),
               manifest = file.path(.repair_dir, "results", "fits", "primary_cache_manifest.csv"))
  i <- 1L
  while (i <= length(args)) {
    token <- args[i]
    if (!startsWith(token, "--")) stop(sprintf("unexpected argument '%s'", token), call. = FALSE)
    keyval <- sub("^--", "", token)
    if (grepl("=", keyval, fixed = TRUE)) {
      bits <- strsplit(keyval, "=", fixed = TRUE)[[1]]
      key <- bits[1]; value <- paste(bits[-1], collapse = "="); i <- i + 1L
    } else {
      key <- keyval
      if (i == length(args)) stop(sprintf("missing value for --%s", key), call. = FALSE)
      value <- args[i + 1L]; i <- i + 2L
    }
    if (!key %in% c("country", "workers", "brazil", "mexico", "outdir", "manifest")) {
      stop(sprintf("unknown option --%s", key), call. = FALSE)
    }
    opts[[key]] <- value
  }
  opts$country <- match.arg(opts$country, c("all", "BR", "MX"))
  opts$workers <- suppressWarnings(as.integer(opts$workers))
  if (length(opts$workers) != 1L || is.na(opts$workers) || opts$workers < 1L ||
      opts$workers > .dlnm_constants$min_workers) {
    stop("--workers must be an integer from 1 through 8", call. = FALSE)
  }
  opts$brazil <- normalizePath(opts$brazil, winslash = "/", mustWork = TRUE)
  opts$mexico <- normalizePath(opts$mexico, winslash = "/", mustWork = TRUE)
  opts$manifest <- path.expand(opts$manifest)
  opts$outdir <- path.expand(opts$outdir)
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
  opts$outdir <- normalizePath(opts$outdir, winslash = "/", mustWork = TRUE)
  opts$manifest <- normalizePath(opts$manifest, winslash = "/", mustWork = FALSE)
  opts
}

read_calendar <- function(path, country) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)
  required <- c("unit", "date", "tmean", "E3", .dlnm_constants$primary_outcomes)
  if (country == "BR") required <- c(required, "population")
  col_classes <- c(unit = "character", date = "character")
  if (grepl("\\.gz$", path, ignore.case = TRUE) && .Platform$OS.type == "windows") {
    # data.table's compressed-file reader requires optional R.utils on Windows;
    # use base gzfile so the repair runner remains dependency-free there.
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    numeric_columns <- setdiff(required, c("unit", "date"))
    classes <- c(unit = "character", date = "character",
      stats::setNames(rep("numeric", length(numeric_columns)), numeric_columns))
    d <- data.table::as.data.table(utils::read.csv(con, stringsAsFactors = FALSE,
      check.names = FALSE, na.strings = c("", "NA"), colClasses = classes))
  } else {
    input <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      paste("gzip -dc --", shQuote(path))
    } else path
    d <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      data.table::fread(cmd = input, select = required, colClasses = col_classes, showProgress = TRUE)
    } else {
      data.table::fread(input, select = required, colClasses = col_classes, showProgress = TRUE)
    }
  }
  missing <- setdiff(required, names(d))
  if (length(missing)) stop(sprintf("%s calendar missing columns: %s", country,
                                     paste(missing, collapse = ", ")), call. = FALSE)
  d[, date := as.Date(date)]
  if (anyNA(d$date) || anyNA(d$unit) || any(!nzchar(d$unit))) {
    stop(sprintf("%s calendar has invalid dates or units", country), call. = FALSE)
  }
  if (anyDuplicated(d, by = c("unit", "date"))) {
    stop(sprintf("%s calendar has duplicate unit/date rows", country), call. = FALSE)
  }
  for (nm in setdiff(required, c("unit", "date"))) {
    if (!is.numeric(d[[nm]])) stop(sprintf("%s column is not numeric: %s", country, nm), call. = FALSE)
  }
  d
}

# Match primary E3 reference eligibility from the core runner. The primary
# P95/median are preserved, while the DLNM spline is trained on raw daily T on
# those same paired analysis dates.
primary_reference <- function(unit_data, country) {
  keep <- is.finite(unit_data$E3)
  for (nm in .dlnm_constants$primary_outcomes) keep <- keep & is.finite(unit_data[[nm]]) & unit_data[[nm]] >= 0
  if (country == "BR") keep <- keep & is.finite(unit_data$population) & unit_data$population > 0
  e <- unit_data$E3[keep]
  if (length(e) < 4L || length(unique(e)) < 2L) stop("fewer than four common primary E3 dates or fewer than two E3 values", call. = FALSE)
  keep_temp <- keep & is.finite(unit_data$tmean)
  tv <- unit_data$tmean[keep_temp]
  if (length(tv) < 4L || length(unique(tv)) < 2L) stop("fewer than four finite paired raw daily temperatures or fewer than two temperatures", call. = FALSE)
  list(p95 = unname(stats::quantile(e, 0.95, type = 7)),
       median = unname(stats::median(e)),
       basis = new_temperature_basis(tv, df = .dlnm_constants$temp_df),
       n_reference = length(e), n_basis = length(tv),
       ref_dates = unit_data$date[keep], basis_dates = unit_data$date[keep_temp])
}

lag_basis_matrix <- function() {
  b <- splines::ns(0:.dlnm_constants$max_lag,
    knots = .dlnm_constants$lag_knot, intercept = TRUE,
    Boundary.knots = c(0, .dlnm_constants$max_lag))
  if (ncol(b) != .dlnm_constants$lag_df || any(!is.finite(b))) {
    stop("lag spline does not have the prespecified three finite columns", call. = FALSE)
  }
  colnames(b) <- paste0("lag", seq_len(ncol(b)))
  b
}

# Build Btemp(T[t-l]) x Blag(l) for exact t-l calendar dates. Calendar gaps,
# non-finite T, and the first 21 days have incomplete histories and stay NA.
build_crossbasis <- function(unit_data, temp_basis, lag_basis = lag_basis_matrix()) {
  d <- as.Date(unit_data$date)
  if (anyNA(d) || anyDuplicated(d)) stop("crossbasis dates must be unique valid dates within unit", call. = FALSE)
  o <- order(d)
  d <- d[o]; temp <- unit_data$tmean[o]
  calendar <- seq(min(d), max(d), by = "day")
  cal_t <- rep(NA_real_, length(calendar))
  pos <- as.integer(d - min(calendar)) + 1L
  cal_t[pos] <- temp
  ncal <- length(calendar)
  cal_cb <- matrix(0, nrow = ncal, ncol = 9L,
                   dimnames = list(NULL, .cb_names))
  complete <- rep(TRUE, ncal)
  if (ncal <= .dlnm_constants$max_lag) complete[] <- FALSE
  if (ncal > .dlnm_constants$max_lag) for (lag in 0:.dlnm_constants$max_lag) {
    src <- seq_len(ncal) - lag
    in_range <- src >= 1L
    hist_t <- rep(NA_real_, ncal)
    hist_t[in_range] <- cal_t[src[in_range]]
    finite <- in_range & is.finite(hist_t)
    complete <- complete & finite
    if (any(finite)) {
      bt <- .basis_at(temp_basis, hist_t[finite])
      add <- do.call(cbind, lapply(seq_len(ncol(lag_basis)), function(k) {
        bt * lag_basis[lag + 1L, k]
      }))
      cal_cb[finite, ] <- cal_cb[finite, , drop = FALSE] + add
    }
  }
  cal_cb[!complete, ] <- NA_real_
  result <- matrix(NA_real_, nrow = length(o), ncol = ncol(cal_cb), dimnames = list(NULL, .cb_names))
  result[o, ] <- cal_cb[pos, , drop = FALSE]
  complete_result <- rep(FALSE, length(o))
  complete_result[o] <- complete[pos]
  list(features = result, complete_history = complete_result, dates = as.Date(unit_data$date),
       calendar_first = min(calendar), calendar_last = max(calendar),
       rows_missing_history = sum(!complete[pos]))
}

make_stack <- function(unit_data, country, features, min_days = .dlnm_constants$min_days,
                       min_cases = .dlnm_constants$min_cases) {
  d <- unit_data
  d$date <- as.Date(d$date)
  if (nrow(features$features) != nrow(d)) stop("crossbasis rows do not match calendar rows", call. = FALSE)
  eligible <- features$complete_history & is.finite(d$y_K80) & d$y_K80 >= 0 &
    is.finite(d$y_urinary) & d$y_urinary >= 0 & is.finite(d$E3)
  if (country == "BR") eligible <- eligible & is.finite(d$population) & d$population > 0
  n <- sum(eligible)
  events <- c(K80 = sum(d$y_K80[eligible]), urinary = sum(d$y_urinary[eligible]))
  if (n < min_days || any(events < min_cases)) {
    return(list(eligible = eligible, n = n, events = events, stack = NULL,
      reason = sprintf("complete lag-history paired days=%d (required >=%d); endpoint events K80=%.0f urinary=%.0f (required >=%d each)",
                       n, min_days, events[1], events[2], min_cases)))
  }
  dates <- d$date[eligible]
  yy <- c(d$y_K80[eligible], d$y_urinary[eligible])
  st <- data.frame(y = yy,
    disease = factor(rep(c("K80", "urinary"), each = n), levels = c("K80", "urinary")),
    date_actual = rep(dates, 2L),
    dow = factor(rep(as.integer(format(dates, "%u")), 2L), levels = 1:7),
    month = factor(rep(as.integer(format(dates, "%m")), 2L), levels = 1:12),
    time = rep(as.numeric(dates - min(d$date)) / 365.25, 2L),
    off = if (country == "BR") rep(log(d$population[eligible]), 2L) else rep(0, 2L),
    features$features[eligible, , drop = FALSE],
    check.names = FALSE)
  list(eligible = eligible, n = n, events = events, stack = st, reason = NULL)
}

stack_formula <- function(country) {
  rhs <- paste(c(.cb_names, "dow", "month", "splines::ns(time, 12)"), collapse = " + ")
  as.formula(paste0("y ~ disease * (", rhs, ")",
                    if (country == "BR") " + offset(off)" else ""))
}

separate_formula <- function(country, time_term = "splines::ns(time, 12)") {
  rhs <- paste(c(.cb_names, "dow", "month", time_term), collapse = " + ")
  as.formula(paste0("y ~ ", rhs, if (country == "BR") " + offset(off)" else ""))
}

.design_row <- function(fit, template, disease, cb_values) {
  nd <- template[1L, , drop = FALSE]
  nd$disease <- factor(disease, levels = levels(template$disease))
  for (j in seq_along(.cb_names)) nd[[.cb_names[j]]] <- cb_values[j]
  tt <- stats::delete.response(stats::terms(fit))
  x <- stats::model.matrix(tt, data = nd, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  x[, names(stats::coef(fit)), drop = FALSE]
}

interaction_basis_contrasts <- function(fit, stack) {
  z <- rep(0, length(.cb_names)); names(z) <- .cb_names
  L <- matrix(NA_real_, nrow = length(.cb_names), ncol = length(stats::coef(fit)),
              dimnames = list(.cb_names, names(stats::coef(fit))))
  for (j in seq_along(.cb_names)) {
    one <- z; one[j] <- 1
    u1 <- .design_row(fit, stack, "urinary", one)
    u0 <- .design_row(fit, stack, "urinary", z)
    k1 <- .design_row(fit, stack, "K80", one)
    k0 <- .design_row(fit, stack, "K80", z)
    L[j, ] <- as.numeric(u1 - u0 - k1 + k0)
  }
  L
}

contrast_weights <- function(temp_delta, lag_weights) {
  w <- numeric(length(.cb_names))
  for (k in seq_len(3L)) for (j in seq_len(3L)) w[(k - 1L) * 3L + j] <- temp_delta[j] * lag_weights[k]
  w
}

summary_contrast <- function(L, fit, V, label) {
  beta <- stats::coef(fit)
  estimate <- as.numeric(L %*% beta)
  variance <- as.numeric(t(L) %*% V %*% L)
  if (!is.finite(variance) || variance <= 0) stop(sprintf("%s has non-positive or non-finite HAC28 variance", label), call. = FALSE)
  se <- sqrt(variance)
  c(estimate = estimate, variance = variance, se = se, RR = exp(estimate),
    lo = exp(estimate - stats::qnorm(0.975) * se),
    hi = exp(estimate + stats::qnorm(0.975) * se))
}

validate_fit_equivalence <- function(fit, stack, country, tol = 1e-6) {
  predvars <- attr(stats::terms(fit), "predvars")
  is_time_ns <- vapply(predvars, function(x) is.call(x) &&
    grepl("ns\\(time,", paste(deparse(x), collapse = "")), logical(1))
  if (sum(is_time_ns) != 1L) stop("cannot recover the stacked fit's ns(time, 12) definition for the separate-GLM check", call. = FALSE)
  time_term <- paste(deparse(predvars[[which(is_time_ns)]], width.cutoff = 500L), collapse = "")
  fsep <- separate_formula(country, time_term = time_term)
  n <- nrow(stack) / 2L
  idx <- list(K80 = seq_len(n), urinary = n + seq_len(n))
  max_diff <- 0
  for (disease in names(idx)) {
    dd <- stack[idx[[disease]], , drop = FALSE]
    dd$y <- stack$y[idx[[disease]]]
    separate <- .fit_glm_audited(fsep, dd)$fit
    predicted <- stats::fitted(separate)
    actual <- fit$fitted.values[idx[[disease]]]
    diff <- max(abs(predicted - actual) / pmax(1, abs(predicted)))
    max_diff <- max(max_diff, diff)
  }
  if (!is.finite(max_diff) || max_diff > tol) {
    stop(sprintf("stacked model differs from equivalent separate endpoint GLMs (max relative fitted-mean difference %.3g; tolerance %.3g)",
                 max_diff, tol), call. = FALSE)
  }
  max_diff
}

fit_one_unit <- function(unit_data, country, unit, manifest_ref = NULL) {
  ref <- primary_reference(unit_data, country)
  features <- build_crossbasis(unit_data, ref$basis)
  prepared <- make_stack(unit_data, country, features)
  if (is.null(prepared$stack)) {
    return(list(failure = data.frame(country = country, unit = unit, phase = "eligibility",
      reason = prepared$reason, n_complete = prepared$n, events_K80 = prepared$events[1],
      events_urinary = prepared$events[2], stringsAsFactors = FALSE)))
  }
  formula <- stack_formula(country)
  audited <- .fit_glm_audited(formula, prepared$stack)
  fit <- audited$fit
  V <- .score_covariances(fit, audited$X, y = fit$y, mu = fit$fitted.values,
    dates = prepared$stack$date_actual, hac_lags = .dlnm_constants$hac_lag,
    full_covariance = TRUE)
  hac <- V$hac[[paste0("hac", .dlnm_constants$hac_lag)]]
  Lbasis <- interaction_basis_contrasts(fit, prepared$stack)
  dB <- as.numeric(.basis_at(ref$basis, ref$p95) - .basis_at(ref$basis, ref$median))
  cum_w <- contrast_weights(dB, colSums(lag_basis_matrix()))
  Lcum <- as.numeric(cum_w %*% Lbasis)
  cumulative <- summary_contrast(Lcum, fit, hac, "cumulative sustained P95 versus median")
  lag_b <- lag_basis_matrix()
  lag_rows <- lapply(0:.dlnm_constants$max_lag, function(lag) {
    w <- contrast_weights(dB, lag_b[lag + 1L, ])
    Llag <- as.numeric(w %*% Lbasis)
    q <- summary_contrast(Llag, fit, hac, sprintf("lag %d", lag))
    data.frame(country = country, unit = unit, lag = lag,
      estimate = unname(q["estimate"]), variance = unname(q["variance"]), se = unname(q["se"]),
      RR = unname(q["RR"]), lo = unname(q["lo"]), hi = unname(q["hi"]),
      p95_E3 = ref$p95, median_E3 = ref$median,
      stringsAsFactors = FALSE)
  })
  # Verify the reported cumulative contrast variance before accepting a unit fit.
  manual_var <- as.numeric(t(Lcum) %*% hac %*% Lcum)
  if (!isTRUE(all.equal(unname(cumulative["variance"]), manual_var, tolerance = 1e-10))) {
    stop("cumulative contrast variance failed the independent matrix identity check", call. = FALSE)
  }
  joint_failure <- ""
  joint_wald <- joint_p <- NA_real_
  joint_cov <- Lbasis %*% hac %*% t(Lbasis)
  joint_est <- as.numeric(Lbasis %*% stats::coef(fit))
  joint_rank <- qr(joint_cov)$rank
  if (joint_rank == nrow(Lbasis) && all(is.finite(joint_cov))) {
    joint_wald <- tryCatch(as.numeric(t(joint_est) %*% solve(joint_cov, joint_est)),
                           error = function(e) NA_real_)
    if (is.finite(joint_wald) && joint_wald >= 0) joint_p <- stats::pchisq(joint_wald, df = nrow(Lbasis), lower.tail = FALSE)
    else joint_failure <- "joint 9-df interaction Wald covariance inversion or statistic failed"
  } else joint_failure <- sprintf("joint 9-df interaction covariance rank %d/%d", joint_rank, nrow(Lbasis))
  eq_diff <- validate_fit_equivalence(fit, prepared$stack, country)
  ref_match <- NA
  if (!is.null(manifest_ref)) {
    ref_match <- abs(ref$p95 - manifest_ref$p95) <= 1e-10 && abs(ref$median - manifest_ref$median) <= 1e-10
    if (!ref_match) stop("canonical primary reference does not match the primary cache manifest", call. = FALSE)
  }
  audit <- data.frame(country = country, unit = unit, converged = audited$audit$converged,
    boundary = audited$audit$boundary, rank = audited$audit$rank, columns = audited$audit$columns,
    n_complete = prepared$n, events_K80 = prepared$events[1], events_urinary = prepared$events[2],
    n_reference_E3 = ref$n_reference, n_basis_rawT = ref$n_basis,
    n_missing_lag_history = features$rows_missing_history,
    history_first_date = as.character(features$calendar_first),
    history_last_date = as.character(features$calendar_last),
    p95_E3 = ref$p95, median_E3 = ref$median,
    manifest_reference_match = ref_match,
    separate_glm_max_relative_fitted_mean_difference = eq_diff,
    joint_interaction_df = nrow(Lbasis), joint_interaction_rank = joint_rank,
    joint_interaction_wald = joint_wald, joint_interaction_p = joint_p,
    joint_interaction_failure = joint_failure,
    warnings = paste(audited$audit$warnings, collapse = " | "),
    formula = audited$audit$formula, stringsAsFactors = FALSE)
  cumulative_row <- data.frame(country = country, unit = unit,
    measure = "cumulative_sustained_E3_P95_vs_median", lag = NA_integer_,
    estimate = unname(cumulative["estimate"]), variance = unname(cumulative["variance"]),
    se = unname(cumulative["se"]), RR = unname(cumulative["RR"]), lo = unname(cumulative["lo"]),
    hi = unname(cumulative["hi"]), n_complete = prepared$n,
    events_K80 = prepared$events[1], events_urinary = prepared$events[2],
    p95_E3 = ref$p95, median_E3 = ref$median, covariance = "hac28",
    interpretation = "Urinary-versus-K80 ratio of rate ratios for a sustained daily temperature contrast at the primary E3 P95 and median across the full 0-21 day history.",
    stringsAsFactors = FALSE)
  list(cumulative = cumulative_row, lag = do.call(rbind, lag_rows), audit = audit, failure = NULL)
}

empty_result <- function() list(
  cumulative = data.frame(country = character(), unit = character(), measure = character(), lag = integer(),
    estimate = numeric(), variance = numeric(), se = numeric(), RR = numeric(), lo = numeric(), hi = numeric(),
    n_complete = integer(), events_K80 = numeric(), events_urinary = numeric(), p95_E3 = numeric(),
    median_E3 = numeric(), covariance = character(), interpretation = character(), stringsAsFactors = FALSE),
  lag = data.frame(country = character(), unit = character(), lag = integer(), estimate = numeric(),
    variance = numeric(), se = numeric(), RR = numeric(), lo = numeric(), hi = numeric(),
    p95_E3 = numeric(), median_E3 = numeric(), stringsAsFactors = FALSE),
  audit = data.frame(), failure = data.frame(country = character(), unit = character(), phase = character(),
    reason = character(), n_complete = integer(), events_K80 = numeric(), events_urinary = numeric(),
    stringsAsFactors = FALSE))

run_country <- function(calendar, country, opts, manifest, expected_input_md5 = NULL) {
  units <- sort(unique(calendar$unit))
  manifest_country <- NULL
  if (!is.null(manifest) && nrow(manifest)) {
    manifest_country <- manifest[manifest$country == country &
      manifest$scenario == "primary_threeoutcome_E3_separate", , drop = FALSE]
    if (nrow(manifest_country)) {
      if (any(!manifest_country$unit %in% units)) {
        stop(sprintf("primary cache manifest has %s unit(s) absent from the canonical %s calendar",
          paste(unique(manifest_country$unit[!manifest_country$unit %in% units]), collapse = ", "), country), call. = FALSE)
      }
      if ("input_md5" %in% names(manifest_country) && !is.null(expected_input_md5)) {
        if (anyNA(manifest_country$input_md5) || any(manifest_country$input_md5 != expected_input_md5)) {
          stop(sprintf("primary cache manifest input hash does not match the current canonical %s calendar", country), call. = FALSE)
        }
      }
      dup <- split(manifest_country, manifest_country$unit)
      inconsistent <- vapply(dup, function(z) length(unique(paste(sprintf("%.17g", z$p95),
        sprintf("%.17g", z$median), sep = "|"))) > 1L, logical(1))
      if (any(inconsistent)) stop(sprintf("primary reference cache contains inconsistent unit entries: %s",
        paste(names(inconsistent)[inconsistent], collapse = ", ")), call. = FALSE)
      manifest_country <- manifest_country[!duplicated(manifest_country$unit), , drop = FALSE]
    }
  }
  group_idx <- split(seq_len(nrow(calendar)), factor(calendar$unit, levels = units), drop = TRUE)
  worker <- function(unit) {
    idx <- group_idx[[unit]]
    d <- calendar[idx]
    mref <- NULL
    if (!is.null(manifest_country)) {
      r <- manifest_country[manifest_country$unit == unit, , drop = FALSE]
      if (nrow(r)) mref <- list(p95 = r$p95[1], median = r$median[1])
    }
    tryCatch(fit_one_unit(d, country, unit, manifest_ref = mref), error = function(e) {
      list(cumulative = NULL, lag = NULL, audit = NULL,
           failure = data.frame(country = country, unit = unit, phase = "fit",
             reason = conditionMessage(e), n_complete = NA_integer_, events_K80 = NA_real_,
             events_urinary = NA_real_, stringsAsFactors = FALSE))
    })
  }
  message(sprintf("%s: DLNM timing sensitivity across %d units with %d worker(s)", country, length(units), opts$workers))
  fits <- if (.Platform$OS.type == "windows" || opts$workers == 1L || length(units) == 1L) {
    lapply(units, worker)
  } else {
    parallel::mclapply(units, worker, mc.cores = min(opts$workers, length(units)),
                       mc.preschedule = FALSE, mc.set.seed = FALSE)
  }
  bind_rows <- function(field, empty) {
    xs <- Filter(function(x) !is.null(x) && NROW(x) > 0L, lapply(fits, `[[`, field))
    if (!length(xs)) return(empty)
    as.data.frame(data.table::rbindlist(xs, fill = TRUE, use.names = TRUE), stringsAsFactors = FALSE)
  }
  z <- empty_result()
  res <- list(cumulative = bind_rows("cumulative", z$cumulative),
    lag = bind_rows("lag", z$lag), audit = bind_rows("audit", z$audit),
    failure = bind_rows("failure", z$failure))
  list(results = res, units = length(units), manifest_units = if (is.null(manifest_country)) 0L else nrow(manifest_country))
}

pool_dlnm <- function(cumulative, lag) {
  all <- list()
  add <- function(d, measure, lag_value) {
    if (!nrow(d)) return(invisible(NULL))
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(d$estimate, d$variance, method = method)
      all[[length(all) + 1L]] <<- data.frame(country = d$country[1], measure = measure,
        lag = lag_value, method = method, k = p$k, estimate = p$estimate, se = p$se,
        RR = p$RR, lo = p$RR95CI[1], hi = p$RR95CI[2], I2 = p$I2, tau2 = p$tau2,
        concordance = p$concordance,
        failures = if (!is.null(p$failures$reason)) p$failures$reason else "",
        stringsAsFactors = FALSE)
    }
  }
  for (country in unique(c(cumulative$country, lag$country))) {
    d <- cumulative[cumulative$country == country, , drop = FALSE]
    if (nrow(d)) add(d, "cumulative_sustained_E3_P95_vs_median", NA_integer_)
    lags <- lag[lag$country == country, , drop = FALSE]
    for (ll in sort(unique(lags$lag))) add(lags[lags$lag == ll, , drop = FALSE], "lag_specific_E3_P95_vs_median", ll)
  }
  if (!length(all)) return(data.frame(country = character(), measure = character(), lag = integer(),
    method = character(), k = integer(), estimate = numeric(), se = numeric(), RR = numeric(), lo = numeric(),
    hi = numeric(), I2 = numeric(), tau2 = numeric(), concordance = numeric(), failures = character(),
    stringsAsFactors = FALSE))
  do.call(rbind, all)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path) && unlink(path) != 0L) stop(sprintf("could not replace output '%s'", path), call. = FALSE)
  if (!file.rename(tmp, path)) stop(sprintf("failed to replace output '%s'", path), call. = FALSE)
}

validate_lag_history <- function() {
  base <- as.Date("2001-01-01")
  d <- base + 0:24
  d <- d[-5L]
  fake <- data.frame(date = d, tmean = seq_along(d))
  fake_basis <- new_temperature_basis(c(0, 1, 2, 3, 4, 5), df = 3L)
  fake$temp <- fake$tmean
  fake_data <- fake[c("date", "tmean")]
  cb <- build_crossbasis(fake_data, fake_basis)
  jan25 <- which(fake_data$date == base + 24)
  if (cb$complete_history[jan25] || any(is.finite(cb$features[jan25, ]))) {
    stop("manual date-gap fixture incorrectly treated a non-calendar row as a one-day lag", call. = FALSE)
  }
  const_basis <- new_temperature_basis(rep(10, 8) + seq_len(8) / 10, df = 3L)
  const_value <- 10.35
  const_data <- data.frame(date = base + 0:24, tmean = rep(const_value, 25L))
  const_cb <- build_crossbasis(const_data, const_basis)
  ll <- lag_basis_matrix()
  analytic <- as.numeric(.basis_at(const_basis, const_value)) %o% colSums(ll)
  analytic <- as.vector(analytic)
  observed <- const_cb$features[25L, ]
  if (!isTRUE(all.equal(unname(observed), unname(analytic), tolerance = 1e-10))) {
    stop("constant-temperature crossbasis failed the analytic outer-product check", call. = FALSE)
  }
  TRUE
}

write_metadata <- function(opts, inputs, countries, run_info) {
  if (!requireNamespace("jsonlite", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
    stop("jsonlite and digest are required to write auditable metadata", call. = FALSE)
  }
  file_hash <- function(path) if (file.exists(path)) digest::digest(file = path, algo = "sha256") else NA_character_
  metadata <- list(
    analysis = "prespecified secondary DLNM timing sensitivity",
    generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    country = countries,
    source_code = list(path = .script_path, sha256 = file_hash(.script_path)),
    repair_engine = list(path = file.path(.script_dir, "repair_engine.R"), sha256 = file_hash(file.path(.script_dir, "repair_engine.R"))),
    input_calendars = lapply(inputs, function(p) list(path = p, sha256 = file_hash(p))),
    primary_manifest = list(path = opts$manifest, sha256 = file_hash(opts$manifest),
      validated_rows = run_info$manifest_rows, units_without_manifest_rows = run_info$units_without_manifest),
    model = list(
      temperature = "Natural spline df=3 of raw daily tmean; unit-specific basis trained on common primary E3 paired dates with finite raw tmean.",
      reference = "Exact primary E3 P95 and median (type 7 quantile) calculated on common finite y_K80, y_K81, y_urinary, E3, and positive Brazil population dates; compared as absolute values on the raw daily temperature basis.",
      lag = "Natural spline df=3 with intercept on integer lag days 0:21, one internal knot sqrt(21), boundaries 0 and 21; 22 lag days.",
      crossbasis = "Nine features sum_l Btemp_j(T[t-l])*Blag_k(l); exact daily calendar lookup; incomplete histories excluded after the features are built.",
      formula = "quasi-Poisson y ~ disease*(CB1_1+...+CB3_3+dow+month+ns(time,12)); Brazil offset log(population), Mexico no offset.",
      eligibility = "At least 400 paired complete-history calendar days and at least 40 counts for each endpoint.",
      covariance = "Full calendar-day score covariance with Bartlett HAC lag 28; no quasi-dispersion multiplier.",
      cumulative_contrast = "Urinary-versus-K80 ratio of rate ratios for a sustained daily temperature contrast at the primary E3 P95 and median across the full 0-21 day history.",
      lag_specific = "Pointwise 95% intervals at each lag; separate pooled estimate reported for every lag and country.",
      joint_test = "Nine-degree-of-freedom Wald test of all temperature-basis interaction coefficients; rank/inversion failures are explicit in fit audit.",
      pooling = "DerSimonian-Laird and REML/Hartung-Knapp per country for cumulative and lag-specific contrasts; errors remain explicit.",
      interpretation_limit = "Cumulative contrast is a sustained-temperature association, distinct from a three-day E3 intervention. Early/late negative coefficients do not establish displacement or a biological mechanism.",
      history_rule = "Use canonical complete actual-calendar data only; first 21 calendar days are omitted when prehistory is unavailable; no row-lagging across gaps."),
    validation = list(constant_temperature_outer_product = TRUE, manual_calendar_gap_fixture = TRUE,
      equivalent_separate_glm_fitted_mean_tolerance = 1e-6,
      cumulative_hac_variance_matrix_identity = TRUE),
    runtime = list(R = R.version.string, workers = opts$workers,
      data_table = as.character(utils::packageVersion("data.table")),
      metafor = if (requireNamespace("metafor", quietly = TRUE)) as.character(utils::packageVersion("metafor")) else NA_character_),
    country_runs = run_info$country_runs
  )
  jsonlite::write_json(metadata, file.path(opts$outdir, "dlnm_analysis_metadata.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
}

run_dlnm <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_args(args)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for SHA-256 provenance", call. = FALSE)
  validate_lag_history()
  manifest <- NULL
  if (file.exists(opts$manifest)) {
    manifest <- utils::read.csv(opts$manifest, stringsAsFactors = FALSE,
      colClasses = c(country = "character", unit = "character", scenario = "character"))
    required_manifest <- c("country", "unit", "scenario", "p95", "median")
    if (length(setdiff(required_manifest, names(manifest)))) {
      stop("primary cache manifest is missing expected identity/reference fields", call. = FALSE)
    }
  } else message(sprintf("Primary cache manifest is not present yet; refs will be reconstructed exactly from canonical eligibility: %s", opts$manifest))
  countries <- if (opts$country == "all") c("BR", "MX") else opts$country
  input_paths <- c(BR = opts$brazil, MX = opts$mexico)
  full_results <- list(); run_info <- list(country_runs = list(), manifest_rows = 0L,
    units_without_manifest = 0L)
  for (country in countries) {
    calendar <- read_calendar(input_paths[[country]], country)
    current_input_md5 <- unname(tools::md5sum(input_paths[[country]]))
    result <- run_country(calendar, country, opts, manifest, expected_input_md5 = current_input_md5)
    full_results[[country]] <- result$results
    run_info$country_runs[[country]] <- list(units = result$units, units_with_primary_manifest = result$manifest_units,
      successful = nrow(result$results$cumulative), failed = nrow(result$results$failure))
    run_info$manifest_rows <- run_info$manifest_rows + result$manifest_units
    run_info$units_without_manifest <- run_info$units_without_manifest + result$units - result$manifest_units
    for (field in names(result$results)) {
      atomic_write_csv(result$results[[field]], file.path(opts$outdir,
        sprintf("dlnm_%s_%s.csv", country, field)))
    }
    pooled <- pool_dlnm(result$results$cumulative, result$results$lag)
    atomic_write_csv(pooled, file.path(opts$outdir, sprintf("dlnm_%s_pooled.csv", country)))
    rm(calendar); gc(verbose = FALSE)
  }
  write_metadata(opts, input_paths[countries], countries, run_info)
  message(sprintf("DLNM timing sensitivity completed; outputs: %s", opts$outdir))
  invisible(full_results)
}

if (sys.nframe() == 0L) run_dlnm()
