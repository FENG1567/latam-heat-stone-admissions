# Reusable statistical helpers for the 2026-11-12 repair analyses.
# The caller constructs lagged weather from the complete actual-calendar
# series before any outcome or analysis-period filtering. Fitting functions
# consume that already-constructed exposure and never row-lag weather.

#' Complete-calendar moving average (current day and preceding k - 1 days).
#' Any absent date or non-finite weather value in the window yields NA.
calendar_ma <- function(x, dates, k = 3L) {
  if (length(x) != length(dates)) stop("x and dates must have identical lengths", call. = FALSE)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1L || k != floor(k)) {
    stop("k must be a positive integer", call. = FALSE)
  }
  d <- as.Date(dates)
  if (anyNA(d)) stop("dates contain missing or invalid values", call. = FALSE)
  if (anyDuplicated(d)) stop("dates must be unique within a unit", call. = FALSE)
  if (!is.numeric(x)) stop("x must be numeric", call. = FALSE)
  names_by_date <- as.character(d)
  out <- rep(NA_real_, length(x))
  # Exact lookup by calendar date prevents gaps from being treated as adjacent.
  for (i in seq_along(d)) {
    window_dates <- d[i] - seq.int(0L, as.integer(k) - 1L)
    j <- match(as.character(window_dates), names_by_date)
    if (all(!is.na(j)) && all(is.finite(x[j]))) out[i] <- mean(x[j])
  }
  out
}

#' Natural-spline temperature basis. Preserve and pass this object to refits.
new_temperature_basis <- function(x, df = 3L, knots = NULL, boundary_knots = NULL) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop("x for the temperature basis must be a non-empty finite numeric vector", call. = FALSE)
  }
  if (is.null(boundary_knots)) boundary_knots <- range(x)
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 3L || df != floor(df)) {
    stop("df must be an integer of at least three", call. = FALSE)
  }
  if (length(boundary_knots) != 2L || any(!is.finite(boundary_knots)) ||
      boundary_knots[1] >= boundary_knots[2]) {
    stop("boundary_knots must be two increasing finite values", call. = FALSE)
  }
  b <- if (is.null(knots)) {
    splines::ns(x, df = as.integer(df), Boundary.knots = boundary_knots)
  } else {
    splines::ns(x, knots = knots, Boundary.knots = boundary_knots)
  }
  if (ncol(b) != 3L) stop("the repair models require a three-column temperature spline", call. = FALSE)
  b
}

.basis_at <- function(basis, x) {
  if (!inherits(basis, "ns")) stop("basis must be an ns object from new_temperature_basis()", call. = FALSE)
  b <- stats::predict(basis, newx = x)
  if (ncol(b) != 3L || any(!is.finite(b))) stop("temperature basis prediction failed", call. = FALSE)
  b
}

.validate_dates <- function(data, date_col) {
  if (!date_col %in% names(data)) stop(sprintf("missing date column '%s'", date_col), call. = FALSE)
  d <- as.Date(data[[date_col]])
  if (anyNA(d)) stop("date column contains missing or invalid dates", call. = FALSE)
  if (anyDuplicated(d)) stop("input must contain one row per actual calendar date within unit", call. = FALSE)
  d
}

.weather_targets <- function(data, exposure, basis, p95, median, mrt_bounds, df) {
  if (!exposure %in% names(data)) stop(sprintf("missing exposure column '%s'", exposure), call. = FALSE)
  x <- data[[exposure]]
  if (!is.numeric(x)) stop("exposure must be numeric", call. = FALSE)
  valid <- x[is.finite(x)]
  if (length(valid) < 4L || length(unique(valid)) < 2L) {
    stop("at least four finite exposure observations and two unique values are required", call. = FALSE)
  }
  if (is.null(basis)) basis <- new_temperature_basis(valid, df = df)
  if (is.null(p95)) p95 <- unname(stats::quantile(valid, 0.95, type = 7))
  if (is.null(median)) median <- stats::median(valid)
  if (is.null(mrt_bounds)) mrt_bounds <- unname(stats::quantile(valid, c(0.10, 0.60), type = 7))
  if (!is.numeric(p95) || length(p95) != 1L || !is.finite(p95) ||
      !is.numeric(median) || length(median) != 1L || !is.finite(median)) {
    stop("p95 and median reference overrides must be finite scalars", call. = FALSE)
  }
  if (!is.numeric(mrt_bounds) || length(mrt_bounds) != 2L ||
      any(!is.finite(mrt_bounds)) || mrt_bounds[1] >= mrt_bounds[2]) {
    stop("mrt_bounds must be two increasing finite values", call. = FALSE)
  }
  list(basis = basis, p95 = as.numeric(p95), median = as.numeric(median),
       mrt_bounds = as.numeric(mrt_bounds), exposure = x)
}

.fit_glm_audited <- function(formula, data) {
  warning_text <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(formula, data = data, family = stats::quasipoisson(link = "log"),
                 na.action = stats::na.fail, x = TRUE, y = TRUE, model = TRUE,
                 control = stats::glm.control(maxit = 50, epsilon = 1e-8)),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) stop(sprintf("quasi-Poisson fit failed: %s", conditionMessage(e)), call. = FALSE)
  )
  X <- stats::model.matrix(fit)
  rank <- qr(X)$rank
  aliases <- names(stats::coef(fit))[is.na(stats::coef(fit))]
  audit <- list(converged = isTRUE(fit$converged), boundary = isTRUE(fit$boundary),
                rank = rank, columns = ncol(X), aliased = aliases, warnings = warning_text,
                n = nrow(X), events = sum(fit$y), formula = paste(deparse(formula), collapse = " "))
  if (!audit$converged) stop(sprintf("quasi-Poisson model did not converge (rank %d/%d); warnings: %s",
                                     rank, ncol(X), paste(warning_text, collapse = " | ")), call. = FALSE)
  if (rank != ncol(X) || length(aliases)) {
    stop(sprintf("rank-deficient model (rank %d/%d; aliased: %s)", rank, ncol(X),
                 paste(aliases, collapse = ", ")), call. = FALSE)
  }
  if (any(!is.finite(stats::coef(fit))) || any(!is.finite(fit$fitted.values))) {
    stop("quasi-Poisson model returned non-finite coefficients or means", call. = FALSE)
  }
  list(fit = fit, X = X, audit = audit)
}

# Bartlett HAC meat as rolling cross-products. Zero-filled actual calendar gaps
# preserve calendar lag distances. Padding by L days on both ends makes each
# date appear in L+1 windows, so dividing by L+1 gives exactly weight
# (L + 1 - lag)/(L + 1) to each lag pair, including boundary dates.
.bartlett_hac_meat <- function(U, dates, lag) {
  if (length(lag) != 1L || is.na(lag) || lag < 0L || lag != as.integer(lag)) {
    stop("HAC lag must be a non-negative integer", call. = FALSE)
  }
  lag <- as.integer(lag)
  dates <- as.Date(dates)
  if (nrow(U) != length(dates) || anyNA(dates) || !length(dates)) {
    stop("daily score/date dimensions are invalid", call. = FALSE)
  }
  calendar <- seq(min(dates), max(dates), by = "day")
  full <- matrix(0, nrow = length(calendar), ncol = ncol(U),
                 dimnames = list(as.character(calendar), colnames(U)))
  full[as.integer(dates - min(dates)) + 1L, ] <- U
  padded <- rbind(matrix(0, nrow = lag, ncol = ncol(U)), full,
                  matrix(0, nrow = lag, ncol = ncol(U)))
  cumulative <- matrix(0, nrow = nrow(padded) + 1L, ncol = ncol(U))
  for (j in seq_len(ncol(U))) cumulative[, j] <- c(0, cumsum(padded[, j]))
  starts <- seq_len(nrow(full) + lag) - 1L
  ends <- starts + lag + 1L
  rolling <- cumulative[ends + 1L, , drop = FALSE] - cumulative[starts + 1L, , drop = FALSE]
  crossprod(rolling) / (lag + 1L)
}

# Score sandwich covariances using actual calendar dates. The Bartlett HAC sums
# only exact t versus t-l calendar matches; a missing day is never treated as
# the next adjacent observation. No quasi-dispersion multiplier is applied to
# HC0 or HAC matrices.
.score_covariances <- function(fit, X, y, mu, dates, hac_lags = c(14L, 28L),
                               contrast_vectors = NULL, full_covariance = TRUE) {
  if (length(y) != nrow(X) || length(mu) != nrow(X) || length(dates) != nrow(X)) {
    stop("internal score dimensions do not match the model design", call. = FALSE)
  }
  dates <- as.Date(dates)
  if (anyNA(dates)) stop("score dates contain invalid values", call. = FALSE)
  score <- X * as.numeric(y - mu)
  unique_dates <- sort(unique(dates))
  # One pass over rows; rowsum avoids repeatedly scanning the full score matrix.
  U0 <- rowsum(score, group = as.character(dates), reorder = FALSE)
  ord <- match(as.character(unique_dates), rownames(U0))
  if (anyNA(ord)) stop("daily score aggregation lost a calendar date", call. = FALSE)
  U <- U0[ord, , drop = FALSE]
  rownames(U) <- as.character(unique_dates)
  H <- crossprod(X, X * as.numeric(mu))
  Hinv <- tryCatch(solve(H), error = function(e) {
    stop(sprintf("working Hessian inversion failed: %s", conditionMessage(e)), call. = FALSE)
  })
  sym <- function(z) (z + t(z)) / 2
  quasi <- stats::vcov(fit)
  if (!full_covariance && is.null(contrast_vectors)) {
    stop("contrast_vectors are required when full_covariance = FALSE", call. = FALSE)
  }
  if (full_covariance) {
    hc0 <- sym(Hinv %*% crossprod(U) %*% Hinv)
  } else {
    hc0 <- NULL
  }
  hac <- list()
  projected <- if (full_covariance) NULL else list()
  influences <- list()
  for (lag in hac_lags) {
    if (length(lag) != 1L || is.na(lag) || lag < 0L || lag != as.integer(lag)) {
      stop("HAC lags must be non-negative integers", call. = FALSE)
    }
    if (full_covariance) {
      meat <- .bartlett_hac_meat(U, unique_dates, lag)
      hac[[paste0("hac", lag)]] <- sym(Hinv %*% meat %*% Hinv)
    } else {
      for (nm in names(contrast_vectors)) {
        L <- contrast_vectors[[nm]]
        if (length(L) != ncol(X)) stop(sprintf("contrast '%s' has the wrong design length", nm), call. = FALSE)
        if (is.null(projected[[nm]])) {
          influences[[nm]] <- as.numeric(U %*% (Hinv %*% L))
          influence <- influences[[nm]]
          projected[[nm]] <- list(day_hc0 = sum(influence^2))
        }
        influence <- influences[[nm]]
        value <- as.numeric(.bartlett_hac_meat(matrix(influence, ncol = 1L),
                                                unique_dates, lag)[1L, 1L])
        projected[[nm]][[paste0("hac", lag)]] <- value
      }
    }
  }
  list(quasi = quasi, day_hc0 = hc0, hac = hac, projected = projected,
       X = X, mu = as.numeric(mu), y = as.numeric(y), dates = dates,
       Hinv = Hinv, daily_scores = U, daily_dates = unique_dates,
       row_scores = score)
}

.design_at <- function(fit, template, endpoint = NULL, basis, temp_value) {
  nd <- template[1L, , drop = FALSE]
  if (!is.null(endpoint)) nd$disease <- factor(endpoint, levels = levels(template$disease))
  b <- as.numeric(.basis_at(basis, temp_value))
  nd$B1 <- b[1L]; nd$B2 <- b[2L]; nd$B3 <- b[3L]
  tt <- stats::delete.response(stats::terms(fit))
  x <- stats::model.matrix(tt, data = nd, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  x[, names(stats::coef(fit)), drop = FALSE]
}

.contrast_summary <- function(L, fit, covariances, projected = NULL) {
  beta <- stats::coef(fit)
  est <- as.numeric(L %*% beta)
  out <- lapply(covariances, function(V) {
    v <- as.numeric(t(L) %*% V %*% L)
    if (!is.finite(v) || v <= 0) stop("contrast variance is non-positive or non-finite", call. = FALSE)
    se <- sqrt(v)
    c(estimate = est, variance = v, se = se, RR = exp(est),
      lo = exp(est - stats::qnorm(0.975) * se), hi = exp(est + stats::qnorm(0.975) * se))
  })
  if (!is.null(projected)) for (nm in names(projected)) {
    v <- as.numeric(projected[[nm]])
    if (length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(sprintf("projected contrast variance '%s' is non-positive or non-finite", nm), call. = FALSE)
    }
    se <- sqrt(v)
    out[[nm]] <- c(estimate = est, variance = v, se = se, RR = exp(est),
                   lo = exp(est - stats::qnorm(0.975) * se),
                   hi = exp(est + stats::qnorm(0.975) * se))
  }
  list(estimate = est, contrast_vector = setNames(as.numeric(L), names(beta)), by_covariance = out)
}

#' Single-endpoint P95 versus constrained P10-P60 minimum-risk temperature.
#' Input exposure must already have been built over the full calendar history.
fit_main_effect <- function(data, outcome, exposure = "E", date = "date", offset = NULL,
                            basis = NULL, p95 = NULL, median = NULL, mrt_bounds = NULL,
                            time_col = NULL, df = 3L, hac_lags = c(14L, 28L)) {
  if (!is.data.frame(data) || !outcome %in% names(data)) stop("data/outcome column is invalid", call. = FALSE)
  dates <- .validate_dates(data, date)
  targ <- .weather_targets(data, exposure, basis, p95, median, mrt_bounds, df)
  y0 <- data[[outcome]]
  if (!is.numeric(y0)) stop("outcome must be numeric counts", call. = FALSE)
  if (any(is.finite(y0) & y0 < 0)) stop("outcome contains negative counts", call. = FALSE)
  if (is.null(time_col)) time <- as.numeric(dates - min(dates)) / 365.25
  else {
    if (!time_col %in% names(data)) stop(sprintf("missing time column '%s'", time_col), call. = FALSE)
    time <- data[[time_col]]
    if (!is.numeric(time)) stop("time column must be numeric", call. = FALSE)
  }
  if (is.null(offset)) off <- rep(0, nrow(data)) else {
    if (length(offset) != 1L || !is.character(offset) || !offset %in% names(data)) stop("offset must be NULL or a population column name", call. = FALSE)
    pop <- data[[offset]]
    if (!is.numeric(pop) || any(is.finite(pop) & pop <= 0)) stop("offset population must be positive numeric", call. = FALSE)
    off <- log(pop)
  }
  keep <- is.finite(y0) & is.finite(targ$exposure) & is.finite(time) & is.finite(off)
  if (!any(keep)) stop("no complete rows remain for main-effect fit", call. = FALSE)
  d <- dates[keep]; e <- targ$exposure[keep]; y <- y0[keep]; ti <- time[keep]
  dow <- factor(as.integer(format(d, "%u")), levels = 1:7)
  month <- factor(as.integer(format(d, "%m")), levels = 1:12)
  bb <- .basis_at(targ$basis, e)
  dfm <- data.frame(y = y, B1 = bb[, 1], B2 = bb[, 2], B3 = bb[, 3],
                    dow = dow, month = month, time = ti, off = off[keep], date_actual = d)
  formula <- if (is.null(offset)) y ~ B1 + B2 + B3 + dow + month + splines::ns(time, 12) else
    y ~ B1 + B2 + B3 + dow + month + splines::ns(time, 12) + offset(off)
  audited <- .fit_glm_audited(formula, dfm)
  fit <- audited$fit; X <- audited$X
  # Use a fine deterministic grid and locally refine every detected grid minimum.
  grid <- seq(targ$mrt_bounds[1], targ$mrt_bounds[2], length.out = 2001L)
  Bgrid <- .basis_at(targ$basis, grid)
  beta <- stats::coef(fit)[c("B1", "B2", "B3")]
  curve <- as.numeric(Bgrid %*% beta)
  mins <- which(curve <= c(Inf, curve[-length(curve)]) & curve <= c(curve[-1L], Inf))
  candidates <- c(grid[1], grid[length(grid)], grid[mins])
  if (length(mins)) for (ii in mins) if (ii > 1L && ii < length(grid)) {
    opt <- stats::optimize(function(z) sum(.basis_at(targ$basis, z) * beta),
                           interval = c(grid[ii - 1L], grid[ii + 1L]))
    candidates <- c(candidates, opt$minimum)
  }
  risks <- vapply(candidates, function(z) sum(.basis_at(targ$basis, z) * beta), numeric(1))
  mrt <- candidates[which.min(risks)]
  X95 <- .design_at(fit, dfm, basis = targ$basis, temp_value = targ$p95)
  Xmrt <- .design_at(fit, dfm, basis = targ$basis, temp_value = mrt)
  L <- as.numeric(X95 - Xmrt); names(L) <- colnames(X)
  V <- .score_covariances(fit, X, y = fit$y, mu = fit$fitted.values,
                         dates = dfm$date_actual, hac_lags = hac_lags)
  covs <- c(list(quasi = V$quasi, day_hc0 = V$day_hc0), V$hac)
  contrast <- .contrast_summary(L, fit, covs)
  dB <- as.numeric(.basis_at(targ$basis, targ$p95) - .basis_at(targ$basis, mrt))
  list(outcome = outcome, exposure = exposure, fit = fit, design = X,
       covariance = covs, quasi_vcov = V$quasi, paired_day_hc0 = V$day_hc0,
       hac_vcov = V$hac, contrast = contrast, dB = dB, p95 = targ$p95,
       median = targ$median, mrt = mrt, mrt_bounds = targ$mrt_bounds,
       basis = targ$basis, data = dfm, row_data = dfm,
       score = V, audit = audited$audit,
       prediction = list(grid = grid, link = curve, temperature_range = range(e)))
}

#' Stacked quasi-Poisson disease contrast with paired-day HC0 and calendar HAC.
#' `outcomes` are ordered: each pair is reported as later outcome versus earlier.
fit_stacked_contrast <- function(data, outcomes, exposure = "E", date = "date", offset = NULL,
                                 basis = NULL, p95 = NULL, median = NULL,
                                 nuisance = c("separate", "shared", "yearmonth", "flexible"),
                                 time_col = NULL, df = 3L, hac_lags = c(14L, 28L),
                                 full_covariance = TRUE, holiday_col = "national_holiday") {
  nuisance <- match.arg(nuisance)
  if (!is.logical(full_covariance) || length(full_covariance) != 1L || is.na(full_covariance)) {
    stop("full_covariance must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.data.frame(data) || length(outcomes) < 2L || anyDuplicated(outcomes) ||
      any(!outcomes %in% names(data))) stop("data/outcome columns are invalid", call. = FALSE)
  dates <- .validate_dates(data, date)
  targ <- .weather_targets(data, exposure, basis, p95, median, NULL, df)
  if (is.null(time_col)) time <- as.numeric(dates - min(dates)) / 365.25
  else {
    if (!time_col %in% names(data)) stop(sprintf("missing time column '%s'", time_col), call. = FALSE)
    time <- data[[time_col]]
    if (!is.numeric(time)) stop("time column must be numeric", call. = FALSE)
  }
  if (is.null(offset)) off <- rep(0, nrow(data)) else {
    if (length(offset) != 1L || !is.character(offset) || !offset %in% names(data)) stop("offset must be NULL or a population column name", call. = FALSE)
    pop <- data[[offset]]
    if (!is.numeric(pop) || any(is.finite(pop) & pop <= 0)) stop("offset population must be positive numeric", call. = FALSE)
    off <- log(pop)
  }
  ys <- lapply(outcomes, function(nm) data[[nm]])
  if (any(!vapply(ys, is.numeric, logical(1)))) stop("all outcomes must be numeric counts", call. = FALSE)
  if (any(vapply(ys, function(y) any(is.finite(y) & y < 0), logical(1)))) {
    stop("an outcome contains negative counts", call. = FALSE)
  }
  names(ys) <- outcomes
  holiday <- if (holiday_col %in% names(data)) suppressWarnings(as.numeric(data[[holiday_col]])) else
    rep(NA_real_, nrow(data))
  if (nuisance == "flexible" && (!holiday_col %in% names(data) || any(!is.finite(holiday)) ||
      any(!holiday %in% c(0, 1)))) {
    stop(sprintf("flexible nuisance requires a complete 0/1 '%s' holiday indicator", holiday_col), call. = FALSE)
  }
  keep <- is.finite(targ$exposure) & is.finite(time) & is.finite(off)
  if (nuisance == "flexible") keep <- keep & is.finite(holiday)
  for (yy in ys) keep <- keep & is.finite(yy)
  if (!any(keep)) stop("no complete paired rows remain for stacked fit", call. = FALSE)
  d <- dates[keep]; e <- targ$exposure[keep]; ti <- time[keep]
  dow <- factor(as.integer(format(d, "%u")), levels = 1:7)
  month <- factor(as.integer(format(d, "%m")), levels = 1:12)
  ym <- factor(format(d, "%Y-%m"))
  national_holiday <- if (nuisance == "flexible") holiday[keep] else rep(0, length(d))
  analysis_years <- length(unique(format(d, "%Y")))
  flexible_time_df <- 6L * analysis_years
  bb <- .basis_at(targ$basis, e)
  common <- data.frame(B1 = bb[, 1], B2 = bb[, 2], B3 = bb[, 3],
                       dow = dow, month = month, ym = ym, time = ti, national_holiday = national_holiday,
                       off = off[keep],
                       date_actual = d)
  st <- do.call(rbind, lapply(seq_along(outcomes), function(j) {
    z <- common
    z$y <- ys[[j]][keep]
    z$disease <- factor(outcomes[j], levels = outcomes)
    z
  }))
  rownames(st) <- NULL
  formula <- if (nuisance == "flexible") {
    offset_term <- if (is.null(offset)) "" else " + offset(off)"
    stats::as.formula(sprintf(
      "y ~ disease * (B1 + B2 + B3 + dow + splines::ns(time, %d) + national_holiday)%s",
      flexible_time_df, offset_term
    ), env = environment())
  } else switch(nuisance,
    separate = if (is.null(offset)) y ~ disease * (B1 + B2 + B3 + dow + month + splines::ns(time, 12)) else
      y ~ disease * (B1 + B2 + B3 + dow + month + splines::ns(time, 12)) + offset(off),
    shared = if (is.null(offset)) y ~ disease * (B1 + B2 + B3) + dow + month + splines::ns(time, 12) else
      y ~ disease * (B1 + B2 + B3) + dow + month + splines::ns(time, 12) + offset(off),
    yearmonth = if (is.null(offset)) y ~ disease * (B1 + B2 + B3 + dow + ym) else
      y ~ disease * (B1 + B2 + B3 + dow + ym) + offset(off))
  audited <- .fit_glm_audited(formula, st)
  fit <- audited$fit; X <- audited$X
  exposure_effects <- lapply(outcomes, function(endpoint) {
    Xhi <- .design_at(fit, st, endpoint = endpoint, basis = targ$basis, temp_value = targ$p95)
    Xmed <- .design_at(fit, st, endpoint = endpoint, basis = targ$basis, temp_value = targ$median)
    as.numeric(Xhi - Xmed)
  })
  names(exposure_effects) <- outcomes
  pairs <- utils::combn(seq_along(outcomes), 2L, simplify = FALSE)
  contrast_vectors <- lapply(pairs, function(p) exposure_effects[[p[2]]] - exposure_effects[[p[1]]])
  contrast_labels <- vapply(pairs, function(p) paste0(outcomes[p[2]], "_vs_", outcomes[p[1]]), character(1))
  names(contrast_vectors) <- contrast_labels
  V <- .score_covariances(fit, X, y = fit$y, mu = fit$fitted.values,
                         dates = st$date_actual, hac_lags = hac_lags,
                         contrast_vectors = contrast_vectors, full_covariance = full_covariance)
  covs <- c(list(quasi = V$quasi), if (full_covariance) list(day_hc0 = V$day_hc0), V$hac)
  contrasts <- lapply(pairs, function(p) {
    i <- p[1]; j <- p[2]
    L <- exposure_effects[[j]] - exposure_effects[[i]]
    names(L) <- colnames(X)
    label <- paste0(outcomes[j], "_vs_", outcomes[i])
    cc <- .contrast_summary(L, fit, covs, projected = if (full_covariance) NULL else V$projected[[label]])
    c(list(label = label, numerator = outcomes[j], denominator = outcomes[i],
           dB = as.numeric(.basis_at(targ$basis, targ$p95) - .basis_at(targ$basis, targ$median)),
           coefcontrast = cc$estimate, coefficient_contrast = cc,
           contrast_variance = lapply(cc$by_covariance, `[[`, "variance")),
      list(contrast_vector = L))
  })
  names(contrasts) <- vapply(contrasts, `[[`, character(1), "label")
  dB <- as.numeric(.basis_at(targ$basis, targ$p95) - .basis_at(targ$basis, targ$median))
  coefcontrast <- vapply(contrasts, `[[`, numeric(1), "coefcontrast")
  list(outcomes = outcomes, exposure = exposure, nuisance = nuisance, fit = fit,
       design = X, covariance = covs, joint_covariance = covs,
       quasi_vcov = V$quasi, paired_day_hc0 = V$day_hc0, hac_vcov = V$hac,
       covariance_scope = if (full_covariance) "full" else "contrast_only",
       primary_covariance = "hac14", sensitivity_covariances = c("hac28", "quasi", "day_hc0"),
       dB = dB, coefcontrast = coefcontrast, contrasts = contrasts,
       exposure_effects = exposure_effects, p95 = targ$p95, median = targ$median,
       flexible_time_df = if (nuisance == "flexible") flexible_time_df else NA_integer_,
       national_holiday_days = if (nuisance == "flexible") sum(national_holiday == 1) else NA_integer_,
       basis = targ$basis, data = st, row_data = st, score = V,
       audit = audited$audit, retained_dates = sort(unique(d)))
}

#' Pool log rate ratios. REML/Hartung-Knapp requires metafor; failures are returned.
pool_meta <- function(estimate, variance, method = c("DL", "REML_HK")) {
  method <- match.arg(method)
  if (length(estimate) != length(variance)) stop("estimate and variance lengths differ", call. = FALSE)
  if (!is.numeric(estimate) || !is.numeric(variance)) stop("estimate and variance must be numeric", call. = FALSE)
  ok <- is.finite(estimate) & is.finite(variance) & variance > 0
  failures <- which(!ok)
  k <- sum(ok)
  result <- list(method = method, k = k, I2 = NA_real_, tau2 = NA_real_, estimate = NA_real_,
                 se = NA_real_, RR = NA_real_, RR95CI = c(NA_real_, NA_real_),
                 concordance = if (k) mean(estimate[ok] > 0) else NA_real_,
                 used = which(ok), failed_rows = failures, failure = NULL,
                 failures = list(rows = failures, reason = NULL))
  if (k < 2L) {
    result$failure <- "fewer than two finite estimates with positive variances"
    result$failures$reason <- result$failure
    return(result)
  }
  yi <- estimate[ok]; vi <- variance[ok]
  if (method == "DL") {
    wi <- 1 / vi; ybar <- sum(wi * yi) / sum(wi)
    Q <- sum(wi * (yi - ybar)^2)
    C <- sum(wi) - sum(wi^2) / sum(wi)
    if (!is.finite(C) || C <= 0) {
      result$failure <- "DerSimonian-Laird denominator is non-positive"
      result$failures$reason <- result$failure
      return(result)
    }
    tau2 <- max(0, (Q - (k - 1)) / C)
    wr <- 1 / (vi + tau2)
    est <- sum(wr * yi) / sum(wr); se <- sqrt(1 / sum(wr))
    I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
    crit <- stats::qnorm(0.975)
    lo <- est - crit * se; hi <- est + crit * se
  } else {
    if (!requireNamespace("metafor", quietly = TRUE)) {
      result$failure <- "metafor is unavailable; REML/Hartung-Knapp was not fitted"
      result$failures$reason <- result$failure
      return(result)
    }
    mf <- tryCatch(metafor::rma.uni(yi = yi, vi = vi, method = "REML", test = "knha"),
                   error = function(e) e)
    if (inherits(mf, "error")) {
      result$failure <- paste("metafor REML/Hartung-Knapp failed:", conditionMessage(mf))
      result$failures$reason <- result$failure
      return(result)
    }
    est <- as.numeric(mf$b[1]); se <- as.numeric(mf$se[1]); tau2 <- as.numeric(mf$tau2)
    I2 <- as.numeric(mf$I2)
    lo <- as.numeric(mf$ci.lb); hi <- as.numeric(mf$ci.ub)
  }
  result$estimate <- est; result$se <- se; result$tau2 <- tau2; result$I2 <- I2
  result$RR <- exp(est); result$RR95CI <- exp(c(lo, hi))
  if (length(failures)) result$failures$reason <- "non-finite estimates or non-positive variances excluded"
  result
}
