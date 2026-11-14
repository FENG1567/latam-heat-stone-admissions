#!/usr/bin/env Rscript

# Extract fixed-reference marginal RRs and their direct urinary-vs-K80 RRR
# from the existing successful joint-fit caches. No model is refitted here.
#
# Usage:
#   Rscript extract_joint_marginals.R [project_root] [output_dir]

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1L) args[[1L]] else "/data/admin/fengshuo/latam_stone_repair_20261112"
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(base_dir, "revision_20261113", "results", "marginals_curves")
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
engine_path <- file.path(base_dir, "scripts", "repair_engine.R")
if (!file.exists(engine_path)) stopf("repair_engine.R is missing: %s", engine_path)
source(engine_path, local = FALSE)
if (!requireNamespace("splines", quietly = TRUE)) stop("splines is required for cached curve validation", call. = FALSE)
if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor is required for REML-HK summaries", call. = FALSE)

sha256_file <- function(path) {
  exe <- Sys.which("sha256sum")
  if (!nzchar(exe)) stop("sha256sum is required for cache provenance", call. = FALSE)
  z <- suppressWarnings(system2(exe, shQuote(path), stdout = TRUE, stderr = TRUE))
  status <- attr(z, "status")
  if (!is.null(status) && status != 0L) stopf("sha256sum failed for %s: %s", path, paste(z, collapse = " "))
  tolower(sub("\\s+.*$", "", z[[1L]]))
}

md5_file <- function(path) unname(tools::md5sum(path)[[1L]])
close_enough <- function(x, y, tol = 1e-8) {
  length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= tol * pmax(1, abs(x), abs(y)))
}
optional_character <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) NA_character_ else as.character(x[[1L]])
}
as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
  as.Date(as.character(x))
}

bartlett_hac_matrix <- function(psi, dates, lag) {
  psi <- as.matrix(psi)
  dates <- as_date(dates)
  if (!nrow(psi) || nrow(psi) != length(dates) || anyNA(dates) || anyDuplicated(dates)) {
    stop("daily influence rows and unique dates must align", call. = FALSE)
  }
  if (any(!is.finite(psi))) stop("daily influence contains non-finite values", call. = FALSE)
  ord <- order(dates)
  dates <- dates[ord]
  psi <- psi[ord, , drop = FALSE]
  calendar <- seq(min(dates), max(dates), by = "day")
  full <- matrix(0, nrow = length(calendar), ncol = ncol(psi),
                 dimnames = list(as.character(calendar), colnames(psi)))
  full[as.integer(dates - min(dates)) + 1L, ] <- psi
  out <- crossprod(full)
  if (lag > 0L && nrow(full) > 1L) {
    for (ell in seq_len(min(as.integer(lag), nrow(full) - 1L))) {
      left <- full[(ell + 1L):nrow(full), , drop = FALSE]
      right <- full[seq_len(nrow(full) - ell), , drop = FALSE]
      cross <- crossprod(left, right)
      out <- out + 2 * (1 - ell / (lag + 1)) * (cross + t(cross)) / 2
    }
  }
  (out + t(out)) / 2
}

psd_audit <- function(V) {
  eig <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
  tol <- 1e-8 * max(1, max(abs(eig)))
  list(min_eigenvalue = min(eig), max_eigenvalue = max(eig),
       eigenvalue_tolerance = tol,
       psd_status = if (min(eig) >= -tol) "PSD_within_numeric_tolerance" else "not_PSD")
}

atomic_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stopf("failed to atomically write %s", path)
}

read_selected_manifest <- function(path, country, scenario) {
  if (!file.exists(path)) stopf("manifest does not exist: %s", path)
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = "character", na.strings = c("", "NA"))
  required <- c("country", "unit", "scenario", "cache_path", "run_id", "p95", "median", "n_complete")
  missing <- setdiff(required, names(m))
  if (length(missing)) stopf("manifest %s lacks: %s", path, paste(missing, collapse = ", "))
  m$country <- toupper(m$country)
  m <- m[m$country == country & m$scenario == scenario, , drop = FALSE]
  if (!nrow(m)) stopf("manifest has no rows for %s/%s: %s", country, scenario, path)
  if (anyDuplicated(m$unit)) stopf("manifest repeats units for %s/%s", country, scenario)
  for (nm in c("p95", "median", "n_complete")) m[[nm]] <- as.numeric(m[[nm]])
  m$cache_path <- normalizePath(m$cache_path, winslash = "/", mustWork = TRUE)
  m$.manifest_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  m
}

load_unit <- function(row) {
  path <- row$cache_path
  cache <- readRDS(path)
  actual_sha <- sha256_file(path)
  expected_sha <- row$cache_sha256
  if (!is.null(expected_sha) && !is.na(expected_sha) && nzchar(expected_sha) &&
      !identical(tolower(expected_sha), actual_sha)) stopf("cache SHA-256 mismatch for %s/%s", row$country, row$unit)
  ef <- cache$engine_fit
  if (is.null(ef)) ef <- cache$fit
  if (is.null(ef) || !is.list(ef)) stopf("cache has no engine fit for %s/%s", row$country, row$unit)
  identity <- c(country = toupper(as.character(cache$country)), unit = as.character(cache$unit),
                scenario = as.character(cache$scenario), run_id = as.character(cache$run_id))
  if (!identical(unname(identity["country"]), row$country) ||
      !identical(unname(identity["unit"]), as.character(row$unit)) ||
      !identical(unname(identity["scenario"]), row$scenario) ||
      !identical(unname(identity["run_id"]), as.character(row$run_id))) {
    stopf("cache identity differs from manifest for %s/%s", row$country, row$unit)
  }
  if (!is.finite(cache$p95) || !is.finite(cache$median) ||
      !close_enough(c(cache$p95, cache$median), c(row$p95, row$median), 1e-10)) {
    stopf("cached P95/shared median differ from manifest for %s/%s", row$country, row$unit)
  }
  if (!is.null(row$input_md5) && !is.na(row$input_md5) && nzchar(row$input_md5)) {
    if (is.null(cache$input_md5) || !identical(tolower(as.character(cache$input_md5)), tolower(row$input_md5))) {
      stopf("cache input MD5 metadata mismatch for %s/%s", row$country, row$unit)
    }
  }
  if (!is.null(row$input_sha256) && !is.na(row$input_sha256) && nzchar(row$input_sha256)) {
    hashes <- cache$input_hashes_sha256
    canonical <- if (!is.null(hashes)) hashes[["canonical"]] else NULL
    if (is.null(canonical) || !identical(tolower(as.character(canonical)), tolower(row$input_sha256))) {
      stopf("cache canonical input SHA-256 metadata mismatch for %s/%s", row$country, row$unit)
    }
  }
  if (!is.null(row$engine_sha256) && !is.na(row$engine_sha256) && nzchar(row$engine_sha256)) {
    actual_engine <- cache$engine_hash_sha256
    if (is.null(actual_engine) || !identical(tolower(as.character(actual_engine)), tolower(row$engine_sha256))) {
      stopf("cache engine SHA-256 metadata mismatch for %s/%s", row$country, row$unit)
    }
  }
  beta <- stats::coef(ef$fit)
  if (is.null(names(beta)) || any(!is.finite(beta))) stopf("invalid coefficient vector for %s/%s", row$country, row$unit)
  covariance <- ef$covariance
  if (!identical(ef$covariance_scope, "full") || is.null(covariance$hac14) || is.null(covariance$hac28)) {
    stopf("cache lacks full HAC14/HAC28 joint coefficient covariance for %s/%s", row$country, row$unit)
  }
  normalize_covariance <- function(V, label) {
    if (!is.matrix(V) || nrow(V) != length(beta) || ncol(V) != length(beta)) stopf("%s covariance dimensions disagree", label)
    if (!is.null(rownames(V)) && !is.null(colnames(V))) {
      if (!setequal(rownames(V), names(beta)) || !setequal(colnames(V), names(beta))) stopf("%s covariance names disagree", label)
      V <- V[names(beta), names(beta), drop = FALSE]
    }
    if (any(!is.finite(V))) stopf("%s covariance contains non-finite values", label)
    if (max(abs(V - t(V))) > 1e-8 * max(1, max(abs(V)))) stopf("%s covariance is not symmetric", label)
    V
  }
  V14 <- normalize_covariance(covariance$hac14, "HAC14")
  V28 <- normalize_covariance(covariance$hac28, "HAC28")
  score <- ef$score
  Hinv <- score$Hinv
  U <- score$daily_scores
  dates <- as_date(score$daily_dates)
  if (!is.matrix(Hinv) || !is.matrix(U) || length(dates) != nrow(U) || ncol(U) != length(beta) ||
      nrow(Hinv) != length(beta) || ncol(Hinv) != length(beta) || anyNA(dates) || anyDuplicated(dates)) {
    stopf("score/Hessian/date dimensions invalid for %s/%s", row$country, row$unit)
  }
  if (!is.null(rownames(Hinv)) && !is.null(colnames(Hinv))) Hinv <- Hinv[names(beta), names(beta), drop = FALSE]
  if (!is.null(colnames(U))) U <- U[, names(beta), drop = FALSE]
  if (max(abs(Hinv - t(Hinv))) > 1e-8 * max(1, max(abs(Hinv)))) stopf("Hinv is not symmetric for %s/%s", row$country, row$unit)
  if (any(!is.finite(Hinv)) || any(!is.finite(U))) stopf("non-finite cached score data for %s/%s", row$country, row$unit)

  effects <- ef$exposure_effects
  if (is.null(effects$y_K80) || is.null(effects$y_urinary)) stopf("cache lacks fixed-reference endpoint contrast vectors")
  direct_obj <- ef$contrasts$y_urinary_vs_y_K80
  if (is.null(direct_obj$contrast_vector)) stopf("cache lacks stored direct urinary-vs-K80 contrast vector")
  L <- cbind(K80_P95_vs_shared_median = as.numeric(effects$y_K80),
             urinary_P95_vs_shared_median = as.numeric(effects$y_urinary),
             direct_urinary_vs_K80_RRR = as.numeric(direct_obj$contrast_vector))
  rownames(L) <- names(beta)
  if (nrow(L) != length(beta) || any(!is.finite(L))) stopf("invalid contrast matrix for %s/%s", row$country, row$unit)
  direct_diff <- max(abs(L[, "direct_urinary_vs_K80_RRR"] -
                           (L[, "urinary_P95_vs_shared_median"] - L[, "K80_P95_vs_shared_median"])))
  if (direct_diff > 1e-10) stopf("stored direct contrast vector is not urinary-minus-K80 for %s/%s", row$country, row$unit)
  direct_est <- sum(L[, "direct_urinary_vs_K80_RRR"] * beta)
  stored_direct_est <- as.numeric(direct_obj$coefficient_contrast$by_covariance$hac14[["estimate"]])
  if (!close_enough(direct_est, stored_direct_est, 1e-8)) stopf("stored direct contrast estimate mismatch for %s/%s", row$country, row$unit)

  HinvL <- Hinv %*% L
  psi <- U %*% HinvL
  colnames(psi) <- colnames(L)
  Vlist <- list(`14` = V14, `28` = V28)
  covariance_rows <- list()
  estimate_rows <- list()
  psd_rows <- list()
  for (bw in c(14L, 28L, 90L)) {
    observed <- bartlett_hac_matrix(psi, dates, bw)
    if (bw %in% c(14L, 28L)) {
      predicted <- t(L) %*% Vlist[[as.character(bw)]] %*% L
      if (!close_enough(observed, predicted, 1e-8)) {
        stopf("projected HAC%d covariance does not reproduce full cached covariance for %s/%s", bw, row$country, row$unit)
      }
    }
    psd <- psd_audit(observed)
    psd_rows[[length(psd_rows) + 1L]] <- data.frame(country = row$country, unit = row$unit,
      scenario = row$scenario, bandwidth_days = bw, psd, stringsAsFactors = FALSE)
    for (a in seq_len(ncol(L))) for (b in seq_len(ncol(L))) {
      va <- observed[a, a]; vb <- observed[b, b]
      corr <- if (va > 0 && vb > 0) observed[a, b] / sqrt(va * vb) else NA_real_
      covariance_rows[[length(covariance_rows) + 1L]] <- data.frame(country = row$country, unit = row$unit,
        scenario = row$scenario, bandwidth_days = bw, contrast_a = colnames(L)[a], contrast_b = colnames(L)[b],
        covariance = observed[a, b], correlation = corr,
        min_eigenvalue = psd$min_eigenvalue, psd_status = psd$psd_status, stringsAsFactors = FALSE)
    }
    for (a in seq_len(ncol(L))) {
      variance <- observed[a, a]
      se <- if (variance >= 0) sqrt(variance) else NA_real_
      estimate <- sum(L[, a] * beta)
      estimate_rows[[length(estimate_rows) + 1L]] <- data.frame(country = row$country, unit = row$unit,
        scenario = row$scenario, run_id = row$run_id, cache_path = path, cache_sha256 = actual_sha,
        input_md5 = optional_character(row$input_md5), input_sha256 = optional_character(row$input_sha256), n_complete = row$n_complete,
        p95 = row$p95, shared_median = row$median, delta_p95_minus_median = row$p95 - row$median,
        contrast = colnames(L)[a], bandwidth_days = bw, estimate_logRR = estimate,
        RR = exp(estimate), variance = variance, se = se,
        lower_RR = if (is.finite(se)) exp(estimate - stats::qnorm(0.975) * se) else NA_real_,
        upper_RR = if (is.finite(se)) exp(estimate + stats::qnorm(0.975) * se) else NA_real_,
        direct_vector_identity_max_abs_difference = direct_diff, stringsAsFactors = FALSE)
    }
  }

  curve26 <- data.frame()
  curve26_mismatch <- data.frame()
  if (row$country == "BR") {
    Lcurve <- list()
    source_rows <- list()
    curve_defs <- list(K80_P95_vs_local_median = c(endpoint = "y_K80", contrast = "RR_vs_shared_median"),
                       urinary_P95_vs_local_median = c(endpoint = "y_urinary", contrast = "RR_vs_shared_median"))
    for (nm in names(curve_defs)) {
      ep <- unname(curve_defs[[nm]][["endpoint"]]); con <- unname(curve_defs[[nm]][["contrast"]])
      Lc <- project_curve_contrast(ef, beta, ep, 26.5); names(Lc) <- names(beta)
      hit <- curve_units[curve_units$unit == row$unit & curve_units$endpoint == ep & curve_units$contrast == con, , drop = FALSE]
      if (nrow(hit) != 1L) stopf("curve source row not unique for %s/%s/%s", row$country, row$unit, nm)
      estimate <- sum(Lc * beta); variance <- as.numeric(t(Lc) %*% V14 %*% Lc)
      Lcurve[[nm]] <- Lc; source_rows[[nm]] <- hit
      if (!close_enough(c(estimate, variance), c(as.numeric(hit$log_estimate), as.numeric(hit$variance)), 1e-8)) {
        curve26_mismatch <- rbind(curve26_mismatch, data.frame(country = "BR", unit = row$unit,
          contrast = nm, status = "curve_source_cache_value_mismatch",
          estimate_difference = estimate - as.numeric(hit$log_estimate), variance_difference = variance - as.numeric(hit$variance),
          stringsAsFactors = FALSE))
      }
    }
    Lcurve$direct_urinary_vs_K80 <- Lcurve$urinary_P95_vs_local_median - Lcurve$K80_P95_vs_local_median
    direct_hit <- curve_units[curve_units$unit == row$unit & curve_units$endpoint == "urinary_vs_K80" & curve_units$contrast == "direct_RRR", , drop = FALSE]
    if (nrow(direct_hit) != 1L) stopf("direct curve source row not unique for BR/%s", row$unit)
    Lcurve_mat <- do.call(cbind, Lcurve)
    Vcurve <- t(Lcurve_mat) %*% V14 %*% Lcurve_mat
    curve_direct_est <- sum(Lcurve$direct_urinary_vs_K80 * beta)
    curve_direct_var <- as.numeric(Vcurve[3L, 3L])
    if (!close_enough(c(curve_direct_est, curve_direct_var), c(as.numeric(direct_hit$log_estimate), as.numeric(direct_hit$variance)), 1e-8)) {
      curve26_mismatch <- rbind(curve26_mismatch, data.frame(country = "BR", unit = row$unit,
        contrast = "direct_urinary_vs_K80", status = "curve_source_cache_value_mismatch",
        estimate_difference = curve_direct_est - as.numeric(direct_hit$log_estimate),
        variance_difference = curve_direct_var - as.numeric(direct_hit$variance), stringsAsFactors = FALSE))
    }
    psd_curve <- psd_audit(Vcurve)
    curve26 <- data.frame(country = "BR", unit = row$unit, temperature = 26.5,
      shared_median = row$median, absolute_distance_from_local_median = abs(26.5 - row$median),
      locally_supported = as.logical(source_rows[[1L]]$locally_supported),
      K80_logRR = sum(Lcurve$K80_P95_vs_local_median * beta),
      urinary_logRR = sum(Lcurve$urinary_P95_vs_local_median * beta),
      direct_logRR = curve_direct_est, direct_RRR = exp(curve_direct_est),
      K80_variance = as.numeric(Vcurve[1L, 1L]), urinary_variance = as.numeric(Vcurve[2L, 2L]),
      direct_variance = curve_direct_var, covariance_K80_urinary = as.numeric(Vcurve[1L, 2L]),
      covariance_K80_direct = as.numeric(Vcurve[1L, 3L]), covariance_urinary_direct = as.numeric(Vcurve[2L, 3L]),
      psd_curve, stringsAsFactors = FALSE)
  }

  list(row = row, dates = dates, psi = psi, estimates = do.call(rbind, estimate_rows),
       covariance = do.call(rbind, covariance_rows), psd = do.call(rbind, psd_rows),
       curve26 = curve26, curve26_mismatch = curve26_mismatch,
       cache_check = data.frame(country = row$country, unit = row$unit, scenario = row$scenario,
         run_id = row$run_id, cache_path = path, cache_sha256 = actual_sha,
         cache_hash_verified_against_manifest = !is.null(expected_sha) && !is.na(expected_sha) && nzchar(expected_sha),
         manifest_input_md5 = optional_character(row$input_md5), manifest_input_sha256 = optional_character(row$input_sha256),
         n_complete = row$n_complete, p95 = row$p95, shared_median = row$median,
         direct_vector_identity_max_abs_difference = direct_diff,
         direct_estimate = direct_est, cached_direct_estimate = stored_direct_est,
         direct_estimate_abs_difference = abs(direct_est - stored_direct_est), stringsAsFactors = FALSE))
}

dl_fit <- function(yi, vi) {
  if (length(yi) < 2L || length(yi) != length(vi) || any(!is.finite(yi)) || any(!is.finite(vi)) || any(vi <= 0)) {
    stop("DL requires all finite estimates and strictly positive variances", call. = FALSE)
  }
  wf <- 1 / vi
  mu <- sum(wf * yi) / sum(wf)
  Q <- sum(wf * (yi - mu)^2)
  C <- sum(wf) - sum(wf^2) / sum(wf)
  if (!is.finite(C) || C <= 0) stop("DL denominator is non-positive", call. = FALSE)
  tau2 <- max(0, (Q - (length(yi) - 1L)) / C)
  wr <- 1 / (vi + tau2)
  a <- wr / sum(wr)
  estimate <- sum(a * yi)
  variance <- 1 / sum(wr)
  if (!is.finite(variance) || variance <= 0) stop("DL variance is invalid", call. = FALSE)
  if (!close_enough(sum(a^2 * (vi + tau2)), variance, 1e-10)) stop("DL variance identity failed", call. = FALSE)
  list(estimate = estimate, variance = variance, se = sqrt(variance), tau2 = tau2,
       weights = a, sum_weight_sq = sum(a^2), Q = Q, C = C)
}

weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  if (!length(x)) return(NA_real_)
  x[which(cumsum(w) >= p * sum(w))[1L]]
}

target_summary <- function(targets, weights = NULL, weighting, contrast) {
  vars <- c("p95", "shared_median", "delta_p95_minus_median")
  z <- lapply(vars, function(v) {
    x <- targets[[v]]
    if (is.null(weights)) weights <- rep(1, length(x))
    data.frame(weighting = weighting, contrast = contrast, variable = v, n = sum(is.finite(x)),
      min = min(x, na.rm = TRUE), q05 = weighted_quantile(x, weights, .05),
      q25 = weighted_quantile(x, weights, .25), median = weighted_quantile(x, weights, .5),
      q75 = weighted_quantile(x, weights, .75), q95 = weighted_quantile(x, weights, .95),
      max = max(x, na.rm = TRUE), weighted_mean = sum(x * weights, na.rm = TRUE) / sum(weights[is.finite(x)], na.rm = TRUE),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, z)
}

project_curve_contrast <- function(ef, beta, endpoint, temperature) {
  Xhi <- .design_at(ef$fit, ef$data, endpoint = endpoint, basis = ef$basis, temp_value = temperature)
  Xmed <- .design_at(ef$fit, ef$data, endpoint = endpoint, basis = ef$basis, temp_value = ef$median)
  as.numeric(Xhi - Xmed)
}

make_attempt_audit <- function(path, country, scenario, sample_rows) {
  a <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = "character", na.strings = c("", "NA"))
  if (!all(c("country", "unit", "scenario") %in% names(a))) stopf("fit audit schema invalid: %s", path)
  a$country <- toupper(a$country)
  a <- a[a$country == country & a$scenario == scenario, , drop = FALSE]
  key <- paste(sample_rows$country, sample_rows$unit, sample_rows$scenario, sep = "\r")
  akey <- paste(a$country, a$unit, a$scenario, sep = "\r")
  cache_key <- paste(sample_rows$country, sample_rows$unit, sample_rows$scenario, sep = "\r")
  a$cache_manifest_match <- akey %in% cache_key
  a$fit_attempt_status <- ifelse(tolower(a$converged %||% "") == "true", "converged",
    ifelse(tolower(a$converged %||% "") == "false", "failed_or_not_converged", "status_unreported"))
  a$source_audit_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  a$source_audit_sha256 <- sha256_file(path)
  a
}

br_manifest_path <- file.path(base_dir, "results", "fits", "primary_cache_manifest.csv")
mx_manifest_path <- file.path(base_dir, "results", "mx_updated", "fits", "fit_cache_manifest.csv")
br_scenario <- "primary_threeoutcome_E3_separate"
mx_scenario <- "primary_complete_source_2015_2022"
br <- read_selected_manifest(br_manifest_path, "BR", br_scenario)
mx <- read_selected_manifest(mx_manifest_path, "MX", mx_scenario)
manifest_cols <- union(names(br), names(mx))
for (nm in setdiff(manifest_cols, names(br))) br[[nm]] <- NA_character_
for (nm in setdiff(manifest_cols, names(mx))) mx[[nm]] <- NA_character_
br <- br[, manifest_cols, drop = FALSE]
mx <- mx[, manifest_cols, drop = FALSE]
selected_rows <- rbind(br, mx)
if (anyDuplicated(paste(selected_rows$country, selected_rows$unit, selected_rows$scenario, sep = "\r"))) {
  stop("selected sample contains duplicate country/unit/scenario", call. = FALSE)
}

expected_counts <- c(BR = 450L, MX = 32L)
curve_unit_path <- file.path(base_dir, "results", "downstream_curves_burden_20261112T210050", "curves_unit.csv")
if (!file.exists(curve_unit_path)) stopf("curve unit source is missing: %s", curve_unit_path)
curve_units <- utils::read.csv(curve_unit_path, stringsAsFactors = FALSE,
                               colClasses = c(unit = "character"), check.names = FALSE)
curve_units <- curve_units[curve_units$country == "BR" & abs(as.numeric(curve_units$temperature) - 26.5) < 1e-10, , drop = FALSE]
if (!setequal(unique(curve_units$unit), br$unit) || nrow(curve_units) != length(br$unit) * 3L) {
  stop("26.5 C curve source does not exactly match the 450 BR primary units and three contrast rows per unit", call. = FALSE)
}
projection_signature <- "fixed-reference-joint-three-contrast-actual-calendar-HAC14-28-90-v1-20261113"
checkpoint_path <- file.path(out_dir, "marginals_projection_checkpoint.rds")
checkpoint_signature <- list(projection_signature = projection_signature,
  engine_sha256 = sha256_file(engine_path),
  br_manifest_sha256 = sha256_file(br_manifest_path), mx_manifest_sha256 = sha256_file(mx_manifest_path),
  curve_unit_sha256 = sha256_file(curve_unit_path),
  selected_unit_keys = paste(selected_rows$country, selected_rows$unit, selected_rows$scenario, sep = "\r"))
cache_results <- list()
failure_rows <- list()
processed_keys <- character()
if (file.exists(checkpoint_path)) {
  checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(e) e)
  if (!inherits(checkpoint, "error") && identical(checkpoint$signature, checkpoint_signature)) {
    cache_results <- checkpoint$cache_results %||% list()
    failure_rows <- checkpoint$failure_rows %||% list()
    processed_keys <- checkpoint$processed_keys %||% names(cache_results)
    # Retry prior cache failures after a software fix; preserve successful projections.
    retry_keys <- vapply(failure_rows, function(z) paste(z$country, z$unit, z$scenario, sep = "\r"), character(1))
    if (length(retry_keys)) {
      prior_failure_path <- file.path(out_dir, "projection_software_retry_audit.csv")
      prior_failure_table <- do.call(rbind, failure_rows)
      if (file.exists(prior_failure_path)) {
        prior_failure_table <- rbind(utils::read.csv(prior_failure_path, stringsAsFactors = FALSE), prior_failure_table)
      }
      atomic_write_csv(prior_failure_table, prior_failure_path)
      failure_rows <- list()
      processed_keys <- setdiff(processed_keys, retry_keys)
    }
    message(sprintf("resuming %d compact projections from checkpoint", length(cache_results)))
  } else message("ignoring projection checkpoint with incompatible signature or unreadable contents")
}
save_projection_checkpoint <- function() {
  checkpoint <- list(signature = checkpoint_signature, processed_keys = processed_keys,
    cache_results = cache_results, failure_rows = failure_rows,
    completed_count = length(processed_keys), validated_count = length(cache_results),
    failure_count = length(failure_rows), saved_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  tmp <- tempfile(pattern = ".marginals_projection_checkpoint.", tmpdir = out_dir, fileext = ".tmp")
  saveRDS(checkpoint, tmp, compress = "gzip")
  if (file.exists(checkpoint_path)) unlink(checkpoint_path)
  if (!file.rename(tmp, checkpoint_path)) stopf("failed to atomically write projection checkpoint")
  invisible(checkpoint_path)
}
processed_since_checkpoint <- 0L
consecutive_error_count <- 0L
last_error_reason <- NULL
for (i in seq_len(nrow(selected_rows))) {
  row <- selected_rows[i, , drop = FALSE]
  key <- paste(row$country, row$unit, row$scenario, sep = "\r")
  if (key %in% processed_keys) next
  result <- tryCatch(load_unit(row), error = function(e) e)
  if (inherits(result, "error")) {
    reason <- conditionMessage(result)
    failure_rows[[length(failure_rows) + 1L]] <- data.frame(country = row$country, unit = row$unit,
      scenario = row$scenario, cache_path = row$cache_path, status = "cache_or_contrast_validation_failed",
      reason = reason, stringsAsFactors = FALSE)
    consecutive_error_count <- if (identical(last_error_reason, reason)) consecutive_error_count + 1L else 1L
    last_error_reason <- reason
  } else {
    cache_results[[key]] <- result
    consecutive_error_count <- 0L
    last_error_reason <- NULL
  }
  processed_keys <- c(processed_keys, key)
  processed_since_checkpoint <- processed_since_checkpoint + 1L
  if (consecutive_error_count >= 2L) {
    save_projection_checkpoint()
    stopf("two consecutive cache projection failures share one reason; checkpoint preserved: %s", last_error_reason)
  }
  if (processed_since_checkpoint >= 25L || i == nrow(selected_rows)) {
    save_projection_checkpoint()
    processed_since_checkpoint <- 0L
  }
  if (i %% 25L == 0L || i == nrow(selected_rows)) message(sprintf("validated %d/%d full joint caches", i, nrow(selected_rows)))
  if (i %% 10L == 0L) invisible(gc(FALSE))
}

cache_checks <- if (length(cache_results)) do.call(rbind, lapply(cache_results, `[[`, "cache_check")) else data.frame()
unit_estimates <- if (length(cache_results)) do.call(rbind, lapply(cache_results, `[[`, "estimates")) else data.frame()
unit_covariance <- if (length(cache_results)) do.call(rbind, lapply(cache_results, `[[`, "covariance")) else data.frame()
unit_psd <- if (length(cache_results)) do.call(rbind, lapply(cache_results, `[[`, "psd")) else data.frame()
if (length(failure_rows)) failure_audit <- do.call(rbind, failure_rows) else failure_audit <- data.frame()
curve26_validation <- if (length(cache_results)) do.call(rbind, lapply(cache_results, `[[`, "curve26")) else data.frame()
curve26_failures <- if (length(cache_results)) do.call(rbind, lapply(cache_results, `[[`, "curve26_mismatch")) else data.frame()
required_curve26_columns <- c("country", "unit", "temperature", "shared_median", "direct_RRR", "direct_variance", "min_eigenvalue", "psd_status")
if (nrow(curve26_validation) > 0L && !all(required_curve26_columns %in% names(curve26_validation))) {
  stopf("26.5 C curve validation table has an incomplete schema: %s", paste(setdiff(required_curve26_columns, names(curve26_validation)), collapse = ", "))
}

targets <- selected_rows[, c("country", "unit", "scenario", "run_id", "cache_path", "p95", "median", "n_complete"), drop = FALSE]
targets$shared_median <- targets$median
targets$delta_p95_minus_median <- targets$p95 - targets$shared_median
targets$cache_sha256 <- if (nrow(cache_checks)) cache_checks$cache_sha256[match(paste(targets$country, targets$unit, targets$scenario), paste(cache_checks$country, cache_checks$unit, cache_checks$scenario))] else NA_character_

attempts_br <- make_attempt_audit(file.path(base_dir, "results", "core", "primary_threeoutcome_E3_separate_BR_fit_audit.csv"), "BR", br_scenario, selected_rows)
attempts_mx <- make_attempt_audit(file.path(base_dir, "results", "mx_updated", "primary_complete_source_2015_2022_MX_fit_audit.csv"), "MX", mx_scenario, selected_rows)
fit_attempts <- rbind(attempts_br, attempts_mx)
if (nrow(fit_attempts)) {
  fit_attempts$cache_validation_status <- ifelse(fit_attempts$cache_manifest_match, "successful_joint_cache_present", "no_successful_joint_cache")
}

country_summaries <- list()
sync_rows <- list()
pooled_cov_rows <- list()
weight_rows <- list()
target_summaries <- list()
pooling_status <- list()
rrr_validation_rows <- list()
contrast_names <- c("K80_P95_vs_shared_median", "urinary_P95_vs_shared_median", "direct_urinary_vs_K80_RRR")
expected_rrr <- c(BR = 1.09806036145693, MX = 1.08339060731185)
saved_summary_paths <- c(BR = file.path(base_dir, "results", "synchronous_meta", "synchronous_meta_summary.csv"),
                         MX = file.path(base_dir, "results", "synchronous_meta", "mx_updated_primary", "synchronous_meta_summary.csv"))
for (country in c("BR", "MX")) {
  sc <- if (country == "BR") br_scenario else mx_scenario
  n_loaded <- sum(vapply(cache_results, function(z) z$row$country == country, logical(1)))
  expected_n <- unname(expected_counts[[country]])
  attempts <- fit_attempts[fit_attempts$country == country, , drop = FALSE]
  attempt_failed <- sum(attempts$fit_attempt_status != "converged")
  valid_country <- n_loaded == expected_n && attempt_failed == 0L
  if (!valid_country) {
    pooling_status[[country]] <- data.frame(country = country, scenario = sc, expected_joint_units = expected_n,
      loaded_joint_units = n_loaded, fit_attempts = nrow(attempts), failed_or_unreported_fit_attempts = attempt_failed,
      pooling_status = "not_pooled_incomplete_or_failed_sample", stringsAsFactors = FALSE)
    next
  }
  per_country <- unit_estimates[unit_estimates$country == country & unit_estimates$bandwidth_days == 14L, , drop = FALSE]
  targets_country <- targets[targets$country == country, , drop = FALSE]
  by_contrast <- vector("list", length(contrast_names)); names(by_contrast) <- contrast_names
  for (contrast in contrast_names) {
    u <- per_country[per_country$contrast == contrast, , drop = FALSE]
    u <- u[match(targets_country$unit, u$unit), , drop = FALSE]
    if (nrow(u) != expected_n || anyNA(u$variance)) stopf("incomplete unit contrast vectors for %s/%s", country, contrast)
    yi <- u$estimate_logRR; vi <- u$variance
    dl <- tryCatch(dl_fit(yi, vi), error = function(e) e)
    if (inherits(dl, "error")) {
      pooling_status[[country]] <- data.frame(country = country, scenario = sc, expected_joint_units = expected_n,
        loaded_joint_units = n_loaded, fit_attempts = nrow(attempts), failed_or_unreported_fit_attempts = attempt_failed,
        pooling_status = paste0("DL_failed:", conditionMessage(dl)), stringsAsFactors = FALSE)
      break
    }
    u$dl_weight_hac14 <- dl$weights
    by_contrast[[contrast]] <- list(unit = u, dl = dl)
    endpoint_name <- switch(contrast,
      K80_P95_vs_shared_median = "y_K80",
      urinary_P95_vs_shared_median = "y_urinary",
      direct_urinary_vs_K80_RRR = "urinary_vs_K80")
    effect_contrast <- if (contrast == "direct_urinary_vs_K80_RRR") "direct_RRR" else "P95_vs_shared_median"
    I2_dl <- if (dl$Q > 0) max(0, (dl$Q - (length(yi) - 1L)) / dl$Q) * 100 else 0
    concordance <- mean(yi > 0)
    weight_rows[[length(weight_rows) + 1L]] <- data.frame(country = country, scenario = sc, contrast = contrast,
      unit = u$unit, estimate_logRR = yi, variance_hac14 = vi, tau2_DL = dl$tau2,
      raw_weight_DL = 1 / (vi + dl$tau2), normalized_weight_DL = dl$weights,
      normalized_weight_sq = dl$weights^2, stringsAsFactors = FALSE)
    reml <- tryCatch(metafor::rma.uni(yi = yi, vi = vi, method = "REML", test = "knha"), error = function(e) e)
    zcrit <- stats::qnorm(.975)
    country_summaries[[length(country_summaries) + 1L]] <- data.frame(country = country, scenario = sc,
      endpoint = endpoint_name, contrast = effect_contrast, method = "DL", covariance = "unit_HAC14",
      estimator = "DL", uncertainty = "conventional_independent_unit_DL", inference = "conventional_independent_unit_DL",
      k = length(yi), estimate_logRR = dl$estimate, RR = exp(dl$estimate), se = dl$se,
      lower_logRR = dl$estimate - zcrit * dl$se, upper_logRR = dl$estimate + zcrit * dl$se,
      lower_RR = exp(dl$estimate - zcrit * dl$se), upper_RR = exp(dl$estimate + zcrit * dl$se),
      lo = exp(dl$estimate - zcrit * dl$se), hi = exp(dl$estimate + zcrit * dl$se),
      I2 = I2_dl, concordance = concordance, tau2 = dl$tau2, sum_normalized_weight_sq = dl$sum_weight_sq,
      maximum_normalized_weight = max(dl$weights), effective_unit_count = 1 / sum(dl$weights^2),
      failure = "", stringsAsFactors = FALSE)
    if (inherits(reml, "error")) {
      country_summaries[[length(country_summaries) + 1L]] <- data.frame(country = country, scenario = sc,
        endpoint = endpoint_name, contrast = effect_contrast, method = "REML_HK", covariance = "unit_HAC14",
        estimator = "REML_HK", uncertainty = "ordinary_independent_unit_REML_HK", inference = "ordinary_independent_unit_REML_HK",
        k = length(yi), estimate_logRR = NA_real_, RR = NA_real_, se = NA_real_, lower_logRR = NA_real_, upper_logRR = NA_real_,
        lower_RR = NA_real_, upper_RR = NA_real_, lo = NA_real_, hi = NA_real_, I2 = NA_real_, concordance = concordance,
        tau2 = NA_real_, sum_normalized_weight_sq = NA_real_,
        maximum_normalized_weight = NA_real_, effective_unit_count = NA_real_, failure = conditionMessage(reml), stringsAsFactors = FALSE)
    } else {
      reml_w <- 1 / (vi + as.numeric(reml$tau2))
      reml_w <- reml_w / sum(reml_w)
      country_summaries[[length(country_summaries) + 1L]] <- data.frame(country = country, scenario = sc,
        endpoint = endpoint_name, contrast = effect_contrast, method = "REML_HK", covariance = "unit_HAC14",
        estimator = "REML_HK", uncertainty = "ordinary_independent_unit_REML_HK", inference = "ordinary_independent_unit_REML_HK",
        k = length(yi), estimate_logRR = as.numeric(reml$b[1L]), RR = exp(as.numeric(reml$b[1L])), se = as.numeric(reml$se),
        lower_logRR = as.numeric(reml$ci.lb), upper_logRR = as.numeric(reml$ci.ub),
        lower_RR = exp(as.numeric(reml$ci.lb)), upper_RR = exp(as.numeric(reml$ci.ub)),
        lo = exp(as.numeric(reml$ci.lb)), hi = exp(as.numeric(reml$ci.ub)), I2 = as.numeric(reml$I2), concordance = concordance,
        tau2 = as.numeric(reml$tau2), sum_normalized_weight_sq = sum(reml_w^2),
        maximum_normalized_weight = max(reml_w), effective_unit_count = 1 / sum(reml_w^2), failure = "", stringsAsFactors = FALSE)
    }
  }
  if (!is.null(pooling_status[[country]])) next
  start_date <- min(vapply(cache_results, function(z) if (z$row$country == country) as.numeric(min(z$dates)) else Inf, numeric(1)), na.rm = TRUE)
  end_date <- max(vapply(cache_results, function(z) if (z$row$country == country) as.numeric(max(z$dates)) else -Inf, numeric(1)), na.rm = TRUE)
  calendar <- seq(as.Date(start_date, origin = "1970-01-01"), as.Date(end_date, origin = "1970-01-01"), by = "day")
  country_psi <- matrix(0, nrow = length(calendar), ncol = length(contrast_names),
                        dimnames = list(as.character(calendar), contrast_names))
  for (contrast in contrast_names) {
    w <- by_contrast[[contrast]]$unit[, c("unit", "dl_weight_hac14"), drop = FALSE]
    for (i in seq_len(nrow(w))) {
      key <- paste(country, w$unit[i], sc, sep = "\r")
      z <- cache_results[[key]]
      idx <- match(z$dates, calendar)
      j <- match(contrast, contrast_names)
      country_psi[idx, j] <- country_psi[idx, j] + z$psi[, contrast] * w$dl_weight_hac14[i]
    }
  }
  for (bw in c(14L, 28L, 90L)) {
    Vs <- bartlett_hac_matrix(country_psi, calendar, bw)
    psd <- psd_audit(Vs)
    for (a in seq_along(contrast_names)) for (b in seq_along(contrast_names)) {
      pooled_cov_rows[[length(pooled_cov_rows) + 1L]] <- data.frame(country = country, scenario = sc,
        bandwidth_days = bw, contrast_a = contrast_names[a], contrast_b = contrast_names[b],
        synchronized_covariance_conditional_on_DL_weights_tau = Vs[a, b],
        synchronized_correlation_conditional_on_DL_weights_tau = if (Vs[a, a] > 0 && Vs[b, b] > 0) Vs[a, b] / sqrt(Vs[a, a] * Vs[b, b]) else NA_real_,
        min_eigenvalue = psd$min_eigenvalue, psd_status = psd$psd_status,
        tau_component_for_diagonal_only = if (a == b) by_contrast[[contrast_names[a]]]$dl$tau2 * by_contrast[[contrast_names[a]]]$dl$sum_weight_sq else NA_real_,
        stringsAsFactors = FALSE)
    }
    for (contrast in contrast_names) {
      j <- match(contrast, contrast_names)
      dl <- by_contrast[[contrast]]$dl
      endpoint_name <- switch(contrast,
        K80_P95_vs_shared_median = "y_K80",
        urinary_P95_vs_shared_median = "y_urinary",
        direct_urinary_vs_K80_RRR = "urinary_vs_K80")
      effect_contrast <- if (contrast == "direct_urinary_vs_K80_RRR") "direct_RRR" else "P95_vs_shared_median"
      yi_current <- by_contrast[[contrast]]$unit$estimate_logRR
      I2_current <- if (dl$Q > 0) max(0, (dl$Q - (length(yi_current) - 1L)) / dl$Q) * 100 else 0
      concordance_current <- mean(yi_current > 0)
      Vsync <- Vs[j, j]
      Vtau <- dl$tau2 * dl$sum_weight_sq
      if (!is.finite(Vsync) || Vsync <= 0 || !is.finite(Vtau) || Vtau < 0) stopf("invalid synchronized variance for %s/%s/%d", country, contrast, bw)
      for (add_tau in c(FALSE, TRUE)) {
        variance <- Vsync + if (add_tau) Vtau else 0
        se <- sqrt(variance)
        point <- dl$estimate
        sync_rows[[length(sync_rows) + 1L]] <- data.frame(country = country, scenario = sc,
          endpoint = endpoint_name, contrast = effect_contrast, method = "DL",
          covariance = paste0("calendar_synchronous_HAC", bw),
          base_DL_weight_variance = "unit_HAC14", bandwidth_days = bw,
          inference = if (add_tau) "calendar_synchronous_plus_DL_tau_component" else "calendar_synchronous_conditional_on_DL_weights_tau",
          reporting_role = if (bw == 14L && add_tau) "primary" else "sensitivity",
          k = expected_n, estimate_logRR = point, RR = exp(point), tau2_DL = dl$tau2,
          sum_normalized_weight_sq = dl$sum_weight_sq, V_sync = Vsync, V_tau_component = if (add_tau) Vtau else 0,
          variance = variance, se = se, lower_logRR = point - stats::qnorm(.975) * se,
          upper_logRR = point + stats::qnorm(.975) * se,
          lower_RR = exp(point - stats::qnorm(.975) * se), upper_RR = exp(point + stats::qnorm(.975) * se),
          lo = exp(point - stats::qnorm(.975) * se), hi = exp(point + stats::qnorm(.975) * se),
          I2 = I2_current, concordance = concordance_current,
          I2_source = "DL_between_unit_heterogeneity_from_HAC14_unit_variances; not re-estimated for synchronized covariance",
          maximum_normalized_weight = max(dl$weights), effective_unit_count = 1 / dl$sum_weight_sq,
          uncertainty_limit = "weights and tau conditioned; synchronized Bartlett HAC captures same-calendar score dependence, not unrestricted spatial dependence",
          failure = "", stringsAsFactors = FALSE)
      }
    }
  }
  if (country == "BR") {
    for (contrast in contrast_names) {
      sumrow <- target_summary(targets_country,
        weights = if (contrast == "all_units_equal") NULL else by_contrast[[contrast]]$dl$weights,
        weighting = if (contrast == "all_units_equal") "equal_unit" else "DL_HAC14_inverse_variance",
        contrast = contrast)
      sumrow$country <- country; sumrow$scenario <- sc
      target_summaries[[length(target_summaries) + 1L]] <- sumrow
    }
    sumrow <- target_summary(targets_country, weights = NULL, weighting = "equal_unit", contrast = "all_units_equal")
    sumrow$country <- country; sumrow$scenario <- sc
    target_summaries[[length(target_summaries) + 1L]] <- sumrow
  } else {
    for (contrast in contrast_names) {
      sumrow <- target_summary(targets_country,
        weights = by_contrast[[contrast]]$dl$weights,
        weighting = "DL_HAC14_inverse_variance", contrast = contrast)
      sumrow$country <- country; sumrow$scenario <- sc
      target_summaries[[length(target_summaries) + 1L]] <- sumrow
    }
    sumrow <- target_summary(targets_country, weights = NULL, weighting = "equal_unit", contrast = "all_units_equal")
    sumrow$country <- country; sumrow$scenario <- sc
    target_summaries[[length(target_summaries) + 1L]] <- sumrow
  }
  if (nrow(fit_attempts[fit_attempts$country == country & fit_attempts$fit_attempt_status != "converged", , drop = FALSE]) == 0L) {
    saved <- saved_summary_paths[[country]]
    saved_rr <- NA_real_; saved_status <- "saved_primary_RRR_not_found"
    if (file.exists(saved)) {
      old <- utils::read.csv(saved, stringsAsFactors = FALSE)
      keep <- old$contrast == "y_urinary_vs_y_K80" & old$base_variance_for_DL_weights == "unit_HAC14" & old$inference == "calendar_synchronous_conditional_on_DL_weights_tau"
      if ("scenario" %in% names(old)) keep <- keep & old$scenario == sc
      old <- old[keep, , drop = FALSE]
      if (nrow(old)) {
        saved_rr <- as.numeric(old$RR[1L])
        our_rr <- by_contrast[["direct_urinary_vs_K80_RRR"]]$dl$estimate |> exp()
        match_expected <- close_enough(c(our_rr, saved_rr), c(expected_rrr[[country]], expected_rrr[[country]]), 1e-12)
        saved_status <- if (match_expected) "matched_saved_primary_RRR_and_expected_value" else "saved_primary_RRR_mismatch"
        checkrow <- data.frame(country = country, scenario = sc, contrast = "direct_urinary_vs_K80_RRR",
          saved_summary_path = normalizePath(saved, winslash = "/", mustWork = TRUE), saved_primary_RRR = saved_rr,
          expected_primary_RRR = expected_rrr[[country]], extracted_primary_RRR = our_rr,
          absolute_difference_from_saved = abs(our_rr - saved_rr),
          absolute_difference_from_expected = abs(our_rr - expected_rrr[[country]]),
          tolerance = 1e-12, validation_status = saved_status, stringsAsFactors = FALSE)
        rrr_validation_rows[[length(rrr_validation_rows) + 1L]] <- checkrow
      }
    }
    pooling_status[[country]] <- data.frame(country = country, scenario = sc, expected_joint_units = expected_n,
      loaded_joint_units = n_loaded, fit_attempts = nrow(attempts), failed_or_unreported_fit_attempts = attempt_failed,
      direct_RRR_validation_status = saved_status, pooling_status = if (saved_status == "matched_saved_primary_RRR_and_expected_value") "accepted" else "rejected_direct_RRR_validation",
      stringsAsFactors = FALSE)
  }
}

if (length(country_summaries)) atomic_write_csv(do.call(rbind, country_summaries), file.path(out_dir, "marginals_country_summary.csv"))
if (length(sync_rows)) atomic_write_csv(do.call(rbind, sync_rows), file.path(out_dir, "marginals_synchronized_uncertainty.csv"))
if (length(pooled_cov_rows)) atomic_write_csv(do.call(rbind, pooled_cov_rows), file.path(out_dir, "marginals_pooled_synchronized_covariance.csv"))
if (length(weight_rows)) atomic_write_csv(do.call(rbind, weight_rows), file.path(out_dir, "marginals_unit_DL_weights.csv"))
if (length(target_summaries)) atomic_write_csv(do.call(rbind, target_summaries), file.path(out_dir, "unit_temperature_target_distributions.csv"))
if (nrow(cache_checks)) atomic_write_csv(cache_checks, file.path(out_dir, "successful_joint_samples.csv"))
if (nrow(unit_estimates)) atomic_write_csv(unit_estimates, file.path(out_dir, "marginals_unit_estimates.csv"))
if (nrow(unit_covariance)) atomic_write_csv(unit_covariance, file.path(out_dir, "marginals_unit_full_contrast_covariance.csv"))
if (nrow(unit_psd)) atomic_write_csv(unit_psd, file.path(out_dir, "marginals_unit_covariance_psd_audit.csv"))
if (nrow(fit_attempts)) atomic_write_csv(fit_attempts, file.path(out_dir, "fit_attempts_preserved.csv"))
if (nrow(failure_audit)) atomic_write_csv(failure_audit, file.path(out_dir, "cache_or_contrast_failure_audit.csv"))
if (nrow(curve26_validation) > 0L) atomic_write_csv(curve26_validation, file.path(out_dir, "curve_26_5_cache_covariance_validation.csv"))
if (nrow(curve26_failures) > 0L) atomic_write_csv(curve26_failures, file.path(out_dir, "curve_26_5_cache_mismatch_audit.csv"))
if (length(pooling_status)) atomic_write_csv(do.call(rbind, pooling_status), file.path(out_dir, "country_sample_acceptance.csv"))
if (length(rrr_validation_rows)) atomic_write_csv(do.call(rbind, rrr_validation_rows), file.path(out_dir, "saved_primary_RRR_validation.csv"))

manifest_files <- c(br_manifest_path, mx_manifest_path,
  file.path(base_dir, "results", "core", "primary_threeoutcome_E3_separate_BR_fit_audit.csv"),
  file.path(base_dir, "results", "mx_updated", "primary_complete_source_2015_2022_MX_fit_audit.csv"),
  curve_unit_path, saved_summary_paths)
manifest_files <- unique(manifest_files[file.exists(manifest_files)])
input_sources <- data.frame(path = normalizePath(manifest_files, winslash = "/", mustWork = TRUE),
  md5 = vapply(manifest_files, md5_file, character(1)), sha256 = vapply(manifest_files, sha256_file, character(1)),
  stringsAsFactors = FALSE)
atomic_write_csv(input_sources, file.path(out_dir, "source_file_hashes.csv"))

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), winslash = "/", mustWork = TRUE)
psd_failure <- nrow(unit_psd) > 0L && any(unit_psd$psd_status == "not_PSD")
curve_psd_failure <- nrow(curve26_validation) > 0L && any(curve26_validation$psd_status == "not_PSD")
overall_status <- if (nrow(failure_audit) == 0L && !psd_failure && !curve_psd_failure &&
  nrow(curve26_failures) == 0L && nrow(curve26_validation) == expected_counts[["BR"]] &&
  length(pooling_status) == 2L && all(vapply(pooling_status, function(x) "accepted" %in% x$pooling_status, logical(1)))) "accepted" else "rejected_or_incomplete"
run_manifest <- list(run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  script_path = script_path, script_sha256 = sha256_file(script_path),
  engine_path = normalizePath(engine_path, winslash = "/", mustWork = TRUE), engine_sha256 = sha256_file(engine_path),
  projection_signature = projection_signature,
  projection_checkpoint_path = normalizePath(checkpoint_path, winslash = "/", mustWork = TRUE),
  projection_checkpoint_sha256 = sha256_file(checkpoint_path),
  prior_software_retry_audit_path = if (file.exists(file.path(out_dir, "projection_software_retry_audit.csv"))) normalizePath(file.path(out_dir, "projection_software_retry_audit.csv"), winslash = "/", mustWork = TRUE) else NULL,
  prior_software_retry_audit_sha256 = if (file.exists(file.path(out_dir, "projection_software_retry_audit.csv"))) sha256_file(file.path(out_dir, "projection_software_retry_audit.csv")) else NULL,
  project_base = base_dir, output_dir = out_dir, primary_samples = list(BR = br_scenario, MX = mx_scenario),
  expected_successful_units = expected_counts, observed_validated_units = c(BR = sum(vapply(cache_results, function(z) z$row$country == "BR", logical(1))),
                                                                             MX = sum(vapply(cache_results, function(z) z$row$country == "MX", logical(1)))),
  saved_primary_RRR_validation_targets = expected_rrr, status = overall_status,
  inference = "DL point and unit weights use per-unit HAC14 contrast variances; synchronized Bartlett HAC14/28/90 uses those fixed weights; DL tau augmentation adds tau2 * sum(weights^2) on each marginal variance. Off-diagonal between-unit tau covariance is not estimated.",
  covariance = "Full 3x3 covariance among K80 P95/shared-median, urinary P95/shared-median, and direct urinary-vs-K80 contrasts projected from cached full joint coefficient covariance (HAC14/28) or score-based actual-calendar Bartlett covariance (HAC90). Direct RRR is the saved joint contrast, checked equal to urinary-minus-K80 contrast vectors; it is never computed by dividing pooled marginal RRs.",
  variance_policy = "No variance floor or clipping applied; negative/non-finite variances or failed PSD checks are preserved in audit outputs and prevent acceptance.",
  countries = if (length(pooling_status)) pooling_status else list(), R_version = R.version.string)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(run_manifest, file.path(out_dir, "run_manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, na = "null", null = "null", digits = NA)
} else stop("jsonlite is required for the provenance manifest", call. = FALSE)
message(sprintf("Fixed-reference joint marginal extraction status: %s; outputs: %s", overall_status, out_dir))
if (overall_status != "accepted") quit(status = 2L)
