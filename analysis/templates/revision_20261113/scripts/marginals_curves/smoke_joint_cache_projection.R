#!/usr/bin/env Rscript

# Fresh-process smoke of one BR and one complete-source MX joint cache,
# including full HAC14/HAC28 projection and the BR 26.5 C curve helper.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args)) args[[1L]] else "/data/admin/fengshuo/latam_stone_repair_20261112"
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
if (!requireNamespace("splines", quietly = TRUE)) stop("splines namespace unavailable", call. = FALSE)
source(file.path(base_dir, "scripts", "repair_engine.R"), local = FALSE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
close_enough <- function(x, y, tol = 1e-8) length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
  all(abs(x - y) <= tol * pmax(1, abs(x), abs(y)))
sha256_file <- function(path) {
  z <- suppressWarnings(system2(Sys.which("sha256sum"), shQuote(path), stdout = TRUE, stderr = TRUE))
  if (!length(z)) stopf("cannot checksum %s", path)
  tolower(sub("\\s+.*$", "", z[[1L]]))
}
hac_matrix <- function(psi, dates, lag) {
  dates <- as.Date(dates); ord <- order(dates); dates <- dates[ord]; psi <- psi[ord, , drop = FALSE]
  cal <- seq(min(dates), max(dates), by = "day")
  full <- matrix(0, length(cal), ncol(psi))
  full[as.integer(dates - min(dates)) + 1L, ] <- psi
  V <- crossprod(full)
  for (ell in seq_len(min(lag, nrow(full) - 1L))) {
    cross <- crossprod(full[(ell + 1L):nrow(full), , drop = FALSE], full[seq_len(nrow(full) - ell), , drop = FALSE])
    V <- V + (1 - ell / (lag + 1)) * (cross + t(cross))
  }
  (V + t(V)) / 2
}

br_manifest <- utils::read.csv(file.path(base_dir, "results", "fits", "primary_cache_manifest.csv"),
  stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
mx_manifest <- utils::read.csv(file.path(base_dir, "results", "mx_updated", "fits", "fit_cache_manifest.csv"),
  stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
sample_cols <- union(names(br_manifest), names(mx_manifest))
br_one <- br_manifest[br_manifest$country == "BR" & br_manifest$unit == "11001" & br_manifest$scenario == "primary_threeoutcome_E3_separate", , drop = FALSE]
mx_one <- mx_manifest[mx_manifest$country == "MX" & mx_manifest$unit == "01" & mx_manifest$scenario == "primary_complete_source_2015_2022", , drop = FALSE]
for (nm in setdiff(sample_cols, names(br_one))) br_one[[nm]] <- NA_character_
for (nm in setdiff(sample_cols, names(mx_one))) mx_one[[nm]] <- NA_character_
samples <- rbind(br_one[, sample_cols, drop = FALSE], mx_one[, sample_cols, drop = FALSE])
if (nrow(samples) != 2L) stop("expected BR 11001 and MX 01 primary caches", call. = FALSE)
curve_path <- file.path(base_dir, "results", "downstream_curves_burden_20261112T210050", "curves_unit.csv")
curve_rows <- utils::read.csv(curve_path, stringsAsFactors = FALSE, colClasses = c(unit = "character"), check.names = FALSE)
for (i in seq_len(nrow(samples))) {
  m <- samples[i, , drop = FALSE]
  x <- readRDS(m$cache_path)
  ef <- x$engine_fit
  if (is.null(ef)) ef <- x$fit
  if (is.null(ef) || !identical(toupper(as.character(x$country)), as.character(m$country)) ||
      !identical(as.character(x$unit), as.character(m$unit)) || !identical(as.character(x$scenario), as.character(m$scenario)) ||
      !close_enough(c(x$p95, x$median), c(as.numeric(m$p95), as.numeric(m$median)), 1e-10)) {
    stopf("cache identity/target smoke failed for %s/%s", m$country, m$unit)
  }
  actual_cache_sha <- sha256_file(m$cache_path)
  if (!is.null(m$cache_sha256) && !is.na(m$cache_sha256) && nzchar(m$cache_sha256) &&
      !identical(tolower(m$cache_sha256), actual_cache_sha)) stopf("cache checksum smoke failed for %s/%s", m$country, m$unit)
  beta <- stats::coef(ef$fit)
  covariance <- ef$covariance
  if (!identical(ef$covariance_scope, "full") || is.null(covariance$hac14) || is.null(covariance$hac28)) {
    stopf("full covariance absent for %s/%s", m$country, m$unit)
  }
  L <- cbind(K80 = ef$exposure_effects$y_K80,
    urinary = ef$exposure_effects$y_urinary,
    direct = ef$contrasts$y_urinary_vs_y_K80$contrast_vector)
  rownames(L) <- names(beta)
  if (!close_enough(L[, "direct"], L[, "urinary"] - L[, "K80"], 1e-10)) stopf("direct vector identity failed for %s/%s", m$country, m$unit)
  Hinv <- ef$score$Hinv
  U <- ef$score$daily_scores
  dates <- as.Date(ef$score$daily_dates)
  if (!is.null(rownames(Hinv))) Hinv <- Hinv[names(beta), names(beta), drop = FALSE]
  if (!is.null(colnames(U))) U <- U[, names(beta), drop = FALSE]
  psi <- U %*% Hinv %*% L
  V14 <- hac_matrix(psi, dates, 14L)
  V28 <- hac_matrix(psi, dates, 28L)
  V90 <- hac_matrix(psi, dates, 90L)
  residual14 <- max(abs(V14 - t(L) %*% covariance$hac14 %*% L))
  residual28 <- max(abs(V28 - t(L) %*% covariance$hac28 %*% L))
  if (residual14 > 1e-8 || residual28 > 1e-8) stopf("HAC projection mismatch for %s/%s", m$country, m$unit)
  eig90 <- eigen((V90 + t(V90)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig90) < -1e-8 * max(1, max(abs(eig90)))) stopf("HAC90 projected covariance not PSD for %s/%s", m$country, m$unit)
  if (m$country == "BR") {
    make_curve <- function(endpoint) {
      xhi <- .design_at(ef$fit, ef$data, endpoint = endpoint, basis = ef$basis, temp_value = 26.5)
      xmed <- .design_at(ef$fit, ef$data, endpoint = endpoint, basis = ef$basis, temp_value = x$median)
      as.numeric(xhi - xmed)
    }
    Lk <- make_curve("y_K80"); Lu <- make_curve("y_urinary"); Ld <- Lu - Lk
    curve_est <- sum(Ld * beta); curve_var <- as.numeric(t(Ld) %*% covariance$hac14 %*% Ld)
    z <- curve_rows[curve_rows$unit == m$unit & curve_rows$temperature == 26.5 &
      curve_rows$endpoint == "urinary_vs_K80" & curve_rows$contrast == "direct_RRR", , drop = FALSE]
    if (nrow(z) != 1L || !close_enough(c(curve_est, curve_var), c(z$log_estimate, z$variance), 1e-8)) {
      stopf("26.5 C curve projection mismatch for BR/%s", m$unit)
    }
  }
  cat(sprintf("PASS %s/%s direct-vector-max=%.3g HAC14=%.3g HAC28=%.3g%s\n", m$country, m$unit,
    max(abs(L[, "direct"] - L[, "urinary"] + L[, "K80"])), residual14, residual28,
    if (m$country == "BR") " curve26=matched" else ""))
}
cat("FRESH_PROCESS_SMOKE_PASS\n")
