#!/usr/bin/env Rscript

# Re-pool existing unit-level primary curves on identical support with DL and
# profiled REML-Hartung-Knapp. This script never fits/re-fits an outcome model.
#
# Usage:
#   Rscript audit_curve_pooling.R [project_root] [output_dir]

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1L) args[[1L]] else "/data/admin/fengshuo/latam_stone_repair_20261112"
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(base_dir, "revision_20261113", "results", "marginals_curves")
smoke_mode <- length(args) >= 3L && identical(args[[3L]], "--smoke26.5")
if (length(args) > 3L || (length(args) == 3L && !smoke_mode)) stop("optional third argument must be --smoke26.5", call. = FALSE)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("metafor", quietly = TRUE)) {
  stop("data.table and metafor are required", call. = FALSE)
}
engine_path <- file.path(base_dir, "scripts", "repair_engine.R")
if (!file.exists(engine_path)) stopf("repair_engine.R is missing: %s", engine_path)
if (!requireNamespace("splines", quietly = TRUE)) stop("splines namespace is required for cached DLNM basis prediction", call. = FALSE)
source(engine_path, local = FALSE)

sha256_file <- function(path) {
  exe <- Sys.which("sha256sum")
  if (!nzchar(exe)) stop("sha256sum is required for source provenance", call. = FALSE)
  z <- suppressWarnings(system2(exe, shQuote(path), stdout = TRUE, stderr = TRUE))
  status <- attr(z, "status")
  if (!is.null(status) && status != 0L) stopf("sha256sum failed for %s: %s", path, paste(z, collapse = " "))
  tolower(sub("\\s+.*$", "", z[[1L]]))
}
md5_file <- function(path) unname(tools::md5sum(path)[[1L]])
close_enough <- function(x, y, tol = 1e-9) length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
  all(abs(x - y) <= tol * pmax(1, abs(x), abs(y)))
weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  if (!length(x)) return(NA_real_)
  x[which(cumsum(w) >= p * sum(w))[1L]]
}
atomic_write_csv <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stopf("failed to atomically write %s", path)
}
unit_set_sha256 <- function(x) {
  temp <- tempfile(tmpdir = out_dir)
  on.exit(unlink(temp), add = TRUE)
  writeLines(sort(as.character(x)), temp, useBytes = TRUE)
  sha256_file(temp)
}

profile_reml_hk <- function(y, v) {
  if (length(y) < 2L || length(y) != length(v) || any(!is.finite(y)) || any(!is.finite(v)) || any(v <= 0)) {
    stop("profiled REML requires all finite estimates and strictly positive variances", call. = FALSE)
  }
  objective <- function(tau) {
    w <- 1 / (v + tau)
    mean <- sum(w * y) / sum(w)
    .5 * (sum(log(v + tau)) + log(sum(w)) + sum(w * (y - mean)^2))
  }
  scale <- max(stats::var(y), max((y - stats::median(y))^2), stats::median(v), 1e-8)
  upper <- 10 * scale
  lower <- max(min(v) * 1e-8, 1e-100)
  best <- NA_real_
  best_obj <- Inf
  for (extension in 0:3) {
    grid <- exp(seq(log(lower), log(upper), length.out = 301L))
    values <- vapply(grid, objective, numeric(1))
    candidates <- c(0, grid)
    minima <- which(values <= c(Inf, values[-length(values)]) & values <= c(values[-1L], Inf))
    for (i in minima) {
      if (i > 1L && i < length(grid)) {
        opt <- stats::optimize(function(z) objective(exp(z)), interval = log(grid[c(i - 1L, i + 1L)]), tol = 1e-10)
        candidates <- c(candidates, exp(opt$minimum))
      }
    }
    candidate_obj <- vapply(candidates, objective, numeric(1))
    j <- which.min(candidate_obj)
    best <- candidates[j]; best_obj <- candidate_obj[j]
    if (best < upper * .99) break
    if (extension == 3L) stop("profile REML optimum remains at its upper search bound", call. = FALSE)
    upper <- upper * 100
  }
  fit <- metafor::rma.uni(yi = y, vi = v, tau2 = best, method = "REML", test = "knha")
  list(fit = fit, tau2 = best, objective = best_obj)
}

primary_manifest_path <- file.path(base_dir, "results", "fits", "primary_cache_manifest.csv")
curve_source_dir <- file.path(base_dir, "results", "downstream_curves_burden_20261112T210050")
curve_unit_path <- file.path(curve_source_dir, "curves_unit.csv")
curve_verified_path <- file.path(curve_source_dir, "curves_country_verified.csv")
cache_validation_path <- file.path(out_dir, "curve_26_5_cache_covariance_validation.csv")
required_paths <- c(primary_manifest_path, curve_unit_path, curve_verified_path,
  if (!smoke_mode) cache_validation_path)
missing <- required_paths[!file.exists(required_paths)]
if (length(missing)) stopf("required curve audit inputs missing: %s", paste(missing, collapse = ", "))

manifest <- utils::read.csv(primary_manifest_path, stringsAsFactors = FALSE, check.names = FALSE,
                            colClasses = "character", na.strings = c("", "NA"))
manifest <- manifest[manifest$country == "BR" & manifest$scenario == "primary_threeoutcome_E3_separate", , drop = FALSE]
if (nrow(manifest) != 450L || anyDuplicated(manifest$unit)) stop("expected 450 unique Brazil primary cache units", call. = FALSE)
manifest$shared_median <- as.numeric(manifest$median)
manifest$p95 <- as.numeric(manifest$p95)

curve <- data.table::fread(curve_unit_path, colClasses = c(unit = "character"), na.strings = c("", "NA"), showProgress = TRUE)
verified <- data.table::fread(curve_verified_path, colClasses = c(country = "character"), na.strings = c("", "NA"))
cache26 <- if (file.exists(cache_validation_path)) {
  utils::read.csv(cache_validation_path, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = c(unit = "character"), na.strings = c("", "NA"))
} else data.frame()
if (!smoke_mode) {
  required_cache26 <- c("country", "unit", "temperature", "shared_median", "locally_supported",
    "K80_variance", "urinary_variance", "direct_variance", "min_eigenvalue", "max_eigenvalue", "psd_status")
  missing_cache26 <- setdiff(required_cache26, names(cache26))
  if (length(missing_cache26)) stopf("26.5 C cache covariance table lacks columns: %s", paste(missing_cache26, collapse = ", "))
  if (anyNA(cache26[, required_cache26, drop = FALSE]) ||
      any(!is.finite(as.matrix(cache26[, c("temperature", "shared_median", "K80_variance", "urinary_variance", "direct_variance", "min_eigenvalue", "max_eigenvalue"), drop = FALSE])))) {
    stop("26.5 C cache covariance table contains missing/non-finite required values", call. = FALSE)
  }
  if (nrow(cache26) != 450L || anyDuplicated(cache26$unit) ||
      !setequal(as.character(cache26$unit), as.character(manifest$unit)) || any(cache26$country != "BR") ||
      any(abs(as.numeric(cache26$temperature) - 26.5) > 1e-10)) {
    stop("26.5 C cache covariance table does not contain exactly one row for each BR primary unit", call. = FALSE)
  }
  mismatch_path <- file.path(out_dir, "curve_26_5_cache_mismatch_audit.csv")
  if (file.exists(mismatch_path) && nrow(utils::read.csv(mismatch_path, stringsAsFactors = FALSE)) > 0L) {
    stop("26.5 C cache predictions differ from saved unit curve rows", call. = FALSE)
  }
}
curve <- curve[country == "BR"]
if (!nrow(curve)) stop("no Brazil unit curves in source file", call. = FALSE)
if (anyDuplicated(curve[, .(country, unit, temperature, endpoint, contrast)])) stop("duplicate unit/temperature/endpoint curve rows", call. = FALSE)
source_units <- sort(unique(as.character(curve$unit)))
manifest_units <- sort(as.character(manifest$unit))
unit_set_match <- identical(source_units, manifest_units)
if (!unit_set_match) stop("unit curve source and BR primary manifest do not have identical unit support", call. = FALSE)
curve[, temperature := as.numeric(temperature)]
curve[, log_estimate := as.numeric(log_estimate)]
curve[, variance := as.numeric(variance)]
curve[, locally_supported := as.logical(locally_supported)]
grid <- sort(unique(curve$temperature))
expected_grid <- seq(14, 30, by = .5)
if (smoke_mode) {
  curve <- curve[abs(temperature - 26.5) < 1e-10 & endpoint == "urinary_vs_K80" & contrast == "direct_RRR"]
  grid <- sort(unique(curve$temperature))
  if (nrow(curve) != 450L || length(grid) != 1L || grid[1L] != 26.5) stop("26.5 C direct-RRR smoke input is incomplete", call. = FALSE)
} else if (!identical(as.numeric(grid), as.numeric(expected_grid))) stop("curve temperatures must cover exactly 14-30 C in 0.5 C steps", call. = FALSE)
if (any(curve$temperature > 30)) stop("curve source contains points beyond 30 C", call. = FALSE)

target_lookup <- data.table::as.data.table(manifest[, c("unit", "shared_median", "p95"), drop = FALSE])
curve <- merge(curve, target_lookup, by = "unit", all.x = TRUE, sort = FALSE)
if (anyNA(curve$shared_median)) stop("curve unit lacks primary-cache target median", call. = FALSE)
curve[, absolute_distance_from_local_median := abs(temperature - shared_median)]

groups <- unique(curve[, .(country, temperature, endpoint, contrast)])
pooled_rows <- list()
audit_rows <- list()
support_rows <- list()
weight_26 <- NULL
failure_rows <- list()
for (j in seq_len(nrow(groups))) {
  g <- groups[j]
  z_all <- curve[country == g$country & temperature == g$temperature & endpoint == g$endpoint & contrast == g$contrast]
  if (anyDuplicated(z_all$unit)) stopf("duplicate unit within %s/%s/%s/%s", g$country, g$temperature, g$endpoint, g$contrast)
  z_supported <- z_all[locally_supported == TRUE]
  for (scope in c("all_successful", "support_restricted")) {
    z <- if (scope == "all_successful") z_all else z_supported
    if (!nrow(z)) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(country = g$country, temperature = g$temperature,
        endpoint = g$endpoint, contrast = g$contrast, scope = scope, failure = "zero units in predefined pooling scope", stringsAsFactors = FALSE)
      next
    }
    n_bad_y <- sum(!is.finite(z$log_estimate))
    n_bad_v <- sum(!is.finite(z$variance))
    n_nonpositive_v <- sum(is.finite(z$variance) & z$variance <= 0)
    if (n_bad_y || n_bad_v || n_nonpositive_v) {
      bad <- z[!is.finite(log_estimate) | !is.finite(variance) | variance <= 0]
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(country = g$country, temperature = g$temperature,
        endpoint = g$endpoint, contrast = g$contrast, scope = scope,
        failure = sprintf("invalid source unit records retained: yi=%d, vi=%d, vi_nonpositive=%d", n_bad_y, n_bad_v, n_nonpositive_v),
        stringsAsFactors = FALSE)
      if (nrow(bad)) {
        badout <- as.data.frame(bad[, .(country, unit, temperature, endpoint, contrast, log_estimate, variance, locally_supported)])
        badout$scope <- scope
        failure_rows[[length(failure_rows) + 1L]] <- badout
      }
      next
    }
    pool_dl <- pool_meta(z$log_estimate, z$variance, method = "DL")
    reml <- tryCatch(profile_reml_hk(z$log_estimate, z$variance), error = function(e) e)
    if (inherits(reml, "error")) {
      failure_rows[[length(failure_rows) + 1L]] <- data.frame(country = g$country, temperature = g$temperature,
        endpoint = g$endpoint, contrast = g$contrast, scope = scope,
        failure = paste("profiled REML-HK failed:", conditionMessage(reml)), stringsAsFactors = FALSE)
      next
    }
    dl_tau <- pool_dl$tau2
    dl_w <- 1 / (z$variance + dl_tau); dl_w <- dl_w / sum(dl_w)
    reml_tau <- reml$tau2
    reml_w <- 1 / (z$variance + reml_tau); reml_w <- reml_w / sum(reml_w)
    old_contrast <- paste(g$contrast, scope, sep = ":")
    old_dl <- verified[country == g$country & abs(temperature - g$temperature) < 1e-10 & endpoint == g$endpoint & contrast == old_contrast & method == "DL"]
    old_reml <- verified[country == g$country & abs(temperature - g$temperature) < 1e-10 & endpoint == g$endpoint & contrast == old_contrast & method == "REML_profile_HK"]
    if (nrow(old_dl) != 1L || nrow(old_reml) != 1L) {
      old_dl_status <- "prior_verified_curve_row_missing_or_nonunique"
      old_reml_status <- "prior_profiled_REML_row_missing_or_nonunique"
      diff_dl <- diff_reml <- NA_real_
    } else {
      diff_dl <- max(abs(c(pool_dl$RR - old_dl$RR, pool_dl$RR95CI[1L] - old_dl$lo, pool_dl$RR95CI[2L] - old_dl$hi,
                           pool_dl$tau2 - old_dl$tau2)))
      diff_reml <- max(abs(c(exp(as.numeric(reml$fit$b[1L])) - old_reml$RR,
                             exp(reml$fit$ci.lb) - old_reml$lo, exp(reml$fit$ci.ub) - old_reml$hi,
                             reml_tau - old_reml$tau2)))
      old_dl_status <- if (diff_dl <= 1e-9) "matched_prior_DL_within_1e-9" else "prior_DL_mismatch"
      old_reml_status <- if (diff_reml <= 1e-9) "matched_prior_profile_REML_within_1e-9" else "prior_profile_REML_mismatch"
    }
    near <- z$absolute_distance_from_local_median
    common <- data.frame(country = g$country, temperature = g$temperature, endpoint = g$endpoint,
      contrast = g$contrast, curve_scope = scope, k = nrow(z),
      support_unit_ids_sha256 = unit_set_sha256(z$unit),
      n_invalid_yi = n_bad_y, n_invalid_variance = n_bad_v, n_nonpositive_variance = n_nonpositive_v,
      minimum_unit_variance = min(z$variance), q01_unit_variance = as.numeric(stats::quantile(z$variance, .01, names = FALSE)),
      median_unit_variance = stats::median(z$variance), maximum_unit_variance = max(z$variance),
      maximum_to_minimum_variance_ratio = max(z$variance) / min(z$variance),
      variance_floor_applied = FALSE,
      absolute_distance_from_local_median_min = min(near),
      absolute_distance_from_local_median_q25 = as.numeric(stats::quantile(near, .25, names = FALSE)),
      absolute_distance_from_local_median_median = stats::median(near),
      absolute_distance_from_local_median_q75 = as.numeric(stats::quantile(near, .75, names = FALSE)),
      absolute_distance_from_local_median_max = max(near),
      units_within_1C_of_local_median = sum(near <= 1),
      unweighted_fraction_within_1C_of_local_median = mean(near <= 1),
      DL_weighted_median_absolute_distance_from_local_median = weighted_quantile(near, dl_w, .5),
      DL_weighted_fraction_within_1C_of_local_median = sum(dl_w[near <= 1]),
      REML_profile_HK_weighted_median_absolute_distance_from_local_median = weighted_quantile(near, reml_w, .5),
      REML_profile_HK_weighted_fraction_within_1C_of_local_median = sum(reml_w[near <= 1]),
      DL_max_normalized_weight = max(dl_w), DL_effective_unit_count = 1 / sum(dl_w^2),
      DL_sum_weight_sq = sum(dl_w^2),
      REML_profile_HK_max_normalized_weight = max(reml_w),
      REML_profile_HK_effective_unit_count = 1 / sum(reml_w^2),
      REML_profile_HK_sum_weight_sq = sum(reml_w^2),
      DL_prior_reproduction_status = old_dl_status, DL_max_abs_difference_from_prior = diff_dl,
      REML_prior_reproduction_status = old_reml_status, REML_max_abs_difference_from_prior = diff_reml,
      stringsAsFactors = FALSE)
    audit_rows[[length(audit_rows) + 1L]] <- common
    support_rows[[length(support_rows) + 1L]] <- data.frame(country = g$country, temperature = g$temperature,
      endpoint = g$endpoint, contrast = g$contrast, curve_scope = scope, k = nrow(z),
      support_unit_ids_sha256 = common$support_unit_ids_sha256, stringsAsFactors = FALSE)
    for (method in c("DL", "REML_profile_HK")) {
      if (method == "DL") {
        fit <- pool_dl
        rr <- fit$RR; lo <- fit$RR95CI[1L]; hi <- fit$RR95CI[2L]
        le <- fit$estimate; se <- fit$se; tau2 <- fit$tau2; I2 <- fit$I2; concordance <- fit$concordance
        max_weight <- max(dl_w); effective_n <- 1 / sum(dl_w^2)
      } else {
        fit <- reml$fit
        le <- as.numeric(fit$b[1L]); se <- as.numeric(fit$se); rr <- exp(le)
        lo <- exp(as.numeric(fit$ci.lb)); hi <- exp(as.numeric(fit$ci.ub)); tau2 <- reml_tau
        I2 <- as.numeric(fit$I2); concordance <- mean(z$log_estimate > 0)
        max_weight <- max(reml_w); effective_n <- 1 / sum(reml_w^2)
      }
      pooled_rows[[length(pooled_rows) + 1L]] <- data.frame(country = g$country,
        endpoint = g$endpoint, contrast = old_contrast, method = method, covariance = "unit_HAC14",
        inference = if (method == "DL") "independent_unit_DL" else "independent_unit_REML_HK",
        k = nrow(z), estimate_logRR = le, se = se, RR = rr, lo = lo, hi = hi,
        I2 = I2, concordance = concordance, tau2 = tau2,
        maximum_normalized_weight = max_weight, effective_unit_count = effective_n,
        temperature = g$temperature, curve_scope = scope,
        prior_reproduction_status = if (method == "DL") old_dl_status else old_reml_status,
        max_abs_difference_from_prior = if (method == "DL") diff_dl else diff_reml,
        stringsAsFactors = FALSE)
    }
    if (abs(g$temperature - 26.5) < 1e-10 && g$endpoint == "urinary_vs_K80" && g$contrast == "direct_RRR") {
      z_out <- as.data.frame(z)
      z_out$DL_tau2 <- dl_tau; z_out$DL_normalized_weight <- dl_w
      z_out$DL_weight_sq <- dl_w^2
      z_out$REML_profile_HK_tau2 <- reml_tau; z_out$REML_profile_HK_normalized_weight <- reml_w
      z_out$REML_profile_HK_weight_sq <- reml_w^2
      z_out$curve_scope <- scope
      if (is.null(weight_26)) weight_26 <- z_out else weight_26 <- rbind(weight_26, z_out)
    }
  }
  if (j %% 50L == 0L || j == nrow(groups)) message(sprintf("pooled %d/%d curve-temperature groups", j, nrow(groups)))
}

if (length(pooled_rows)) atomic_write_csv(do.call(rbind, pooled_rows), file.path(out_dir, "curves_pooled_revision.csv"))
if (length(audit_rows)) atomic_write_csv(do.call(rbind, audit_rows), file.path(out_dir, "curve_pooling_weight_variance_audit.csv"))
if (length(support_rows)) atomic_write_csv(do.call(rbind, support_rows), file.path(out_dir, "curve_support_identity_audit.csv"))
if (!is.null(weight_26)) atomic_write_csv(weight_26, file.path(out_dir, "curve_26_5_direct_unit_weights.csv"))
if (length(failure_rows)) atomic_write_csv(do.call(data.table::rbindlist, c(failure_rows, list(fill = TRUE))), file.path(out_dir, "curve_pooling_failure_audit.csv"))

if (nrow(cache26) > 0L) {
  psd_summary <- data.frame(n_units = nrow(cache26), n_not_psd = sum(cache26$psd_status == "not_PSD", na.rm = TRUE),
    n_nonpositive_component_variances = sum(cache26$K80_variance <= 0 | cache26$urinary_variance <= 0 | cache26$direct_variance <= 0, na.rm = TRUE),
    minimum_projected_eigenvalue = min(cache26$min_eigenvalue, na.rm = TRUE),
    maximum_absolute_projected_eigenvalue = max(abs(cache26$max_eigenvalue), na.rm = TRUE),
    n_curve_cache_value_mismatches = if (file.exists(file.path(out_dir, "curve_26_5_cache_mismatch_audit.csv"))) nrow(utils::read.csv(file.path(out_dir, "curve_26_5_cache_mismatch_audit.csv"))) else 0L,
    validation_status = if (nrow(cache26) == 450L && all(!is.na(cache26$psd_status)) &&
      all(cache26$psd_status == "PSD_within_numeric_tolerance") &&
      all(cache26$K80_variance > 0 & cache26$urinary_variance > 0 & cache26$direct_variance > 0)) "all_450_projected_covariances_PSD_and_positive" else "projected_covariance_or_variance_issue",
    stringsAsFactors = FALSE)
} else psd_summary <- data.frame(n_units = 0L, n_not_psd = 0L, n_nonpositive_component_variances = 0L,
  minimum_projected_eigenvalue = NA_real_, maximum_absolute_projected_eigenvalue = NA_real_,
  n_curve_cache_value_mismatches = NA_integer_, validation_status = "not_run_in_26_5_smoke_mode", stringsAsFactors = FALSE)
atomic_write_csv(psd_summary, file.path(out_dir, "curve_cache_covariance_psd_summary.csv"))

example <- if (length(pooled_rows)) do.call(rbind, pooled_rows) else data.frame()
if (nrow(example)) example <- example[example$country == "BR" & example$temperature == 26.5 &
  example$endpoint == "urinary_vs_K80" & example$curve_scope == "support_restricted", , drop = FALSE]
audit_table <- if (length(audit_rows)) do.call(rbind, audit_rows) else data.frame()
example_audit <- if (nrow(audit_table)) audit_table[audit_table$country == "BR" & audit_table$temperature == 26.5 &
  audit_table$endpoint == "urinary_vs_K80" & audit_table$contrast == "direct_RRR" &
  audit_table$curve_scope == "support_restricted", , drop = FALSE] else data.frame()
example_check <- if (nrow(example) == 2L && nrow(example_audit) == 1L) {
  dl <- example[example$method == "DL", , drop = FALSE]
  reml <- example[example$method == "REML_profile_HK", , drop = FALSE]
  a <- example_audit[1L, , drop = FALSE]
  data.frame(country = "BR", temperature = 26.5, endpoint = "urinary_vs_K80", curve_scope = "support_restricted",
    expected_k = 322L, observed_k = dl$k,
    DL_RRR = dl$RR, DL_lo = dl$lo, DL_hi = dl$hi,
    REML_profile_HK_RRR = reml$RR, REML_profile_HK_lo = reml$lo, REML_profile_HK_hi = reml$hi,
    DL_max_normalized_weight = dl$maximum_normalized_weight, DL_effective_unit_count = dl$effective_unit_count,
    REML_max_normalized_weight = reml$maximum_normalized_weight, REML_effective_unit_count = reml$effective_unit_count,
    unweighted_median_abs_median_distance_C = a$absolute_distance_from_local_median_median,
    max_abs_median_distance_C = a$absolute_distance_from_local_median_max,
    DL_weighted_median_abs_median_distance_C = a$DL_weighted_median_absolute_distance_from_local_median,
    REML_weighted_median_abs_median_distance_C = a$REML_profile_HK_weighted_median_absolute_distance_from_local_median,
    DL_prior_reproduction_status = dl$prior_reproduction_status,
    REML_prior_reproduction_status = reml$prior_reproduction_status,
    exact_support_identity_hash = a$support_unit_ids_sha256,
    expected_value_check = if (dl$k == 322L && close_enough(c(dl$RR, reml$RR), c(1.00838604411303, 1.02965774870279), 1e-9)) "matched" else "mismatch",
    stringsAsFactors = FALSE)
} else data.frame(country = "BR", temperature = 26.5, endpoint = "urinary_vs_K80", curve_scope = "support_restricted",
  expected_k = 322L, observed_k = NA_integer_, expected_value_check = "missing_or_nonunique_curve_rows", stringsAsFactors = FALSE)
atomic_write_csv(example_check, file.path(out_dir, "curve_26_5_identical_support_example.csv"))

grid_complete <- if (smoke_mode) length(grid) == 1L && grid[1L] == 26.5 else length(grid) == 33L && min(grid) == 14 && max(grid) == 30
prior_mismatch <- length(audit_rows) > 0L && any(vapply(audit_rows, function(x) {
  x$DL_prior_reproduction_status != "matched_prior_DL_within_1e-9" ||
    x$REML_prior_reproduction_status != "matched_prior_profile_REML_within_1e-9"
}, logical(1)))
has_failures <- length(failure_rows) > 0L
status <- if (unit_set_match && grid_complete && !prior_mismatch && !has_failures &&
  (smoke_mode || psd_summary$validation_status == "all_450_projected_covariances_PSD_and_positive") &&
  example_check$expected_value_check == "matched") "accepted" else "rejected_or_incomplete"

source_files <- c(primary_manifest_path, curve_unit_path, curve_verified_path, if (file.exists(cache_validation_path)) cache_validation_path)
source_hashes <- data.frame(path = normalizePath(source_files, winslash = "/", mustWork = TRUE),
  md5 = vapply(source_files, md5_file, character(1)), sha256 = vapply(source_files, sha256_file, character(1)),
  stringsAsFactors = FALSE)
atomic_write_csv(source_hashes, file.path(out_dir, "curve_source_file_hashes.csv"))
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]), winslash = "/", mustWork = TRUE)
run_manifest <- list(run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  script_path = script_path, script_sha256 = sha256_file(script_path),
  engine_path = normalizePath(engine_path, winslash = "/", mustWork = TRUE), engine_sha256 = sha256_file(engine_path),
  output_dir = out_dir, curve_unit_source = normalizePath(curve_unit_path, winslash = "/", mustWork = TRUE),
  exact_source_unit_set_matches_BR_primary_manifest = unit_set_match,
  temperature_grid = grid, temperature_grid_check_passed = grid_complete,
  point_through_30C_verified = !smoke_mode && grid_complete, smoke_mode = smoke_mode,
  method = "DL and profiled REML-HK on identical per-point, per-scope unit sets; all supported points retained",
  variance_policy = "No unit variance floor, clipping, or significance-based unit removal; invalid source rows are retained in the failure audit and prevent acceptance.",
  diagnostic_definitions = list(maximum_normalized_weight = "max of 1/(unit variance + estimator tau2), normalized to sum 1",
    effective_unit_count = "1/sum(normalized_weight^2)",
    near_median_distance = "absolute temperature minus each included unit's cached shared median E3 target; both weighted diagnostics and raw unit values retained at 26.5C"),
  example_26_5C = example_check, cache_covariance_psd = psd_summary,
  status = status, n_curve_groups = nrow(groups), n_pooling_scopes = 2L,
  R_version = R.version.string)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for run_manifest.json", call. = FALSE)
jsonlite::write_json(run_manifest, file.path(out_dir, "curve_run_manifest.json"), auto_unbox = TRUE,
                     pretty = TRUE, na = "null", null = "null", digits = NA)
message(sprintf("Curve pooling audit status: %s; outputs: %s", status, out_dir))
if (status != "accepted") quit(status = 2L)
