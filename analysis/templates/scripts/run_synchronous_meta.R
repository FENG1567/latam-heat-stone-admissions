#!/usr/bin/env Rscript
# Calendar-synchronous uncertainty sensitivity for pooled urinary-versus-K80
# log rate-ratio contrasts. This script consumes full primary fit caches; it
# does not refit any model.
#
# Supported cache formats:
#   * Frozen primary_cache_manifest.csv: manifest cache_path/input_md5/
#     protocol_signature; RDS fields $engine_fit, $input_md5 and
#     $protocol_signature.
#   * Updated-MX fit_cache_manifest.csv: manifest cache_path/cache_sha256/
#     input_sha256/engine_sha256; RDS fields $fit, $input_hashes_sha256 and
#     $engine_hash_sha256. This updated manifest can contain primary and
#     controlled-fill scenarios; the default selects the primary scenario.
#
# Example:
# Rscript run_synchronous_meta.R --manifest results/fits/primary_cache_manifest.csv \
#   --country both --outdir results/synchronous_meta

.arg_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.arg_file)) stop("cannot resolve script path from --file", call. = FALSE)
.script_path <- normalizePath(sub("^--file=", "", .arg_file[1L]), winslash = "/", mustWork = TRUE)
.script_dir <- dirname(.script_path)
.repair_dir <- normalizePath(file.path(.script_dir, ".."), winslash = "/", mustWork = TRUE)
.engine_path <- file.path(.script_dir, "repair_engine.R")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

parse_args <- function(args) {
  out <- list(manifest = NULL, country = "both", outdir = NULL,
              scenario = NULL, self_test = FALSE, help = FALSE)
  i <- 1L
  while (i <= length(args)) {
    token <- args[i]
    if (token %in% c("--self-test", "--help", "-h")) {
      if (token == "--self-test") out$self_test <- TRUE else out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(token, "--")) stopf("unexpected argument '%s'", token)
    kv <- substring(token, 3L)
    if (grepl("=", kv, fixed = TRUE)) {
      bits <- strsplit(kv, "=", fixed = TRUE)[[1L]]
      key <- bits[1L]; value <- paste(bits[-1L], collapse = "="); i <- i + 1L
    } else {
      key <- kv
      if (i == length(args)) stopf("missing value for --%s", key)
      value <- args[i + 1L]; i <- i + 2L
    }
    if (!key %in% c("manifest", "country", "outdir", "scenario")) stopf("unknown option --%s", key)
    out[[key]] <- value
  }
  if (out$help) return(out)
  if (!out$self_test && is.null(out$manifest)) stopf("--manifest PATH is required")
  out$country <- toupper(out$country)
  if (!out$country %in% c("BR", "MX", "BOTH")) stopf("--country must be BR, MX, or both")
  if (!is.null(out$manifest)) {
    if (!file.exists(out$manifest)) stopf("manifest does not exist: %s", out$manifest)
    out$manifest <- normalizePath(out$manifest, winslash = "/", mustWork = TRUE)
  }
  if (is.null(out$outdir)) out$outdir <- file.path(.repair_dir, "results", "synchronous_meta")
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", out$outdir)) out$outdir <- file.path(.repair_dir, out$outdir)
  if (!out$self_test) {
    dir.create(out$outdir, recursive = TRUE, showWarnings = FALSE)
    out$outdir <- normalizePath(out$outdir, winslash = "/", mustWork = TRUE)
  }
  out
}

show_help <- function() {
  cat(paste(
    "Usage: Rscript run_synchronous_meta.R --manifest PATH [--country BR|MX|both]",
    "       [--scenario NAME] [--outdir PATH]",
    "       Rscript run_synchronous_meta.R --self-test",
    "",
    "The frozen manifest uses $engine_fit and protocol_signature/input_md5 metadata.",
    "The updated-MX manifest uses $fit and cache_sha256/input_sha256/engine_sha256 metadata.",
    "For a mixed-scenario updated-MX manifest, the default chooses primary_complete_source_2015_2022;
    use --scenario to select another exact scenario.",
    sep = "\n"))
}

resolve_column <- function(d, candidates, required = TRUE, label = candidates[1L]) {
  hit <- candidates[candidates %in% names(d)]
  if (length(hit)) return(hit[1L])
  if (required) stopf("manifest lacks %s (accepted columns: %s)", label, paste(candidates, collapse = ", "))
  NULL
}

column_values <- function(d, candidates, required = TRUE, label = candidates[1L]) {
  nm <- resolve_column(d, candidates, required, label)
  if (is.null(nm)) return(NULL)
  as.character(d[[nm]])
}

as_date <- function(x, label = "dates") {
  if (inherits(x, "Date")) ans <- x
  else if (is.numeric(x)) ans <- as.Date(x, origin = "1970-01-01")
  else ans <- as.Date(as.character(x))
  if (anyNA(ans)) stopf("%s contain invalid values", label)
  ans
}

bartlett_hac_scalar <- function(psi, dates, lag) {
  if (!is.numeric(psi) || length(psi) != length(dates) || !length(psi)) {
    stopf("influence/date vectors must have equal positive lengths")
  }
  if (length(lag) != 1L || !is.finite(lag) || lag < 0 || lag != as.integer(lag)) {
    stopf("lag must be one non-negative integer")
  }
  dates <- as_date(dates, "HAC dates")
  if (anyDuplicated(dates)) stopf("HAC dates must be unique within an influence series")
  if (any(!is.finite(psi))) stopf("influence vector contains non-finite values")
  ord <- order(dates)
  dates <- dates[ord]; psi <- psi[ord]
  calendar <- seq(min(dates), max(dates), by = "day")
  full <- numeric(length(calendar))
  full[as.integer(dates - min(dates)) + 1L] <- psi
  out <- sum(full^2)
  lag <- as.integer(lag)
  if (lag > 0L && length(full) > 1L) {
    for (ell in seq_len(min(lag, length(full) - 1L))) {
      cross <- sum(full[(ell + 1L):length(full)] * full[seq_len(length(full) - ell)])
      out <- out + 2 * (1 - ell / (lag + 1)) * cross
    }
  }
  if (!is.finite(out) || out < -1e-8) stopf("Bartlett HAC variance is invalid")
  max(0, out)
}

check_close <- function(a, b, tol = 1e-8, label = "values") {
  if (length(a) != length(b) || any(!is.finite(a)) || any(!is.finite(b)) ||
      any(abs(a - b) > tol * pmax(1, abs(a), abs(b)))) {
    stopf("%s do not agree within relative tolerance %.3g", label, tol)
  }
  invisible(TRUE)
}

run_self_tests <- function() {
  # Two observations separated by 14 actual calendar days receive Bartlett
  # lag-14 weight 1/15 even when the intervening dates were not observed.
  d0 <- as.Date("2020-01-01")
  v14 <- bartlett_hac_scalar(c(1, 1), c(d0, d0 + 14), lag = 14L)
  check_close(v14, 2 + 2 / 15, label = "actual-calendar lag-14 fixture")
  v15 <- bartlett_hac_scalar(c(1, 1), c(d0, d0 + 15), lag = 14L)
  check_close(v15, 2, label = "actual-calendar lag-15 fixture")
  perm <- bartlett_hac_scalar(c(1, 1), c(d0 + 14, d0), lag = 14L)
  check_close(perm, v14, label = "date/value permutation fixture")
  check_close(bartlett_hac_scalar(2 * c(1, -0.3, 0.4), d0 + 0:2, 14L),
              4 * bartlett_hac_scalar(c(1, -0.3, 0.4), d0 + 0:2, 14L),
              label = "doubled-influence variance fixture")
  one_unit <- bartlett_hac_scalar(c(1, -0.3, 0.4), d0 + 0:2, 14L)
  duplicate_mean <- (c(1, -0.3, 0.4) + c(1, -0.3, 0.4)) / 2
  check_close(bartlett_hac_scalar(duplicate_mean, d0 + 0:2, 14L), one_unit,
              label = "normalized identical duplicated-unit variance fixture")

  if (!file.exists(.engine_path)) stopf("repair_engine.R is missing for projection self-test")
  source(.engine_path, local = environment())
  U <- matrix(c(1, 0.5, -0.2, 0.3, 0.8, -0.4, 1.1, -0.6), ncol = 2L,
              dimnames = list(NULL, c("b1", "b2")))
  dates <- as.Date("2021-01-01") + c(0, 1, 4, 5)
  Hinv <- matrix(c(0.8, 0.1, 0.1, 0.6), nrow = 2L)
  L <- c(b1 = 0.7, b2 = -0.2)
  influence <- as.numeric(crossprod(L, Hinv %*% t(U)))
  engine_influence <- as.numeric(U %*% (Hinv %*% L))
  check_close(influence, engine_influence, label = "influence matrix orientation")
  stored <- as.numeric(t(L) %*% Hinv %*% .bartlett_hac_meat(U, dates, 14L) %*% Hinv %*% L)
  reconstructed <- bartlett_hac_scalar(influence, dates, 14L)
  check_close(reconstructed, stored, label = "synthetic projected HAC reconstruction")
  message("Synchronous meta self-tests passed: calendar gaps, lag kernel, permutation, scaling, and influence orientation.")
  invisible(TRUE)
}

sha256_file <- function(path) {
  exe <- Sys.which("sha256sum")
  if (nzchar(exe)) {
    z <- suppressWarnings(system2(exe, shQuote(path), stdout = TRUE, stderr = TRUE))
    status <- attr(z, "status")
    if (!is.null(status) && status != 0L) stopf("sha256sum failed for %s: %s", path, paste(z, collapse = " "))
    return(tolower(sub("\\s+.*$", "", z[1L])))
  }
  exe <- Sys.which("shasum")
  if (nzchar(exe)) {
    z <- suppressWarnings(system2(exe, c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE))
    status <- attr(z, "status")
    if (!is.null(status) && status != 0L) stopf("shasum failed for %s: %s", path, paste(z, collapse = " "))
    return(tolower(sub("\\s+.*$", "", z[1L])))
  }
  stopf("a sha256sum or shasum executable is required to verify this cache")
}

md5_file <- function(path) unname(tools::md5sum(path)[1L])

is_absolute_path <- function(path) grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
resolve_cache_path <- function(path, manifest_path) {
  if (is.na(path) || !nzchar(path)) stopf("manifest contains an empty cache path")
  if (!is_absolute_path(path)) path <- file.path(dirname(manifest_path), path)
  if (!file.exists(path)) stopf("cache file does not exist: %s", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

pick_scenario <- function(rows, requested = NULL) {
  scen <- unique(rows$scenario)
  if (!is.null(requested)) {
    if (!requested %in% scen) stopf("requested scenario '%s' is absent; found: %s", requested, paste(scen, collapse = ", "))
    return(requested)
  }
  priority <- c("primary_threeoutcome_E3_separate", "primary_complete_source_2015_2022")
  hit <- priority[priority %in% scen]
  if (length(hit)) return(hit[1L])
  primary <- scen[grepl("^primary", scen)]
  if (length(primary) == 1L) return(primary)
  if (!length(primary)) stopf("manifest has no primary scenario")
  stopf("manifest has multiple primary scenarios; select one with --scenario: %s", paste(primary, collapse = ", "))
}

read_manifest <- function(path, country_request, requested_scenario = NULL) {
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                              colClasses = "character", na.strings = c("", "NA"))
  if (!nrow(manifest)) stopf("manifest is empty: %s", path)
  cn <- resolve_column(manifest, c("country"), TRUE, "country")
  un <- resolve_column(manifest, c("unit", "state", "geography"), TRUE, "unit identifier")
  sc <- resolve_column(manifest, c("scenario", "analysis", "cohort"), TRUE, "scenario")
  cp <- resolve_column(manifest, c("cache_path", "rds_path", "path"), TRUE, "cache path")
  manifest$country <- toupper(as.character(manifest[[cn]]))
  manifest$unit <- as.character(manifest[[un]])
  manifest$scenario <- as.character(manifest[[sc]])
  manifest$cache_path <- as.character(manifest[[cp]])
  wanted <- if (country_request == "BOTH") c("BR", "MX") else country_request
  absent <- setdiff(wanted, unique(manifest$country))
  if (length(absent)) stopf("requested country/countries absent from manifest: %s", paste(absent, collapse = ", "))
  selected <- list()
  for (cty in wanted) {
    rows <- manifest[manifest$country == cty, , drop = FALSE]
    chosen <- pick_scenario(rows, requested_scenario)
    rows <- rows[rows$scenario == chosen, , drop = FALSE]
    if (anyDuplicated(rows$unit)) stopf("manifest repeats unit IDs for %s/%s", cty, chosen)
    rows$cache_path <- vapply(rows$cache_path, resolve_cache_path, character(1), manifest_path = path)
    rows$.cache_sha256 <- column_values(rows, c("cache_sha256", "rds_sha256", "sha256"), FALSE)
    rows$.cache_md5 <- column_values(rows, c("cache_md5", "rds_md5"), FALSE)
    rows$.input_md5 <- column_values(rows, c("input_md5"), FALSE)
    rows$.protocol_signature <- column_values(rows, c("protocol_signature"), FALSE)
    rows$.input_sha256 <- column_values(rows, c("input_sha256", "input_hash"), FALSE)
    rows$.engine_sha256 <- column_values(rows, c("engine_sha256", "engine_hash_sha256"), FALSE)
    selected[[cty]] <- list(scenario = chosen, rows = rows)
  }
  selected
}

first_nonempty <- function(x, candidates, fallback = NULL) {
  for (nm in candidates) if (!is.null(x[[nm]]) && length(x[[nm]]) && !all(is.na(x[[nm]]))) return(x[[nm]])
  fallback
}

residual_acf_rows <- function(engine_fit, country, unit, scenario) {
  # Residual ACF is calculated against actual date-minus-lag matches. A missing
  # date removes that pair; rows are never compressed into adjacent positions.
  d <- engine_fit$row_data
  s <- engine_fit$score
  if (is.null(d) || !"disease" %in% names(d) || is.null(s$y) || is.null(s$mu) || is.null(s$dates)) {
    stopf("cache lacks stacked row_data disease or score y/mu/dates required for residual ACF")
  }
  if (nrow(d) != length(s$y) || length(s$y) != length(s$mu) || length(s$y) != length(s$dates)) {
    stopf("residual row_data/score dimensions disagree")
  }
  disease <- as.character(d$disease)
  dates <- as_date(s$dates, "residual score dates")
  y <- as.numeric(s$y); mu <- as.numeric(s$mu)
  if (any(!is.finite(y)) || any(!is.finite(mu)) || any(mu <= 0)) stopf("invalid y or fitted means in cached residuals")
  pearson <- (y - mu) / sqrt(mu)
  endpoints <- c("y_K80", "y_K81", "y_urinary")
  lags <- c(1L, 7L, 14L, 28L, 90L)
  out <- list()
  for (endpoint in endpoints) {
    idx <- which(disease == endpoint)
    if (!length(idx)) stopf("cached stacked model lacks residual endpoint %s", endpoint)
    ep_dates <- dates[idx]
    if (anyDuplicated(ep_dates)) stopf("residual endpoint %s has duplicate actual dates", endpoint)
    r <- pearson[idx]
    date_key <- as.character(ep_dates)
    for (lag in lags) {
      prior <- match(as.character(ep_dates - lag), date_key)
      keep <- !is.na(prior)
      rho <- if (sum(keep) >= 3L) suppressWarnings(stats::cor(r[keep], r[prior[keep]])) else NA_real_
      out[[length(out) + 1L]] <- data.frame(country = country, scenario = scenario, unit = unit,
        endpoint = endpoint, lag_calendar_days = lag, pairwise_Pearson_residual_ACF = rho,
        n_actual_calendar_pairs = sum(keep), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

load_projected_unit <- function(row, country, scenario) {
  path <- row$cache_path
  if (!is.null(row$.cache_sha256) && !is.na(row$.cache_sha256) && nzchar(row$.cache_sha256)) {
    actual <- sha256_file(path)
    if (!identical(tolower(actual), tolower(row$.cache_sha256))) stopf("cache SHA-256 mismatch: %s", path)
    checksum_method <- "sha256_verified"
  } else if (!is.null(row$.cache_md5) && !is.na(row$.cache_md5) && nzchar(row$.cache_md5)) {
    actual <- md5_file(path)
    if (!identical(tolower(actual), tolower(row$.cache_md5))) stopf("cache MD5 mismatch: %s", path)
    checksum_method <- "md5_verified"
  } else {
    checksum_method <- "manifest_has_no_cache_checksum"
  }
  cache <- readRDS(path)
  if (!is.list(cache)) stopf("cache is not an R list: %s", path)
  engine_fit <- first_nonempty(cache, c("engine_fit", "fit"))
  if (is.null(engine_fit) || !is.list(engine_fit)) stopf("cache lacks $engine_fit or $fit: %s", path)
  cached_country <- toupper(as.character(first_nonempty(cache, c("country"), country))[1L])
  cached_unit <- as.character(first_nonempty(cache, c("unit"), row$unit)[1L])
  cached_scenario <- as.character(first_nonempty(cache, c("scenario"), scenario)[1L])
  if (!identical(cached_country, country) || !identical(cached_unit, as.character(row$unit)) ||
      !identical(cached_scenario, scenario)) {
    stopf("cache identity mismatch for %s/%s/%s (cache says %s/%s/%s)",
          country, row$unit, scenario, cached_country, cached_unit, cached_scenario)
  }
  if (!is.null(row$.input_md5) && !is.na(row$.input_md5) && nzchar(row$.input_md5)) {
    cached_hash <- first_nonempty(cache, c("input_md5"))
    if (is.null(cached_hash) || !identical(tolower(as.character(cached_hash[1L])), tolower(row$.input_md5))) {
      stopf("input MD5 metadata mismatch for %s/%s", country, row$unit)
    }
  }
  if (!is.null(row$.protocol_signature) && !is.na(row$.protocol_signature) && nzchar(row$.protocol_signature)) {
    cached_sig <- first_nonempty(cache, c("protocol_signature"))
    if (is.null(cached_sig) || !identical(as.character(cached_sig[1L]), row$.protocol_signature)) {
      stopf("protocol signature metadata mismatch for %s/%s", country, row$unit)
    }
  }
  if (!is.null(row$.input_sha256) && !is.na(row$.input_sha256) && nzchar(row$.input_sha256)) {
    cached_inputs <- first_nonempty(cache, c("input_hashes_sha256"))
    cached_canonical <- if (!is.null(cached_inputs) && "canonical" %in% names(cached_inputs)) cached_inputs[["canonical"]] else NULL
    if (is.null(cached_canonical) || !identical(tolower(as.character(cached_canonical)), tolower(row$.input_sha256))) {
      stopf("input SHA-256 metadata mismatch for %s/%s", country, row$unit)
    }
  }
  if (!is.null(row$.engine_sha256) && !is.na(row$.engine_sha256) && nzchar(row$.engine_sha256)) {
    cached_engine <- first_nonempty(cache, c("engine_hash_sha256"))
    if (is.null(cached_engine) || !identical(tolower(as.character(cached_engine[1L])), tolower(row$.engine_sha256))) {
      stopf("engine SHA-256 metadata mismatch for %s/%s", country, row$unit)
    }
  }

  contrasts <- engine_fit$contrasts
  contrast_name <- "y_urinary_vs_y_K80"
  if (is.null(contrasts) || is.null(contrasts[[contrast_name]])) stopf("cache lacks contrast %s", contrast_name)
  contrast <- contrasts[[contrast_name]]
  L <- first_nonempty(contrast, c("contrast_vector", "L", "coefficient_contrast_vector"))
  if (is.null(L)) stopf("cache contrast lacks a coefficient contrast vector: %s", path)
  beta <- tryCatch(stats::coef(engine_fit$fit), error = function(e) NULL)
  V <- first_nonempty(engine_fit, c("covariance"))
  V <- if (is.list(V)) V[["hac14"]] else NULL
  score <- engine_fit$score
  if (is.null(beta) || is.null(names(beta)) || is.null(V) || is.null(score)) stopf("cache lacks coefficients, HAC covariance, or score data: %s", path)
  Hinv <- score$Hinv
  U <- score$daily_scores
  dates <- score$daily_dates
  if (is.null(Hinv) || is.null(U) || is.null(dates)) stopf("cache score object lacks Hinv, daily_scores, or daily_dates: %s", path)
  if (!is.matrix(V) || nrow(V) != length(beta) || ncol(V) != length(beta) ||
      !is.matrix(Hinv) || nrow(Hinv) != length(beta) || ncol(Hinv) != length(beta)) {
    stopf("coefficient/covariance/Hinv dimensions disagree in %s", path)
  }
  if (!is.matrix(U) || ncol(U) != length(beta) || nrow(U) != length(dates)) stopf("score/date dimensions disagree in %s", path)
  if (!is.null(rownames(V)) && !is.null(colnames(V))) {
    if (!setequal(rownames(V), names(beta)) || !setequal(colnames(V), names(beta))) stopf("HAC14 covariance names do not match coefficients in %s", path)
    V <- V[names(beta), names(beta), drop = FALSE]
  }
  if (!is.null(rownames(Hinv)) && !is.null(colnames(Hinv))) {
    if (!setequal(rownames(Hinv), names(beta)) || !setequal(colnames(Hinv), names(beta))) stopf("Hinv names do not match coefficients in %s", path)
    Hinv <- Hinv[names(beta), names(beta), drop = FALSE]
  }
  if (!is.null(colnames(U))) {
    if (!setequal(colnames(U), names(beta))) stopf("daily score columns do not match coefficients in %s", path)
    U <- U[, names(beta), drop = FALSE]
  }
  if (!is.null(rownames(U))) {
    if (!setequal(rownames(U), as.character(dates))) stopf("daily score row names do not match score dates in %s", path)
    U <- U[match(as.character(dates), rownames(U)), , drop = FALSE]
  }
  if (max(abs(Hinv - t(Hinv))) > 1e-8 * max(1, max(abs(Hinv)))) stopf("Hinv is not symmetric within tolerance in %s", path)
  if (!is.null(names(L))) {
    if (!setequal(names(L), names(beta))) stopf("contrast and coefficient names differ in %s", path)
    L <- as.numeric(L[names(beta)])
  } else {
    L <- as.numeric(L)
    if (length(L) != length(beta)) stopf("contrast and coefficient vector lengths differ in %s", path)
  }
  names(L) <- names(beta)
  if (any(!is.finite(L)) || any(!is.finite(beta)) || any(!is.finite(V)) || any(!is.finite(Hinv)) || any(!is.finite(U))) {
    stopf("cache has non-finite coefficient, covariance, or score values: %s", path)
  }
  dates <- as_date(dates, "daily score dates")
  if (anyDuplicated(dates)) stopf("daily score dates repeat in %s", path)

  estimate <- sum(L * beta)
  stored_variance <- as.numeric(t(L) %*% V %*% L)
  if (!is.finite(stored_variance) || stored_variance <= 0) stopf("stored HAC14 contrast variance is invalid in %s", path)
  stored_q <- first_nonempty(contrast$coefficient_contrast$by_covariance$hac14, c("estimate"))
  stored_v <- first_nonempty(contrast$coefficient_contrast$by_covariance$hac14, c("variance"))
  if (!is.null(stored_q)) check_close(estimate, as.numeric(stored_q), 1e-8, "cache coefficient contrast estimate")
  if (!is.null(stored_v)) check_close(stored_variance, as.numeric(stored_v), 1e-8, "cache coefficient contrast covariance")

  # Explicit influence orientation is L' Hinv U_t. Hinv must be symmetric,
  # making this identical to the engine's U_t' Hinv L projection.
  influence <- as.numeric(crossprod(L, Hinv %*% t(U)))
  engine_orientation <- as.numeric(U %*% (Hinv %*% L))
  check_close(influence, engine_orientation, 1e-10, "influence orientation identity")
  reconstructed_v <- bartlett_hac_scalar(influence, dates, 14L)
  check_close(reconstructed_v, stored_variance, 1e-8, "unit stored/reconstructed HAC14 variance")
  variance_hac90 <- bartlett_hac_scalar(influence, dates, 90L)
  residual_acf <- residual_acf_rows(engine_fit, country, as.character(row$unit), scenario)
  list(country = country, unit = as.character(row$unit), scenario = scenario,
       estimate = estimate, variance = stored_variance, reconstructed_variance = reconstructed_v,
       variance_abs_difference = abs(stored_variance - reconstructed_v), variance_hac90 = variance_hac90,
       residual_acf = residual_acf, influence = influence,
       dates = dates, n_score_dates = length(dates), score_start = min(dates), score_end = max(dates),
       checksum_method = checksum_method, cache_path = path,
       cache_input_hash = first_nonempty(cache, c("input_md5", "input_hashes_sha256")))
}

dl_weights <- function(yi, vi) {
  if (length(yi) != length(vi) || length(yi) < 2L || any(!is.finite(yi)) || any(!is.finite(vi)) || any(vi <= 0)) {
    stopf("DL pooling needs at least two finite estimates and positive variances")
  }
  w_fixed <- 1 / vi
  y_fixed <- sum(w_fixed * yi) / sum(w_fixed)
  Q <- sum(w_fixed * (yi - y_fixed)^2)
  C <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
  if (!is.finite(C) || C <= 0) stopf("DL denominator is non-positive")
  tau2 <- max(0, (Q - (length(yi) - 1)) / C)
  w_raw <- 1 / (vi + tau2)
  a <- w_raw / sum(w_raw)
  theta <- sum(a * yi)
  dl_var <- 1 / sum(w_raw)
  check_close(sum(a^2 * (vi + tau2)), dl_var, 1e-10,
              "DL normalized-weight variance identity")
  list(tau2 = tau2, raw_weight = w_raw, weight = a, estimate = theta,
       variance = dl_var, se = sqrt(dl_var), Q = Q, C = C)
}

find_pooled_csv <- function(manifest_path, scenario, country) {
  results_dir <- dirname(dirname(manifest_path))
  stem <- scenario
  candidates <- c(file.path(results_dir, "core", paste0(stem, "_", country, "_pooled.csv")),
                  file.path(results_dir, paste0(stem, "_", country, "_pooled.csv")))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[1L] else NULL
}

compare_primary_pool <- function(manifest_path, scenario, country, dl) {
  path <- find_pooled_csv(manifest_path, scenario, country)
  if (is.null(path)) return(list(status = "not_available", path = NA_character_, estimate_difference = NA_real_, se_difference = NA_real_, tau2_difference = NA_real_))
  # The preferred-MX export previously recycled three-element dB vectors into
  # three rows per unit. Never use its pooled CSV unless its per-unit contrast
  # slice has exactly one row for every cache-derived state.
  if (identical(country, "MX") && identical(scenario, "primary_complete_source_2015_2022")) {
    per_unit_path <- sub("_pooled\\.csv$", "_per_unit.csv", path)
    if (!file.exists(per_unit_path)) return(list(status = "updated_MX_per_unit_export_not_available_skip_pooled_check", path = path,
      estimate_difference = NA_real_, se_difference = NA_real_, tau2_difference = NA_real_))
    per_unit <- utils::read.csv(per_unit_path, stringsAsFactors = FALSE)
    if (!all(c("unit", "contrast", "covariance") %in% names(per_unit))) return(list(status = "updated_MX_per_unit_export_schema_invalid_skip_pooled_check", path = path,
      estimate_difference = NA_real_, se_difference = NA_real_, tau2_difference = NA_real_))
    unit_rows <- per_unit[per_unit$contrast == "y_urinary_vs_y_K80" & per_unit$covariance == "hac14", , drop = FALSE]
    if (nrow(unit_rows) != length(dl$weight) || anyDuplicated(as.character(unit_rows$unit))) {
      return(list(status = "updated_MX_per_unit_export_has_duplicate_or_missing_states_skip_pooled_check", path = path,
        estimate_difference = NA_real_, se_difference = NA_real_, tau2_difference = NA_real_))
    }
  }
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  needed <- c("contrast", "covariance", "method", "estimate", "se", "tau2")
  if (!all(needed %in% names(tab))) return(list(status = "pooled_file_schema_unrecognized", path = path, estimate_difference = NA_real_, se_difference = NA_real_, tau2_difference = NA_real_))
  keep <- tab$contrast == "y_urinary_vs_y_K80" & tab$covariance == "hac14" & tab$method == "DL"
  if ("scenario" %in% names(tab)) keep <- keep & tab$scenario == scenario
  z <- tab[keep, , drop = FALSE]
  if (nrow(z) != 1L) return(list(status = "matching_DL_row_not_unique", path = path, estimate_difference = NA_real_, se_difference = NA_real_, tau2_difference = NA_real_))
  diff <- c(estimate = dl$estimate - z$estimate[1L], se = dl$se - z$se[1L], tau2 = dl$tau2 - z$tau2[1L])
  check_close(c(dl$estimate, dl$se, dl$tau2), c(z$estimate[1L], z$se[1L], z$tau2[1L]), 1e-8,
              sprintf("cache-derived DL pool versus existing %s pooled CSV", country))
  list(status = "matched_within_1e-8", path = path, estimate_difference = unname(diff["estimate"]),
       se_difference = unname(diff["se"]), tau2_difference = unname(diff["tau2"]))
}

summarize_residual_acf <- function(x) {
  groups <- split(x, interaction(x$country, x$scenario, x$endpoint, x$lag_calendar_days,
                                 drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    v <- z$pairwise_Pearson_residual_ACF[is.finite(z$pairwise_Pearson_residual_ACF)]
    qs <- if (length(v)) as.numeric(stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE, type = 7)) else rep(NA_real_, 3L)
    data.frame(country = z$country[1L], scenario = z$scenario[1L], endpoint = z$endpoint[1L],
      lag_calendar_days = z$lag_calendar_days[1L], units_total = length(unique(z$unit)),
      units_with_finite_ACF = length(v), median_ACF = qs[2L], Q1_ACF = qs[1L], Q3_ACF = qs[3L],
      IQR_ACF = qs[3L] - qs[1L], median_actual_calendar_pairs = stats::median(z$n_actual_calendar_pairs),
      stringsAsFactors = FALSE)
  }))
}

pool_hac90 <- function(yi, vi, country, scenario) {
  dl <- dl_weights(yi, vi)
  zcrit <- stats::qnorm(0.975)
  rows <- list(data.frame(country = country, scenario = scenario, contrast = "y_urinary_vs_y_K80",
    unit_variance = "HAC90", method = "DL", k = length(yi), estimate_logRR = dl$estimate,
    RR = exp(dl$estimate), se = dl$se, lo_logRR = dl$estimate - zcrit * dl$se,
    hi_logRR = dl$estimate + zcrit * dl$se, lo_RR = exp(dl$estimate - zcrit * dl$se),
    hi_RR = exp(dl$estimate + zcrit * dl$se), tau2 = dl$tau2, failure = "", stringsAsFactors = FALSE))
  if (requireNamespace("metafor", quietly = TRUE)) {
    hk <- tryCatch(metafor::rma.uni(yi = yi, vi = vi, method = "REML", test = "knha"), error = function(e) e)
    if (inherits(hk, "error")) {
      rows[[2L]] <- data.frame(country = country, scenario = scenario, contrast = "y_urinary_vs_y_K80",
        unit_variance = "HAC90", method = "REML_HK", k = length(yi), estimate_logRR = NA_real_,
        RR = NA_real_, se = NA_real_, lo_logRR = NA_real_, hi_logRR = NA_real_, lo_RR = NA_real_, hi_RR = NA_real_,
        tau2 = NA_real_, failure = conditionMessage(hk), stringsAsFactors = FALSE)
    } else rows[[2L]] <- data.frame(country = country, scenario = scenario, contrast = "y_urinary_vs_y_K80",
      unit_variance = "HAC90", method = "REML_HK", k = length(yi), estimate_logRR = as.numeric(hk$b[1L]),
      RR = exp(as.numeric(hk$b[1L])), se = as.numeric(hk$se), lo_logRR = as.numeric(hk$ci.lb),
      hi_logRR = as.numeric(hk$ci.ub), lo_RR = exp(as.numeric(hk$ci.lb)), hi_RR = exp(as.numeric(hk$ci.ub)),
      tau2 = as.numeric(hk$tau2), failure = "", stringsAsFactors = FALSE)
  } else {
    rows[[2L]] <- data.frame(country = country, scenario = scenario, contrast = "y_urinary_vs_y_K80",
      unit_variance = "HAC90", method = "REML_HK", k = length(yi), estimate_logRR = NA_real_,
      RR = NA_real_, se = NA_real_, lo_logRR = NA_real_, hi_logRR = NA_real_, lo_RR = NA_real_, hi_RR = NA_real_,
      tau2 = NA_real_, failure = "metafor unavailable", stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stopf("failed to atomically replace %s", path)
  invisible(path)
}

safe_stem <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

process_country <- function(country, selected, manifest_path, outdir) {
  scenario <- selected$scenario
  rows <- selected$rows
  message(sprintf("%s %s: reading %d validated cache rows", country, scenario, nrow(rows)))
  if (identical(country, "MX") && identical(scenario, "primary_complete_source_2015_2022") &&
      (nrow(rows) != 32L || length(unique(rows$unit)) != 32L)) {
    stopf("preferred Mexico primary must contain exactly one manifest/cache row for each of 32 unique states")
  }
  units <- vector("list", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    units[[i]] <- load_projected_unit(rows[i, , drop = FALSE], country, scenario)
    if (i %% 25L == 0L || i == nrow(rows)) message(sprintf("%s %s: projected %d/%d unit caches", country, scenario, i, nrow(rows)))
    if (i %% 10L == 0L) invisible(gc(FALSE))
  }
  yi <- vapply(units, `[[`, numeric(1), "estimate")
  vi <- vapply(units, `[[`, numeric(1), "variance")
  dl <- dl_weights(yi, vi)
  # Check the requested independent-unit DL identity directly. This is the
  # same standard error reported by the existing DL pooling implementation.
  independent_identity <- sum(dl$weight^2 * (vi + dl$tau2))
  check_close(independent_identity, dl$se^2, 1e-8,
              "independent DL pooled variance identity sum(a_i^2(v_i+tau2)) = DL SE^2")

  start_day <- as.Date(min(vapply(units, function(z) as.numeric(z$score_start), numeric(1))), origin = "1970-01-01")
  end_day <- as.Date(max(vapply(units, function(z) as.numeric(z$score_end), numeric(1))), origin = "1970-01-01")
  calendar <- seq(start_day, end_day, by = "day")
  country_psi <- numeric(length(calendar))
  nonzero_units <- integer(length(calendar))
  per_unit <- vector("list", length(units))
  for (i in seq_along(units)) {
    z <- units[[i]]
    idx <- as.integer(z$dates - calendar[1L]) + 1L
    country_psi[idx] <- country_psi[idx] + dl$weight[i] * z$influence
    nonzero_units[idx] <- nonzero_units[idx] + as.integer(z$influence != 0)
    per_unit[[i]] <- data.frame(country = country, scenario = scenario, unit = z$unit,
      estimate_logRR = z$estimate, variance_hac14_stored = z$variance,
      variance_hac14_reconstructed = z$reconstructed_variance,
      reconstruction_abs_difference = z$variance_abs_difference,
      variance_hac90 = z$variance_hac90,
      tau2_DL = dl$tau2, normalized_DL_weight = dl$weight[i],
      normalized_weight_sq_times_within_between = dl$weight[i]^2 * (vi[i] + dl$tau2),
      n_score_dates = z$n_score_dates, score_start = as.character(z$score_start),
      score_end = as.character(z$score_end), checksum_method = z$checksum_method,
      stringsAsFactors = FALSE)
  }
  per_unit <- do.call(rbind, per_unit)
  point <- dl$estimate
  v_tau <- dl$tau2 * sum(dl$weight^2)
  zcrit <- stats::qnorm(0.975)
  summaries <- list()
  for (bw in c(14L, 28L, 90L)) {
    vs <- bartlett_hac_scalar(country_psi, calendar, bw)
    for (augmented in c(FALSE, TRUE)) {
      vv <- vs + if (augmented) v_tau else 0
      se <- sqrt(vv)
      method <- if (augmented) "calendar_synchronous_plus_DL_random_effect" else "calendar_synchronous_conditional_on_DL_weights_tau"
      summaries[[length(summaries) + 1L]] <- data.frame(country = country, scenario = scenario,
        contrast = "y_urinary_vs_y_K80", base_variance_for_DL_weights = "unit_HAC14",
        synchronized_Bartlett_bandwidth_days = bw, inference = method, k = length(units),
        estimate_logRR = point, RR = exp(point), tau2_DL = dl$tau2,
        sum_normalized_weight_sq = sum(dl$weight^2), V_sync = vs,
        V_tau_component = if (augmented) v_tau else 0, variance = vv, se = se,
        lo_logRR = point - zcrit * se, hi_logRR = point + zcrit * se,
        lo_RR = exp(point - zcrit * se), hi_RR = exp(point + zcrit * se),
        independent_DL_variance = dl$se^2,
        independent_DL_identity_variance = independent_identity,
        stringsAsFactors = FALSE)
    }
  }
  summaries <- do.call(rbind, summaries)
  long_pool <- pool_hac90(yi, vapply(units, `[[`, numeric(1), "variance_hac90"), country, scenario)
  comparison <- compare_primary_pool(manifest_path, scenario, country, dl)
  summaries$primary_DL_point_status <- comparison$status
  summaries$primary_DL_pooled_csv <- comparison$path
  summaries$primary_DL_estimate_difference <- comparison$estimate_difference
  summaries$primary_DL_se_difference <- comparison$se_difference
  summaries$primary_DL_tau2_difference <- comparison$tau2_difference
  summaries$interpretation_limit <- paste(
    "Normal intervals condition on the observed DL unit weights and tau-squared.",
    "Calendar-synchronous HAC captures observed cross-region same-date score dependence,",
    "not unrestricted spatial dependence or all uncertainty in the random-effects structure.")

  stem <- safe_stem(paste(country, scenario, sep = "_"))
  atomic_write_csv(per_unit, file.path(outdir, paste0(stem, "_unit_weights_influence_checks.csv")))
  acf_unit <- do.call(rbind, lapply(units, `[[`, "residual_acf"))
  acf_country <- summarize_residual_acf(acf_unit)
  daily <- data.frame(country = country, scenario = scenario, date = as.character(calendar),
                      weighted_country_influence = country_psi,
                      units_with_nonzero_influence = nonzero_units,
                      stringsAsFactors = FALSE)
  atomic_write_csv(daily, file.path(outdir, paste0(stem, "_country_daily_influence.csv")))
  atomic_write_csv(summaries, file.path(outdir, paste0(stem, "_summary.csv")))
  atomic_write_csv(acf_unit, file.path(outdir, paste0(stem, "_unit_residual_acf.csv")))
  atomic_write_csv(acf_country, file.path(outdir, paste0(stem, "_residual_acf_country_summary.csv")))
  atomic_write_csv(long_pool, file.path(outdir, paste0(stem, "_hac90_DL_REML_HK_sensitivity.csv")))
  list(summary = summaries, per_unit = per_unit, daily = daily, dl = dl,
       checks = list(k = length(units), variance_reconstruction_max_abs_difference = max(per_unit$reconstruction_abs_difference),
         independent_identity_abs_difference = abs(independent_identity - dl$se^2),
         primary_pool_comparison = comparison, cache_checksum_verified = sum(grepl("verified$", per_unit$checksum_method)),
         residual_acf_units = length(unique(acf_unit$unit)), hac90_pool_rows = nrow(long_pool),
         manifest_rows = nrow(rows), scenario = scenario))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  if (opts$help) { show_help(); return(invisible(NULL)) }
  if (opts$self_test) { run_self_tests(); return(invisible(NULL)) }
  if (!file.exists(.engine_path)) stopf("repair_engine.R is missing: %s", .engine_path)
  source(.engine_path, local = FALSE)
  selected <- read_manifest(opts$manifest, opts$country, opts$scenario)
  countries <- names(selected)
  results <- lapply(countries, function(cty) process_country(cty, selected[[cty]], opts$manifest, opts$outdir))
  names(results) <- countries
  all_summary <- do.call(rbind, lapply(results, `[[`, "summary"))
  atomic_write_csv(all_summary, file.path(opts$outdir, "synchronous_meta_summary.csv"))
  manifest_hash <- md5_file(opts$manifest)
  out_manifest <- list(
    run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    script = normalizePath(.script_path, winslash = "/", mustWork = TRUE),
    script_md5 = md5_file(.script_path), repair_engine_md5 = md5_file(.engine_path),
    input_manifest = opts$manifest, input_manifest_md5 = manifest_hash,
    requested_country = opts$country, requested_scenario = opts$scenario,
    contrast = "y_urinary_vs_y_K80", estimate = "DL normalized inverse-(HAC14 variance + DL tau^2) weighted mean of unit logRR contrasts",
    influence = "psi_it = L_i' Hinv_i U_it; Hinv symmetry and influence orientation checked; Bartlett HAC uses actual daily dates with zero-filled gaps, no 1/T and no dispersion multiplier",
    uncertainty = "95% normal intervals conditional on estimated DL weights and tau^2; synchronized variance captures observed calendar-synchronous score dependence, not unrestricted spatial dependence; augmented variance adds tau^2 * sum(a_i^2) under independent unit random effects conditional on tau^2/weights",
    bandwidths = c(14L, 28L, 90L), weight_variance_source = "unit HAC14 (held fixed for synchronized HAC14/28/90 comparisons)",
    self_test_passed = TRUE, country_checks = lapply(results, `[[`, "checks"),
    R_version = R.version.string)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stopf("jsonlite is required to write run_manifest.json")
  jsonlite::write_json(out_manifest, file.path(opts$outdir, "run_manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, na = "null", null = "null", digits = NA)
  message(sprintf("Synchronous meta uncertainty completed; outputs: %s", opts$outdir))
  invisible(results)
}

if (identical(environment(), globalenv())) {
  # Run the synthetic invariants for every production invocation so output is
  # never emitted unless the tested calendar/HAC projection identities hold.
  .cli <- commandArgs(trailingOnly = TRUE)
  if (!any(.cli %in% c("--help", "-h", "--self-test"))) run_self_tests()
  main()
}
