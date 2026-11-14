#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(splines)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Pass <base-root> <temporal-root> <main-output-dir>", call. = FALSE)
}
base_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
temporal_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
main_outdir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)
setDTthreads(1L)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
engine_path <- file.path(temporal_root, "scripts/temporal/repair_engine_temporal.R")
source(engine_path, local = .GlobalEnv)

read_unit_csv <- function(path, select = NULL) {
  if (!file.exists(path)) stop(sprintf("Missing input: %s", path), call. = FALSE)
  header <- names(fread(path, nrows = 0L))
  text_columns <- intersect(c("country", "unit", "scenario", "nuisance", "contrast", "covariance", "status"), header)
  z <- fread(path, select = select, colClasses = list(character = text_columns))
  if ("unit" %in% names(z)) z[, unit := as.character(unit)]
  z
}

main_path <- file.path(main_outdir, "main_per_unit.csv")
main <- read_unit_csv(main_path)
br_old_path <- file.path(base_root, "results/core/primary_threeoutcome_E3_separate_BR_per_unit.csv")
mx_old_path <- file.path(base_root, "results/mx_updated/primary_complete_source_2015_2022_MX_per_unit.csv")
br_old <- read_unit_csv(br_old_path)
mx_old <- read_unit_csv(mx_old_path)

specs <- list(
  list(country = "BR", tag = "br_flexible", nuisance = "flexible", old = br_old,
       old_scenario = "primary_threeoutcome_E3_separate"),
  list(country = "MX", tag = "mx_yearmonth", nuisance = "yearmonth", old = mx_old,
       old_scenario = "primary_complete_source_2015_2022"),
  list(country = "MX", tag = "mx_flexible", nuisance = "flexible", old = mx_old,
       old_scenario = "primary_complete_source_2015_2022")
)
rows <- list()
set_summaries <- list()
for (spec in specs) {
  new <- main[country == spec$country & nuisance == spec$nuisance &
    contrast == "y_urinary_vs_y_K80" & covariance == "hac14" & status == "ok"]
  old <- spec$old[country == spec$country & scenario == spec$old_scenario & nuisance == "separate" &
    contrast == "y_urinary_vs_y_K80" & covariance == "hac14" & is.finite(estimate) &
    is.finite(variance) & variance > 0]
  if (anyDuplicated(new$unit) || anyDuplicated(old$unit)) {
    stop(sprintf("duplicate per-unit rows in common-pool inputs for %s", spec$tag), call. = FALSE)
  }
  common <- sort(intersect(new$unit, old$unit))
  if (length(common) != nrow(new)) stop(sprintf("new success set is not fully covered by old baseline for %s", spec$tag), call. = FALSE)
  set_summaries[[spec$tag]] <- list(country = spec$country, tag = spec$tag,
    new_success_units = uniqueN(new$unit), baseline_finite_units = uniqueN(old$unit),
    common_success_units = length(common), units_only_new = setdiff(new$unit, old$unit),
    units_only_baseline = setdiff(old$unit, common))
  for (role in c("old_baseline", "new_model")) {
    z <- if (role == "old_baseline") old[unit %in% common] else new[unit %in% common]
    for (method in c("DL", "REML_HK")) {
      q <- pool_meta(z$estimate, z$variance, method = method)
      rows[[length(rows) + 1L]] <- data.table(country = spec$country, tag = spec$tag,
        role = role, baseline_scenario = if (role == "old_baseline") spec$old_scenario else NA_character_,
        model_nuisance = spec$nuisance, covariance = "hac14", pooling = method,
        common_success_units = length(common), k = q$k, estimate = q$estimate, se = q$se,
        RR = q$RR, lo = q$RR95CI[1L], hi = q$RR95CI[2L], I2 = q$I2, tau2 = q$tau2,
        concordance = q$concordance,
        failure = if (is.null(q$failure)) NA_character_ else q$failure)
    }
  }
}
pooled <- rbindlist(rows, fill = TRUE)
out_path <- file.path(main_outdir, "main_common_success_pooled.csv")
fwrite(pooled, out_path, na = "")

full_pool_sources <- c(
  BR = file.path(base_root, "results/core/primary_threeoutcome_E3_separate_BR_pooled.csv"),
  MX = file.path(base_root, "results/mx_updated/primary_complete_source_2015_2022_MX_pooled.csv")
)
full_pool_copies <- file.path(main_outdir, c("original_primary_full_pooled_BR.csv", "original_primary_full_pooled_MX.csv"))
if (!all(file.copy(unname(full_pool_sources), full_pool_copies, overwrite = FALSE))) {
  stop("failed to preserve original full primary pooled tables", call. = FALSE)
}
sha256_file <- function(path) {
  ans <- suppressWarnings(system2(Sys.which("sha256sum"), shQuote(path), stdout = TRUE, stderr = TRUE))
  status <- attr(ans, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("sha256sum failed for %s", path), call. = FALSE)
  sub("[[:space:]].*$", "", ans[[1L]])
}
hash_inputs <- c(main_path, br_old_path, mx_old_path, engine_path)
hash_outputs <- c(out_path, full_pool_copies)
audit <- list(created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  method = "Re-pool saved per-unit primary estimates and new M04 estimates on each exact new-success intersection; no unit refits.",
  selected_contrast = "y_urinary_vs_y_K80", covariance = "hac14", pooling = c("DL", "REML_HK"),
  set_summaries = set_summaries,
  original_full_pool_sources = as.list(full_pool_sources),
  input_sha256 = as.list(stats::setNames(vapply(hash_inputs, sha256_file, character(1)), basename(hash_inputs))),
  output_sha256 = as.list(stats::setNames(vapply(hash_outputs, sha256_file, character(1)), basename(hash_outputs)))
)
audit_path <- file.path(main_outdir, "main_common_success_audit.json")
writeLines(as.character(toJSON(audit, auto_unbox = TRUE, null = "null", na = "null",
  digits = NA, pretty = TRUE)), audit_path, useBytes = TRUE)
cat(sprintf("Wrote %d paired common-success pool rows across %d model variants.\n", nrow(pooled), length(specs)))
