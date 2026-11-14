#!/usr/bin/env Rscript

# Merge only per-unit outputs from the four Brazil round-robin batches and
# recompute national pools. Batch-level pooled summaries are never inputs.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

raw <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  i <- match(name, raw)
  if (is.na(i) || i == length(raw)) default else raw[[i + 1L]]
}
root <- normalizePath(option_value("--root"), winslash = "/", mustWork = TRUE)
outdir <- normalizePath(option_value("--outdir", file.path(root, "results/clinical_full_BR")),
                        winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(root, "scripts"))) stop("--root must be the repair directory", call. = FALSE)
merge_script_path <- normalizePath(file.path(root, "scripts/merge_brazil_clinical_batches.R"),
                                   winslash = "/", mustWork = TRUE)
merge_script_sha256_output <- system2("sha256sum", args = shQuote(merge_script_path),
                                      stdout = TRUE, stderr = TRUE)
if (!is.null(attr(merge_script_sha256_output, "status")) &&
    attr(merge_script_sha256_output, "status") != 0L) {
  stop("Unable to calculate the merger-script SHA-256", call. = FALSE)
}
merge_script_sha256 <- sub("^([[:xdigit:]]{64}).*$", "\\1", merge_script_sha256_output[[1L]])
if (!grepl("^[[:xdigit:]]{64}$", merge_script_sha256)) {
  stop("Unexpected sha256sum output for the merger script", call. = FALSE)
}

setDTthreads(1L)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
source(file.path(root, "scripts/repair_engine.R"), local = .GlobalEnv)

batch_count <- 4L
batch_dirs <- file.path(root, "results/clinical_full_BR", paste0("batch_", seq_len(batch_count)))
if (!all(dir.exists(batch_dirs))) stop("All four Brazil batch directories must exist", call. = FALSE)

sensitivity_plan_path <- normalizePath(file.path(root, "admin/clinical_postresult_sensitivity_plan.md"),
                                       winslash = "/", mustWork = TRUE)
sensitivity_plan_sha256_output <- system2("sha256sum", args = shQuote(sensitivity_plan_path),
                                          stdout = TRUE, stderr = TRUE)
if (!is.null(attr(sensitivity_plan_sha256_output, "status")) &&
    attr(sensitivity_plan_sha256_output, "status") != 0L) {
  stop("Unable to calculate the frozen sensitivity-plan SHA-256", call. = FALSE)
}
sensitivity_plan_sha256 <- sub("^([[:xdigit:]]{64}).*$", "\\1", sensitivity_plan_sha256_output[[1L]])
if (!grepl("^[[:xdigit:]]{64}$", sensitivity_plan_sha256)) {
  stop("Unexpected sha256sum output for the frozen sensitivity plan", call. = FALSE)
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing required per-unit output: %s", path), call. = FALSE)
  z <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  setDT(z)
  if ("unit" %in% names(z)) set(z, j = "unit", value = as.character(z[["unit"]]))
  z
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  fwrite(x, temporary, na = "")
  if (!file.rename(temporary, path)) stop(sprintf("Unable to atomically write %s", path), call. = FALSE)
  invisible(path)
}

write_json_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
                                             digits = NA, pretty = TRUE)), temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop(sprintf("Unable to atomically write %s", path), call. = FALSE)
  invisible(path)
}

read_batch_file <- function(filename) {
  lapply(seq_len(batch_count), function(i) {
    x <- read_csv_safe(file.path(batch_dirs[[i]], filename))
    x[, batch_index := as.integer(i)]
    x
  })
}

combine_batches <- function(filename) rbindlist(read_batch_file(filename), fill = TRUE, use.names = TRUE)

manifest_paths <- file.path(batch_dirs, "cache", "manifest.json")
manifests <- lapply(manifest_paths, jsonlite::fromJSON, simplifyVector = FALSE)
if (!all(vapply(manifests, function(x) identical(x$reference_md5, manifests[[1L]]$reference_md5), logical(1)))) {
  stop("Brazil batch manifests do not share the same analysis-reference MD5", call. = FALSE)
}
if (!all(vapply(manifests, function(x) identical(x$source_hashes, manifests[[1L]]$source_hashes), logical(1)))) {
  stop("Brazil batch manifests do not share identical source hashes", call. = FALSE)
}
weight_md5 <- vapply(manifests, function(x) x$weight_signatures$brazil$md5, character(1))
if (length(unique(weight_md5)) != 1L) stop("Brazil batch runtime weight-signature MD5s differ", call. = FALSE)
runner_path <- file.path(root, "scripts/run_clinical_validated.R")
runner_path <- normalizePath(runner_path, winslash = "/", mustWork = TRUE)
if (!identical(unname(tools::md5sum(runner_path)), unname(manifests[[1L]]$source_hashes[[runner_path]]))) {
  stop("Current validated runner MD5 does not match the migrated batch manifests", call. = FALSE)
}

calendar_path <- file.path(root, "data/br_calendar.csv.gz")
calendar_units <- sort(unique(as.character(fread(cmd = paste("gzip -dc --", shQuote(calendar_path)), select = "unit")[["unit"]])))
calendar_units <- calendar_units[!is.na(calendar_units) & nzchar(calendar_units)]
expected_units <- lapply(seq_len(batch_count), function(i) {
  calendar_units[((seq_along(calendar_units) - 1L) %% batch_count) + 1L == i]
})
expected_sizes <- lengths(expected_units)
if (!identical(as.integer(expected_sizes), c(113L, 113L, 112L, 112L))) {
  stop(sprintf("Unexpected four-way Brazil batch sizes: %s", paste(expected_sizes, collapse = ",")), call. = FALSE)
}
unit_to_batch <- rep(seq_len(batch_count), times = expected_sizes)
names(unit_to_batch) <- unlist(expected_units, use.names = FALSE)

check_group_partitions <- function(x, group_keys, table_name) {
  if (!nrow(x)) stop(sprintf("%s is empty", table_name), call. = FALSE)
  missing_batch <- setdiff(unique(x$unit), names(unit_to_batch))
  if (length(missing_batch)) stop(sprintf("%s contains units outside the Brazil calendar: %s", table_name,
                                           paste(missing_batch, collapse = ",")), call. = FALSE)
  if (any(as.integer(unit_to_batch[x$unit]) != x$batch_index)) {
    stop(sprintf("%s contains a unit assigned to the wrong round-robin batch", table_name), call. = FALSE)
  }
  groups <- unique(x[, ..group_keys])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    keep <- rep(TRUE, nrow(x))
    for (key in group_keys) keep <- keep & !is.na(x[[key]]) & x[[key]] == g[[key]][[1L]]
    z <- x[keep]
    if (uniqueN(z$unit) != length(calendar_units) || !setequal(z$unit, calendar_units)) {
      stop(sprintf("%s group %s is missing/duplicating Brazil units", table_name,
                   paste(unlist(g, use.names = FALSE), collapse = "/")), call. = FALSE)
    }
    for (batch in seq_len(batch_count)) {
      observed <- z[batch_index == batch, unique(unit)]
      if (!setequal(observed, expected_units[[batch]])) {
        stop(sprintf("%s group %s does not match batch %d's declared unit partition", table_name,
                     paste(unlist(g, use.names = FALSE), collapse = "/"), batch), call. = FALSE)
      }
    }
  }
  invisible(list(groups = nrow(groups), rows = nrow(x)))
}

assert_unique <- function(x, keys, label) {
  if (nrow(x) && uniqueN(x, by = keys) != nrow(x)) {
    stop(sprintf("Duplicate per-unit keys in %s: %s", label, paste(keys, collapse = ",")), call. = FALSE)
  }
  invisible(TRUE)
}

assert_expected_groups <- function(x, expected_groups, group_keys, table_name) {
  observed <- unique(x[, ..group_keys])
  expected <- unique(expected_groups[, ..group_keys])
  missing <- data.table::fsetdiff(expected, observed)
  extra <- data.table::fsetdiff(observed, expected)
  if (nrow(missing) || nrow(extra)) {
    stop(sprintf("%s planned-group mismatch; missing=%s; extra=%s", table_name,
                 paste(capture.output(print(missing)), collapse = " "),
                 paste(capture.output(print(extra)), collapse = " ")), call. = FALSE)
  }
  invisible(TRUE)
}

check_success_failure_partitions <- function(success, failures, group_keys, unique_keys, expected_groups, table_name) {
  partition_keys <- c(group_keys, "unit", "batch_index")
  combined <- rbindlist(list(
    success[, ..partition_keys],
    failures[, ..partition_keys]
  ), use.names = TRUE, fill = FALSE)
  assert_unique(combined, unique_keys, paste(table_name, "success/failure union"))
  assert_expected_groups(combined, expected_groups, group_keys, table_name)
  check_group_partitions(combined, group_keys, paste(table_name, "success/failure union"))
  combined
}

endpoint <- combine_batches("clinical_brazil_endpoint_per_unit.csv")
endpoint_keys <- c("country", "scenario", "nuisance", "endpoint", "unit")
assert_unique(endpoint, endpoint_keys, "endpoint per-unit table")

interaction <- combine_batches("clinical_brazil_interaction_per_unit.csv")
interaction_keys <- c("country", "scenario", "nuisance", "mode", "unit")
assert_unique(interaction, interaction_keys, "interaction per-unit table")

endpoint_failures <- combine_batches("clinical_brazil_endpoint_failures.csv")
interaction_failures <- combine_batches("clinical_brazil_interaction_failures.csv")
endpoint_all <- rbindlist(list(endpoint, endpoint_failures), use.names = TRUE, fill = TRUE)
interaction_all <- rbindlist(list(interaction, interaction_failures), use.names = TRUE, fill = TRUE)
case_mix <- combine_batches("clinical_brazil_case_mix_standardized_per_unit.csv")
case_mix_strata <- combine_batches("clinical_brazil_case_mix_strata_per_unit.csv")
case_mix_ageband <- combine_batches("clinical_brazil_case_mix_ageband_per_unit.csv")
case_mix_failures <- combine_batches("clinical_brazil_case_mix_failures.csv")
case_mix_availability_batches <- combine_batches("clinical_brazil_case_mix_availability.csv")

assert_unique(case_mix, c("country", "scenario", "unit", "weight_set"), "case-mix standardized per-unit table")
assert_unique(case_mix_strata, c("country", "scenario", "unit", "disease", "stratum"), "case-mix strata table")
assert_unique(case_mix_ageband, c("country", "scenario", "unit", "age_group"), "case-mix age-band table")
assert_unique(case_mix_failures, c("country", "scenario", "unit"), "case-mix failures table")
assert_unique(endpoint_failures, c("country", "scenario", "nuisance", "endpoint", "unit"), "endpoint failures table")
assert_unique(interaction_failures, c("country", "scenario", "nuisance", "mode", "unit"), "interaction failures table")
if (any(endpoint$status != "ok", na.rm = TRUE) || any(endpoint_failures$status == "ok", na.rm = TRUE)) {
  stop("Endpoint success/failure tables contain unexpected status labels", call. = FALSE)
}
if (any(interaction$status != "ok", na.rm = TRUE) || any(interaction_failures$status == "ok", na.rm = TRUE)) {
  stop("Interaction success/failure tables contain unexpected status labels", call. = FALSE)
}

scenario_names <- c("matched_source_2014_2020", "latebilling_reconstruction_2014_2021")
nuisance_names <- c("shared", "separate")
endpoint_names <- c("strict_N20_only_vs_K80", "confirmed_N20_N21_vs_K80", "composition_N23_vs_K80",
                    "N20_vs_K80_0_acute_cholecystitis_plus_K80_3_cholangitis",
                    "combined_original_N20_N23_vs_K80")
expected_endpoint_groups <- CJ(country = "brazil", scenario = scenario_names, nuisance = nuisance_names,
                               endpoint = endpoint_names, unique = TRUE)
expected_interaction_groups <- CJ(country = "brazil", scenario = scenario_names, nuisance = nuisance_names,
                                  mode = c("admission_type", "sex"), unique = TRUE)
expected_case_mix_groups <- CJ(country = "brazil", scenario = scenario_names, unique = TRUE)
endpoint_partition <- check_success_failure_partitions(
  endpoint, endpoint_failures, c("country", "scenario", "nuisance", "endpoint"), endpoint_keys,
  expected_endpoint_groups,
  "endpoint per-unit table"
)
interaction_partition <- check_success_failure_partitions(
  interaction, interaction_failures, c("country", "scenario", "nuisance", "mode"), interaction_keys,
  expected_interaction_groups,
  "interaction per-unit table"
)
case_mix_units <- unique(case_mix[, .(country, scenario, unit, batch_index)])
case_mix_partition <- check_success_failure_partitions(
  case_mix_units, case_mix_failures, c("country", "scenario"), c("country", "scenario", "unit"),
  expected_case_mix_groups,
  "case-mix standardized table"
)

subset_group <- function(x, g, keys) {
  keep <- rep(TRUE, nrow(x))
  for (key in keys) keep <- keep & !is.na(x[[key]]) & x[[key]] == g[[key]][[1L]]
  x[keep]
}

pool_endpoint <- function(x, sensitivity = FALSE) {
  covariance_methods <- c("hac14", "hac28", "day_hc0", "quasi")
  group_keys <- c("country", "scenario", "scenario_label", "nuisance", "endpoint", "endpoint_label")
  if (sensitivity) x <- x[nuisance == "separate"]
  groups <- unique(x[, ..group_keys])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z_all <- subset_group(x, g, c("country", "scenario", "nuisance", "endpoint"))
    event_count_classified <- is.finite(z_all$events_denominator) & is.finite(z_all$events_numerator)
    event_count_pass <- event_count_classified & z_all$events_denominator >= 30 & z_all$events_numerator >= 30
    if (sensitivity) {
      z_count <- z_all[event_count_pass]
      z <- z_count[status == "ok"]
    } else {
      z <- z_all[status == "ok"]
      z_count <- z
    }
    for (method in covariance_methods) {
      estimates <- z[[paste0("logRRR_", method)]]
      variances <- z[[paste0("var_", method)]]
      for (pool_method in c("DL", "REML_HK")) {
        p <- pool_meta(estimates, variances, method = pool_method)
        row <- data.table::as.data.table(as.list(g))
        row[, `:=`(covariance = method, pooling = pool_method,
                   status = if (is.null(p$failure)) "ok" else "partial_or_failed",
                   units_total = nrow(z_all),
                   units_passing_event_count = if (sensitivity) sum(event_count_pass) else NA_integer_,
                   units_below_event_count = if (sensitivity) sum(event_count_classified & !event_count_pass) else NA_integer_,
                   units_event_count_unclassified = sum(!event_count_classified),
                   units_fit = nrow(z), units_fit_failed = sum(z_all$status != "ok"),
                   units_count_eligible_but_fit_failed = if (sensitivity) sum(event_count_pass & z_all$status != "ok") else NA_integer_,
                   k = p$k, logRRR = p$estimate, se = p$se,
                   RRR = p$RR, RR_lo = p$RR95CI[1L], RR_hi = p$RR95CI[2L],
                   I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
                   min_events_per_disease = if (sensitivity) 30L else NA_integer_,
                   failure = if (is.null(p$failure)) NA_character_ else p$failure)]
        rows[[length(rows) + 1L]] <- row
      }
    }
  }
  rbindlist(rows, fill = TRUE, use.names = TRUE)
}

pool_interaction <- function(x) {
  groups <- unique(x[, .(country, scenario, scenario_label, nuisance, mode)])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z_all <- subset_group(x, g, c("country", "scenario", "nuisance", "mode"))
    z <- z_all[status == "ok"]
    for (pool_method in c("DL", "REML_HK")) {
      p <- pool_meta(z$logRRR_difference, z$variance_hac14, method = pool_method)
      row <- data.table::as.data.table(as.list(g))
      row[, `:=`(pooling = pool_method, units_total = nrow(z_all), units_fit = nrow(z),
                 units_fit_failed = sum(z_all$status != "ok"), k = p$k,
                 logRRR_difference = p$estimate, se = p$se, RR_ratio = p$RR,
                 RR_lo = p$RR95CI[1L], RR_hi = p$RR95CI[2L], I2 = p$I2,
                 tau2 = p$tau2, concordance = p$concordance,
                 median_joint_wald_chisq3 = median(z$joint_wald_chisq3, na.rm = TRUE),
                 joint_wald_p_lt_0_05_units = sum(z$joint_wald_p < 0.05, na.rm = TRUE),
                 failure = if (is.null(p$failure)) NA_character_ else p$failure)]
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (!length(rows)) return(data.table(country = character(), scenario = character(), scenario_label = character(),
                                       nuisance = character(), mode = character(), pooling = character(),
                                       units_fit = integer(), k = integer(), logRRR_difference = numeric(),
                                       se = numeric(), RR_ratio = numeric(), RR_lo = numeric(), RR_hi = numeric()))
  rbindlist(rows, fill = TRUE, use.names = TRUE)
}

pool_case_mix <- function(x) {
  groups <- unique(x[, .(scenario, scenario_label, weight_set)])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- subset_group(x, g, c("scenario", "weight_set"))
    failures <- subset_group(case_mix_failures, g, c("scenario"))
    for (pool_method in c("DL", "REML_HK")) {
      p <- pool_meta(z$logRRR, z$variance, method = pool_method)
      row <- data.table::as.data.table(as.list(g))
      row[, `:=`(country = "brazil", pooling = pool_method, units_total = length(calendar_units),
                 units_fit = uniqueN(z$unit), units_fit_failed = nrow(failures), k = p$k,
                 logRRR = p$estimate, se = p$se, RRR = p$RR, RR_lo = p$RR95CI[1L],
                 RR_hi = p$RR95CI[2L], I2 = p$I2, tau2 = p$tau2,
                 failure = if (is.null(p$failure)) NA_character_ else p$failure)]
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (!length(rows)) return(data.table(country = character(), scenario = character(), scenario_label = character(),
                                       weight_set = character(), pooling = character(), units_fit = integer(),
                                       k = integer(), logRRR = numeric(), se = numeric(), RRR = numeric(),
                                       RR_lo = numeric(), RR_hi = numeric()))
  rbindlist(rows, fill = TRUE, use.names = TRUE)
}

pool_ageband <- function(x) {
  groups <- unique(x[, .(scenario, scenario_label, age_group)])
  rows <- list()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- subset_group(x, g, c("scenario", "age_group"))
    for (pool_method in c("DL", "REML_HK")) {
      p <- pool_meta(z$logRRR, z$variance, method = pool_method)
      row <- data.table::as.data.table(as.list(g))
      row[, `:=`(country = "brazil", pooling = pool_method, units_fit = nrow(z), k = p$k,
                 logRRR = p$estimate, se = p$se, RRR = p$RR, RR_lo = p$RR95CI[1L],
                 RR_hi = p$RR95CI[2L], I2 = p$I2, tau2 = p$tau2,
                 failure = if (is.null(p$failure)) NA_character_ else p$failure)]
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (!length(rows)) return(data.table(country = character(), scenario = character(), scenario_label = character(),
                                       age_group = character(), pooling = character(), units_fit = integer(),
                                       k = integer(), logRRR = numeric(), se = numeric(), RRR = numeric(),
                                       RR_lo = numeric(), RR_hi = numeric()))
  rbindlist(rows, fill = TRUE, use.names = TRUE)
}

endpoint_pooled <- pool_endpoint(endpoint_all)
endpoint_min30_eligibility <- endpoint_all[nuisance == "separate", {
  event_count_classified <- is.finite(events_denominator) & is.finite(events_numerator)
  event_count_pass <- event_count_classified & events_denominator >= 30 & events_numerator >= 30
  event_count_below <- event_count_classified & !event_count_pass
  .(
    units_total = .N,
    units_event_count_classified = sum(event_count_classified),
    units_passing_event_count = sum(event_count_pass),
    units_excluded_event_count = sum(event_count_below),
    units_event_count_unclassified = sum(!event_count_classified),
    units_fit_passing_event_count = sum(event_count_pass & status == "ok"),
    units_count_eligible_but_fit_failed = sum(event_count_pass & status != "ok"),
    units_fit_failed_total = sum(status != "ok")
  )
}, by = .(country, scenario, scenario_label, endpoint, endpoint_label)]
endpoint_min30 <- endpoint_all[nuisance == "separate" & status == "ok" &
                             is.finite(events_denominator) & is.finite(events_numerator) &
                             events_denominator >= 30 & events_numerator >= 30]
endpoint_pooled_min30 <- pool_endpoint(endpoint_all, sensitivity = TRUE)
interaction_pooled <- pool_interaction(interaction_all)
case_mix_pooled <- pool_case_mix(case_mix)
case_mix_ageband_pooled <- pool_ageband(case_mix_ageband)

availability <- case_mix_availability_batches[, .(
  units_in_calendar_window = sum(units_in_calendar_window),
  units_fit_successfully = sum(units_fit_successfully),
  units_failed_or_unavailable = sum(units_failed_or_unavailable)
), by = .(country, scenario, scenario_label)]
availability[, failure_reasons := vapply(seq_len(.N), function(i) {
  g <- as.data.table(as.list(.BY))
  fails <- subset_group(case_mix_failures, g, c("country", "scenario"))
  if (!nrow(fails)) return("[]")
  counts <- fails[, .N, by = "failure"]
  values <- stats::setNames(as.list(counts$N), counts$failure)
  as.character(jsonlite::toJSON(values, auto_unbox = TRUE, null = "null", na = "null"))
}, character(1)), by = .(country, scenario, scenario_label)]

if (!identical(as.integer(availability$units_in_calendar_window), rep.int(length(calendar_units), nrow(availability)))) {
  stop("Merged Brazil case-mix availability does not cover the full calendar units", call. = FALSE)
}

outputs <- list(
  clinical_brazil_endpoint_per_unit.csv = endpoint,
  clinical_brazil_endpoint_failures.csv = endpoint_failures,
  clinical_brazil_endpoint_pooled.csv = endpoint_pooled,
  clinical_brazil_endpoint_per_unit_min30_events.csv = endpoint_min30,
  clinical_brazil_endpoint_min30_events_eligibility.csv = endpoint_min30_eligibility,
  clinical_brazil_endpoint_pooled_min30_events.csv = endpoint_pooled_min30,
  clinical_brazil_interaction_per_unit.csv = interaction,
  clinical_brazil_interaction_failures.csv = interaction_failures,
  clinical_brazil_interaction_coverage.csv = interaction_pooled,
  clinical_brazil_case_mix_standardized_per_unit.csv = case_mix,
  clinical_brazil_case_mix_strata_per_unit.csv = case_mix_strata,
  clinical_brazil_case_mix_ageband_per_unit.csv = case_mix_ageband,
  clinical_brazil_case_mix_failures.csv = case_mix_failures,
  clinical_brazil_case_mix_availability.csv = availability,
  clinical_brazil_case_mix_pooled.csv = case_mix_pooled,
  clinical_brazil_case_mix_ageband_pooled.csv = case_mix_ageband_pooled
)
for (nm in names(outputs)) atomic_fwrite(outputs[[nm]], file.path(outdir, nm))

weights_tables <- lapply(batch_dirs, function(d) read_csv_safe(file.path(d, "clinical_brazil_standardization_weights.csv")))
eligibility_tables <- lapply(batch_dirs, function(d) read_csv_safe(file.path(d, "clinical_brazil_case_mix_eligibility.csv")))
if (!all(vapply(weights_tables[-1L], function(x) identical(x, weights_tables[[1L]]), logical(1))) ||
    !all(vapply(eligibility_tables[-1L], function(x) identical(x, eligibility_tables[[1L]]), logical(1)))) {
  stop("Brazil full-cohort standardization weights/eligibility differ across batches", call. = FALSE)
}
atomic_fwrite(weights_tables[[1L]], file.path(outdir, "clinical_brazil_standardization_weights.csv"))
atomic_fwrite(eligibility_tables[[1L]], file.path(outdir, "clinical_brazil_case_mix_eligibility.csv"))

output_paths <- c(file.path(outdir, names(outputs)),
                  file.path(outdir, c("clinical_brazil_standardization_weights.csv",
                                      "clinical_brazil_case_mix_eligibility.csv")))
output_md5 <- tools::md5sum(output_paths)
failure_counts <- list(
  endpoint = endpoint_failures[, .N, by = .(scenario, nuisance, endpoint, failure)],
  interaction = interaction_failures[, .N, by = .(scenario, nuisance, mode, failure)],
  case_mix = case_mix_failures[, .N, by = .(scenario, failure)]
)
min30_summary <- endpoint_min30_eligibility
audit <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  merge_script_sha256 = merge_script_sha256,
  runner_md5 = unname(tools::md5sum(runner_path)),
  source_hashes = manifests[[1L]]$source_hashes,
  reference_md5 = manifests[[1L]]$reference_md5,
  configuration_md5_by_batch = stats::setNames(lapply(manifests, `[[`, "config_md5"), paste0("batch_", seq_len(batch_count))),
  cache_signature_md5_by_batch = stats::setNames(lapply(manifests, `[[`, "signature_md5"), paste0("batch_", seq_len(batch_count))),
  runtime_weight_signature_md5 = weight_md5[1L],
  round_robin_unit_batch_sizes = expected_sizes,
  round_robin_partition = lapply(expected_units, as.list),
  table_checks = list(
    endpoint = list(keys = endpoint_keys, groups_checked = uniqueN(endpoint_partition, by = c("country", "scenario", "nuisance", "endpoint")),
    successful_rows = nrow(endpoint), failure_rows = nrow(endpoint_failures),
                    fit_failure_rows = nrow(endpoint_failures),
                    success_failure_union_rows = nrow(endpoint_partition), unique_keys = TRUE,
                    success_failure_disjoint = TRUE, all_planned_groups_present = TRUE),
    interaction = list(keys = interaction_keys, groups_checked = uniqueN(interaction_partition, by = c("country", "scenario", "nuisance", "mode")),
                       successful_rows = nrow(interaction), failure_rows = nrow(interaction_failures),
                       fit_failure_rows = nrow(interaction_failures),
                       success_failure_union_rows = nrow(interaction_partition), unique_keys = TRUE,
                       success_failure_disjoint = TRUE, all_planned_groups_present = TRUE),
    case_mix_success_units = nrow(case_mix_units), case_mix_failures = nrow(case_mix_failures),
    case_mix_success_failure_union_units = nrow(case_mix_partition),
    case_mix_standardized_rows = nrow(case_mix), case_mix_stratum_rows = nrow(case_mix_strata),
    case_mix_ageband_rows = nrow(case_mix_ageband), full_calendar_units = length(calendar_units)
  ),
  failure_counts = failure_counts,
  min30_events_sensitivity = list(
    definition = "Separate nuisance model only; count eligibility uses all success and failure rows with observed K80 denominator and endpoint N-code numerator each >=30. Missing event counts are unclassified. Pools use only status=ok estimates; no refitting or effect/p-value selection.",
    frozen_plan_path = sensitivity_plan_path,
    frozen_plan_sha256 = sensitivity_plan_sha256,
    eligibility_table = "clinical_brazil_endpoint_min30_events_eligibility.csv",
    eligible_per_unit_table = "clinical_brazil_endpoint_per_unit_min30_events.csv",
    counts_by_scenario_endpoint = min30_summary,
    sensitivity_pool_rows = nrow(endpoint_pooled_min30),
    full_cohort_pool_preserved = TRUE
  ),
  model_scope = list(
    endpoint_scenarios = sort(unique(endpoint$scenario)),
    endpoint_nuisance = sort(unique(endpoint$nuisance)),
    endpoints = sort(unique(endpoint$endpoint)),
    endpoint_covariances = c("hac14", "hac28", "day_hc0", "quasi"),
    endpoint_pooling = c("DL", "REML_HK"),
    interaction_modes = sort(unique(interaction$mode)),
    interaction_pooling = c("DL", "REML_HK"),
    case_mix_weight_sets = sort(unique(case_mix$weight_set)),
    case_mix_pooling = c("DL", "REML_HK")
  ),
  output_md5 = stats::setNames(as.list(unname(output_md5)), names(output_md5)),
  note = "National pools were recomputed exclusively from merged per-unit estimates; batch-level pooled summaries were not used."
)
write_json_atomic(audit, file.path(outdir, "clinical_brazil_batch_merge_audit.json"))
cat(sprintf("Merged %d Brazil calendar units from four batches. Endpoint rows=%d, interaction rows=%d, case-mix successes=%d, min30 endpoint pool rows=%d.\n",
            length(calendar_units), nrow(endpoint), nrow(interaction), nrow(case_mix), nrow(endpoint_pooled_min30)))
