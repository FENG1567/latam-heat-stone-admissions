#!/usr/bin/env Rscript

# Extract within-admission/sex RRRs from existing joint fits and re-pool them.
# No unit models are refitted. The formal interaction remains the separately
# pooled per-unit high-minus-low contrast in clinical_brazil_interaction_coverage.csv.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

raw <- commandArgs(trailingOnly = TRUE)
root_i <- match("--root", raw)
if (is.na(root_i) || root_i == length(raw)) stop("Pass --root <repair-directory>", call. = FALSE)
root <- normalizePath(raw[[root_i + 1L]], winslash = "/", mustWork = TRUE)
setDTthreads(1L)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
source(file.path(root, "scripts/repair_engine.R"), local = .GlobalEnv)

outdir <- file.path(root, "results/clinical_full_BR")
per_unit_path <- file.path(outdir, "clinical_brazil_interaction_per_unit.csv")
failures_path <- file.path(outdir, "clinical_brazil_interaction_failures.csv")
per_unit_output <- file.path(outdir, "clinical_brazil_interaction_group_per_unit.csv")
pooled_output <- file.path(outdir, "clinical_brazil_interaction_group_pooled.csv")
audit_output <- file.path(outdir, "clinical_brazil_interaction_group_pool_audit.json")

read_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing input: %s", path), call. = FALSE)
  z <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  setDT(z)
  if ("unit" %in% names(z)) set(z, j = "unit", value = as.character(z[["unit"]]))
  z
}

atomic_fwrite <- function(x, path) {
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  fwrite(x, temporary, na = "")
  if (!file.rename(temporary, path)) stop(sprintf("Unable to atomically write %s", path), call. = FALSE)
  invisible(path)
}

group_order <- list(admission_type = c("elective", "urgent"), sex = c("female", "male"))
group_contrast <- function(mode, group) {
  if (mode == "admission_type" && group == "elective") return("N20-only versus K80 among elective admissions (CAR_INT=01)")
  if (mode == "admission_type" && group == "urgent") return("N20-only versus K80 among urgent admissions (CAR_INT=02)")
  if (mode == "sex" && group == "female") return("N20-only versus K80 among female admissions")
  if (mode == "sex" && group == "male") return("N20-only versus K80 among male admissions")
  stop(sprintf("Unexpected interaction group %s / %s", mode, group), call. = FALSE)
}

read_group_estimates <- function(txt, mode) {
  z <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  wanted <- group_order[[mode]]
  if (is.null(wanted) || length(z) != length(wanted)) {
    stop(sprintf("Unexpected group count for interaction mode %s", mode), call. = FALSE)
  }
  out <- vector("list", length(z))
  names(out) <- wanted
  value_names <- c("group", "logRRR", "variance_hac14", "se_hac14", "RRR_hac14", "RR_lo_hac14", "RR_hi_hac14")
  for (i in seq_along(z)) {
    values <- as.character(unlist(z[[i]], use.names = FALSE))
    if (length(values) != length(value_names)) stop("Unexpected group-effects JSON field count", call. = FALSE)
    names(values) <- value_names
    group <- values[["group"]]
    if (!group %in% wanted || !is.null(out[[group]])) stop("Unexpected/duplicate group name in group-effects JSON", call. = FALSE)
    numeric_values <- as.list(as.numeric(values[-1L]))
    names(numeric_values) <- value_names[-1L]
    out[[group]] <- list(group = group, values = numeric_values)
  }
  if (!all(wanted %in% names(out)) || any(vapply(out[wanted], is.null, logical(1)))) {
    stop(sprintf("Group-effects JSON is missing a required group for %s", mode), call. = FALSE)
  }
  out[wanted]
}

per_unit <- read_csv(per_unit_path)
failures <- read_csv(failures_path)
if (any(per_unit$status != "ok", na.rm = TRUE) || any(failures$status == "ok", na.rm = TRUE)) {
  stop("Interaction success/failure tables have unexpected status labels", call. = FALSE)
}

unit_rows <- list()
for (i in seq_len(nrow(per_unit))) {
  source_row <- per_unit[i]
  mode <- source_row[["mode"]]
  group_rows <- read_group_estimates(source_row[["group_effects_json"]], mode)
  for (j in seq_along(group_rows)) {
    group_name <- names(group_rows)[[j]]
    values <- group_rows[[j]][["values"]]
    unit_rows[[length(unit_rows) + 1L]] <- data.table(
      country = source_row[["country"]], scenario = source_row[["scenario"]],
      scenario_label = source_row[["scenario_label"]], nuisance = source_row[["nuisance"]],
      mode = mode, group = group_name, group_order = as.integer(j),
      contrast = group_contrast(mode, group_name), unit = source_row[["unit"]],
      joint_fit_status = "ok", status = "ok",
      n_days = source_row[["n_days"]], events_K80_low = source_row[["events_K80_low"]],
      events_N20_low = source_row[["events_N20_low"]], events_K80_high = source_row[["events_K80_high"]],
      events_N20_high = source_row[["events_N20_high"]], p95 = source_row[["p95"]],
      median = source_row[["median"]], exposure_days = source_row[["exposure_days"]],
      logRRR = values[["logRRR"]], variance_hac14 = values[["variance_hac14"]],
      se_hac14 = values[["se_hac14"]], RRR_hac14 = values[["RRR_hac14"]],
      RR_lo_hac14 = values[["RR_lo_hac14"]], RR_hi_hac14 = values[["RR_hi_hac14"]],
      joint_logRRR_difference = source_row[["logRRR_difference"]],
      joint_RR_ratio_hac14 = source_row[["RR_ratio_hac14"]],
      joint_RR_ratio_lo = source_row[["RR_ratio_lo"]], joint_RR_ratio_hi = source_row[["RR_ratio_hi"]],
      joint_wald_chisq3 = source_row[["joint_wald_chisq3"]], joint_wald_p = source_row[["joint_wald_p"]],
      failure = NA_character_, batch_index = source_row[["batch_index"]]
    )
  }
}

for (i in seq_len(nrow(failures))) {
  source_row <- failures[i]
  mode <- source_row[["mode"]]
  wanted <- group_order[[mode]]
  if (is.null(wanted)) stop(sprintf("Unexpected interaction mode in failures: %s", mode), call. = FALSE)
  for (j in seq_along(wanted)) {
    group_name <- wanted[[j]]
    unit_rows[[length(unit_rows) + 1L]] <- data.table(
      country = source_row[["country"]], scenario = source_row[["scenario"]],
      scenario_label = source_row[["scenario_label"]], nuisance = source_row[["nuisance"]],
      mode = mode, group = group_name, group_order = as.integer(j),
      contrast = group_contrast(mode, group_name), unit = source_row[["unit"]],
      joint_fit_status = "failed", status = "joint_fit_failed",
      n_days = source_row[["n_days"]], events_K80_low = NA_real_, events_N20_low = NA_real_,
      events_K80_high = NA_real_, events_N20_high = NA_real_, p95 = NA_real_, median = NA_real_,
      exposure_days = NA_real_, logRRR = NA_real_, variance_hac14 = NA_real_, se_hac14 = NA_real_,
      RRR_hac14 = NA_real_, RR_lo_hac14 = NA_real_, RR_hi_hac14 = NA_real_,
      joint_logRRR_difference = NA_real_, joint_RR_ratio_hac14 = NA_real_,
      joint_RR_ratio_lo = NA_real_, joint_RR_ratio_hi = NA_real_, joint_wald_chisq3 = NA_real_,
      joint_wald_p = NA_real_, failure = source_row[["failure"]], batch_index = source_row[["batch_index"]]
    )
  }
}

group_per_unit <- rbindlist(unit_rows, use.names = TRUE, fill = TRUE)
group_keys <- c("country", "scenario", "nuisance", "mode", "group", "unit")
if (uniqueN(group_per_unit, by = group_keys) != nrow(group_per_unit)) {
  stop("Duplicate interaction-group per-unit rows", call. = FALSE)
}
coverage <- group_per_unit[, .(
  units_total = uniqueN(unit), units_joint_fit_successful = uniqueN(unit[joint_fit_status == "ok"]),
  units_joint_fit_failed = uniqueN(unit[joint_fit_status != "ok"]),
  group_rows = .N, successful_group_rows = sum(status == "ok"), failed_group_rows = sum(status != "ok")
), by = .(country, scenario, nuisance, mode, group)]
if (any(coverage$units_total != 450L) || any(coverage$group_rows != 450L)) {
  stop("Interaction-group output does not preserve all 450 units per group", call. = FALSE)
}

pool_groups <- unique(group_per_unit[, .(country, scenario, scenario_label, nuisance, mode, group, contrast)])
pooled_rows <- list()
for (i in seq_len(nrow(pool_groups))) {
  g <- pool_groups[i]
  z <- group_per_unit[country == g$country & scenario == g$scenario & nuisance == g$nuisance &
                        mode == g$mode & group == g$group]
  successful <- z[status == "ok"]
  n_fit <- uniqueN(successful$unit)
  n_failed <- uniqueN(z[status != "ok", unit])
  for (pool_method in c("DL", "REML_HK")) {
    p <- pool_meta(successful$logRRR, successful$variance_hac14, method = pool_method)
    row <- as.data.table(as.list(g))
    row[, `:=`(covariance = "hac14", pooling = pool_method,
               units_total = uniqueN(z$unit), units_joint_fit_successful = n_fit,
               units_joint_fit_failed = n_failed, units_fit = n_fit, k = p$k,
               logRRR = p$estimate, se = p$se, RRR = p$RR,
               RR_lo = p$RR95CI[1L], RR_hi = p$RR95CI[2L], I2 = p$I2, tau2 = p$tau2,
               concordance = p$concordance,
               status = if (is.null(p$failure)) "ok" else "partial_or_failed",
               failure = if (is.null(p$failure)) NA_character_ else p$failure)]
    pooled_rows[[length(pooled_rows) + 1L]] <- row
  }
}
group_pooled <- rbindlist(pooled_rows, fill = TRUE, use.names = TRUE)
if (nrow(group_pooled) != 32L) stop(sprintf("Expected 32 pooled group rows, observed %d", nrow(group_pooled)), call. = FALSE)

atomic_fwrite(group_per_unit, per_unit_output)
atomic_fwrite(group_pooled, pooled_output)
output_hashes <- tools::md5sum(c(per_unit_output, pooled_output))
audit <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  method = "Existing joint-fit group_effects_json extracted and pooled by country/unit; no unit model refits.",
  source_per_unit_md5 = unname(tools::md5sum(per_unit_path)),
  source_failure_md5 = unname(tools::md5sum(failures_path)),
  interaction_mode_groups = list(admission_type = c("elective", "urgent"), sex = c("female", "male")),
  group_effect_json_order = c("group", "estimate", "variance", "se", "RRR", "lo", "hi"),
  per_unit_rows = nrow(group_per_unit),
  successful_joint_unit_rows = nrow(per_unit),
  failed_joint_unit_rows = nrow(failures),
  group_coverage = coverage,
  pooled_rows = nrow(group_pooled),
  pooling = c("DL", "REML_HK"),
  preferred_selector = list(
    country = "brazil", scenario = "latebilling_reconstruction_2014_2021", nuisance = "separate",
    mode = "admission_type", groups = c("elective", "urgent"), covariance = "hac14", pooling = "DL",
    selected_unit_set = "units with status=ok in the four-endpoint joint admission_type model",
    selected_unit_count = uniqueN(per_unit[scenario == "latebilling_reconstruction_2014_2021" &
                                             nuisance == "separate" & mode == "admission_type", unit])
  ),
  formal_interaction_note = "The urgent-minus-elective formal interaction remains the separately pooled per-unit logRRR_difference from clinical_brazil_interaction_coverage.csv. It is not derived from these pooled elective and urgent RRRs.",
  output_md5 = stats::setNames(as.list(unname(output_hashes)), names(output_hashes))
)
temporary_audit <- tempfile(pattern = ".clinical_brazil_interaction_group_pool_audit.tmp-", tmpdir = dirname(audit_output))
on.exit(unlink(temporary_audit), add = TRUE)
writeLines(as.character(jsonlite::toJSON(audit, auto_unbox = TRUE, null = "null", na = "null",
                                          digits = NA, pretty = TRUE)), temporary_audit, useBytes = TRUE)
if (!file.rename(temporary_audit, audit_output)) stop("Unable to atomically write interaction group audit", call. = FALSE)
cat(sprintf("Extracted %d group-specific unit estimates and pooled %d country rows from %d joint-fit success/failure units.\n",
            nrow(group_per_unit), nrow(group_pooled), nrow(per_unit) + nrow(failures)))
