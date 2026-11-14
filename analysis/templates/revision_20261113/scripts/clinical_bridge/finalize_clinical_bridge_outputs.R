#!/usr/bin/env Rscript

# Finalize audit and freeze-candidate artifacts from completed bridge CSVs.
# This script never loads inputs for fitting and never refits a model.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
raw <- commandArgs(trailingOnly = TRUE)
if (length(raw) != 6L) {
  stop("Usage: finalize_clinical_bridge_outputs.R <base-root> <output-dir> <production-log> <exit-code-file> <source-script> <pid>", call. = FALSE)
}
base_root <- normalizePath(raw[[1L]], winslash = "/", mustWork = TRUE)
outdir <- normalizePath(raw[[2L]], winslash = "/", mustWork = TRUE)
production_log <- normalizePath(raw[[3L]], winslash = "/", mustWork = TRUE)
exit_file <- normalizePath(raw[[4L]], winslash = "/", mustWork = TRUE)
source_script <- normalizePath(raw[[5L]], winslash = "/", mustWork = TRUE)
production_pid <- raw[[6L]]

write_csv_once <- function(x, path) {
  if (file.exists(path)) {
    tmp <- tempfile(pattern = "clinical-bridge-finalizer-check-", tmpdir = dirname(path), fileext = ".csv")
    on.exit(unlink(tmp), add = TRUE)
    fwrite(x, tmp, na = "")
    if (!identical(unname(tools::md5sum(tmp)), unname(tools::md5sum(path)))) {
      stop(sprintf("refusing to overwrite existing nonidentical audit artifact: %s", path), call. = FALSE)
    }
    return(invisible(path))
  }
  fwrite(x, path, na = "")
  invisible(path)
}
write_json_once <- function(x, path) {
  if (file.exists(path)) stop(sprintf("refusing to overwrite existing JSON artifact: %s", path), call. = FALSE)
  write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
  invisible(path)
}

required_files <- c(
  "clinical_bridge_unit_membership.csv", "clinical_bridge_matched_unit_fit_audit.csv",
  "clinical_bridge_matched_k_audit.csv", "clinical_bridge_standardization_weights.csv",
  "clinical_bridge_preferred_weight_sensitivity_weights.csv", "clinical_bridge_targets_audit.csv",
  "clinical_bridge_case_eligibility.csv", "clinical_bridge_case_composition.csv",
  "clinical_bridge_adult_unstandardized_pooled.csv", "clinical_bridge_standardized_pooled.csv",
  "clinical_bridge_pool_summary.csv", "mexico_strict_main_RRR_unit_membership.csv",
  "mexico_strict_phenotypes_per_unit.csv", "mexico_strict_phenotypes_failures.csv",
  "mexico_strict_phenotypes_pooled.csv", "mexico_strict_eligibility_summary.csv",
  "mexico_strict_broad_same_unit_reference_check.csv"
)
missing <- required_files[!file.exists(file.path(outdir, required_files))]
if (length(missing)) stop(sprintf("saved production outputs are incomplete: %s", paste(missing, collapse = ", ")), call. = FALSE)

membership <- fread(file.path(outdir, "clinical_bridge_unit_membership.csv"), colClasses = list(character = "unit"))
matched <- fread(file.path(outdir, "clinical_bridge_matched_unit_fit_audit.csv"), colClasses = list(character = "unit"))
k_audit <- fread(file.path(outdir, "clinical_bridge_matched_k_audit.csv"))
weights <- fread(file.path(outdir, "clinical_bridge_standardization_weights.csv"))
preferred_weights <- fread(file.path(outdir, "clinical_bridge_preferred_weight_sensitivity_weights.csv"))
target_audit <- fread(file.path(outdir, "clinical_bridge_targets_audit.csv"), colClasses = list(character = "unit"))
main_membership <- fread(file.path(outdir, "mexico_strict_main_RRR_unit_membership.csv"), colClasses = list(character = "unit"))
strict_per_unit <- fread(file.path(outdir, "mexico_strict_phenotypes_per_unit.csv"), colClasses = list(character = "unit"))
strict_eligibility <- fread(file.path(outdir, "mexico_strict_eligibility_summary.csv"))
reproduction_summary <- fread(file.path(outdir, "reproduction_check_attempt2/mexico_strict_reproduction_summary.csv"))
pool_path <- file.path(outdir, "clinical_bridge_pool_summary.csv")
pool <- fread(pool_path)
pool[, upstream_source_row_ordinal := .I]

membership_summary <- membership[, .(
  clinical_units = .N,
  standardization_success_units = sum(standardization_status == "success", na.rm = TRUE),
  standardization_failed_units = sum(standardization_status == "failed", na.rm = TRUE)
), by = .(country, scenario)]
fit_summary <- matched[, .(
  adult_fit_attempted = sum(adult_fit_attempted, na.rm = TRUE),
  adult_fit_successes = sum(adult_status == "ok", na.rm = TRUE),
  adult_hac14_pool_eligible = sum(hac14_pool_eligible, na.rm = TRUE),
  common_bridge_success_units = sum(common_bridge_success, na.rm = TRUE),
  common_bridge_failed_or_ineligible = sum(!common_bridge_success, na.rm = TRUE)
), by = .(country, scenario)]
weight_summary <- weights[, .(
  primary_weight_strata = .N,
  primary_weight_sum = sum(weight),
  primary_case_count_total = sum(pooled_K80_N20_N23_cases),
  primary_weights_reconstruct_from_counts = all(abs(weight - pooled_K80_N20_N23_cases / sum(pooled_K80_N20_N23_cases)) < 1e-12),
  primary_weight_source = paste(unique(weight_source_scenario), collapse = ";")
), by = .(country)]
preferred_summary <- preferred_weights[, .(
  preferred_weight_strata = .N,
  preferred_weight_sum = sum(weight),
  preferred_case_count_total = sum(pooled_K80_N20_N23_cases),
  preferred_weights_reconstruct_from_counts = all(abs(weight - pooled_K80_N20_N23_cases / sum(pooled_K80_N20_N23_cases)) < 1e-12),
  preferred_weight_source = paste(unique(weight_source_scenario), collapse = ";")
), by = .(country)]
m06_audit <- Reduce(function(x, y) merge(x, y, by = c("country", "scenario"), all = TRUE, sort = FALSE),
                    list(membership_summary, fit_summary, k_audit))
m06_audit <- merge(m06_audit, weight_summary, by = "country", all.x = TRUE, sort = FALSE)
m06_audit <- merge(m06_audit, preferred_summary, by = "country", all.x = TRUE, sort = FALSE)
m06_audit[, common_success_set_matches_standardization_success := {
  success <- membership[country == .BY$country & scenario == .BY$scenario & standardization_status == "success", unit]
  common <- matched[country == .BY$country & scenario == .BY$scenario & common_bridge_success == TRUE, unit]
  setequal(success, common)
}, by = .(country, scenario)]
m06_audit[, all_adult_successes_hac14_eligible := adult_fit_successes == adult_hac14_pool_eligible]
m06_audit[, all_success_units_used_in_primary_pool := primary_pool_k_DL == common_bridge_success_units]
m06_audit[, primary_weights_sum_to_one := abs(primary_weight_sum - 1) < 1e-12]
m06_audit[, preferred_weights_sum_to_one := is.na(preferred_weight_sum) | abs(preferred_weight_sum - 1) < 1e-12]

expected_success <- c(brazil = 208L, mexico = 16L)
for (country_name in names(expected_success)) {
  row <- m06_audit[country == country_name]
  if (nrow(row) != 1L || row$standardization_success_units != expected_success[[country_name]] ||
      row$adult_fit_attempted != expected_success[[country_name]] ||
      row$adult_fit_successes != expected_success[[country_name]] ||
      row$adult_hac14_pool_eligible != expected_success[[country_name]] ||
      row$common_bridge_success_units != expected_success[[country_name]] ||
      !row$common_success_set_matches_standardization_success ||
      !row$exact_same_set_pooling || !row$all_adult_successes_hac14_eligible ||
      !row$all_success_units_used_in_primary_pool || !row$primary_weights_reconstruct_from_counts ||
      !row$primary_weights_sum_to_one) {
    stop(sprintf("M06 same-sample or primary-weight audit failed for %s", country_name), call. = FALSE)
  }
}
if (m06_audit[country == "brazil", standardization_failed_units] != 242L ||
    m06_audit[country == "mexico", standardization_failed_units] != 16L ||
    m06_audit[country == "brazil", primary_weight_source] != "matched_source_2014_2020" ||
    m06_audit[country == "mexico", primary_weight_source] != "mx_complete_source_2015_2022" ||
    m06_audit[country == "brazil", preferred_weight_source] != "latebilling_reconstruction_2014_2021" ||
    !m06_audit[country == "brazil", preferred_weights_reconstruct_from_counts] ||
    !m06_audit[country == "brazil", preferred_weights_sum_to_one]) {
  stop("M06 failure counts or weight source/normalization checks failed", call. = FALSE)
}

m06_units <- merge(
  membership,
  matched[, .(country, scenario, unit, adult_status, adult_failure, hac14_pool_eligible,
              adult_fit_attempted, common_bridge_success)],
  by = c("country", "scenario", "unit"), all.x = TRUE, sort = FALSE
)
m06_units[, analysis_origin := "new M06 clinical bridge production"]
write_csv_once(m06_audit, file.path(outdir, "clinical_bridge_m06_set_weight_audit.csv"))
write_csv_once(m06_units, file.path(outdir, "clinical_bridge_m06_same_sample_membership.csv"))

if (nrow(main_membership) != 32L || any(!main_membership$primary_RRR_and_clinical_unit_match)) {
  stop("M07 main RRR/clinical target membership does not exactly match all 32 states", call. = FALSE)
}
coverage <- strict_per_unit[, .(
  clinical_units = uniqueN(unit),
  fits_attempted = .N,
  fits_successful = sum(status == "ok"),
  fits_failed = sum(status != "ok"),
  primary_RRR_units = uniqueN(unit[primary_RRR_unit_included == TRUE])
), by = .(endpoint)]
coverage <- merge(
  coverage,
  strict_eligibility[, .(endpoint, primary_minimum_days = primary_RRR_minimum_complete_days,
                         primary_minimum_events_each_endpoint = primary_RRR_minimum_events_each_endpoint,
                         clinical_eligibility = eligibility_rule)],
  by = "endpoint", all.x = TRUE, sort = FALSE
)
expected_endpoint_success <- c(strict_N20_only_vs_K80 = 32L, confirmed_N20_N21_vs_K80 = 32L,
                               composition_N23_vs_K80 = 31L, combined_original_N20_N23_vs_K80 = 32L)
if (nrow(coverage) != 4L || any(coverage$clinical_units != 32L) || any(coverage$fits_attempted != 32L) ||
    any(coverage$fits_successful != expected_endpoint_success[coverage$endpoint]) ||
    any(coverage$fits_failed != 32L - expected_endpoint_success[coverage$endpoint]) ||
    any(reproduction_summary$statuses_exact != 32L) ||
    any(reproduction_summary$event_counts_exact != 32L) ||
    any(reproduction_summary$hac14_hac28_estimates_exact != reproduction_summary$bridge_successes)) {
  stop("M07 endpoint coverage or saved-fit reproduction audit failed", call. = FALSE)
}
write_csv_once(coverage, file.path(outdir, "clinical_bridge_m07_coverage_audit.csv"))
write_csv_once(reproduction_summary, file.path(outdir, "clinical_bridge_m07_reproduction_audit.csv"))

spec <- data.table(
  tag = c("br_adult_matched", "mx_adult_matched", "br_standardized_matched", "mx_standardized_matched",
          "mx_strict_N20", "mx_confirmed_N20_N21", "mx_N23"),
  country = c("brazil", "mexico", "brazil", "mexico", "mexico", "mexico", "mexico"),
  scenario = c("latebilling_reconstruction_2014_2021", "mx_complete_source_2015_2022",
               "latebilling_reconstruction_2014_2021", "mx_complete_source_2015_2022",
               "mx_complete_source_2015_2022", "mx_complete_source_2015_2022", "mx_complete_source_2015_2022"),
  analysis_scope = c("adult_unstandardized_matched_standardization_success_units",
                     "adult_unstandardized_matched_standardization_success_units",
                     "observed_case_mix_standardized", "observed_case_mix_standardized",
                     "all_age_complete_source_same_primary_RRR_units",
                     "all_age_complete_source_same_primary_RRR_units",
                     "all_age_complete_source_same_primary_RRR_units"),
  endpoint = c("broad_N20_N23_vs_K80", "broad_N20_N23_vs_K80",
               "broad_N20_N23_vs_K80", "broad_N20_N23_vs_K80",
               "strict_N20_only_vs_K80", "confirmed_N20_N21_vs_K80", "composition_N23_vs_K80")
)
display <- merge(spec, pool, by = c("country", "scenario", "analysis_scope", "endpoint"),
                 all.x = TRUE, sort = FALSE)
display <- display[pooling == "DL" & covariance == "hac14" & status == "ok"]
if (nrow(display) != 7L || anyDuplicated(display$tag) ||
    !setequal(display$tag, spec$tag) ||
    any(display$units_fit != display$k) ||
    any(display$units_total != c(208L, 16L, 208L, 16L, 32L, 32L, 32L))) {
  stop("freeze/display row filters did not select exactly the seven expected DL/HAC14 rows", call. = FALSE)
}
display[, method := pooling]
display[, upstream_source := source]
display[, upstream_source_file := normalizePath(pool_path, winslash = "/", mustWork = TRUE)]
display[, upstream_source_row_ordinal := as.integer(upstream_source_row_ordinal)]
display[, upstream_source_filter := paste0("country=", country, ";scenario=", scenario,
                                             ";analysis_scope=", analysis_scope, ";endpoint=", endpoint,
                                             ";pooling=DL;covariance=hac14;status=ok")]
display[, upstream_field_map := "RR=RR;lo=lo;hi=hi;k=k;I2=I2;concordance=concordance;method=pooling;covariance=covariance"]
display[, analysis_origin := fcase(
  tag %in% c("br_adult_matched", "mx_adult_matched"),
  "new M06 bridge fit on exact standardization-success/common HAC14-eligible unit set",
  tag %in% c("br_standardized_matched", "mx_standardized_matched"),
  "validated runner per-unit standardized estimates re-pooled with observed case-mix weights",
  default = "M07 bridge fit exactly reproduced saved validated complete-source estimates"
)]
display <- display[, .(country, tag, RR, lo, hi, k, I2, concordance, method, covariance,
                       scenario, analysis_scope, endpoint, units_total, units_fit, units_failed, status,
                       analysis_origin, upstream_source, upstream_source_file, upstream_source_row_ordinal,
                       upstream_source_filter, upstream_field_map)]
write_csv_once(display, file.path(outdir, "clinical_bridge_display_rows.csv"))

exit_code <- suppressWarnings(as.integer(readLines(exit_file, warn = FALSE)[1L]))
audit_paths <- c(
  source_script, production_log, file.path(base_root, "results/mx_updated/primary_complete_source_2015_2022_MX_per_unit.csv"),
  file.path(base_root, "results/mx_updated/run_manifest.json"),
  file.path(base_root, "results/clinical_full_MX_complete/clinical_mexico_endpoint_per_unit.csv"),
  file.path(base_root, "results/clinical_full_MX_complete/clinical_mexico_endpoint_failures.csv"),
  file.path(outdir, required_files),
  file.path(outdir, "reproduction_check_attempt2/mexico_strict_reproduction_summary.csv")
)
audit_paths <- unique(audit_paths[file.exists(audit_paths)])
source_hashes <- as.list(stats::setNames(unname(tools::md5sum(audit_paths)), audit_paths))
audit <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  task = "Sensitivity M06 same-sample case-mix bridge and M07 Mexico strict phenotype results",
  audit_mode = "reconstructed from completed production outputs; no model refits",
  output_directory = outdir,
  production_run = list(pid = production_pid, exit_code = exit_code,
                        log = production_log,
                        terminal_error = "object 'preferred_weight_table' not found in the original final-audit block; all M06/M07 numeric result CSVs had been written first.",
                        outputs_completed_before_terminal_error = TRUE),
  m06_same_sample_and_weights = list(
    audit_file = file.path(outdir, "clinical_bridge_m06_set_weight_audit.csv"),
    detailed_membership_file = file.path(outdir, "clinical_bridge_m06_same_sample_membership.csv"),
    country_summary = lapply(seq_len(nrow(m06_audit)), function(i) as.list(m06_audit[i])),
    primary_weight_sources = list(brazil = "matched_source_2014_2020", mexico = "mx_complete_source_2015_2022"),
    preferred_brazil_sensitivity_source = "latebilling_reconstruction_2014_2021"
  ),
  m07_strict_phenotypes = list(
    main_target_units = nrow(main_membership),
    main_target_units_exactly_matching_calendar = sum(main_membership$primary_RRR_and_clinical_unit_match),
    coverage_file = file.path(outdir, "clinical_bridge_m07_coverage_audit.csv"),
    reproduction_file = file.path(outdir, "clinical_bridge_m07_reproduction_audit.csv"),
    endpoint_summary = lapply(seq_len(nrow(coverage)), function(i) as.list(coverage[i])),
    analysis_origin = "M07 bridge refit exactly reproduced saved validated complete-source per-unit results"
  ),
  freeze_display_rows = list(
    file = file.path(outdir, "clinical_bridge_display_rows.csv"),
    tags = display$tag,
    method = "DL",
    covariance = "HAC14",
    row_count = nrow(display)
  ),
  interpretation = "Aggregate administrative diagnosis associations; pooled outputs retain all clinical units and distinguish fit failures from primary RRR eligibility.",
  software = list(R = R.version.string, data_table = as.character(packageVersion("data.table")),
                  jsonlite = as.character(packageVersion("jsonlite")),
                  libPaths = .libPaths()),
  source_files_md5 = source_hashes
)
write_json_once(audit, file.path(outdir, "clinical_bridge_analysis_audit.json"))

manifest_files <- c(
  file.path(outdir, "clinical_bridge_analysis_audit.json"),
  file.path(outdir, "clinical_bridge_display_rows.csv"),
  file.path(outdir, "clinical_bridge_m06_set_weight_audit.csv"),
  file.path(outdir, "clinical_bridge_m06_same_sample_membership.csv"),
  file.path(outdir, "clinical_bridge_m07_coverage_audit.csv"),
  file.path(outdir, "clinical_bridge_m07_reproduction_audit.csv"),
  file.path(outdir, "clinical_bridge_pool_summary.csv"),
  file.path(outdir, "mexico_strict_phenotypes_per_unit.csv"),
  file.path(outdir, "mexico_strict_phenotypes_failures.csv"),
  file.path(outdir, "mexico_strict_phenotypes_pooled.csv"),
  file.path(outdir, "mexico_strict_eligibility_summary.csv"),
  file.path(outdir, "mexico_strict_broad_same_unit_reference_check.csv")
)
manifest_files <- manifest_files[file.exists(manifest_files)]
write_json_once(list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = "finalized_from_saved_outputs_without_refitting",
  production_pid = production_pid, production_exit_code = exit_code,
  production_log = production_log, finalizer_source = source_script,
  display_rows = file.path(outdir, "clinical_bridge_display_rows.csv"),
  display_tags = display$tag,
  artifact_md5 = as.list(stats::setNames(unname(tools::md5sum(manifest_files)), manifest_files))
), file.path(outdir, "clinical_bridge_finalization_manifest.json"))
cat("CLINICAL_BRIDGE_FINALIZE_OK\n")

