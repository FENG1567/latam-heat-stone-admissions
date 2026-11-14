#!/usr/bin/env Rscript

# Lightweight M07 target/membership smoke using the completed M06 checkpoint
# and saved validated complete-source endpoint fits. The production bridge
# script refits the M07 endpoints after root approval.

raw <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  idx <- match(name, raw)
  if (is.na(idx) || idx == length(raw)) default else raw[[idx + 1L]]
}
root <- normalizePath(arg_value("--root"), winslash = "/", mustWork = TRUE)
base_root <- root
local_rlib <- file.path(base_root, "Rlib")
if (dir.exists(local_rlib)) .libPaths(c(local_rlib, .libPaths()))
suppressPackageStartupMessages(library(data.table))
smoke_dir <- normalizePath(arg_value("--smoke-dir"), winslash = "/", mustWork = TRUE)
message("R library paths: ", paste(.libPaths(), collapse = .Platform$path.sep))
message(paste(capture.output(sessionInfo()), collapse = "\n"))

canonical_unit <- function(x) {
  z <- trimws(as.character(x))
  numeric_code <- grepl("^\\d+$", z)
  z[numeric_code] <- sub("^0+(?=\\d)", "", z[numeric_code], perl = TRUE)
  z[!nzchar(z)] <- "unknown"
  z
}

target_path <- file.path(smoke_dir, "clinical_bridge_targets_audit.csv")
membership_path <- file.path(smoke_dir, "clinical_bridge_unit_membership.csv")
target <- fread(target_path, colClasses = list(character = "unit"))[
  country == "mexico" & scenario == "mx_complete_source_2015_2022"
]
target[, unit := canonical_unit(unit)]
if (nrow(target) != 32L || anyDuplicated(target$unit)) stop("M06 checkpoint lacks the exact 32-state complete-source target set", call. = FALSE)

main_path <- file.path(root, "results/mx_updated/primary_complete_source_2015_2022_MX_per_unit.csv")
manifest_path <- file.path(root, "results/mx_updated/run_manifest.json")
main <- fread(main_path, colClasses = list(character = "unit"))[contrast == "y_urinary_vs_y_K80"]
main[, unit := canonical_unit(unit)]
main_unit <- main[, .(
  main_rrr_days = unique(n_complete)[1L],
  main_rrr_target_days = unique(target_n_complete)[1L],
  main_rrr_p95 = unique(p95)[1L],
  main_rrr_median = unique(median)[1L]
), by = unit]
if (anyDuplicated(main_unit$unit)) stop("updated main RRR has duplicate units after canonicalization", call. = FALSE)

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
scenario <- manifest$scenarios$primary_complete_source_2015_2022
metadata <- scenario$metadata
min_main_days <- as.integer(metadata$minimum_complete_dates)
min_main_events <- as.integer(metadata$minimum_events_each_endpoint)
if (as.integer(scenario$units_expected) != 32L || as.integer(scenario$units_fit) != 32L ||
    as.integer(scenario$units_failed) != 0L || min_main_days != 250L || min_main_events != 40L ||
    !setequal(main_unit$unit, target$unit)) {
  stop("updated primary RRR cohort/eligibility does not match the frozen 32-state manifest", call. = FALSE)
}

joined <- merge(
  main_unit,
  target[, .(unit, clinical_calendar_days = eligible_days,
             clinical_calendar_p95 = p95, clinical_calendar_median = median)],
  by = "unit", all = TRUE, sort = FALSE
)
joined[, unit_set_equal := !is.na(clinical_calendar_days) & !is.na(main_rrr_days)]
joined[, dates_match := main_rrr_days == clinical_calendar_days &
         main_rrr_target_days == clinical_calendar_days]
joined[, targets_match := abs(main_rrr_p95 - clinical_calendar_p95) < 1e-10 &
         abs(main_rrr_median - clinical_calendar_median) < 1e-10]
if (nrow(joined) != 32L || any(!joined$unit_set_equal) ||
    any(!joined$dates_match) || any(!joined$targets_match)) {
  stop("updated main RRR unit/date/target basis does not exactly match the complete clinical calendar")
}

member <- fread(membership_path, colClasses = list(character = "unit"))[
  country == "mexico" & scenario == "mx_complete_source_2015_2022"
]
if (nrow(member) != 32L || sum(member$standardization_status == "success") != 16L ||
    sum(member$standardization_status == "failed") != 16L) {
  stop("M06 checkpoint membership does not preserve the 16-success/16-failure standardization set")
}
joined <- merge(joined, member[, .(unit, standardization_status)], by = "unit", all.x = TRUE, sort = FALSE)
joined[, primary_RRR_minimum_complete_days := min_main_days]
joined[, primary_RRR_minimum_events_each_endpoint := min_main_events]
fwrite(joined, file.path(smoke_dir, "mexico_strict_main_RRR_target_join_smoke.csv"), na = "")

endpoints <- c("strict_N20_only_vs_K80", "confirmed_N20_N21_vs_K80",
               "composition_N23_vs_K80", "combined_original_N20_N23_vs_K80")
per_unit_path <- file.path(root, "results/clinical_full_MX_complete/clinical_mexico_endpoint_per_unit.csv")
failures_path <- file.path(root, "results/clinical_full_MX_complete/clinical_mexico_endpoint_failures.csv")
old <- fread(per_unit_path, colClasses = list(character = "unit"))[
  scenario == "mx_complete_source_2015_2022" & nuisance == "separate" & endpoint %in% endpoints
]
old[, unit := canonical_unit(unit)]
old_joined <- merge(old, target[, .(unit, clinical_calendar_days = eligible_days,
                                     clinical_calendar_p95 = p95, clinical_calendar_median = median)],
                    by = "unit", all.x = TRUE, sort = FALSE)
old_joined[, exact_target_match := n_days == clinical_calendar_days &
             abs(p95 - clinical_calendar_p95) < 1e-10 &
             abs(median - clinical_calendar_median) < 1e-10]
if (nrow(old_joined) && any(!old_joined$exact_target_match)) {
  stop("saved complete-source endpoint fits do not match the full clinical calendar targets")
}
failure_rows <- fread(failures_path, colClasses = list(character = "unit"))[
  scenario == "mx_complete_source_2015_2022" & nuisance == "separate" & endpoint %in% endpoints
]
endpoint_summary <- rbindlist(lapply(endpoints, function(id) {
  ok_n <- uniqueN(old[endpoint == id, unit])
  fail_n <- uniqueN(failure_rows[endpoint == id, unit])
  data.table(
    endpoint = id, clinical_units = 32L, saved_fit_successes = ok_n,
    saved_fit_failures = fail_n, main_RRR_units = 32L,
    main_RRR_minimum_days = min_main_days, main_RRR_minimum_events_each_endpoint = min_main_events,
    clinical_eligibility = ">=400 eligible days and positive K80/numerator counts"
  )
}), fill = TRUE)
fwrite(endpoint_summary, file.path(smoke_dir, "mexico_strict_existing_fit_coverage_smoke.csv"), na = "")

audit <- list(
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  mode = "M07 target/membership smoke; reuses completed M06 checkpoint",
  m06_smoke_directory = smoke_dir,
  main_rrr_per_unit = main_path,
  main_rrr_manifest = manifest_path,
  main_rrr_minimum_complete_days = min_main_days,
  main_rrr_minimum_events_each_endpoint = min_main_events,
  m07_clinical_endpoint_rule = "At least 400 complete eligible exposure days and positive K80/numerator counts; no extra 30-event cutoff.",
  exact_main_unit_date_target_matches = nrow(joined),
  standardization_success_units = sum(member$standardization_status == "success"),
  standardization_failure_units = sum(member$standardization_status == "failed"),
  saved_complete_source_endpoint_fits_checked = nrow(old_joined),
  saved_endpoint_fits_exact_target_matches = sum(old_joined$exact_target_match, na.rm = TRUE),
  endpoint_coverage = endpoint_summary,
  claim_boundary = "Uses aggregate state-level summaries and previously validated complete-source endpoint fits for this smoke check; full M07 refit remains in clinical_bridge.R and is gated on root approval."
)
jsonlite::write_json(audit, file.path(smoke_dir, "mexico_strict_main_RRR_target_join_smoke.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
cat("M07_TARGET_JOIN_SMOKE_OK\n")
