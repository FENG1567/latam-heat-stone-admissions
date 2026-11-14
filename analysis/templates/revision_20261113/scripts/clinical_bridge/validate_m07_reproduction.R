#!/usr/bin/env Rscript

# Compare completed M07 bridge fits with saved validated complete-source fits.
# This reads aggregate state-level outputs and does not refit models.
suppressPackageStartupMessages(library(data.table))
raw <- commandArgs(trailingOnly = TRUE)
if (length(raw) != 4L) {
  stop("Usage: validate_m07_reproduction.R <saved-per-unit.csv> <saved-failures.csv> <bridge-per-unit.csv> <outdir>", call. = FALSE)
}
saved_path <- raw[[1L]]
failure_path <- raw[[2L]]
bridge_path <- raw[[3L]]
outdir <- raw[[4L]]

canonical_unit <- function(x) {
  z <- trimws(as.character(x))
  numeric_code <- grepl("^\\d+$", z)
  z[numeric_code] <- sub("^0+(?=\\d)", "", z[numeric_code], perl = TRUE)
  z[!nzchar(z)] <- "unknown"
  z
}

endpoints <- c("strict_N20_only_vs_K80", "confirmed_N20_N21_vs_K80",
               "composition_N23_vs_K80", "combined_original_N20_N23_vs_K80")
saved <- fread(saved_path, colClasses = list(character = "unit"))[
  scenario == "mx_complete_source_2015_2022" & nuisance == "separate" & endpoint %in% endpoints
]
saved_failures <- fread(failure_path, colClasses = list(character = "unit"))[
  scenario == "mx_complete_source_2015_2022" & nuisance == "separate" & endpoint %in% endpoints
]
bridge <- fread(bridge_path, colClasses = list(character = "unit"))[endpoint %in% endpoints]
saved[, unit := canonical_unit(unit)]
saved_failures[, unit := canonical_unit(unit)]
bridge[, unit := canonical_unit(unit)]
saved_all <- rbindlist(list(saved, saved_failures), fill = TRUE, use.names = TRUE)

required <- c("status", "n_days", "events_denominator", "events_numerator", "p95", "median",
              "formula", "logRRR_hac14", "var_hac14", "se_hac14", "RR_hac14", "lo_hac14", "hi_hac14",
              "logRRR_hac28", "var_hac28", "se_hac28", "RR_hac28", "lo_hac28", "hi_hac28")
missing_saved <- setdiff(required, names(saved_all))
missing_bridge <- setdiff(required, names(bridge))
if (length(missing_saved) || length(missing_bridge)) {
  stop(sprintf("required comparison fields missing (saved: %s; bridge: %s)",
               paste(missing_saved, collapse = ","), paste(missing_bridge, collapse = ",")), call. = FALSE)
}
if (anyDuplicated(saved_all[, .(endpoint, unit)]) || anyDuplicated(bridge[, .(endpoint, unit)])) {
  stop("saved or bridge M07 output has duplicate endpoint/unit rows", call. = FALSE)
}

setnames(saved_all, required, paste0("saved_", required))
setnames(bridge, required, paste0("bridge_", required))
joined <- merge(bridge, saved_all, by = c("endpoint", "unit"), all = TRUE, sort = FALSE)
joined[, date_match := !is.na(bridge_n_days) & bridge_n_days == saved_n_days & bridge_n_days == 2920L]
joined[, event_counts_match := !is.na(bridge_events_denominator) &
         bridge_events_denominator == saved_events_denominator &
         bridge_events_numerator == saved_events_numerator &
         bridge_events_denominator > 0 & bridge_events_numerator > 0]
joined[, target_fields_available := is.finite(saved_p95) & is.finite(saved_median)]
joined[, target_match := target_fields_available & abs(bridge_p95 - saved_p95) < 1e-10 &
         abs(bridge_median - saved_median) < 1e-10]
joined[, status_match := !is.na(bridge_status) & !is.na(saved_status) & bridge_status == saved_status]
joined[, formula_match := bridge_status != "ok" |
         gsub("\\s+", " ", trimws(bridge_formula)) == gsub("\\s+", " ", trimws(saved_formula))]
hac_cols <- c("logRRR_hac14", "var_hac14", "se_hac14", "RR_hac14", "lo_hac14", "hi_hac14",
              "logRRR_hac28", "var_hac28", "se_hac28", "RR_hac28", "lo_hac28", "hi_hac28")
reduce_and <- function(a, b) a & b
joined[, estimates_complete := bridge_status != "ok" |
         (Reduce(reduce_and, lapply(hac_cols, function(nm) is.finite(get(paste0("bridge_", nm))))) &
          Reduce(reduce_and, lapply(hac_cols, function(nm) is.finite(get(paste0("saved_", nm))))))]
joined[, estimates_match := bridge_status != "ok" |
         Reduce(reduce_and, lapply(hac_cols, function(nm) abs(get(paste0("bridge_", nm)) - get(paste0("saved_", nm))) < 1e-8))]
joined[, analysis_origin := "current bridge fits reproduced existing validated complete-source results"]

summary <- joined[, .(
  units = uniqueN(unit), saved_rows = sum(!is.na(saved_status)), bridge_rows = sum(!is.na(bridge_status)),
  saved_successes = sum(saved_status == "ok", na.rm = TRUE), bridge_successes = sum(bridge_status == "ok", na.rm = TRUE),
  saved_failures = sum(saved_status != "ok", na.rm = TRUE), bridge_failures = sum(bridge_status != "ok", na.rm = TRUE),
  dates_exact = sum(date_match, na.rm = TRUE), event_counts_exact = sum(event_counts_match, na.rm = TRUE),
  saved_target_rows = sum(target_fields_available, na.rm = TRUE),
  targets_exact_on_saved_rows = sum(target_match, na.rm = TRUE), statuses_exact = sum(status_match, na.rm = TRUE),
  formulas_exact_on_success = sum(formula_match & bridge_status == "ok", na.rm = TRUE),
  complete_hac14_hac28_rows = sum(estimates_complete & bridge_status == "ok", na.rm = TRUE),
  hac14_hac28_estimates_exact = sum(estimates_match & bridge_status == "ok", na.rm = TRUE),
  analysis_origin = unique(analysis_origin)
), by = .(endpoint)]

fwrite(joined, file.path(outdir, "mexico_strict_reproduction_vs_saved.csv"), na = "")
fwrite(summary, file.path(outdir, "mexico_strict_reproduction_summary.csv"), na = "")

if (nrow(joined) != 128L || nrow(summary) != length(endpoints) || any(summary$units != 32L) ||
    any(summary$saved_rows != 32L) || any(summary$bridge_rows != 32L) ||
    any(summary$dates_exact != 32L) || any(summary$event_counts_exact != 32L) ||
    any(summary$saved_target_rows != summary$saved_successes) ||
    any(summary$targets_exact_on_saved_rows != summary$saved_successes) || any(summary$statuses_exact != 32L) ||
    any(summary$formulas_exact_on_success != summary$bridge_successes) ||
    any(summary$complete_hac14_hac28_rows != summary$bridge_successes) ||
    any(summary$hac14_hac28_estimates_exact != summary$bridge_successes)) {
  stop("M07 bridge outputs do not fully reproduce the saved validated results", call. = FALSE)
}
jsonlite::write_json(list(
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  mode = "read-only aggregate comparison; no model refits",
  analysis_origin = "current bridge fits reproduced existing validated complete-source results",
  rows_compared = nrow(joined), endpoints = summary,
  fields_checked = list(days = "n_days", events = c("events_denominator", "events_numerator"),
                        targets = c("p95", "median"), estimates = hac_cols),
  claim_boundary = "Confirms same-unit saved-versus-bridge numerical reproduction for complete-source Mexico strict outcomes. Saved failure rows lack P95/median fields; the separate M07 target smoke confirms all 32 unit targets against the clinical calendar."
), file.path(outdir, "mexico_strict_reproduction_summary.json"),
auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
cat("M07_REPRODUCTION_CHECK_OK\n")

