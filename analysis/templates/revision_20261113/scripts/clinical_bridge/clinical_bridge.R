#!/usr/bin/env Rscript

# Sensitivity M06 same-sample case-mix bridge and M07 Mexico strict-code
# replication. Reads aggregate inputs and existing validated runner outputs;
# writes only to the revision output directory passed with --outdir.

raw_args <- commandArgs(trailingOnly = TRUE)
root_pos <- match("--root", raw_args)
root_cli <- if (!is.na(root_pos) && root_pos < length(raw_args)) raw_args[[root_pos + 1L]] else NA_character_
if (is.na(root_cli)) stop("Pass --root pointing to the validated repair workspace", call. = FALSE)
base_root <- root_cli
local_rlib <- file.path(base_root, "Rlib")
if (dir.exists(local_rlib)) .libPaths(c(local_rlib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
message("R library paths: ", paste(.libPaths(), collapse = .Platform$path.sep))
message(paste(capture.output(sessionInfo()), collapse = "\n"))

runner_path <- file.path(root_cli, "scripts", "run_clinical_validated.R")
if (!file.exists(runner_path)) stop(sprintf("validated runner not found: %s", runner_path), call. = FALSE)
runner_expr <- parse(file = runner_path)
runner_tail <- paste(deparse(runner_expr[[length(runner_expr)]]), collapse = "")
if (!identical(runner_tail, "main()")) stop("validated runner no longer has the expected final main() call", call. = FALSE)
# Reuse the validated runner's own setup and functions but omit only main(),
# so no old result or fit-cache directory is written.
for (i in seq_len(length(runner_expr) - 1L)) eval(runner_expr[[i]], envir = .GlobalEnv)
if (!identical(args$country, "both") || !identical(args$mx_cohort, "complete")) {
  stop("Run with --country both --mx-cohort complete", call. = FALSE)
}
setDTthreads(1L)

smoke <- identical(Sys.getenv("CLINICAL_BRIDGE_SMOKE", unset = "0"), "1")
smoke_unit_cap <- if (smoke) 2L else Inf
analysis_outdir <- args$outdir
dir.create(analysis_outdir, recursive = TRUE, showWarnings = FALSE)

hash_file_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

hash_object_md5 <- function(x) {
  path <- tempfile(pattern = "clinical-bridge-object-")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3L)
  hash_file_md5(path)
}

eligible_calendar_for <- function(calendar, scenario) {
  keep <- as.integer(format(calendar$date, "%Y")) %in% as.integer(scenario$analysis_years) &
    is.finite(calendar$E3_actual)
  if (!is.null(scenario$offset_col)) {
    keep <- keep & is.finite(calendar[[scenario$offset_col]]) & calendar[[scenario$offset_col]] > 0
  }
  out <- calendar[keep, .(unit, date, E3 = E3_actual, population)]
  setorder(out, unit, date)
  out
}

load_standardization <- function(country, scenario_name, expected_success, clinical_units) {
  result_dir <- if (country == "brazil") {
    file.path(args$root, "results", "clinical_full_BR")
  } else {
    file.path(args$root, "results", "clinical_full_MX_complete")
  }
  prefix <- if (country == "brazil") "clinical_brazil" else "clinical_mexico"
  standardized <- fread(file.path(result_dir, paste0(prefix, "_case_mix_standardized_per_unit.csv")))
  failures <- fread(file.path(result_dir, paste0(prefix, "_case_mix_failures.csv")), fill = TRUE)
  strata <- fread(file.path(result_dir, paste0(prefix, "_case_mix_strata_per_unit.csv")), fill = TRUE)
  pooled <- fread(file.path(result_dir, paste0(prefix, "_case_mix_pooled.csv")), fill = TRUE)
  weights <- fread(file.path(result_dir, paste0(prefix, "_standardization_weights.csv")), fill = TRUE)

  observed <- standardized[scenario == scenario_name & weight_set == "pooled_observed_case_mix"]
  if (anyDuplicated(observed$unit)) stop(sprintf("duplicate observed-standardization rows for %s", country), call. = FALSE)
  success_units <- sort(unique(as.character(observed$unit)))
  failure_rows <- failures[scenario == scenario_name]
  if (anyDuplicated(failure_rows$unit)) stop(sprintf("duplicate standardization failure units for %s", country), call. = FALSE)
  failed_units <- sort(unique(as.character(failure_rows$unit)))
  if (length(success_units) != expected_success) {
    stop(sprintf("%s observed standardization success count is %d; expected %d", country, length(success_units), expected_success), call. = FALSE)
  }
  if (!setequal(c(success_units, failed_units), clinical_units) || length(intersect(success_units, failed_units))) {
    stop(sprintf("%s standardization success/failure rows do not partition all clinical units", country), call. = FALSE)
  }
  if (nrow(weights) != 6L) stop(sprintf("%s standardization weights do not contain six age-sex strata", country), call. = FALSE)

  membership <- data.table(country = country, scenario = scenario_name, unit = clinical_units)
  membership[, standardization_status := fifelse(unit %in% success_units, "success", "failed")]
  membership[, standardization_failure := failure_rows$failure[match(unit, failure_rows$unit)]]
  list(
    directory = result_dir, standardized = standardized, failures = failure_rows,
    strata = strata[scenario == scenario_name], pooled = pooled[scenario == scenario_name],
    weights = weights, success_units = success_units, failed_units = failed_units,
    membership = membership,
    source_files = c(
      file.path(result_dir, paste0(prefix, "_case_mix_standardized_per_unit.csv")),
      file.path(result_dir, paste0(prefix, "_case_mix_failures.csv")),
      file.path(result_dir, paste0(prefix, "_case_mix_strata_per_unit.csv")),
      file.path(result_dir, paste0(prefix, "_case_mix_pooled.csv")),
      file.path(result_dir, paste0(prefix, "_standardization_weights.csv"))
    )
  )
}

make_target_audit <- function(ctx, standardization) {
  rows <- lapply(ctx$units, function(unit_id) {
    cal <- ctx$eligible_calendar[unit == unit_id]
    target <- unit_weather_target(cal)
    ref <- standardization$strata[unit == unit_id]
    ref_p95 <- if (nrow(ref)) unique(as.numeric(ref$p95)) else numeric()
    ref_median <- if (nrow(ref)) unique(as.numeric(ref$median)) else numeric()
    ref_days <- if (nrow(ref)) unique(as.integer(ref$n_days)) else integer()
    is_success <- unit_id %in% standardization$success_units
    target_matches <- if (!is_success) NA else {
      length(ref_p95) == 1L && length(ref_median) == 1L && length(ref_days) == 1L &&
        abs(target$p95 - ref_p95) < 1e-10 && abs(target$median - ref_median) < 1e-10 &&
        length(cal$date) == ref_days
    }
    if (is_success && !isTRUE(target_matches)) {
      stop(sprintf("calendar/target mismatch against stored standardization fit for %s unit %s", ctx$country, unit_id), call. = FALSE)
    }
    data.table(
      country = ctx$country, scenario = ctx$scenario$name, unit = unit_id,
      eligible_days = nrow(cal), first_date = min(cal$date), last_date = max(cal$date),
      unique_exposure_days = uniqueN(cal$E3), p95 = target$p95, median = target$median,
      basis_df = 3L, basis_md5 = hash_object_md5(target$basis),
      date_vector_md5 = hash_object_md5(as.character(cal$date)),
      matches_stored_standardization_target = target_matches
    )
  })
  rbindlist(rows, fill = TRUE)
}

case_composition <- function(ctx, membership) {
  x <- copy(ctx$panel$selected)
  x <- merge(x, membership[, .(unit, standardization_status)], by = "unit", all.x = TRUE, sort = FALSE)
  x[, country := ctx$country]
  x[, scenario := ctx$scenario$name]
  x[, age_group_raw := trimws(as.character(age_group))]
  x[is.na(age_group_raw) | !nzchar(age_group_raw), age_group_raw := "unknown"]
  x[, sex_group_raw := tolower(trimws(as.character(sex_group)))]
  x[is.na(sex_group_raw) | !nzchar(sex_group_raw), sex_group_raw := "unknown"]
  x[, age_sex_eligibility := fcase(
    age_group_raw %in% c("18-44", "45-64", "65plus") & sex_group_raw %in% c("male", "female"), "adult_known_sex_in_standardization_weights",
    age_group_raw %in% c("18-44", "45-64", "65plus"), "adult_unknown_sex_excluded_from_weights",
    age_group_raw %in% c("<18", "under18", "0-17"), "under18_excluded_from_weights",
    age_group_raw %in% c("unknown", ""), "unknown_age_excluded_from_weights",
    default = "other_age_label_review"
  )]
  x[, diagnosis_group := fcase(
    dx3 == "K80", "K80",
    dx3 %in% c("N20", "N21", "N22", "N23"), "N20-N23",
    default = NA_character_
  )]
  x <- x[!is.na(diagnosis_group)]
  exact <- x[, .(admission_records = sum(n, na.rm = TRUE)),
             by = .(country, scenario, standardization_status,
                    age_sex_eligibility, age_group = age_group_raw, sex_group = sex_group_raw,
                    diagnosis_group, dx3)]
  exact[, pct_within_unit_status_and_diagnosis := admission_records / sum(admission_records),
        by = .(country, scenario, standardization_status, diagnosis_group, dx3)]
  list(
    exact = exact,
    audit = x[, .(aggregate_rows = .N, admission_records = sum(n, na.rm = TRUE)),
              by = .(country, scenario, standardization_status,
                     age_sex_eligibility, diagnosis_group)]
  )
}

build_daily_outcomes <- function(ctx, unit_id, adult_only = FALSE) {
  selected <- ctx$panel$selected[unit == unit_id]
  if (adult_only) selected <- selected[age_group %in% c("18-44", "45-64", "65plus")]
  daily <- selected[, .(
    y_K80_all = sum(y_K80_all, na.rm = TRUE),
    y_N20_only = sum(y_N20_only, na.rm = TRUE),
    y_N20_N21 = sum(y_N20_N21, na.rm = TRUE),
    y_N23 = sum(y_N23, na.rm = TRUE),
    y_N20_N23 = sum(y_N20_N23, na.rm = TRUE)
  ), by = .(date)]
  cal <- copy(ctx$eligible_calendar[unit == unit_id, .(date, E3, population)])
  if (!nrow(cal)) stop(sprintf("unit %s has no eligible calendar dates", unit_id), call. = FALSE)
  out <- merge(cal, daily, by = "date", all.x = TRUE, sort = FALSE)
  for (nm in c("y_K80_all", "y_N20_only", "y_N20_N21", "y_N23", "y_N20_N23")) {
    set(out, which(is.na(out[[nm]])), nm, 0)
  }
  setorder(out, date)
  out
}

fit_pair_unit <- function(ctx, unit_id, numerator, endpoint, sample_scope,
                          membership_status, target, adult_only = FALSE) {
  denom <- "y_K80_all"
  label <- paste0(numerator, "_vs_", denom)
  d <- NULL
  tryCatch({
    d <- build_daily_outcomes(ctx, unit_id, adult_only = adult_only)
    events_denominator <- sum(d[[denom]], na.rm = TRUE)
    events_numerator <- sum(d[[numerator]], na.rm = TRUE)
    if (nrow(d) < 400L) stop(sprintf("only %d eligible exposure days; 400 required", nrow(d)), call. = FALSE)
    if (uniqueN(d$E3) < 2L) stop("fewer than two distinct eligible E3 values", call. = FALSE)
    if (events_denominator <= 0 || events_numerator <= 0) {
      stop(sprintf("nonpositive total admissions: K80=%s; numerator=%s", events_denominator, events_numerator), call. = FALSE)
    }
    fit <- fit_stacked_contrast(
      data = d, outcomes = c(denom, numerator), exposure = "E3", date = "date",
      offset = ctx$scenario$offset_col, basis = target$basis, p95 = target$p95,
      median = target$median, nuisance = "separate", hac_lags = c(14L, 28L),
      full_covariance = FALSE
    )
    cc <- fit$contrasts[[label]]$coefficient_contrast
    if (is.null(cc)) stop(sprintf("expected joint contrast %s not returned", label), call. = FALSE)
    vals <- covariance_fields(cc)
    row <- data.table(
      country = ctx$country, scenario = ctx$scenario$name, unit = unit_id,
      endpoint = endpoint, sample_scope = sample_scope,
      standardization_status = membership_status, status = "ok",
      n_days = nrow(d), events_denominator = events_denominator,
      events_numerator = events_numerator, p95 = target$p95, median = target$median,
      exposure_days = target$exposure_days, formula = fit$audit$formula,
      warnings = as.character(toJSON(fit$audit$warnings, auto_unbox = TRUE))
    )
    for (nm in names(vals)) set(row, j = nm, value = vals[[nm]])
    row
  }, error = function(e) {
    data.table(
      country = ctx$country, scenario = ctx$scenario$name, unit = unit_id,
      endpoint = endpoint, sample_scope = sample_scope,
      standardization_status = membership_status, status = "failed",
      n_days = if (is.null(d)) NA_integer_ else nrow(d),
      events_denominator = if (is.null(d)) NA_real_ else sum(d[[denom]], na.rm = TRUE),
      events_numerator = if (is.null(d)) NA_real_ else sum(d[[numerator]], na.rm = TRUE),
      p95 = target$p95, median = target$median, exposure_days = target$exposure_days,
      formula = NA_character_, warnings = NA_character_, failure = conditionMessage(e)
    )
  })
}

pool_per_unit <- function(per_unit, analysis_scope) {
  out <- list()
  groups <- unique(per_unit[, .(country, scenario, endpoint)])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    rows <- per_unit[country == g[["country"]] & scenario == g[["scenario"]] & endpoint == g[["endpoint"]]]
    for (cov in c("hac14", "hac28", "day_hc0", "quasi")) {
      estimate_col <- paste0("logRRR_", cov)
      variance_col <- paste0("var_", cov)
      if (!estimate_col %in% names(rows) || !variance_col %in% names(rows)) next
      for (method in c("DL", "REML_HK")) {
        p <- pool_meta(rows[[estimate_col]], rows[[variance_col]], method = method)
        out[[length(out) + 1L]] <- data.table(
          country = g[["country"]], scenario = g[["scenario"]], analysis_scope = analysis_scope,
          endpoint = g[["endpoint"]], covariance = cov, pooling = method,
          status = if (is.null(p$failure)) "ok" else "failed",
          units_total = nrow(rows), units_fit = sum(rows$status == "ok"),
          units_failed = sum(rows$status != "ok"), k = p$k,
          logRRR = p$estimate, se = p$se, RR = p$RR,
          lo = p$RR95CI[1L], hi = p$RR95CI[2L],
          I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
          failure = if (is.null(p$failure)) NA_character_ else p$failure,
          source = "clinical_bridge_refit"
        )
      }
    }
  }
  rbindlist(out, fill = TRUE)
}

pool_standardized <- function(std, country, scenario_name) {
  selected <- std$standardized[scenario == scenario_name &
                                 weight_set %in% c("pooled_observed_case_mix", "uniform_1_over_6")]
  out <- list()
  for (w in unique(selected$weight_set)) {
    rows <- selected[weight_set == w]
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(rows$logRRR, rows$variance, method = method)
      original <- std$pooled[weight_set == w & pooling == method]
      if (nrow(original) != 1L) stop(sprintf("missing/duplicate stored standardized pool: %s %s", country, w), call. = FALSE)
      if (p$k != as.integer(original$k) ||
          abs(p$estimate - as.numeric(original$logRRR)) > 1e-10 ||
          abs(p$RR95CI[1L] - as.numeric(original$RR_lo)) > 1e-10 ||
          abs(p$RR95CI[2L] - as.numeric(original$RR_hi)) > 1e-10) {
        stop(sprintf("recomputed standardized pool differs from validated output: %s %s %s", country, w, method), call. = FALSE)
      }
      estimand <- if (w == "pooled_observed_case_mix") "observed_case_mix_standardized" else "uniform_case_mix_standardized"
      out[[length(out) + 1L]] <- data.table(
        country = country, scenario = scenario_name, analysis_scope = estimand,
        endpoint = "broad_N20_N23_vs_K80", covariance = "hac14", pooling = method,
        status = "ok", units_total = nrow(rows), units_fit = nrow(rows), units_failed = 0L,
        k = p$k, logRRR = p$estimate, se = p$se, RR = p$RR,
        lo = p$RR95CI[1L], hi = p$RR95CI[2L], I2 = p$I2, tau2 = p$tau2,
        concordance = p$concordance, failure = NA_character_,
        source = "validated_runner_per_unit_repool"
      )
    }
  }
  rbindlist(out, fill = TRUE)
}

pool_standardized_subset <- function(std, country, scenario_name, unit_set, suffix) {
  selected <- std$standardized[
    scenario == scenario_name & unit %in% unit_set &
      weight_set %in% c("pooled_observed_case_mix", "uniform_1_over_6")
  ]
  out <- list()
  for (w in c("pooled_observed_case_mix", "uniform_1_over_6")) {
    rows <- selected[weight_set == w]
    if (nrow(rows) != length(unit_set) || anyDuplicated(rows$unit)) {
      stop(sprintf("standardized common-set membership mismatch for %s %s", country, w), call. = FALSE)
    }
    for (method in c("DL", "REML_HK")) {
      p <- pool_meta(rows$logRRR, rows$variance, method = method)
      estimand <- if (w == "pooled_observed_case_mix") "observed_case_mix_standardized" else "uniform_case_mix_standardized"
      out[[length(out) + 1L]] <- data.table(
        country = country, scenario = scenario_name,
        analysis_scope = paste0(estimand, "_", suffix),
        endpoint = "broad_N20_N23_vs_K80", covariance = "hac14", pooling = method,
        status = if (is.null(p$failure)) "ok" else "failed",
        units_total = length(unit_set), units_fit = nrow(rows), units_failed = length(unit_set) - nrow(rows),
        k = p$k, logRRR = p$estimate, se = p$se, RR = p$RR,
        lo = p$RR95CI[1L], hi = p$RR95CI[2L], I2 = p$I2, tau2 = p$tau2,
        concordance = p$concordance,
        failure = if (is.null(p$failure)) NA_character_ else p$failure,
        source = "validated_runner_per_unit_repool_common_success_set"
      )
    }
  }
  rbindlist(out, fill = TRUE)
}

load_all_age_baseline <- function(ctx, standardization) {
  prefix <- if (ctx$country == "brazil") "clinical_brazil" else "clinical_mexico"
  result_dir <- if (ctx$country == "brazil") file.path(args$root, "results", "clinical_full_BR") else file.path(args$root, "results", "clinical_full_MX_complete")
  per_unit_path <- file.path(result_dir, paste0(prefix, "_endpoint_per_unit.csv"))
  pooled_path <- file.path(result_dir, paste0(prefix, "_endpoint_pooled.csv"))
  failure_path <- file.path(result_dir, paste0(prefix, "_endpoint_failures.csv"))
  if (!file.exists(per_unit_path) || !file.exists(pooled_path)) stop(sprintf("validated all-age endpoint outputs missing for %s", ctx$country), call. = FALSE)
  per_unit <- fread(per_unit_path, colClasses = list(character = "unit"), fill = TRUE)[
    scenario == ctx$scenario$name & nuisance == "separate" & endpoint == "combined_original_N20_N23_vs_K80"
  ]
  pooled <- fread(pooled_path, fill = TRUE)[
    scenario == ctx$scenario$name & nuisance == "separate" & endpoint == "combined_original_N20_N23_vs_K80"
  ]
  failures <- if (file.exists(failure_path)) fread(failure_path, colClasses = list(character = "unit"), fill = TRUE)[
    scenario == ctx$scenario$name & nuisance == "separate" & endpoint == "combined_original_N20_N23_vs_K80"
  ] else data.table()
  if (!nrow(pooled)) stop(sprintf("validated all-age broad pooled baseline missing for %s", ctx$country), call. = FALSE)
  membership <- copy(standardization$membership)
  setnames(membership, "standardization_status", "unit_standardization_status")
  baseline_membership <- merge(membership, per_unit[, .(unit, endpoint_status = status)], by = "unit", all.x = TRUE, sort = FALSE)
  if (nrow(failures)) baseline_membership <- merge(
    baseline_membership, failures[, .(unit, endpoint_failure = failure)], by = "unit", all.x = TRUE, sort = FALSE
  )
  baseline_membership[is.na(endpoint_status), endpoint_status := "failed_or_absent_from_success_file"]
  pooled[, analysis_scope := "all_age_unstandardized_all_clinical_units"]
  pooled[, source := "validated_runner_existing_pooled"]
  pooled[, analysis_age_scope := "all ages, including unknown age/sex"]
  if (!"units_total" %in% names(pooled)) pooled[, units_total := nrow(baseline_membership)]
  if (!"units_failed" %in% names(pooled)) pooled[, units_failed := pmax(0L, as.integer(units_total) - as.integer(units_fit))]
  if (!"RR" %in% names(pooled)) setnames(pooled, "RRR", "RR")
  if (!"lo" %in% names(pooled)) setnames(pooled, "RR_lo", "lo")
  if (!"hi" %in% names(pooled)) setnames(pooled, "RR_hi", "hi")
  list(per_unit = per_unit, pooled = pooled, membership = baseline_membership,
       source_files = c(per_unit_path, pooled_path, failure_path))
}

new_contexts <- function() {
  br_calendar <- read_calendar(args$calendar)
  br_clinical <- read_clinical(args$brazil_aggregate)
  br_scenarios <- source_scenarios(include_latebilling = TRUE)
  br_idx <- which(vapply(br_scenarios, function(z) z$name == "latebilling_reconstruction_2014_2021", logical(1)))[1L]
  br_scenario <- br_scenarios[[br_idx]]
  br_panel <- build_daily_panel(br_clinical, br_calendar$data, br_scenario)
  br_units <- sort(unique(br_panel$panel$unit))
  br_ctx <- list(
    country = "brazil", scenario = br_scenario, panel = br_panel, clinical = br_clinical,
    calendar = br_calendar$data, calendar_audit = br_calendar$audit,
    eligible_calendar = eligible_calendar_for(br_calendar$data, br_scenario), units = br_units
  )

  mx_calendar <- read_calendar(args$mx_calendar)
  mx_clinical <- read_mexico_clinical(args$mexico_aggregate, cohort = "complete")
  mx_scenario <- mexico_primary_scenario(cohort = "complete")
  mx_panel <- build_daily_panel(mx_clinical, mx_calendar$data, mx_scenario)
  mx_units <- sort(unique(mx_panel$panel$unit))
  mx_ctx <- list(
    country = "mexico", scenario = mx_scenario, panel = mx_panel, clinical = mx_clinical,
    calendar = mx_calendar$data, calendar_audit = mx_calendar$audit,
    eligible_calendar = eligible_calendar_for(mx_calendar$data, mx_scenario), units = mx_units
  )
  list(brazil = br_ctx, mexico = mx_ctx)
}

contexts <- new_contexts()
standardizations <- list(
  brazil = load_standardization("brazil", "latebilling_reconstruction_2014_2021", 208L, contexts$brazil$units),
  mexico = load_standardization("mexico", "mx_complete_source_2015_2022", 16L, contexts$mexico$units)
)
for (country in names(contexts)) {
  contexts[[country]]$target_audit <- make_target_audit(contexts[[country]], standardizations[[country]])
}

membership_all <- rbindlist(lapply(names(contexts), function(country) standardizations[[country]]$membership), fill = TRUE)
fwrite(membership_all, file.path(analysis_outdir, "clinical_bridge_unit_membership.csv"), na = "")

composition_rows <- list()
composition_audit <- list()
for (country in names(contexts)) {
  z <- case_composition(contexts[[country]], standardizations[[country]]$membership)
  composition_rows[[country]] <- z$exact
  composition_audit[[country]] <- z$audit
}
fwrite(rbindlist(composition_rows, fill = TRUE), file.path(analysis_outdir, "clinical_bridge_case_composition.csv"), na = "")
fwrite(rbindlist(composition_audit, fill = TRUE), file.path(analysis_outdir, "clinical_bridge_case_eligibility.csv"), na = "")
target_audit <- rbindlist(lapply(contexts, function(ctx) ctx$target_audit), fill = TRUE)
fwrite(target_audit, file.path(analysis_outdir, "clinical_bridge_targets_audit.csv"), na = "")

weight_rows <- list()
for (country in names(standardizations)) {
  z <- copy(standardizations[[country]]$weights)
  set(z, j = "country", value = country)
  weight_source <- if (country == "brazil") "matched_source_2014_2020" else "mx_complete_source_2015_2022"
  z[, weight_source_scenario := weight_source]
  weight_rows[[country]] <- z
}
fwrite(rbindlist(weight_rows, fill = TRUE), file.path(analysis_outdir, "clinical_bridge_standardization_weights.csv"), na = "")

# The preferred Brazil late-billing cohort has a different adult case-mix
# distribution from the matched-source cohort used for the previously frozen
# weights. Reuse the saved late-billing six-stratum fits and their full HAC14
# covariance to calculate this limited reweighting sensitivity without refit.
preferred_weight_rows <- list()
preferred_reweight_rows <- list()
br_ctx <- contexts$brazil
br_std <- standardizations$brazil
preferred_info <- case_mix_inputs(br_ctx$clinical, list(br_ctx$panel), br_ctx$calendar)
preferred_weights <- preferred_info$weights
preferred_strata <- preferred_info$strata
preferred_case_counts <- preferred_info$weight_counts$cases[
  match(names(preferred_weights), preferred_info$weight_counts$stratum)
]
if (anyNA(preferred_case_counts) ||
    any(abs(as.numeric(preferred_weights) - preferred_case_counts / sum(preferred_case_counts)) > 1e-12)) {
  stop("preferred-cohort case weights do not equal counts in explicit standardized-stratum order", call. = FALSE)
}
preferred_weights_table <- data.table(
  country = "brazil", weight_source_scenario = br_ctx$scenario$name,
  stratum = names(preferred_weights), pooled_K80_N20_N23_cases = as.numeric(preferred_case_counts),
  weight = as.numeric(preferred_weights)
)
preferred_weight_rows[["brazil"]] <- preferred_weights_table
for (unit_id in br_std$success_units) {
  cache_candidates <- file.path(
    args$root, "results", "clinical_full_BR", paste0("batch_", 1:4),
    "cache", "standardized", "brazil", br_ctx$scenario$name, paste0(unit_id, ".rds")
  )
  cache_path <- cache_candidates[file.exists(cache_candidates)]
  if (length(cache_path) != 1L) stop(sprintf("expected one saved Brazil standardization fit cache for unit %s", unit_id), call. = FALSE)
  fit_cache <- readRDS(cache_path)
  # Older validated caches persist fitted stratum contrasts and their joint
  # covariance, but not a beta vector. Reconstruct it from the exact stored
  # stratum estimates and assert covariance order before use.
  if (!identical(fit_cache$status, "ok") || as.character(fit_cache$unit) != unit_id ||
      !identical(fit_cache$scenario, br_ctx$scenario$name) ||
      is.null(fit_cache$stratum_estimates) || is.null(fit_cache$V12)) {
    stop(sprintf("saved Brazil fit cache lacks its successful stratum estimates or V12 for unit %s", unit_id), call. = FALSE)
  }
  beta_names <- unlist(lapply(c("K80", "N20_N23"), function(disease) paste(disease, preferred_strata, sep = "__")))
  beta <- stats::setNames(rep(NA_real_, length(beta_names)), beta_names)
  strata_estimates <- as.data.table(fit_cache$stratum_estimates)
  for (disease in c("K80", "N20_N23")) for (stratum in preferred_strata) {
    disease_id <- disease
    stratum_id <- stratum
    coefficient_name <- paste(disease_id, stratum_id, sep = "__")
    estimate_row <- strata_estimates[disease == disease_id & stratum == stratum_id]
    if (nrow(estimate_row) != 1L) stop(sprintf("missing/duplicate cached coefficient %s for unit %s", coefficient_name, unit_id), call. = FALSE)
    beta[[coefficient_name]] <- as.numeric(estimate_row$logRR)
  }
  if (any(!is.finite(beta)) || !identical(rownames(fit_cache$V12), beta_names) ||
      !identical(colnames(fit_cache$V12), beta_names)) {
    stop(sprintf("cached V12 dimensions/order do not match reconstructed 12-coefficient vector for unit %s", unit_id), call. = FALSE)
  }
  original_weights <- stats::setNames(as.numeric(br_std$weights$weight), as.character(br_std$weights$stratum))
  uniform_weights <- stats::setNames(rep(1 / length(preferred_strata), length(preferred_strata)), preferred_strata)
  reproduced_observed <- summarize_standardized(beta, fit_cache$V12, original_weights,
                                                 preferred_strata, "pooled_observed_case_mix")
  reproduced_uniform <- summarize_standardized(beta, fit_cache$V12, uniform_weights,
                                                preferred_strata, "uniform_1_over_6")
  existing <- br_std$standardized[
    scenario == br_ctx$scenario$name & unit == unit_id & weight_set == "pooled_observed_case_mix"
  ]
  existing_uniform <- br_std$standardized[
    scenario == br_ctx$scenario$name & unit == unit_id & weight_set == "uniform_1_over_6"
  ]
  if (nrow(existing) != 1L ||
      nrow(existing_uniform) != 1L ||
      abs(as.numeric(existing$logRRR) - as.numeric(reproduced_observed$logRRR)) > 1e-10 ||
      abs(as.numeric(existing$variance) - as.numeric(reproduced_observed$variance)) > 1e-10 ||
      abs(as.numeric(existing_uniform$logRRR) - as.numeric(reproduced_uniform$logRRR)) > 1e-10 ||
      abs(as.numeric(existing_uniform$variance) - as.numeric(reproduced_uniform$variance)) > 1e-10) {
    stop(sprintf("reconstructed cache coefficients/V12 do not reproduce stored observed and uniform standardization for unit %s", unit_id), call. = FALSE)
  }
  reweighted <- summarize_standardized(
    beta, fit_cache$V12, preferred_weights, preferred_strata,
    "latebilling_preferred_adult_known_sex_case_mix"
  )
  set(reweighted, j = "country", value = "brazil")
  set(reweighted, j = "unit", value = unit_id)
  set(reweighted, j = "scenario", value = br_ctx$scenario$name)
  set(reweighted, j = "scenario_label", value = unique(br_std$standardized[scenario == br_ctx$scenario$name, scenario_label]))
  set(reweighted, j = "n_days", value = fit_cache$n_days)
  set(reweighted, j = "p95", value = fit_cache$p95)
  set(reweighted, j = "median", value = fit_cache$median)
  set(reweighted, j = "weight_source_scenario", value = br_ctx$scenario$name)
  set(reweighted, j = "source", value = "saved_joint_stratum_fits_and_HAC14_covariance")
  preferred_reweight_rows[[length(preferred_reweight_rows) + 1L]] <- reweighted
}
preferred_reweighted <- rbindlist(preferred_reweight_rows, fill = TRUE)
preferred_weight_pooled <- list()
for (method in c("DL", "REML_HK")) {
  p <- pool_meta(preferred_reweighted$logRRR, preferred_reweighted$variance, method = method)
  preferred_weight_pooled[[method]] <- data.table(
    country = "brazil", scenario = br_ctx$scenario$name,
    analysis_scope = "latebilling_preferred_case_mix_reweighted",
    endpoint = "broad_N20_N23_vs_K80", pooling = method, covariance = "hac14",
    status = if (is.null(p$failure)) "ok" else "failed",
    units_total = nrow(preferred_reweighted), units_fit = nrow(preferred_reweighted),
    units_failed = 0L, k = p$k, logRRR = p$estimate, se = p$se,
    RR = p$RR, lo = p$RR95CI[1L], hi = p$RR95CI[2L],
    I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
    failure = if (is.null(p$failure)) NA_character_ else p$failure,
    source = "saved_joint_stratum_fits_and_HAC14_covariance"
  )
}
preferred_weight_pooled <- rbindlist(preferred_weight_pooled, fill = TRUE)
fwrite(rbindlist(preferred_weight_rows, fill = TRUE),
       file.path(analysis_outdir, "clinical_bridge_preferred_weight_sensitivity_weights.csv"), na = "")
fwrite(preferred_reweighted,
       file.path(analysis_outdir, "clinical_bridge_preferred_weight_reweighted_per_unit.csv"), na = "")
fwrite(preferred_weight_pooled,
       file.path(analysis_outdir, "clinical_bridge_preferred_weight_reweighted_pooled.csv"), na = "")

adult_perunit_rows <- list()
for (country in names(contexts)) {
  ctx <- contexts[[country]]
  std <- standardizations[[country]]
  fit_units <- ctx$units
  if (smoke) fit_units <- head(std$success_units, smoke_unit_cap)
  for (unit_id in fit_units) {
    target <- unit_weather_target(ctx$eligible_calendar[unit == unit_id])
    member_status <- std$membership[unit == unit_id, standardization_status]
    adult_perunit_rows[[length(adult_perunit_rows) + 1L]] <- fit_pair_unit(
      ctx, unit_id, numerator = "y_N20_N23", endpoint = "broad_N20_N23_vs_K80",
      sample_scope = "age_known_adult_sex_unrestricted", membership_status = member_status,
      target = target, adult_only = TRUE
    )
  }
}
adult_perunit <- rbindlist(adult_perunit_rows, fill = TRUE)
if (!"failure" %in% names(adult_perunit)) adult_perunit[, failure := NA_character_]
adult_perunit[is.na(failure), failure := NA_character_]
fwrite(adult_perunit, file.path(analysis_outdir, "clinical_bridge_adult_unstandardized_per_unit.csv"), na = "")
adult_pool_rows <- list()
matched_membership_rows <- list()
matched_k_audit_rows <- list()
common_success_pool_rows <- list()
common_success_unit_rows <- list()
for (country in names(contexts)) {
  rows <- adult_perunit[adult_perunit[["country"]] == country]
  adult_pool_rows[[length(adult_pool_rows) + 1L]] <- pool_per_unit(rows, "adult_unstandardized_fitted_clinical_units")
  std <- standardizations[[country]]
  attempted_units <- if (smoke) head(std$success_units, smoke_unit_cap) else std$success_units
  matched <- rows[unit %in% attempted_units]
  if (nrow(matched) != length(attempted_units) || anyDuplicated(matched$unit)) {
    stop(sprintf("adult bridge did not attempt each matched standardization unit exactly once for %s", country), call. = FALSE)
  }
  matched_pool <- pool_per_unit(matched, "adult_unstandardized_matched_standardization_success_units")
  adult_pool_rows[[length(adult_pool_rows) + 1L]] <- matched_pool
  matched_hac14 <- matched[
    status == "ok" & is.finite(logRRR_hac14) & is.finite(var_hac14) & var_hac14 > 0
  ]
  common_units <- sort(unique(as.character(matched_hac14$unit)))
  if (nrow(matched_hac14) != length(common_units)) stop(sprintf("duplicate adult primary estimates for %s", country), call. = FALSE)
  fit_attempted <- length(attempted_units)
  expected_attempted <- if (smoke) min(smoke_unit_cap, length(std$success_units)) else length(std$success_units)
  matched_membership <- copy(std$membership)
  adult_attempt_rows <- matched[, .(unit, adult_status = status, adult_failure = if ("failure" %in% names(matched)) failure else NA_character_,
                                    hac14_pool_eligible = status == "ok" & is.finite(logRRR_hac14) &
                                      is.finite(var_hac14) & var_hac14 > 0)]
  matched_membership <- merge(matched_membership, adult_attempt_rows, by = "unit", all.x = TRUE, sort = FALSE)
  matched_membership[, adult_fit_attempted := unit %in% attempted_units]
  matched_membership[is.na(adult_fit_attempted), adult_fit_attempted := FALSE]
  matched_membership[, common_bridge_success := unit %in% common_units]
  matched_membership[adult_fit_attempted == FALSE, adult_status := if (smoke) "not_attempted_smoke" else "not_attempted"]
  matched_membership[adult_fit_attempted == FALSE, hac14_pool_eligible := FALSE]
  if (any(matched_membership$standardization_status == "success" & matched_membership$adult_fit_attempted == FALSE) && !smoke) {
    stop(sprintf("some standardization-success units were not attempted in the full adult bridge for %s", country), call. = FALSE)
  }
  matched_membership_rows[[country]] <- matched_membership

  expected_full_success <- length(std$success_units)
  matched_k <- matched_pool[covariance == "hac14"]
  matched_k_audit_rows[[country]] <- data.table(
    country = country, scenario = unique(std$membership$scenario),
    expected_standardization_success_units = expected_full_success,
    adult_bridge_units_attempted = fit_attempted, adult_bridge_status_ok = sum(matched$status == "ok"),
    adult_bridge_hac14_pool_eligible = length(common_units),
    smoke_mode = smoke,
    primary_pool_k_DL = matched_k[pooling == "DL", k],
    primary_pool_k_REML_HK = matched_k[pooling == "REML_HK", k],
    exact_same_set_pooling = if (smoke) all(matched_k$k == expected_attempted) else all(matched_k$k == expected_full_success),
    failures_or_missing_hac14 = fit_attempted - length(common_units)
  )
  if (smoke && (fit_attempted != expected_attempted || length(common_units) != expected_attempted ||
                any(matched_k$k != expected_attempted))) {
    stop(sprintf("smoke adult bridge did not pool all %d attempted units for %s", expected_attempted, country), call. = FALSE)
  }
  if (!smoke && (fit_attempted != expected_full_success || length(common_units) != expected_full_success ||
                 any(matched_k$k != expected_full_success))) {
    if (!length(common_units)) stop(sprintf("no common adult/standardization-success units remain for %s", country), call. = FALSE)
    # Retain the failed bridge attempt above, then report standardized and
    # unstandardized pools on the identical common successful-unit set.
    common_success_pool_rows[[country]] <- rbindlist(list(
      pool_per_unit(matched_hac14[unit %in% common_units], "adult_unstandardized_common_success_set"),
      pool_standardized_subset(std, country, unique(std$membership$scenario), common_units, "common_success_set")
    ), fill = TRUE)
    common_success_unit_rows[[country]] <- matched_membership[unit %in% std$success_units]
  }
}
adult_pooled <- rbindlist(adult_pool_rows, fill = TRUE)
fwrite(adult_pooled, file.path(analysis_outdir, "clinical_bridge_adult_unstandardized_pooled.csv"), na = "")
matched_membership_audit <- rbindlist(matched_membership_rows, fill = TRUE)
fwrite(matched_membership_audit, file.path(analysis_outdir, "clinical_bridge_matched_unit_fit_audit.csv"), na = "")
matched_k_audit <- rbindlist(matched_k_audit_rows, fill = TRUE)
fwrite(matched_k_audit, file.path(analysis_outdir, "clinical_bridge_matched_k_audit.csv"), na = "")
if (length(common_success_pool_rows)) {
  fwrite(rbindlist(common_success_pool_rows, fill = TRUE),
         file.path(analysis_outdir, "clinical_bridge_common_success_pooled.csv"), na = "")
  fwrite(rbindlist(common_success_unit_rows, fill = TRUE),
         file.path(analysis_outdir, "clinical_bridge_common_success_units.csv"), na = "")
}

standardized_pooled <- rbindlist(lapply(names(contexts), function(country) {
  pool_standardized(
    standardizations[[country]], country,
    if (country == "brazil") "latebilling_reconstruction_2014_2021" else "mx_complete_source_2015_2022"
  )
}), fill = TRUE)
fwrite(standardized_pooled, file.path(analysis_outdir, "clinical_bridge_standardized_pooled.csv"), na = "")

all_age_baselines <- lapply(names(contexts), function(country) load_all_age_baseline(contexts[[country]], standardizations[[country]]))
names(all_age_baselines) <- names(contexts)
all_age_pool <- rbindlist(lapply(all_age_baselines, function(x) x$pooled), fill = TRUE)
fwrite(all_age_pool, file.path(analysis_outdir, "clinical_bridge_all_age_validated_baseline.csv"), na = "")
fwrite(rbindlist(lapply(all_age_baselines, function(x) x$membership), fill = TRUE),
       file.path(analysis_outdir, "clinical_bridge_all_age_baseline_membership.csv"), na = "")

mx <- contexts$mexico
mexico_fit_units <- mx$units
if (smoke) mexico_fit_units <- head(standardizations$mexico$success_units, smoke_unit_cap)
mx_endpoint_ids <- c("strict_N20_only_vs_K80", "confirmed_N20_N21_vs_K80",
                     "composition_N23_vs_K80", "combined_original_N20_N23_vs_K80")
mx_numerators <- c(
  strict_N20_only_vs_K80 = "y_N20_only",
  confirmed_N20_N21_vs_K80 = "y_N20_N21",
  composition_N23_vs_K80 = "y_N23",
  combined_original_N20_N23_vs_K80 = "y_N20_N23"
)
mexico_strict_rows <- list()
for (unit_id in mexico_fit_units) {
  target <- unit_weather_target(mx$eligible_calendar[unit == unit_id])
  member_status <- standardizations$mexico$membership[unit == unit_id, standardization_status]
  for (endpoint in mx_endpoint_ids) {
    mexico_strict_rows[[length(mexico_strict_rows) + 1L]] <- fit_pair_unit(
      mx, unit_id, numerator = unname(mx_numerators[[endpoint]]), endpoint = endpoint,
      sample_scope = "all_age_complete_source_2015_2022",
      membership_status = member_status, target = target, adult_only = FALSE
    )
  }
}
mexico_strict <- rbindlist(mexico_strict_rows, fill = TRUE)

main_rrr_path <- file.path(args$root, "results", "mx_updated", "primary_complete_source_2015_2022_MX_per_unit.csv")
main_manifest_path <- file.path(args$root, "results", "mx_updated", "run_manifest.json")
if (!file.exists(main_rrr_path)) stop("frozen primary Mexico RRR per-unit file is missing", call. = FALSE)
if (!file.exists(main_manifest_path)) stop("updated Mexico primary run manifest is missing", call. = FALSE)
main_rrr <- fread(main_rrr_path, colClasses = list(character = "unit"))[
  contrast == "y_urinary_vs_y_K80"
]
main_rrr[, unit := canonical_unit(unit)]
main_unit_audit <- main_rrr[, .(
  primary_RRR_unit_included = TRUE,
  primary_RRR_dates = unique(n_complete)[1L],
  primary_RRR_target_n_complete = unique(target_n_complete)[1L],
  primary_RRR_p95 = unique(p95)[1L],
  primary_RRR_median = unique(median)[1L]
), by = unit]
if (anyDuplicated(main_unit_audit$unit)) stop("primary Mexico RRR membership contains duplicate canonical units", call. = FALSE)
primary_manifest <- jsonlite::fromJSON(main_manifest_path, simplifyVector = FALSE)
primary_scenario_manifest <- primary_manifest$scenarios$primary_complete_source_2015_2022
primary_metadata <- primary_scenario_manifest$metadata
primary_minimum_days <- as.integer(primary_metadata$minimum_complete_dates)
primary_minimum_events <- as.integer(primary_metadata$minimum_events_each_endpoint)
if (as.integer(primary_scenario_manifest$units_fit) != length(mx$units) ||
    as.integer(primary_scenario_manifest$units_expected) != length(mx$units) ||
    as.integer(primary_scenario_manifest$units_failed) != 0L ||
    !setequal(main_unit_audit$unit, mx$units)) {
  stop("updated complete-source primary Mexico RRR unit set does not match the 32-state clinical cohort", call. = FALSE)
}
main_unit_audit <- merge(
  main_unit_audit,
  mx$target_audit[, .(unit, clinical_calendar_days = eligible_days,
                      clinical_calendar_p95 = p95, clinical_calendar_median = median)],
  by = "unit", all.x = TRUE, sort = FALSE
)
main_unit_audit[, clinical_full_calendar_unit := unit %in% mx$units]
main_unit_audit[, primary_RRR_and_clinical_unit_match :=
                  clinical_full_calendar_unit & primary_RRR_dates == clinical_calendar_days &
                  primary_RRR_target_n_complete == clinical_calendar_days &
                  abs(primary_RRR_p95 - clinical_calendar_p95) < 1e-10 &
                  abs(primary_RRR_median - clinical_calendar_median) < 1e-10]
if (any(!main_unit_audit$primary_RRR_and_clinical_unit_match)) {
  stop("updated complete-source primary RRR unit/date/target basis does not exactly match clinical full-calendar basis", call. = FALSE)
}
fwrite(main_unit_audit, file.path(analysis_outdir, "mexico_strict_main_RRR_unit_membership.csv"), na = "")
mexico_strict[, unit := canonical_unit(unit)]
mexico_strict <- merge(mexico_strict, main_unit_audit, by = "unit", all.x = TRUE, sort = FALSE)
mexico_strict[is.na(primary_RRR_unit_included), primary_RRR_unit_included := FALSE]
mexico_strict[, main_RRR_target_matches_clinical_full_calendar :=
                abs(primary_RRR_p95 - clinical_calendar_p95) < 1e-10 &
                abs(primary_RRR_median - clinical_calendar_median) < 1e-10]

# Independently compare the new all-age broad Mexico model against the saved
# complete-source validated runner per-unit estimates on the same states.
mx_old <- all_age_baselines$mexico$per_unit[
  endpoint == "combined_original_N20_N23_vs_K80" & nuisance == "separate"
]
mx_new_broad <- mexico_strict[endpoint == "combined_original_N20_N23_vs_K80"]
mx_reproduction <- merge(
  mx_new_broad, mx_old[, .(unit, old_status = status, old_events_K80 = events_denominator,
                           old_events_N20_N23 = events_numerator, old_p95 = p95, old_median = median,
                           old_logRRR_hac14 = logRRR_hac14, old_var_hac14 = var_hac14,
                           old_logRRR_hac28 = logRRR_hac28, old_var_hac28 = var_hac28)],
  by = "unit", all.x = TRUE, suffixes = c("_new", "_old"), sort = FALSE
)
mx_reproduction[, logRRR_hac14_abs_difference := abs(logRRR_hac14 - old_logRRR_hac14)]
mx_reproduction[, var_hac14_abs_difference := abs(var_hac14 - old_var_hac14)]
mx_reproduction[, logRRR_hac28_abs_difference := abs(logRRR_hac28 - old_logRRR_hac28)]
mx_reproduction[, var_hac28_abs_difference := abs(var_hac28 - old_var_hac28)]
mx_reproduction[, event_counts_match := events_denominator == old_events_K80 & events_numerator == old_events_N20_N23]
mx_reproduction[, target_match := abs(p95 - old_p95) < 1e-10 & abs(median - old_median) < 1e-10]
mx_reproduction[, estimate_match := status == "ok" & old_status == "ok" &
                  logRRR_hac14_abs_difference < 1e-8 & var_hac14_abs_difference < 1e-8 &
                  logRRR_hac28_abs_difference < 1e-8 & var_hac28_abs_difference < 1e-8]
fwrite(mx_reproduction, file.path(analysis_outdir, "mexico_strict_broad_same_unit_reference_check.csv"), na = "")
if (any((mx_reproduction$status == "ok" & !mx_reproduction$estimate_match) |
        (mx_reproduction$status != "ok" & mx_reproduction$old_status == "ok") |
        !mx_reproduction$event_counts_match | !mx_reproduction$target_match, na.rm = TRUE)) {
  stop("new Mexico broad fit failed same-unit reproduction against validated complete-source outputs", call. = FALSE)
}

fwrite(mexico_strict, file.path(analysis_outdir, "mexico_strict_phenotypes_per_unit.csv"), na = "")
fwrite(mexico_strict[status != "ok"], file.path(analysis_outdir, "mexico_strict_phenotypes_failures.csv"), na = "")
mexico_strict_pooled_full <- pool_per_unit(mexico_strict, "all_age_complete_source_all_clinical_units")
mexico_strict_main_units <- mexico_strict[primary_RRR_unit_included == TRUE]
mexico_strict_pooled_main <- pool_per_unit(mexico_strict_main_units, "all_age_complete_source_same_primary_RRR_units")
mexico_strict_pooled <- rbindlist(list(mexico_strict_pooled_full, mexico_strict_pooled_main), fill = TRUE)
fwrite(mexico_strict_pooled, file.path(analysis_outdir, "mexico_strict_phenotypes_pooled.csv"), na = "")

mexico_eligibility <- mexico_strict[, .(
  units_total_in_clinical_calendar = uniqueN(mx$units),
  units_primary_RRR = uniqueN(unit[primary_RRR_unit_included == TRUE]),
  primary_RRR_minimum_complete_days = primary_minimum_days,
  primary_RRR_minimum_events_each_endpoint = primary_minimum_events,
  units_attempted = uniqueN(unit),
  fits_attempted = .N,
  fits_ok = sum(status == "ok"),
  fits_failed = sum(status != "ok"),
  units_with_successful_fit = uniqueN(unit[status == "ok"]),
  eligibility_rule = "validated clinical endpoint rule: >=400 complete eligible exposure days, >=2 E3 values, K80 and numerator totals >0; no additional 30-event cutoff"
), by = .(endpoint)]
fwrite(mexico_eligibility, file.path(analysis_outdir, "mexico_strict_eligibility_summary.csv"), na = "")

pool_summary <- rbindlist(list(standardized_pooled, adult_pooled, all_age_pool, preferred_weight_pooled,
                               mexico_strict_pooled_full, mexico_strict_pooled_main), fill = TRUE)
if (nrow(pool_summary)) {
  wanted <- c("country", "scenario", "analysis_scope", "endpoint", "pooling", "covariance",
              "k", "RR", "lo", "hi", "I2", "concordance", "logRRR", "se", "tau2",
              "units_total", "units_fit", "units_failed", "status", "failure", "source")
  setcolorder(pool_summary, c(intersect(wanted, names(pool_summary)), setdiff(names(pool_summary), wanted)))
}
fwrite(pool_summary, file.path(analysis_outdir, "clinical_bridge_pool_summary.csv"), na = "")

cohort_code_counts <- lapply(contexts, function(ctx) ctx$panel$selected[
  , .(admission_records = sum(n, na.rm = TRUE)), by = .(dx3)][order(dx3)
])
names(cohort_code_counts) <- names(contexts)
target_summary <- target_audit[, .(
  units = .N, days_min = min(eligible_days), days_max = max(eligible_days),
  days_unique = uniqueN(eligible_days), p95_min = min(p95), p95_max = max(p95),
  median_min = min(median), median_max = max(median),
  standardized_target_matches = sum(matches_stored_standardization_target %in% TRUE, na.rm = TRUE)
), by = .(country, scenario)]

hash_paths <- unique(c(
  runner_path, args$engine, args$analysis_plan, args$brazil_aggregate, args$calendar,
  args$mexico_aggregate, args$mx_calendar,
  unlist(lapply(standardizations, function(x) x$source_files)),
  unlist(lapply(all_age_baselines, function(x) x$source_files)),
  main_rrr_path, main_manifest_path
))
preferred_weight_table <- rbindlist(preferred_weight_rows, fill = TRUE)
audit <- list(
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  task = "Sensitivity M06 same-sample case-mix bridges and M07 Mexico strict-phenotype replication",
  smoke = smoke, smoke_unit_cap_per_country = if (smoke) smoke_unit_cap else NULL,
  output_directory = normalizePath(analysis_outdir, winslash = "/", mustWork = TRUE),
  source_files_md5 = as.list(stats::setNames(vapply(hash_paths, hash_file_md5, character(1)), hash_paths)),
  input_cohort_audits = lapply(contexts, function(ctx) list(source_panel = ctx$panel$audit, calendar = ctx$calendar_audit)),
  standardization_success_and_failure = lapply(standardizations, function(std) list(
    scenario = unique(std$membership$scenario), units_in_calendar = nrow(std$membership),
    successful_units = length(std$success_units), failed_or_unavailable_units = length(std$failed_units),
    success_units = std$success_units, failed_units = std$failed_units,
    failure_reasons = as.list(table(std$failures$failure)), stored_case_mix_weights = std$weights,
    weight_source_scenario = if (unique(std$membership$country) == "brazil") "matched_source_2014_2020" else "mx_complete_source_2015_2022"
  )),
  matched_bridge_set_audit = matched_k_audit,
  matched_bridge_failures_by_country = lapply(unique(matched_membership_audit$country), function(country_id) {
    z <- matched_membership_audit[country == country_id]
    z[, .(standardization_success_units = sum(standardization_status == "success"),
          adult_fit_attempted = sum(adult_fit_attempted),
          adult_fit_failures = sum(adult_fit_attempted & adult_status != "ok"),
          missing_or_nonpositive_hac14 = sum(adult_fit_attempted & !hac14_pool_eligible),
          common_success_units = sum(common_bridge_success)), by = .(country)]
  }),
  case_mix_weight_metadata = list(
    validated_runner_weights = rbindlist(weight_rows, fill = TRUE),
    validated_runner_weight_sources = list(brazil = "matched_source_2014_2020", mexico = "mx_complete_source_2015_2022"),
    preferred_brazil_latebilling_weights = preferred_weight_table
  ),
  target_summary = target_summary,
  age_sex_composition = rbindlist(composition_audit, fill = TRUE),
  diagnosis_code_counts = cohort_code_counts,
  same_unit_validation = list(
    mexico_primary_RRR_unit_count = nrow(main_unit_audit),
    mexico_primary_RRR_units_matching_complete_clinical_units = sum(main_unit_audit$clinical_full_calendar_unit),
    mexico_primary_RRR_exact_date_and_target_matches = sum(main_unit_audit$primary_RRR_and_clinical_unit_match),
    mexico_primary_RRR_minimum_complete_days = primary_minimum_days,
    mexico_primary_RRR_minimum_events_each_endpoint = primary_minimum_events,
    mexico_new_broad_success_rows_reproduced = sum(mx_reproduction$estimate_match, na.rm = TRUE),
    mexico_new_broad_reference_rows_checked = nrow(mx_reproduction)
  ),
  analysis_definitions = list(
    adult_bridge = "Age-known adults (18-44, 45-64, 65plus); sex unrestricted, so unknown sex remains in the unstandardized adult comparison. The matched bridge uses the exact standardization-success units, unit-day calendars, E3 P95/shared-median targets, and df=3 natural-spline basis.",
    standardization_weights = "Observed and uniform pools reproduce the existing validated runner's successful-unit per-unit outputs. Brazil late-billing fits reuse the runner's frozen matched-source 2014-2020 adult known-sex weights; this source is retained in weight metadata.",
    preferred_weight_sensitivity = "For Brazil late-billing, recomputed pooled adult known-sex case weights from the preferred 2014-2021 source/2014-2020 admission cohort and reweighted the already saved successful six-stratum fits using their joint HAC14 covariance; no refit and no change to the frozen primary weights.",
    all_age_baseline = "Existing validated runner full-cohort all-age broad N20-N23 versus K80 pooled output, retained as the original all-clinical-units baseline.",
    mexico_strict_phenotypes = "All-age complete-source Mexico; K80 versus N20-only, N20-N21, N23, and original broad N20-N23; all outcomes come from the complete-source aggregate and full-period state calendar.",
    existing_endpoint_eligibility = "At least 400 complete eligible exposure days, at least two distinct E3 values, and positive K80 and numerator totals; no additional event-count threshold. Age-sex standardization separately retains its existing 400-day and >=30-cases-per-required-stratum criteria. The Mexico strict output also reports its overlap with the frozen 32-state primary RRR sample and retains the original broad same-unit reference check.",
    frozen_primary_Mexico_RRR_eligibility = paste0(
      "Updated primary_complete_source_2015_2022 run_manifest rule: at least ",
      primary_minimum_days, " complete dates and at least ", primary_minimum_events,
      " events for each main endpoint; the frozen primary run retained all 32 states. M07 retains the validated clinical endpoint rule (>=400 days and positive events) without altering either threshold."
    ),
    model = "Validated engine fit_stacked_contrast; separate disease-specific day-of-week, month, and 12-df time-spline nuisance terms; quasi-Poisson, P95 versus shared median; HAC14/HAC28 and other covariance outputs when available; DL and REML-HK pooling.",
    interpretation = "Aggregate administrative diagnosis associations and case-mix bridges; no individual-risk or causal mechanism inference."
  ),
  software = list(R = R.version.string, data_table = as.character(packageVersion("data.table")),
                  jsonlite = as.character(packageVersion("jsonlite")))
)
jsonlite::write_json(audit, file.path(analysis_outdir, "clinical_bridge_analysis_audit.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
log_event("clinical_bridge_complete", smoke = smoke, output = analysis_outdir,
          adult_bridge_rows = nrow(adult_perunit), mexico_strict_rows = nrow(mexico_strict),
          pool_rows = nrow(pool_summary), audit = file.path(analysis_outdir, "clinical_bridge_analysis_audit.json"))
