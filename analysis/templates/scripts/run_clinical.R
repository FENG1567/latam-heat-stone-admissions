#!/usr/bin/env Rscript

# Staged Brazil and Mexico clinical endpoint analyses for the 2026-11-12 repair.
# Outcomes are SIH-RD admission rows, not emergency-department encounters.
# This runner consumes the aggregate extractor output only; it never handles
# N_AIH values or patient-level data.

# Prefer the repair-local library when run on the isolated Linux workspace;
# the installed Windows R library remains the fallback on the desktop.
file_arg_early <- grep("^--file=", commandArgs(FALSE), value = TRUE)
runner_script_path <- if (length(file_arg_early)) {
  normalizePath(sub("^--file=", "", file_arg_early[[1L]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("work/repair_20261112/scripts/run_clinical.R", winslash = "/", mustWork = FALSE)
}
if (length(file_arg_early)) {
  script_path_early <- sub("^--file=", "", file_arg_early[[1L]])
  repair_root_early <- dirname(dirname(normalizePath(script_path_early, winslash = "/", mustWork = FALSE)))
  repair_rlib <- file.path(repair_root_early, "Rlib")
  if (dir.exists(repair_rlib)) .libPaths(c(repair_rlib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

parse_cli <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "work/repair_20261112/scripts/run_clinical.R"
  script_file <- normalizePath(script_file, winslash = "/", mustWork = FALSE)
  default_root <- dirname(dirname(script_file))
  out <- list(
    root = default_root,
    country = "brazil",
    mx_cohort = "baseline",
    stage = "endpoint",
    include_latebilling = TRUE,
    force = FALSE,
    max_units = Inf,
    batch_index = 1L,
    batch_count = 1L,
    brazil_aggregate = NULL,
    mexico_aggregate = NULL,
    calendar = NULL,
    mx_calendar = NULL,
    engine = NULL,
    analysis_plan = NULL,
    outdir = NULL,
    io_cache_compatibility = NULL
  )
  raw <- commandArgs(trailingOnly = TRUE)
  if (any(raw %in% c("--help", "-h"))) {
    cat(paste(
      "Usage: Rscript run_clinical.R [options]",
      "  --root PATH                    repair workspace root",
      "  --stage endpoint|interaction|standardized|all",
      "  --country brazil|mexico|both",
      "  --mx-cohort baseline|complete  Mexico source cohort (default baseline)",
      "  --brazil-aggregate PATH        extractor aggregate CSV.GZ",
      "  --mexico-aggregate PATH        Mexico clinical long CSV.GZ",
      "  --calendar PATH                canonical Brazil calendar CSV.GZ",
      "  --mx-calendar PATH             canonical Mexico calendar CSV.GZ",
      "  --engine PATH                  repair_engine.R",
      "  --analysis-plan PATH           analysis_plan.json",
      "  --outdir PATH                  output directory (default ROOT/results/clinical)",
      "  --io-cache-compatibility PATH JSON authorization for an I/O-only cache reader fix",
      "  --max-units N                  optional deterministic unit cap for smoke runs",
      "  --batch-index N                deterministic round-robin batch number (1-based)",
      "  --batch-count N                total round-robin batches (default 1)",
      "  --force                        ignore existing per-unit cache rows",
      sep = "\n"
    ), "\n")
    quit(status = 0L)
  }
  i <- 1L
  while (i <= length(raw)) {
    token <- raw[[i]]
    if (token == "--force") {
      out$force <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(token, "--")) stop(sprintf("Unexpected argument: %s", token), call. = FALSE)
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
      key <- bits[[1L]]; value <- paste(bits[-1L], collapse = "=")
    } else {
      key <- sub("^--", "", token)
      if (i == length(raw) || startsWith(raw[[i + 1L]], "--")) {
        stop(sprintf("Option --%s needs a value", key), call. = FALSE)
      }
      value <- raw[[i + 1L]]
      i <- i + 1L
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    if (!key %in% names(out)) stop(sprintf("Unknown option --%s", key), call. = FALSE)
    if (key == "force") stop("Use --force without a value", call. = FALSE)
    if (key == "include_latebilling") out[[key]] <- tolower(value) %in% c("1", "true", "yes")
    else if (key %in% c("max_units", "batch_index", "batch_count")) out[[key]] <- as.numeric(value)
    else out[[key]] <- value
    i <- i + 1L
  }
  if (!out$stage %in% c("endpoint", "interaction", "standardized", "all")) {
    stop("--stage must be endpoint, interaction, standardized, or all", call. = FALSE)
  }
  if (!out$country %in% c("brazil", "mexico", "both")) stop("--country must be brazil, mexico, or both", call. = FALSE)
  if (!out$mx_cohort %in% c("baseline", "complete")) stop("--mx-cohort must be baseline or complete", call. = FALSE)
  if (!is.finite(out$max_units) && !is.infinite(out$max_units)) stop("--max-units must be positive", call. = FALSE)
  if (is.finite(out$max_units) && out$max_units < 1) stop("--max-units must be at least 1", call. = FALSE)
  if (!is.finite(out$batch_index) || !is.finite(out$batch_count) || out$batch_index < 1 || out$batch_count < 1 ||
      out$batch_index != floor(out$batch_index) || out$batch_count != floor(out$batch_count) ||
      out$batch_index > .Machine$integer.max || out$batch_count > .Machine$integer.max ||
      out$batch_index > out$batch_count) {
    stop("--batch-index and --batch-count must be positive integers with index <= count", call. = FALSE)
  }
  out$batch_index <- as.integer(out$batch_index)
  out$batch_count <- as.integer(out$batch_count)
  out$root <- normalizePath(out$root, winslash = "/", mustWork = FALSE)
  if (is.null(out$brazil_aggregate)) out$brazil_aggregate <- file.path(out$root, "results/brazil_clinical_aggregate_all_admission_years_by_source_year.csv.gz")
  if (is.null(out$mexico_aggregate)) {
    mx_filename <- if (out$mx_cohort == "complete") "mx_complete_source_long_2015_2022.csv.gz" else "mx_clinical_long.csv.gz"
    out$mexico_aggregate <- file.path(out$root, "data", mx_filename)
  }
  if (is.null(out$calendar)) out$calendar <- file.path(out$root, "data/br_calendar.csv.gz")
  if (is.null(out$mx_calendar)) {
    mx_calendar_filename <- if (out$mx_cohort == "complete") "mx_updated_calendar.csv.gz" else "mx_calendar.csv.gz"
    out$mx_calendar <- file.path(out$root, "data", mx_calendar_filename)
  }
  if (is.null(out$engine)) out$engine <- file.path(out$root, "scripts/repair_engine.R")
  if (is.null(out$analysis_plan)) out$analysis_plan <- file.path(out$root, "admin/analysis_plan.json")
  if (is.null(out$outdir)) {
    outdir_name <- if (out$mx_cohort == "complete") "clinical_mx_complete" else "clinical"
    out$outdir <- file.path(out$root, "results", outdir_name)
  }
  for (nm in c("brazil_aggregate", "mexico_aggregate", "calendar", "mx_calendar", "engine", "analysis_plan", "outdir")) {
    out[[nm]] <- normalizePath(out[[nm]], winslash = "/", mustWork = FALSE)
  }
  if (!is.null(out$io_cache_compatibility)) {
    out$io_cache_compatibility <- normalizePath(out$io_cache_compatibility, winslash = "/", mustWork = FALSE)
  }
  out
}

args <- parse_cli()
if (!file.exists(args$engine)) stop(sprintf("repair engine not found: %s", args$engine), call. = FALSE)
if (!file.exists(args$analysis_plan)) stop(sprintf("analysis plan not found: %s", args$analysis_plan), call. = FALSE)
if (args$country %in% c("brazil", "both") && !file.exists(args$calendar)) stop(sprintf("Brazil calendar not found: %s", args$calendar), call. = FALSE)
if (args$country %in% c("brazil", "both") && !file.exists(args$brazil_aggregate)) stop(sprintf("Brazil clinical aggregate not found: %s", args$brazil_aggregate), call. = FALSE)
if (args$country %in% c("mexico", "both") && !file.exists(args$mx_calendar)) stop(sprintf("Mexico calendar not found: %s", args$mx_calendar), call. = FALSE)
if (args$country %in% c("mexico", "both") && args$mx_cohort == "complete" && !file.exists(args$mexico_aggregate)) {
  stop(sprintf("Complete Mexico clinical source not found: %s", args$mexico_aggregate), call. = FALSE)
}
if (args$country %in% c("mexico", "both") && args$mx_cohort == "baseline" && !file.exists(args$mexico_aggregate) &&
    !length(list.files(dirname(args$mexico_aggregate), pattern = "^mx_clinical_year_[0-9]{4}\\.csv\\.gz$"))) {
  stop(sprintf("Mexico clinical aggregate or per-year caches not found beside: %s", args$mexico_aggregate), call. = FALSE)
}
source(args$engine, local = .GlobalEnv)

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "derived"), recursive = TRUE, showWarnings = FALSE)

log_event <- function(event, ...) {
  payload <- c(list(event = event, time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), list(...))
  cat(as.character(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", na = "null")), "\n", sep = "")
  flush.console()
}

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  z <- trimws(as.character(x))
  out <- as.Date(z, format = "%Y-%m-%d")
  miss <- is.na(out) & grepl("^\\d{8}$", z)
  if (any(miss)) out[miss] <- as.Date(z[miss], format = "%Y%m%d")
  out
}

canonical_unit <- function(x) {
  z <- trimws(as.character(x))
  numeric_code <- grepl("^\\d+$", z)
  z[numeric_code] <- sub("^0+(?=\\d)", "", z[numeric_code], perl = TRUE)
  z[!nzchar(z)] <- "unknown"
  z
}

normalize_dx4 <- function(x) {
  z <- toupper(gsub("[^A-Z0-9]", "", trimws(as.character(x))))
  substr(z, 1L, 4L)
}

map_sex <- function(x) {
  code <- suppressWarnings(as.integer(trimws(as.character(x))))
  fifelse(code == 1L, "male", fifelse(code == 3L, "female", "unknown"))
}

md5_json <- function(x) {
  encoded <- as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = FALSE))
  temporary <- tempfile(pattern = "clinical-json-md5-")
  on.exit(unlink(temporary), add = TRUE)
  writeBin(charToRaw(enc2utf8(encoded)), temporary)
  unname(tools::md5sum(temporary))
}

cache_source_paths <- function() {
  paths <- c(args$engine, args$analysis_plan, runner_script_path)
  if (args$country %in% c("brazil", "both")) paths <- c(paths, args$brazil_aggregate, args$calendar)
  if (args$country %in% c("mexico", "both")) {
    paths <- c(paths, args$mx_calendar)
    if (args$mx_cohort == "complete") {
      paths <- c(paths, args$mexico_aggregate)
    } else {
      cache_files <- list.files(dirname(args$mexico_aggregate), pattern = "^mx_clinical_year_[0-9]{4}\\.csv\\.gz$", full.names = TRUE)
      cache_years <- suppressWarnings(as.integer(sub("^.*mx_clinical_year_([0-9]{4})\\.csv\\.gz$", "\\1", cache_files)))
      needed_years <- c(2015:2018, 2020:2023)
      if (all(needed_years %in% cache_years)) {
        paths <- c(paths, cache_files[cache_years %in% c(2015:2019, 2020:2023)])
      } else if (file.exists(args$mexico_aggregate)) {
        paths <- c(paths, args$mexico_aggregate)
      } else {
        stop("Mexico source files are insufficient to establish the cache fingerprint", call. = FALSE)
      }
    }
  }
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = TRUE))
  paths
}

initialize_cache_manifest <- function(analysis_plan, weight_signatures) {
  paths <- cache_source_paths()
  source_hashes <- tools::md5sum(paths)
  source_hashes <- stats::setNames(as.list(unname(source_hashes)), names(source_hashes))
  config <- list(
    root = args$root, outdir = args$outdir, country = args$country, mx_cohort = args$mx_cohort, stage = args$stage,
    batch_index = args$batch_index, batch_count = args$batch_count,
    include_latebilling = args$include_latebilling,
    max_units = if (is.finite(args$max_units)) as.integer(args$max_units) else "all",
    inputs = list(brazil_aggregate = args$brazil_aggregate, mexico_aggregate = args$mexico_aggregate,
                  brazil_calendar = args$calendar, mexico_calendar = args$mx_calendar),
    R = R.version.string,
    data_table = as.character(utils::packageVersion("data.table")),
    jsonlite = as.character(utils::packageVersion("jsonlite"))
  )
  reference_md5 <- md5_json(analysis_plan$reference)
  config_md5 <- md5_json(config)
  signature_core <- list(source_hashes = source_hashes, reference_md5 = reference_md5,
                         config_md5 = config_md5, weight_signatures = weight_signatures)
  manifest <- c(signature_core, list(signature_md5 = md5_json(signature_core), config = config,
                                     analysis_plan_reference = analysis_plan$reference))
  manifest_path <- file.path(args$outdir, "cache", "manifest.json")
  write_json_atomic <- function(value, path, prefix) {
    temporary <- tempfile(pattern = prefix, tmpdir = dirname(path))
    on.exit(unlink(temporary), add = TRUE)
    writeLines(as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null", na = "null",
                                              digits = NA, pretty = TRUE)),
               temporary, useBytes = TRUE)
    atomic_replace_file(temporary, path)
  }
  sha256_file <- function(path) {
    command <- Sys.which("sha256sum")
    command_args <- if (nzchar(command)) shQuote(path) else character()
    if (!nzchar(command)) {
      command <- Sys.which("shasum")
      if (nzchar(command)) command_args <- c("-a", "256", shQuote(path))
    }
    if (!nzchar(command)) stop("sha256sum or shasum is required to audit cache compatibility", call. = FALSE)
    output <- suppressWarnings(system2(command, args = command_args, stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L || !length(output)) {
      stop("Unable to calculate SHA-256 for I/O cache compatibility document", call. = FALSE)
    }
    sub("[[:space:]].*$", "", output[[1L]])
  }
  if (file.exists(manifest_path)) {
    existing <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(existing)) {
      stop(sprintf("Output cache fingerprint mismatch at %s; choose a fresh output directory", manifest_path), call. = FALSE)
    }
    if (identical(existing$signature_md5, manifest$signature_md5) && !is.null(args$io_cache_compatibility)) {
      audit_path <- file.path(dirname(manifest_path), "manifest_io_compatibility_audit.json")
      if (file.exists(audit_path)) {
        pending_audit <- tryCatch(jsonlite::fromJSON(audit_path, simplifyVector = FALSE), error = function(e) NULL)
        if (!is.null(pending_audit) && identical(pending_audit$status, "verified_pending_manifest_commit") &&
            identical(pending_audit$new$signature_md5, manifest$signature_md5)) {
          pending_audit$status <- "manifest_committed"
          pending_audit$committed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
          write_json_atomic(pending_audit, audit_path, ".manifest-compat-audit.tmp-")
          log_event("clinical_cache_io_compatibility_audit_finalized", manifest = manifest_path,
                    audit = audit_path, signature_md5 = manifest$signature_md5)
        }
      }
    }
    if (!identical(existing$signature_md5, manifest$signature_md5)) {
      if (is.null(args$io_cache_compatibility)) {
        stop(sprintf("Output cache fingerprint mismatch at %s; choose a fresh output directory", manifest_path), call. = FALSE)
      }
      if (!file.exists(args$io_cache_compatibility)) {
        stop(sprintf("I/O cache compatibility document not found: %s", args$io_cache_compatibility), call. = FALSE)
      }
      compatibility <- tryCatch(
        jsonlite::fromJSON(args$io_cache_compatibility, simplifyVector = FALSE),
        error = function(e) stop(sprintf("Unable to read I/O cache compatibility document: %s", conditionMessage(e)), call. = FALSE)
      )
      required_compatibility_fields <- c("old_runner_md5", "new_runner_md5")
      if (!all(required_compatibility_fields %in% names(compatibility))) {
        stop("I/O cache compatibility document must specify old_runner_md5 and new_runner_md5", call. = FALSE)
      }

      hash_map <- function(x, description) {
        if (is.null(x) || !length(x)) return(stats::setNames(character(), character()))
        values <- unlist(x, recursive = TRUE, use.names = TRUE)
        if (is.null(names(values)) || any(!nzchar(names(values)))) {
          stop(sprintf("Invalid %s in cache manifest", description), call. = FALSE)
        }
        values <- stats::setNames(as.character(values), names(values))
        values[order(names(values))]
      }
      weight_md5_map <- function(x, description) {
        if (is.null(x) || !length(x)) return(stats::setNames(character(), character()))
        countries <- names(x)
        if (is.null(countries) || any(!nzchar(countries))) {
          stop(sprintf("Invalid %s in cache manifest", description), call. = FALSE)
        }
        values <- vapply(countries, function(country) {
          item <- x[[country]]
          if (is.null(item$md5) || length(item$md5) != 1L) {
            stop(sprintf("Missing weight-signature MD5 for %s in cache manifest", country), call. = FALSE)
          }
          as.character(item$md5)
        }, character(1))
        values[order(names(values))]
      }
      old_source_hashes <- hash_map(existing$source_hashes, "source_hashes")
      new_source_hashes <- hash_map(manifest$source_hashes, "source_hashes")
      runner_path <- normalizePath(runner_script_path, winslash = "/", mustWork = TRUE)
      if (!runner_path %in% names(old_source_hashes) || !runner_path %in% names(new_source_hashes)) {
        stop(sprintf("Runner path is not present in both cache manifests: %s", runner_path), call. = FALSE)
      }
      old_runner_md5 <- unname(old_source_hashes[[runner_path]])
      new_runner_md5 <- unname(new_source_hashes[[runner_path]])
      if (!identical(old_runner_md5, as.character(compatibility$old_runner_md5)) ||
          !identical(new_runner_md5, as.character(compatibility$new_runner_md5)) ||
          !identical(new_runner_md5, unname(tools::md5sum(runner_path)))) {
        stop("I/O cache compatibility document does not pin the existing and current runner hashes", call. = FALSE)
      }
      old_other_sources <- old_source_hashes[names(old_source_hashes) != runner_path]
      new_other_sources <- new_source_hashes[names(new_source_hashes) != runner_path]
      if (!identical(old_other_sources, new_other_sources)) {
        stop("I/O cache compatibility rejected: non-runner source hashes differ", call. = FALSE)
      }
      if (!identical(existing$config_md5, manifest$config_md5)) {
        stop("I/O cache compatibility rejected: configuration hashes differ", call. = FALSE)
      }
      if (!identical(existing$reference_md5, manifest$reference_md5)) {
        stop("I/O cache compatibility rejected: analysis reference hashes differ", call. = FALSE)
      }
      old_weight_md5 <- weight_md5_map(existing$weight_signatures, "weight_signatures")
      new_weight_md5 <- weight_md5_map(manifest$weight_signatures, "weight_signatures")
      if (!identical(old_weight_md5, new_weight_md5)) {
        stop("I/O cache compatibility rejected: case-mix weight-signature hashes differ", call. = FALSE)
      }

      backup_path <- file.path(dirname(manifest_path), "manifest.pre_io_fix.json")
      if (file.exists(backup_path)) {
        if (!identical(unname(tools::md5sum(backup_path)), unname(tools::md5sum(manifest_path)))) {
          stop(sprintf("Existing pre-I/O-fix manifest differs from current manifest: %s", backup_path), call. = FALSE)
        }
      } else if (!file.copy(manifest_path, backup_path, overwrite = FALSE, copy.mode = TRUE)) {
        stop(sprintf("Unable to preserve original cache manifest at %s", backup_path), call. = FALSE)
      }

      compatibility_sha256 <- sha256_file(args$io_cache_compatibility)
      audit <- list(
        status = "verified_pending_manifest_commit",
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        compatibility_document = args$io_cache_compatibility,
        compatibility_document_sha256 = compatibility_sha256,
        runner_path = runner_path,
        old = list(signature_md5 = existing$signature_md5, source_hashes = as.list(old_source_hashes),
                   config_md5 = existing$config_md5, reference_md5 = existing$reference_md5,
                   weight_signature_md5 = as.list(old_weight_md5)),
        new = list(signature_md5 = manifest$signature_md5, source_hashes = as.list(new_source_hashes),
                   config_md5 = manifest$config_md5, reference_md5 = manifest$reference_md5,
                   weight_signature_md5 = as.list(new_weight_md5)),
        verified_invariants = list(
          non_runner_source_hashes_identical = TRUE,
          configuration_md5_identical = TRUE,
          reference_md5_identical = TRUE,
          weight_signature_md5_identical = TRUE
        ),
        preserved_manifest = backup_path
      )
      audit_path <- file.path(dirname(manifest_path), "manifest_io_compatibility_audit.json")
      write_json_atomic(audit, audit_path, ".manifest-compat-audit.tmp-")
      write_json_atomic(manifest, manifest_path, ".manifest.json.tmp-")
      audit$status <- "manifest_committed"
      audit$committed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
      write_json_atomic(audit, audit_path, ".manifest-compat-audit.tmp-")
      log_event("clinical_cache_io_compatibility_applied", manifest = manifest_path,
                compatibility_document = args$io_cache_compatibility,
                compatibility_document_sha256 = compatibility_sha256,
                old_signature_md5 = existing$signature_md5, new_signature_md5 = manifest$signature_md5)
    }
  } else {
    cached_files <- list.files(file.path(args$outdir, "cache"), recursive = TRUE, full.names = TRUE, all.files = TRUE)
    cached_files <- cached_files[!dir.exists(cached_files)]
    if (length(cached_files)) {
      stop(sprintf("Unbound cache files already exist under %s; choose a fresh output directory", file.path(args$outdir, "cache")), call. = FALSE)
    }
    temporary <- tempfile(pattern = ".manifest.json.tmp-", tmpdir = file.path(args$outdir, "cache"))
    on.exit(unlink(temporary), add = TRUE)
    writeLines(as.character(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null", na = "null", digits = NA, pretty = TRUE)),
               temporary, useBytes = TRUE)
    atomic_replace_file(temporary, manifest_path)
  }
  manifest
}

atomic_replace_file <- function(temporary, path) {
  if (file.rename(temporary, path)) return(invisible(path))
  if (!file.exists(path)) stop(sprintf("Unable to move temporary file into place: %s", path), call. = FALSE)
  backup <- tempfile(pattern = paste0(".", basename(path), ".backup-"), tmpdir = dirname(path))
  if (!file.rename(path, backup)) stop(sprintf("Unable to replace existing file: %s", path), call. = FALSE)
  if (!file.rename(temporary, path)) {
    file.rename(backup, path)
    stop(sprintf("Unable to move temporary file into place: %s", path), call. = FALSE)
  }
  unlink(backup)
  invisible(path)
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  data.table::fwrite(x, temporary, na = "")
  atomic_replace_file(temporary, path)
}

cache_failure_path <- function(path) sub("\\.csv$", "_failures.csv", path)

cache_read_rows <- function(path) {
  if (!file.exists(path)) return(data.table())
  # Cache rows contain JSON text columns. In the installed data.table version,
  # fread leaves doubled CSV quotes inside these values; each rewrite then
  # escapes them again and grows the cache exponentially. Base read.csv
  # round-trips the quoted JSON correctly and avoids mmap on very large caches.
  z <- utils::read.csv(path, colClasses = c(unit = "character"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  z <- data.table::as.data.table(z)
  if (!all(c("unit", "status") %in% names(z))) {
    stop(sprintf("Invalid per-unit cache schema in %s; use a fresh output directory", path), call. = FALSE)
  }
  z
}

cache_rewrite_without_unit <- function(path, unit_id) {
  if (!file.exists(path)) return(invisible(NULL))
  z <- cache_read_rows(path)
  keep <- z[which(z[["unit"]] != unit_id)]
  if (nrow(keep)) atomic_fwrite(keep, path) else unlink(path)
  invisible(NULL)
}

cache_write_to <- function(row, path) {
  old <- cache_read_rows(path)
  new <- data.table::rbindlist(list(old, row), fill = TRUE, use.names = TRUE)
  unit_id <- as.character(row[["unit"]][1L])
  new <- new[which(new[["unit"]] != unit_id)]
  new <- data.table::rbindlist(list(new, row), fill = TRUE, use.names = TRUE)
  if (nrow(new)) data.table::setorderv(new, "unit")
  atomic_fwrite(new, path)
}

cache_write_row <- function(row, path) {
  row <- data.table::as.data.table(data.table::copy(row))
  if (!all(c("unit", "status") %in% names(row)) || nrow(row) != 1L) {
    stop("Per-unit cache rows must contain one unit and a status", call. = FALSE)
  }
  unit_id <- as.character(row[["unit"]][1L])
  failure_path <- cache_failure_path(path)
  is_failure <- !identical(as.character(row[["status"]][1L]), "ok")
  if (is_failure) {
    cache_write_to(row, failure_path)
    cache_rewrite_without_unit(path, unit_id)
  } else {
    cache_write_to(row, path)
    cache_rewrite_without_unit(failure_path, unit_id)
  }
}

cache_done_units <- function(path) {
  z <- cache_read_rows(path)
  if (!nrow(z)) character() else unique(z[status == "ok", as.character(unit)])
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]", "_", as.character(x))

scenario_identity <- function(scenario) {
  name <- scenario$name
  label <- scenario$label
  if (length(name) != 1L || is.na(name) || !nzchar(as.character(name)) ||
      length(label) != 1L || is.na(label) || !nzchar(as.character(label))) {
    stop("Scenario metadata must contain scalar name and label fields", call. = FALSE)
  }
  list(name = as.character(name), label = as.character(label))
}

select_unit_batch <- function(units, batch_index, batch_count, max_units = Inf) {
  unit_ids <- sort(unique(as.character(units)))
  unit_ids <- unit_ids[!is.na(unit_ids) & nzchar(unit_ids)]
  if (length(unit_ids)) {
    assigned_batch <- ((seq_along(unit_ids) - 1L) %% as.integer(batch_count)) + 1L
    unit_ids <- unit_ids[assigned_batch == as.integer(batch_index)]
  }
  if (is.finite(max_units)) unit_ids <- head(unit_ids, as.integer(max_units))
  unit_ids
}

fread_path <- function(path, ...) {
  is_gzip <- grepl("\\.gz$", path, ignore.case = TRUE)
  gzip <- if (is_gzip && .Platform$OS.type == "unix") Sys.which("gzip") else ""
  if (nzchar(gzip) && !requireNamespace("R.utils", quietly = TRUE)) {
    command <- paste(shQuote(gzip), "-dc --", shQuote(path))
    return(data.table::fread(cmd = command, ...))
  }
  data.table::fread(path, ...)
}

read_calendar <- function(path) {
  need <- c("date", "unit", "tmean", "population", "E3", "y_K80", "y_urinary")
  cal <- fread_path(path, select = need, colClasses = list(character = c("date", "unit")), showProgress = TRUE)
  missing <- setdiff(need, names(cal))
  if (length(missing)) stop(sprintf("Brazil calendar lacks columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  cal[, date := parse_date_safe(date)]
  if (anyNA(cal$date)) stop("Brazil calendar has invalid dates", call. = FALSE)
  cal[, unit := canonical_unit(unit)]
  for (nm in c("tmean", "population", "E3", "y_K80", "y_urinary")) set(cal, j = nm, value = suppressWarnings(as.numeric(cal[[nm]])))
  setorder(cal, unit, date)
  duplicate_dates <- cal[, .N, by = .(unit, date)][N > 1L]
  if (nrow(duplicate_dates)) stop(sprintf("Brazil calendar has %d duplicated unit-date pairs", nrow(duplicate_dates)), call. = FALSE)

  # Audit canonical E3 against a fresh actual-calendar lag0-2 mean before any
  # admission-year or outcome filtering. The input calendar contains one row
  # per unit-day; date-distance checks ensure gaps are not treated as adjacent.
  calendar_ma_fast <- function(x, dates) {
    n <- length(x); out <- rep(NA_real_, n)
    if (n < 3L) return(out)
    i <- 3L:n
    adjacent <- as.integer(dates[i] - dates[i - 1L]) == 1L & as.integer(dates[i - 1L] - dates[i - 2L]) == 1L
    valid <- adjacent & is.finite(x[i]) & is.finite(x[i - 1L]) & is.finite(x[i - 2L])
    out[i[valid]] <- (x[i[valid]] + x[i[valid] - 1L] + x[i[valid] - 2L]) / 3
    out
  }
  cal[, E3_rebuilt := calendar_ma_fast(tmean, date), by = unit]
  compare <- cal[is.finite(E3) & is.finite(E3_rebuilt), abs(E3 - E3_rebuilt)]
  exposure_audit <- list(
    calendar_rows = nrow(cal),
    calendar_units = uniqueN(cal$unit),
    finite_canonical_E3 = sum(is.finite(cal$E3)),
    finite_rebuilt_E3 = sum(is.finite(cal$E3_rebuilt)),
    compared_E3_rows = length(compare),
    E3_mismatch_rows_gt_1e_6 = sum(compare > 1e-6),
    E3_max_absolute_difference = if (length(compare)) max(compare) else NA_real_,
    exposure_used = "E3_rebuilt from tmean over complete actual calendar lag0-2; only finite complete windows retained"
  )
  cal[, E3_actual := E3_rebuilt]
  list(data = cal, audit = exposure_audit)
}

read_clinical <- function(path) {
  required <- c("source_year", "hrid", "admission_date", "diag_princ4", "sexo_raw", "age_group", "car_int_group", "admission_records")
  x <- fread_path(
    path,
    select = c(required, "car_int_raw"),
    colClasses = list(character = c("source_year", "hrid", "admission_date", "diag_princ4", "sexo_raw", "age_group", "car_int_group", "car_int_raw")),
    showProgress = TRUE
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(sprintf("Brazil clinical aggregate lacks columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  x[, source_year := suppressWarnings(as.integer(source_year))]
  x[, date := parse_date_safe(admission_date)]
  x[, date_year := suppressWarnings(as.integer(format(date, "%Y")))]
  x[, unit := canonical_unit(hrid)]
  x[, diag4 := normalize_dx4(diag_princ4)]
  x[, dx3 := substr(diag4, 1L, 3L)]
  x[, sex_group := map_sex(sexo_raw)]
  x[, n := suppressWarnings(as.numeric(admission_records))]
  if (any(!is.finite(x$n) | x$n < 0)) stop("Clinical aggregate contains invalid admission_records", call. = FALSE)
  x
}

read_mexico_clinical <- function(path, cohort = "baseline") {
  required <- c("source_role", "state", "date", "fullcode", "age_group", "sex", "n")
  if (!cohort %in% c("baseline", "complete")) stop("Mexico cohort must be baseline or complete", call. = FALSE)
  if (cohort == "complete") {
    x <- fread_path(path, select = c(required, "source_year"),
                    colClasses = list(character = c("source_role", "state", "date", "fullcode", "age_group", "sex", "source_year")),
                    showProgress = TRUE)
    missing <- setdiff(c(required, "source_year"), names(x))
    if (length(missing)) stop(sprintf("Complete Mexico source lacks columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
    x[, source_year := suppressWarnings(as.integer(source_year))]
  } else {
    cache_files <- list.files(dirname(path), pattern = "^mx_clinical_year_[0-9]{4}\\.csv\\.gz$", full.names = TRUE)
    cache_years <- suppressWarnings(as.integer(sub("^.*mx_clinical_year_([0-9]{4})\\.csv\\.gz$", "\\1", cache_files)))
    needed_years <- c(2015:2018, 2020:2023)
    if (all(needed_years %in% cache_years)) {
      cache_files <- cache_files[cache_years %in% c(2015:2019, 2020:2023)]
      cache_years <- cache_years[cache_years %in% c(2015:2019, 2020:2023)]
      cache_required <- setdiff(required, "source_role")
      pieces <- lapply(seq_along(cache_files), function(i) {
        d <- fread_path(
          cache_files[[i]], select = cache_required,
          colClasses = list(character = c("state", "date", "fullcode", "age_group", "sex")),
          showProgress = TRUE
        )
        missing <- setdiff(cache_required, names(d))
        if (length(missing)) stop(sprintf("Mexico per-year cache %s lacks columns: %s", basename(cache_files[[i]]), paste(missing, collapse = ", ")), call. = FALSE)
        d[, source_year := cache_years[[i]]]
        d[, event_year := suppressWarnings(as.integer(substr(date, 1L, 4L)))]
        d <- d[!is.na(event_year) & event_year >= 2015L & event_year <= 2023L]
        d[, source_role := "primary_panel_year"]
        d[source_year == 2019L, source_role := "recovered_additional_2019"]
        d[source_year == 2023L & event_year == 2023L, source_role := "extension_2023"]
        d[, event_year := NULL]
        d
      })
      x <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
    } else if (file.exists(path)) {
      x <- fread_path(path, colClasses = list(character = c("source_role", "state", "date", "fullcode", "age_group", "sex")), showProgress = TRUE)
      if (!"source_year" %in% names(x)) {
        stop(sprintf(
          "Mexico combined long lacks source_year and per-year caches are incomplete; need caches for %s",
          paste(needed_years, collapse = ",")
        ), call. = FALSE)
      }
      missing <- setdiff(c(required, "source_year"), names(x))
      if (length(missing)) stop(sprintf("Mexico clinical aggregate lacks columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
      x[, source_year := suppressWarnings(as.integer(source_year))]
    } else {
      stop(sprintf("Mexico clinical aggregate or required per-year caches not found beside: %s", path), call. = FALSE)
    }
  }
  x[, source_role := trimws(as.character(source_role))]
  x[, date := parse_date_safe(date)]
  x[, date_year := suppressWarnings(as.integer(format(date, "%Y")))]
  x[, unit := canonical_unit(state)]
  x[, diag4 := normalize_dx4(fullcode)]
  x[, dx3 := substr(diag4, 1L, 3L)]
  x[, age_group := trimws(as.character(age_group))]
  x[, sex_group := tolower(trimws(as.character(sex)))]
  x[!sex_group %in% c("male", "female"), sex_group := "unknown"]
  x[, car_int_raw := "not_available"]
  x[, car_int_group := "not_available"]
  x[, n := suppressWarnings(as.numeric(n))]
  if (any(!is.finite(x$n) | x$n < 0)) stop("Mexico clinical aggregate contains invalid n", call. = FALSE)
  x
}

source_scenarios <- function(include_latebilling = TRUE) {
  out <- list(list(
    country = "brazil", name = "matched_source_2014_2020",
    label = "matched-source replication; source years 2014-2020",
    years = 2014:2020, analysis_years = 2014:2020, offset_col = "population",
    modes = c("admission_type", "sex")
  ))
  if (include_latebilling) {
    out[[length(out) + 1L]] <- list(
      country = "brazil",
      name = "latebilling_reconstruction_2014_2021",
      label = "late-billing reconstruction; source years 2014-2021; source-year 2021 contributes late 2020 admissions",
      years = 2014:2021, analysis_years = 2014:2020, offset_col = "population",
      modes = c("admission_type", "sex")
    )
  }
  out
}

mexico_primary_scenario <- function(cohort = "baseline") {
  if (cohort == "complete") {
    return(list(
      country = "mexico", name = "mx_complete_source_2015_2022",
      label = "Complete source-year cohort 2015-2022; all source roles and admission years 2015-2022.",
      years = 2015:2022, analysis_years = 2015:2022,
      source_role = NULL, offset_col = NULL, modes = "sex"
    ))
  }
  if (cohort != "baseline") stop("Mexico cohort must be baseline or complete", call. = FALSE)
  list(
    country = "mexico", name = "mx_primary_panel_year_exact_baseline",
    label = paste0(
      "primary_panel_year only; source years 2015-2018 and 2020-2023; admission years 2015-2018 and 2020-2022; ",
      "recovered_additional_2019 excluded; source-year 2023 late billings for admission years 2021-2022 retained"
    ),
    years = c(2015:2018, 2020:2023), analysis_years = c(2015:2018, 2020:2022),
    source_role = "primary_panel_year", offset_col = NULL, modes = "sex"
  )
}

build_daily_panel <- function(clinical, calendar, scenario) {
  identity <- scenario_identity(scenario)
  date_keep <- !is.na(clinical$date) & clinical$date_year %in% scenario$analysis_years
  source_keep <- clinical$source_year %in% scenario$years
  target_keep <- clinical$dx3 %in% c("K80", "N20", "N21", "N22", "N23")
  role_keep <- if (is.null(scenario$source_role) || !"source_role" %in% names(clinical)) rep(TRUE, nrow(clinical)) else clinical$source_role == scenario$source_role
  excluded_source_year_records <- clinical[date_keep & !source_keep & target_keep, sum(n)]
  excluded_source_role_records <- if ("source_role" %in% names(clinical) && !is.null(scenario$source_role)) {
    clinical[date_keep & source_keep & target_keep & !role_keep, sum(n)]
  } else 0
  selected <- clinical[date_keep & source_keep & target_keep & role_keep]
  matched_units <- unique(calendar$unit)
  unknown_unit <- selected$unit %in% c("unknown", "")
  matched_unit <- selected$unit %in% matched_units
  dropped_unmapped <- selected[unknown_unit, sum(n)]
  dropped_unknown_unit <- selected[!unknown_unit & !matched_unit, sum(n)]
  selected <- selected[!unknown_unit & matched_unit]
  if (!nrow(selected)) stop(sprintf("No target clinical rows match the Brazil calendar in scenario %s", scenario$name), call. = FALSE)

  selected[, car_code_group := trimws(as.character(car_int_group))]
  selected[, y_K80_all := fifelse(dx3 == "K80", n, 0)]
  selected[, y_K80_acute_mix := fifelse(diag4 %in% c("K800", "K803"), n, 0)]
  selected[, y_N20_only := fifelse(dx3 == "N20", n, 0)]
  selected[, y_N20_N21 := fifelse(dx3 %in% c("N20", "N21"), n, 0)]
  selected[, y_N23 := fifelse(dx3 == "N23", n, 0)]
  selected[, y_N20_N23 := fifelse(dx3 %in% c("N20", "N21", "N22", "N23"), n, 0)]
  selected[, y_K80_elective := fifelse(dx3 == "K80" & car_code_group == "01", n, 0)]
  selected[, y_N20_elective := fifelse(dx3 == "N20" & car_code_group == "01", n, 0)]
  selected[, y_K80_urgent := fifelse(dx3 == "K80" & car_code_group == "02", n, 0)]
  selected[, y_N20_urgent := fifelse(dx3 == "N20" & car_code_group == "02", n, 0)]
  selected[, y_K80_male := fifelse(dx3 == "K80" & sex_group == "male", n, 0)]
  selected[, y_N20_male := fifelse(dx3 == "N20" & sex_group == "male", n, 0)]
  selected[, y_K80_female := fifelse(dx3 == "K80" & sex_group == "female", n, 0)]
  selected[, y_N20_female := fifelse(dx3 == "N20" & sex_group == "female", n, 0)]
  model_outcomes <- c(
    "y_K80_all", "y_K80_acute_mix", "y_N20_only", "y_N20_N21", "y_N23", "y_N20_N23",
    "y_K80_elective", "y_N20_elective", "y_K80_urgent", "y_N20_urgent",
    "y_K80_male", "y_N20_male", "y_K80_female", "y_N20_female"
  )
  daily <- selected[, lapply(.SD, sum, na.rm = TRUE), by = .(unit, date), .SDcols = model_outcomes]

  calendar_keep <- format(calendar$date, "%Y") %in% as.character(scenario$analysis_years) & is.finite(calendar$E3_actual)
  if (!is.null(scenario$offset_col)) calendar_keep <- calendar_keep & is.finite(calendar$population) & calendar$population > 0
  calendar_window <- calendar[
    calendar_keep,
    .(unit, date, E3_actual, population, legacy_y_K80 = y_K80, legacy_y_urinary = y_urinary)
  ]
  setorder(calendar_window, unit, date)
  panel <- merge(calendar_window, daily, by = c("unit", "date"), all.x = TRUE, sort = FALSE)
  for (nm in model_outcomes) set(panel, which(is.na(panel[[nm]])), nm, 0)
  panel[, E3 := E3_actual]
  setorder(panel, unit, date)
  latebill_2020_available <- if (scenario$country == "brazil") clinical[
    source_year == 2021L & date_year == 2020L & dx3 %in% c("K80", "N20", "N21", "N22", "N23"), sum(n)
  ] else 0
  source_2023_original_dates <- selected[source_year == 2023L, sum(n)]
  recovered_2019_excluded <- recovered_2019_included <- extension_2023_excluded <- extension_2023_included <- 0
  if ("source_role" %in% names(clinical)) {
    cohort_selected <- date_keep & source_keep & role_keep
    recovered_rows <- clinical$source_role == "recovered_additional_2019" & target_keep
    extension_rows <- clinical$source_role == "extension_2023" & target_keep
    recovered_2019_excluded <- clinical[recovered_rows & !cohort_selected, sum(n)]
    recovered_2019_included <- clinical[recovered_rows & cohort_selected, sum(n)]
    extension_2023_excluded <- clinical[extension_rows & !cohort_selected, sum(n)]
    extension_2023_included <- clinical[extension_rows & cohort_selected, sum(n)]
  }
  scenario_summary <- list(
    country = scenario$country,
    scenario = scenario$name,
    name = identity$name,
    source_years = as.integer(scenario$years),
    analysis_years = as.integer(scenario$analysis_years),
    scenario_label = scenario$label,
    label = identity$label,
    source_role = scenario$source_role,
    source_role_description = if (is.null(scenario$source_role)) "all available source roles" else scenario$source_role,
    offset_col = scenario$offset_col,
    available_interaction_modes = scenario$modes,
    target_records_excluded_by_source_year = excluded_source_year_records,
    target_records_excluded_by_source_role = excluded_source_role_records,
    recovered_additional_2019_target_records_excluded = recovered_2019_excluded,
    recovered_additional_2019_target_records_included = recovered_2019_included,
    extension_2023_target_records_excluded = extension_2023_excluded,
    extension_2023_target_records_included = extension_2023_included,
    selected_diagnosis_records_before_calendar_match = selected[, sum(n)] + dropped_unmapped + dropped_unknown_unit,
    selected_diagnosis_records_retained = selected[, sum(n)],
    selected_rows_unknown_hrid = dropped_unmapped,
    selected_rows_hrid_not_in_calendar = dropped_unknown_unit,
    source_year_rows = as.list(selected[, .(records = sum(n)), by = source_year][order(source_year)]),
    source_2021_admission_2020_rows_available = latebill_2020_available,
    source_2021_admission_2020_rows_included = if (2021L %in% scenario$years) latebill_2020_available else 0,
    source_2023_rows_included_for_original_admission_dates = source_2023_original_dates,
    panel_units = uniqueN(panel$unit),
    panel_unit_dates = nrow(panel),
    panel_dates = uniqueN(panel$date),
    panel_dates_dropped_for_incomplete_weather_or_population = calendar[
      format(date, "%Y") %in% as.character(scenario$analysis_years) &
        (!is.finite(E3_actual) | (!is.null(scenario$offset_col) & (!is.finite(population) | population <= 0))), .N
    ]
  )
  list(panel = panel, selected = selected, audit = scenario_summary,
       offset_col = scenario$offset_col, modes = scenario$modes)
}

unit_weather_target <- function(d) {
  e <- d$E3[is.finite(d$E3)]
  if (length(e) < 400L || length(unique(e)) < 2L) stop("fewer than 400 complete exposure days or fewer than two E3 values", call. = FALSE)
  list(
    basis = new_temperature_basis(e, df = 3L),
    p95 = unname(stats::quantile(e, 0.95, type = 7)),
    median = stats::median(e),
    exposure_days = length(e)
  )
}

covariance_fields <- function(cc) {
  methods <- c("hac14", "hac28", "day_hc0", "quasi")
  out <- list()
  for (nm in methods) {
    z <- cc$by_covariance[[nm]]
    if (is.null(z)) {
      out[[paste0("logRRR_", nm)]] <- NA_real_
      out[[paste0("var_", nm)]] <- NA_real_
      out[[paste0("se_", nm)]] <- NA_real_
      out[[paste0("RR_", nm)]] <- NA_real_
      out[[paste0("lo_", nm)]] <- NA_real_
      out[[paste0("hi_", nm)]] <- NA_real_
    } else {
      out[[paste0("logRRR_", nm)]] <- unname(z["estimate"])
      out[[paste0("var_", nm)]] <- unname(z["variance"])
      out[[paste0("se_", nm)]] <- unname(z["se"])
      out[[paste0("RR_", nm)]] <- unname(z["RR"])
      out[[paste0("lo_", nm)]] <- unname(z["lo"])
      out[[paste0("hi_", nm)]] <- unname(z["hi"])
    }
  }
  out
}

endpoint_specs <- list(
  strict_N20_only_vs_K80 = list(denominator = "y_K80_all", numerator = "y_N20_only", label = "N20-only versus K80"),
  confirmed_N20_N21_vs_K80 = list(denominator = "y_K80_all", numerator = "y_N20_N21", label = "N20+N21 versus K80"),
  composition_N23_vs_K80 = list(denominator = "y_K80_all", numerator = "y_N23", label = "N23 renal-colic composition versus K80"),
  N20_vs_K80_0_acute_cholecystitis_plus_K80_3_cholangitis = list(
    denominator = "y_K80_acute_mix", numerator = "y_N20_only",
    label = "N20 versus K80.0 cholecystitis plus K80.3 cholangitis"
  ),
  combined_original_N20_N23_vs_K80 = list(denominator = "y_K80_all", numerator = "y_N20_N23", label = "combined original N20-N23 versus K80")
)

fit_endpoint_unit <- function(d, scenario, nuisance, endpoint, spec, unit, offset_col) {
  identity <- scenario_identity(scenario)
  tryCatch({
    target <- unit_weather_target(d)
    event_counts <- c(denominator = sum(d[[spec$denominator]], na.rm = TRUE),
                      numerator = sum(d[[spec$numerator]], na.rm = TRUE))
    if (any(!is.finite(event_counts) | event_counts <= 0)) {
      empty <- names(event_counts)[!is.finite(event_counts) | event_counts <= 0]
      stop(sprintf("endpoint has zero or invalid total admissions: %s", paste(empty, collapse = ",")), call. = FALSE)
    }
    fit <- fit_stacked_contrast(
      data = d,
      outcomes = c(spec$denominator, spec$numerator),
      exposure = "E3", date = "date", offset = offset_col,
      basis = target$basis, p95 = target$p95, median = target$median,
      nuisance = nuisance, hac_lags = c(14L, 28L), full_covariance = FALSE
    )
    pair_label <- paste0(spec$numerator, "_vs_", spec$denominator)
    cc <- fit$contrasts[[pair_label]]$coefficient_contrast
    if (is.null(cc)) stop(sprintf("expected contrast %s not returned", pair_label), call. = FALSE)
    values <- covariance_fields(cc)
    row <- data.table(
      country = scenario$country, scenario = identity$name, scenario_label = identity$label, nuisance = nuisance,
      endpoint = endpoint, endpoint_label = spec$label, unit = unit,
      status = "ok", n_days = nrow(d), events_denominator = sum(d[[spec$denominator]]),
      events_numerator = sum(d[[spec$numerator]]), p95 = target$p95, median = target$median,
      exposure_days = target$exposure_days, formula = fit$audit$formula,
      warnings = as.character(jsonlite::toJSON(fit$audit$warnings, auto_unbox = TRUE))
    )
    for (nm in names(values)) set(row, j = nm, value = values[[nm]])
    row
  }, error = function(e) {
    data.table(country = scenario$country, scenario = identity$name, scenario_label = identity$label, nuisance = nuisance,
               endpoint = endpoint, endpoint_label = spec$label, unit = unit, status = "failed",
               n_days = nrow(d), events_denominator = if (spec$denominator %in% names(d)) sum(d[[spec$denominator]]) else NA_real_,
               events_numerator = if (spec$numerator %in% names(d)) sum(d[[spec$numerator]]) else NA_real_,
               p95 = NA_real_, median = NA_real_, exposure_days = NA_integer_, formula = NA_character_,
               warnings = NA_character_, failure = conditionMessage(e))
  })
}

run_endpoint_stage <- function(panels, units, force = FALSE) {
  all_rows <- list()
  failure_rows <- list()
  combos <- 0L
  for (scenario in panels) {
    scenario_identity(scenario$audit)
    panel <- scenario$panel
    for (nuisance in c("separate", "shared")) {
      for (endpoint in names(endpoint_specs)) {
        spec <- endpoint_specs[[endpoint]]
        path <- file.path(args$outdir, "cache", sprintf("endpoint_%s_%s_%s_%s.csv", scenario$audit$country, scenario$audit$scenario, nuisance, endpoint))
        if (force && file.exists(path)) unlink(path)
        done <- cache_done_units(path)
        for (u in units) {
          if (u %in% done) next
          d <- panel[unit == u]
          row <- fit_endpoint_unit(d, scenario$audit, nuisance, endpoint, spec, u, scenario$offset_col)
          cache_write_row(row, path)
        }
        if (file.exists(path)) all_rows[[length(all_rows) + 1L]] <- cache_read_rows(path)
        failure_path <- cache_failure_path(path)
        if (file.exists(failure_path)) failure_rows[[length(failure_rows) + 1L]] <- cache_read_rows(failure_path)
        combos <- combos + 1L
        log_event("endpoint_model_complete", country = scenario$audit$country, scenario = scenario$audit$scenario, nuisance = nuisance,
                  endpoint = endpoint, units_cached = length(cache_done_units(path)))
      }
    }
  }
  per_unit <- rbindlist(all_rows, fill = TRUE, use.names = TRUE)
  if (is.null(per_unit) || !ncol(per_unit)) {
    per_unit <- data.table(country = character(), scenario = character(), scenario_label = character(), nuisance = character(),
                           endpoint = character(), endpoint_label = character(), unit = character(), status = character())
  }
  if (nrow(per_unit)) {
    key <- c("scenario", "nuisance", "endpoint", "unit")
    setorderv(per_unit, key)
    per_unit <- unique(per_unit, by = key, fromLast = TRUE)
  }
  countries <- if (nrow(per_unit)) paste(sort(unique(per_unit$country)), collapse = "_") else
    paste(sort(unique(vapply(panels, function(x) x$audit$country, character(1)))), collapse = "_")
  fwrite(per_unit, file.path(args$outdir, paste0("clinical_", countries, "_endpoint_per_unit.csv")), na = "")
  summary_rows <- list()
  covariance_methods <- c("hac14", "hac28", "day_hc0", "quasi")
  group_keys <- c("country", "scenario", "scenario_label", "nuisance", "endpoint", "endpoint_label")
  groups <- unique(per_unit[, ..group_keys])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- per_unit[country == g$country & scenario == g$scenario & nuisance == g$nuisance & endpoint == g$endpoint & status == "ok"]
    for (method in covariance_methods) {
      estimates <- z[[paste0("logRRR_", method)]]
      variances <- z[[paste0("var_", method)]]
      for (pool_method in c("DL", "REML_HK")) {
        pooled <- pool_meta(estimates, variances, method = pool_method)
        summary_rows[[length(summary_rows) + 1L]] <- data.table(
          country = g$country, scenario = g$scenario, scenario_label = g$scenario_label, nuisance = g$nuisance,
          endpoint = g$endpoint, endpoint_label = g$endpoint_label, covariance = method,
          pooling = pool_method, status = if (is.null(pooled$failure)) "ok" else "partial_or_failed",
          units_fit = nrow(z), k = pooled$k, logRRR = pooled$estimate, se = pooled$se,
          RRR = pooled$RR, RR_lo = pooled$RR95CI[1L], RR_hi = pooled$RR95CI[2L],
          I2 = pooled$I2, tau2 = pooled$tau2, concordance = pooled$concordance,
          failure = if (is.null(pooled$failure)) NA_character_ else pooled$failure
        )
      }
    }
  }
  pooled_dt <- rbindlist(summary_rows, fill = TRUE)
  if (is.null(pooled_dt) || !ncol(pooled_dt)) {
    pooled_dt <- data.table(country = character(), scenario = character(), scenario_label = character(), nuisance = character(),
                            endpoint = character(), endpoint_label = character(), covariance = character(), pooling = character(),
                            status = character(), units_fit = integer(), k = integer(), logRRR = numeric(), se = numeric(),
                            RRR = numeric(), RR_lo = numeric(), RR_hi = numeric(), I2 = numeric(), tau2 = numeric(),
                            concordance = numeric(), failure = character())
  }
  fwrite(pooled_dt, file.path(args$outdir, paste0("clinical_", countries, "_endpoint_pooled.csv")), na = "")
  failures_dt <- rbindlist(failure_rows, fill = TRUE, use.names = TRUE)
  if (is.null(failures_dt) || !ncol(failures_dt)) {
    failures_dt <- data.table(country = character(), scenario = character(), scenario_label = character(), nuisance = character(),
                              endpoint = character(), endpoint_label = character(), unit = character(), status = character(), failure = character())
  }
  fwrite(failures_dt, file.path(args$outdir, paste0("clinical_", countries, "_endpoint_failures.csv")), na = "")
  list(per_unit = per_unit, pooled = pooled_dt, combinations = combos)
}

design_at_endpoint <- function(fit_result, endpoint, temp) {
  .design_at(fit_result$fit, fit_result$data, endpoint = endpoint,
             basis = fit_result$basis, temp_value = temp)
}

pure_basis_coef <- function(fit_result, endpoint, basis_name) {
  fit <- fit_result$fit
  nd <- as.data.frame(fit_result$data[1L, , drop = FALSE])
  nd$disease <- factor(endpoint, levels = levels(fit_result$data$disease))
  nd$B1 <- 0; nd$B2 <- 0; nd$B3 <- 0
  tt <- stats::delete.response(stats::terms(fit))
  x0 <- stats::model.matrix(tt, data = nd, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  nd[[basis_name]] <- 1
  x1 <- stats::model.matrix(tt, data = nd, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  nm <- names(stats::coef(fit))
  as.numeric(x1[, nm, drop = FALSE] - x0[, nm, drop = FALSE])
}

linear_effect <- function(L, fit, V) {
  beta <- stats::coef(fit)
  L <- as.numeric(L); names(L) <- names(beta)
  est <- as.numeric(L %*% beta)
  variance <- as.numeric(t(L) %*% V %*% L)
  if (!is.finite(variance) || variance < 0) stop("linear contrast has invalid HAC variance", call. = FALSE)
  se <- sqrt(variance)
  c(estimate = est, variance = variance, se = se, RR = exp(est),
    lo = exp(est - qnorm(0.975) * se), hi = exp(est + qnorm(0.975) * se))
}

interaction_results <- function(fit_result, endpoint_names, high_group, low_group, mode) {
  # endpoint_names are K80/N20-only split across elective/urgent or sex.
  target <- list(p95 = fit_result$p95, median = fit_result$median)
  by_group <- list()
  for (g in c(low_group, high_group)) {
    k <- endpoint_names[[paste0("K80_", g)]]
    u <- endpoint_names[[paste0("N20_", g)]]
    Lk <- as.numeric(design_at_endpoint(fit_result, k, target$p95) - design_at_endpoint(fit_result, k, target$median))
    Lu <- as.numeric(design_at_endpoint(fit_result, u, target$p95) - design_at_endpoint(fit_result, u, target$median))
    by_group[[g]] <- Lu - Lk
  }
  L_high_low <- by_group[[high_group]] - by_group[[low_group]]
  V_hac14 <- fit_result$hac_vcov$hac14
  scalar <- linear_effect(L_high_low, fit_result$fit, V_hac14)

  component_vectors <- lapply(c("B1", "B2", "B3"), function(bn) {
    diff_group <- function(g) {
      k <- endpoint_names[[paste0("K80_", g)]]
      u <- endpoint_names[[paste0("N20_", g)]]
      pure_basis_coef(fit_result, u, bn) - pure_basis_coef(fit_result, k, bn)
    }
    diff_group(high_group) - diff_group(low_group)
  })
  Lmat <- do.call(rbind, component_vectors)
  rownames(Lmat) <- c("B1", "B2", "B3")
  colnames(Lmat) <- names(stats::coef(fit_result$fit))
  beta <- stats::coef(fit_result$fit)
  basis_est <- as.numeric(Lmat %*% beta)
  basis_cov <- Lmat %*% V_hac14 %*% t(Lmat)
  wald <- tryCatch(as.numeric(t(basis_est) %*% solve(basis_cov, basis_est)), error = function(e) NA_real_)
  p_wald <- if (is.finite(wald)) stats::pchisq(wald, df = 3L, lower.tail = FALSE) else NA_real_

  group_estimates <- lapply(names(by_group), function(g) {
    z <- linear_effect(by_group[[g]], fit_result$fit, V_hac14)
    c(group = g, estimate = z["estimate"], variance = z["variance"], se = z["se"],
      RRR = z["RR"], lo = z["lo"], hi = z["hi"])
  })
  list(
    mode = mode,
    logRRR_difference = unname(scalar["estimate"]),
    variance_hac14 = unname(scalar["variance"]),
    se_hac14 = unname(scalar["se"]),
    RR_ratio_hac14 = unname(scalar["RR"]),
    RR_ratio_lo = unname(scalar["lo"]), RR_ratio_hi = unname(scalar["hi"]),
    joint_wald_chisq3 = wald, joint_wald_p = p_wald,
    joint_basis_coefficients = as.list(setNames(basis_est, rownames(Lmat))),
    joint_basis_covariance = unclass(basis_cov),
    groups = group_estimates
  )
}

fit_four_endpoint_unit <- function(d, scenario, nuisance, mode, unit, offset_col) {
  identity <- scenario_identity(scenario)
  tryCatch({
    target <- unit_weather_target(d)
    if (mode == "admission_type") {
      outcomes <- c("y_K80_elective", "y_N20_elective", "y_K80_urgent", "y_N20_urgent")
      endpoints <- list(K80_elective = "y_K80_elective", N20_elective = "y_N20_elective",
                        K80_urgent = "y_K80_urgent", N20_urgent = "y_N20_urgent")
      high <- "urgent"; low <- "elective"
    } else {
      outcomes <- c("y_K80_male", "y_N20_male", "y_K80_female", "y_N20_female")
      endpoints <- list(K80_male = "y_K80_male", N20_male = "y_N20_male",
                        K80_female = "y_K80_female", N20_female = "y_N20_female")
      high <- "male"; low <- "female"
    }
    event_counts <- vapply(outcomes, function(nm) sum(d[[nm]], na.rm = TRUE), numeric(1))
    if (any(!is.finite(event_counts) | event_counts <= 0)) {
      empty <- names(event_counts)[!is.finite(event_counts) | event_counts <= 0]
      stop(sprintf("interaction endpoint has zero or invalid total admissions: %s", paste(empty, collapse = ",")), call. = FALSE)
    }
    fit <- fit_stacked_contrast(
      data = d, outcomes = outcomes, exposure = "E3", date = "date", offset = offset_col,
      basis = target$basis, p95 = target$p95, median = target$median, nuisance = nuisance,
      hac_lags = c(14L, 28L), full_covariance = TRUE
    )
    fit$p95 <- target$p95; fit$median <- target$median
    effect <- interaction_results(fit, endpoints, high, low, mode)
    data.table(
      country = scenario$country, scenario = identity$name, scenario_label = identity$label, nuisance = nuisance,
      mode = mode, unit = unit, status = "ok", n_days = nrow(d),
      events_K80_low = sum(d[[endpoints[[paste0("K80_", low)]]]]),
      events_N20_low = sum(d[[endpoints[[paste0("N20_", low)]]]]),
      events_K80_high = sum(d[[endpoints[[paste0("K80_", high)]]]]),
      events_N20_high = sum(d[[endpoints[[paste0("N20_", high)]]]]),
      p95 = target$p95, median = target$median, exposure_days = target$exposure_days,
      logRRR_difference = effect$logRRR_difference, variance_hac14 = effect$variance_hac14,
      se_hac14 = effect$se_hac14, RR_ratio_hac14 = effect$RR_ratio_hac14,
      RR_ratio_lo = effect$RR_ratio_lo, RR_ratio_hi = effect$RR_ratio_hi,
      joint_wald_chisq3 = effect$joint_wald_chisq3, joint_wald_p = effect$joint_wald_p,
      joint_basis_coefficients_json = as.character(toJSON(effect$joint_basis_coefficients, auto_unbox = TRUE)),
      joint_basis_covariance_json = as.character(toJSON(effect$joint_basis_covariance, auto_unbox = TRUE)),
      group_effects_json = as.character(toJSON(effect$groups, auto_unbox = TRUE)),
      formula = fit$audit$formula, warnings = as.character(toJSON(fit$audit$warnings, auto_unbox = TRUE))
    )
  }, error = function(e) {
    data.table(country = scenario$country, scenario = identity$name, scenario_label = identity$label, nuisance = nuisance,
               mode = mode, unit = unit, status = "failed", n_days = nrow(d),
               failure = conditionMessage(e))
  })
}

run_interaction_stage <- function(panels, units, force = FALSE) {
  all_rows <- list()
  failure_rows <- list()
  for (scenario in panels) for (mode in scenario$modes) for (nuisance in c("separate", "shared")) {
    scenario_identity(scenario$audit)
    path <- file.path(args$outdir, "cache", sprintf("interaction_%s_%s_%s_%s.csv", scenario$audit$country, scenario$audit$scenario, nuisance, mode))
    if (force && file.exists(path)) unlink(path)
    done <- cache_done_units(path)
    for (u in units) {
      if (u %in% done) next
      d <- scenario$panel[unit == u]
      row <- fit_four_endpoint_unit(d, scenario$audit, nuisance, mode, u, scenario$offset_col)
      cache_write_row(row, path)
    }
    if (file.exists(path)) all_rows[[length(all_rows) + 1L]] <- cache_read_rows(path)
    failure_path <- cache_failure_path(path)
    if (file.exists(failure_path)) failure_rows[[length(failure_rows) + 1L]] <- cache_read_rows(failure_path)
    log_event("interaction_model_complete", country = scenario$audit$country, scenario = scenario$audit$scenario,
              nuisance = nuisance, mode = mode, units_cached = length(cache_done_units(path)))
  }
  per_unit <- rbindlist(all_rows, fill = TRUE, use.names = TRUE)
  if (is.null(per_unit) || !ncol(per_unit)) {
    per_unit <- data.table(country = character(), scenario = character(), scenario_label = character(), nuisance = character(),
                           mode = character(), unit = character(), status = character())
  }
  if (nrow(per_unit)) per_unit <- unique(per_unit, by = c("scenario", "nuisance", "mode", "unit"), fromLast = TRUE)
  countries <- if (nrow(per_unit)) paste(sort(unique(per_unit$country)), collapse = "_") else
    paste(sort(unique(vapply(panels, function(x) x$audit$country, character(1)))), collapse = "_")
  fwrite(per_unit, file.path(args$outdir, paste0("clinical_", countries, "_interaction_per_unit.csv")), na = "")
  pooled_rows <- list()
  groups <- unique(per_unit[status == "ok", .(country, scenario, scenario_label, nuisance, mode)])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i]
    z <- per_unit[country == g$country & scenario == g$scenario & nuisance == g$nuisance & mode == g$mode & status == "ok"]
    for (pool_method in c("DL", "REML_HK")) {
      p <- pool_meta(z$logRRR_difference, z$variance_hac14, method = pool_method)
      pooled_rows[[length(pooled_rows) + 1L]] <- data.table(
        country = g$country, scenario = g$scenario, scenario_label = g$scenario_label,
        nuisance = g$nuisance, mode = g$mode, pooling = pool_method,
        units_fit = nrow(z), k = p$k, logRRR_difference = p$estimate, se = p$se,
        RR_ratio = p$RR, RR_lo = p$RR95CI[1L], RR_hi = p$RR95CI[2L],
        I2 = p$I2, tau2 = p$tau2, concordance = p$concordance,
        median_joint_wald_chisq3 = median(z$joint_wald_chisq3, na.rm = TRUE),
        joint_wald_p_lt_0_05_units = sum(z$joint_wald_p < 0.05, na.rm = TRUE),
        failure = if (is.null(p$failure)) NA_character_ else p$failure
      )
    }
  }
  pooled <- rbindlist(pooled_rows, fill = TRUE)
  if (is.null(pooled) || !ncol(pooled)) pooled <- data.table(country = character(), scenario = character(), scenario_label = character(),
                                              nuisance = character(), mode = character(), pooling = character(),
                                              units_fit = integer(), k = integer(), logRRR_difference = numeric(),
                                              se = numeric(), RR_ratio = numeric(), RR_lo = numeric(), RR_hi = numeric())
  fwrite(pooled, file.path(args$outdir, paste0("clinical_", countries, "_interaction_coverage.csv")), na = "")
  failures_dt <- rbindlist(failure_rows, fill = TRUE, use.names = TRUE)
  if (is.null(failures_dt) || !ncol(failures_dt)) {
    failures_dt <- data.table(country = character(), scenario = character(), scenario_label = character(), nuisance = character(),
                              mode = character(), unit = character(), status = character(), failure = character())
  }
  fwrite(failures_dt, file.path(args$outdir, paste0("clinical_", countries, "_interaction_failures.csv")), na = "")
  list(per_unit = per_unit, coverage = pooled)
}

build_adult_daily <- function(clinical, scenario, calendar_units) {
  date_keep <- !is.na(clinical$date) & clinical$date_year %in% scenario$analysis_years
  x <- clinical[date_keep & source_year %in% scenario$source_years & unit %in% calendar_units]
  if (!is.null(scenario$source_role) && "source_role" %in% names(x)) x <- x[source_role == scenario$source_role]
  x <- x[age_group %in% c("18-44", "45-64", "65plus") & sex_group %in% c("male", "female")]
  x[, stratum := paste(age_group, sex_group, sep = "_")]
  x[, disease := fcase(dx3 == "K80", "K80", dx3 %in% c("N20", "N21", "N22", "N23"), "N20_N23", default = NA_character_)]
  x <- x[!is.na(disease)]
  x[, .(n = sum(n, na.rm = TRUE)), by = .(unit, date, stratum, disease)]
}

case_mix_inputs <- function(clinical, panels, calendar) {
  strata <- c("18-44_male", "45-64_male", "65plus_male", "18-44_female", "45-64_female", "65plus_female")
  scenario_primary <- panels[[1L]]$audit
  scenario_id <- scenario_identity(scenario_primary)
  scope_keep <- !is.na(clinical$date) & clinical$date_year %in% scenario_primary$analysis_years &
    clinical$source_year %in% scenario_primary$source_years & clinical$dx3 %in% c("K80", "N20", "N21", "N22", "N23")
  if (!is.null(scenario_primary$source_role) && "source_role" %in% names(clinical)) {
    scope_keep <- scope_keep & clinical$source_role == scenario_primary$source_role
  }
  scope <- clinical[scope_keep]
  case_mix_eligibility <- data.table(
    country = scenario_primary$country,
    scenario = scenario_id$name,
    scenario_label = scenario_id$label,
    target_records_in_date_and_source_role_window = sum(scope$n),
    adult_known_sex_records_used_for_weights = scope[age_group %in% c("18-44", "45-64", "65plus") & sex_group %in% c("male", "female"), sum(n)],
    under18_records_excluded = scope[age_group == "<18", sum(n)],
    unknown_age_records_excluded = scope[age_group == "unknown", sum(n)],
    unknown_sex_records_excluded = scope[!sex_group %in% c("male", "female"), sum(n)],
    age_known_sex_unknown_records = scope[age_group %in% c("18-44", "45-64", "65plus") & !sex_group %in% c("male", "female"), sum(n)],
    age_unknown_sex_known_records = scope[age_group == "unknown" & sex_group %in% c("male", "female"), sum(n)],
    country_specific_weights_frozen_before_stratum_models = TRUE,
    standardization_interpretation = "Adult case-mix-standardized association conditional on common pooled K80+N20-N23 case weights; not a population-standardized rate."
  )
  primary_adult <- build_adult_daily(clinical, scenario_primary, unique(calendar$unit))
  eligible_keep <- format(calendar$date, "%Y") %in% as.character(scenario_primary$analysis_years) & is.finite(calendar$E3_actual)
  if (!is.null(scenario_primary$offset_col)) eligible_keep <- eligible_keep & is.finite(calendar$population) & calendar$population > 0
  eligible_calendar <- calendar[eligible_keep]
  primary_adult <- primary_adult[eligible_calendar[, .(unit, date)], on = .(unit, date), nomatch = 0L]
  weight_counts <- primary_adult[, .(cases = sum(n)), by = stratum]
  weights <- stats::setNames(rep(0, length(strata)), strata)
  weights[weight_counts$stratum] <- weight_counts$cases
  if (!sum(weights) || any(weights <= 0)) stop("At least one adult age-sex stratum has no cases; case-mix weights are undefined", call. = FALSE)
  weights <- weights / sum(weights)
  weight_signature_md5 <- md5_json(list(scenario = scenario_id$name, strata = strata,
                                        weights = as.list(weights), definition = "adult common K80+N20-N23 case mix"))
  list(strata = strata, scenario_primary = scenario_primary, case_mix_eligibility = case_mix_eligibility,
       eligible_calendar = eligible_calendar, weight_counts = weight_counts, weights = weights,
       weight_signature_md5 = weight_signature_md5)
}

summarize_standardized <- function(beta, V, weights, strata, label) {
  idxK <- match(paste("K80", strata, sep = "__"), names(beta))
  idxU <- match(paste("N20_N23", strata, sep = "__"), names(beta))
  if (anyNA(c(idxK, idxU))) stop("standardization coefficient order is incomplete", call. = FALSE)
  w <- as.numeric(weights[strata]); names(w) <- strata
  rK <- exp(beta[idxK]); rU <- exp(beta[idxU])
  sumK <- sum(w * rK); sumU <- sum(w * rU)
  grad <- rep(0, length(beta)); grad[idxU] <- w * rU / sumU; grad[idxK] <- -w * rK / sumK
  est <- log(sumU / sumK); variance <- as.numeric(t(grad) %*% V %*% grad)
  if (!is.finite(variance) || variance < 0) stop("case-mix delta-method variance is invalid", call. = FALSE)
  se <- sqrt(variance)
  data.table(weight_set = label, standardized_RR_K80 = sumK, standardized_RR_N20_N23 = sumU,
             logRRR = est, variance = variance, se = se, RRR = exp(est),
             RR_lo = exp(est - qnorm(0.975) * se), RR_hi = exp(est + qnorm(0.975) * se))
}

cache_standardized_result <- function(result, cache_path, cache_signature) {
  result$cache_signature_md5 <- cache_signature
  is_success <- identical(result$status, "ok")
  output_path <- if (is_success) cache_path else sub("\\.rds$", ".failure.rds", cache_path)
  temporary <- tempfile(pattern = paste0(".", basename(output_path), ".tmp-"), tmpdir = dirname(output_path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(result, temporary)
  atomic_replace_file(temporary, output_path)
  counterpart <- if (is_success) sub("\\.rds$", ".failure.rds", cache_path) else cache_path
  if (file.exists(counterpart)) unlink(counterpart)
  result
}

run_standardized_unit <- function(unit, scenario, adult_daily, calendar_unit, weights, strata, cache_dir, offset_col, cache_signature) {
  identity <- scenario_identity(scenario)
  unit_id <- as.character(unit)
  cache_path <- file.path(cache_dir, identity$name, paste0(safe_name(unit_id), ".rds"))
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  if (!args$force && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (!identical(cached$cache_signature_md5, cache_signature)) {
      stop(sprintf("Standardization cache fingerprint mismatch for %s; use a fresh output directory", cache_path), call. = FALSE)
    }
    return(cached)
  }
  fail <- function(reason) list(country = scenario$country, unit = unit_id, scenario = identity$name,
                                scenario_label = identity$label, status = "failed", failure = reason)
  if (nrow(calendar_unit) < 400L) {
    result <- fail(sprintf("only %d eligible calendar days; 400 required", nrow(calendar_unit)))
    return(cache_standardized_result(result, cache_path, cache_signature))
  }
  counts <- adult_daily[unit == unit_id, .(cases = sum(n)), by = .(stratum, disease)]
  full_counts <- merge(CJ(stratum = strata, disease = c("K80", "N20_N23"), unique = TRUE),
                       counts, by = c("stratum", "disease"), all.x = TRUE, sort = FALSE)
  full_counts[is.na(cases), cases := 0]
  low <- full_counts[cases < 30]
  if (nrow(low)) {
    detail <- paste(sprintf("%s/%s=%s", low$stratum, low$disease, low$cases), collapse = ";")
    result <- fail(paste("fewer than 30 cases in required stratum(s):", detail))
    return(cache_standardized_result(result, cache_path, cache_signature))
  }
  tryCatch({
    target <- unit_weather_target(calendar_unit)
    beta <- numeric(12L); names(beta) <- unlist(lapply(c("K80", "N20_N23"), function(d) paste(d, strata, sep = "__")))
    influence <- matrix(0, nrow = nrow(calendar_unit), ncol = length(beta),
                        dimnames = list(NULL, names(beta)))
    stratum_rows <- list()
    k <- 0L
    for (disease in c("K80", "N20_N23")) for (stratum in strata) {
      k <- k + 1L
      stratum_id <- stratum
      disease_id <- disease
      dd <- copy(calendar_unit[, .(date, E3 = E3_actual, population)])
      sub <- adult_daily[unit == unit_id & stratum == stratum_id & disease == disease_id, .(date, n)]
      if (nrow(sub)) setnames(sub, "n", "y") else sub <- data.table(date = as.Date(character()), y = numeric())
      dd <- merge(dd, sub, by = "date", all.x = TRUE, sort = FALSE)
      setorder(dd, date)
      dd[is.na(y), y := 0]
      fit <- fit_main_effect(
        data = dd, outcome = "y", exposure = "E3", date = "date", offset = offset_col,
        basis = target$basis, p95 = target$p95, median = target$median,
        hac_lags = 14L
      )
      X95 <- .design_at(fit$fit, fit$data, basis = fit$basis, temp_value = target$p95)
      Xcommon_median <- .design_at(fit$fit, fit$data, basis = fit$basis, temp_value = target$median)
      L <- as.numeric(X95 - Xcommon_median)
      names(L) <- names(stats::coef(fit$fit))
      effect <- linear_effect(L, fit$fit, fit$hac_vcov$hac14)
      beta[k] <- unname(effect["estimate"])
      sc <- fit$score
      infl <- as.numeric(sc$daily_scores %*% (sc$Hinv %*% L))
      pos <- match(as.character(sc$daily_dates), as.character(calendar_unit$date))
      if (anyNA(pos)) stop("daily influence dates could not be aligned to the common unit calendar", call. = FALSE)
      influence[pos, k] <- infl
      hac <- effect
      stratum_rows[[length(stratum_rows) + 1L]] <- data.table(
        unit = unit_id, scenario = identity$name, disease = disease_id, stratum = stratum_id,
        cases = full_counts[which(full_counts[["stratum"]] == stratum_id & full_counts[["disease"]] == disease_id), cases][1L],
        logRR = unname(hac["estimate"]), variance_hac14 = unname(hac["variance"]),
        RR = unname(hac["RR"]), RR_lo = unname(hac["lo"]), RR_hi = unname(hac["hi"]),
        p95 = target$p95, median = target$median, n_days = nrow(dd)
      )
    }
    V12 <- .bartlett_hac_meat(influence, calendar_unit$date, lag = 14L)
    dimnames(V12) <- list(names(beta), names(beta))
    primary <- summarize_standardized(beta, V12, weights, strata, "pooled_observed_case_mix")
    uniform_weights <- stats::setNames(rep(1 / length(strata), length(strata)), strata)
    uniform <- summarize_standardized(beta, V12, uniform_weights, strata, "uniform_1_over_6")
    age_rows <- list()
    for (age in c("18-44", "45-64", "65plus")) {
      age_strata <- strata[grepl(paste0("^", gsub("\\+", "\\\\+", age)), strata)]
      age_weights <- weights[age_strata] / sum(weights[age_strata])
      tmp <- summarize_standardized(beta, V12, age_weights, age_strata, paste0("age_band_", gsub("-", "_", age)))
      tmp[, age_group := age]
      age_rows[[length(age_rows) + 1L]] <- tmp
    }
    # The per-fit projected HAC14 variance should equal each diagonal element
    # in the date-aligned 12-score covariance, up to numerical tolerance.
    diagonal_audit <- list()
    for (i in seq_along(stratum_rows)) {
      row <- stratum_rows[[i]]
      nm <- paste(row$disease, row$stratum, sep = "__")
      v_fit <- row$variance_hac14
      v_joint <- V12[nm, nm]
      diagonal_audit[[nm]] <- abs(v_fit - v_joint)
    }
    result <- list(
      country = scenario$country, unit = unit_id, scenario = identity$name, scenario_label = identity$label, status = "ok",
      n_days = nrow(calendar_unit), p95 = target$p95, median = target$median,
      weights = as.list(weights), standardized = primary, uniform = uniform,
      age_band = rbindlist(age_rows), stratum_estimates = rbindlist(stratum_rows),
      influence_dates = as.character(calendar_unit$date), influence = influence,
      V12 = V12, diagonal_variance_max_abs_diff = max(unlist(diagonal_audit), na.rm = TRUE),
      interpretation = "Adult case-mix-standardized association conditional on common pooled K80+N20-N23 case weights; not a population-standardized rate."
    )
    cache_standardized_result(result, cache_path, cache_signature)
  }, error = function(e) {
    result <- fail(conditionMessage(e))
    cache_standardized_result(result, cache_path, cache_signature)
  })
}

run_standardized_stage <- function(clinical, panels, calendar, units, force = FALSE, weight_info, cache_signature) {
  strata <- weight_info$strata
  scenario_primary <- weight_info$scenario_primary
  scenario_identity(scenario_primary)
  case_mix_eligibility <- weight_info$case_mix_eligibility
  fwrite(case_mix_eligibility, file.path(args$outdir, paste0("clinical_", scenario_primary$country, "_case_mix_eligibility.csv")), na = "")
  eligible_calendar <- weight_info$eligible_calendar
  weight_counts <- weight_info$weight_counts
  weights <- weight_info$weights
  country_name <- scenario_primary$country
  weight_path <- file.path(args$outdir, paste0("clinical_", country_name, "_standardization_weights.csv"))
  fwrite(data.table(country = country_name, stratum = names(weights), pooled_K80_N20_N23_cases = as.numeric(weight_counts$cases[match(names(weights), weight_counts$stratum)]), weight = as.numeric(weights)),
         weight_path, na = "")
  rds_cache <- file.path(args$outdir, "cache", "standardized", country_name)
  outputs <- list(); failures <- list(); fit_rows <- list(); stratum_rows <- list(); age_rows <- list(); availability_rows <- list()
  for (scenario in panels) {
    scenario_meta <- scenario_identity(scenario$audit)
    adult_daily <- build_adult_daily(clinical, scenario$audit, unique(calendar$unit))
    adult_daily <- adult_daily[eligible_calendar[, .(unit, date)], on = .(unit, date), nomatch = 0L]
    cached_results <- list()
    for (u in units) {
      cu <- eligible_calendar[unit == u]
      result <- run_standardized_unit(u, scenario$audit, adult_daily, cu, weights, strata, rds_cache,
                                      scenario$offset_col, cache_signature)
      cached_results[[length(cached_results) + 1L]] <- result
      if (!identical(result$status, "ok")) {
        failures[[length(failures) + 1L]] <- data.table(country = country_name, unit = u, scenario = scenario_meta$name, status = "failed", failure = result$failure)
      } else {
        standardized_row <- copy(result$standardized)
        standardized_row[, `:=`(country = country_name, unit = u, scenario = scenario_meta$name, scenario_label = scenario_meta$label,
                                n_days = result$n_days, diagonal_variance_max_abs_diff = result$diagonal_variance_max_abs_diff,
                                interpretation = result$interpretation)]
        uniform_row <- copy(result$uniform)
        uniform_row[, `:=`(country = country_name, unit = u, scenario = scenario_meta$name, scenario_label = scenario_meta$label, n_days = result$n_days)]
        fit_rows[[length(fit_rows) + 1L]] <- rbindlist(list(standardized_row, uniform_row), fill = TRUE)
        strata_out <- copy(result$stratum_estimates)
        strata_out[, `:=`(country = country_name, scenario_label = scenario_meta$label)]
        stratum_rows[[length(stratum_rows) + 1L]] <- strata_out
        age_out <- copy(result$age_band); age_out[, `:=`(country = country_name, unit = u, scenario = scenario_meta$name, scenario_label = scenario_meta$label)]
        age_rows[[length(age_rows) + 1L]] <- age_out
      }
      if (length(cached_results) %% 10L == 0L || u == tail(units, 1L)) {
        log_event("standardized_unit_progress", scenario = scenario_meta$name,
                  completed_units = length(cached_results), total_units = length(units),
                  completed_ok = sum(vapply(cached_results, function(x) identical(x$status, "ok"), logical(1))),
                  failed = sum(vapply(cached_results, function(x) !identical(x$status, "ok"), logical(1))))
      }
    }
    outputs[[scenario_meta$name]] <- cached_results
    n_ok <- sum(vapply(cached_results, function(x) identical(x$status, "ok"), logical(1)))
    n_failed <- sum(vapply(cached_results, function(x) !identical(x$status, "ok"), logical(1)))
    failed_results <- cached_results[vapply(cached_results, function(x) !identical(x$status, "ok"), logical(1))]
    failure_reasons <- if (length(failed_results)) {
      as.list(table(vapply(failed_results, function(x) if (is.null(x$failure)) "unspecified" else x$failure, character(1))))
    } else list()
    availability_rows[[length(availability_rows) + 1L]] <- data.table(
      country = country_name, scenario = scenario_meta$name, scenario_label = scenario_meta$label,
      units_in_calendar_window = length(units), units_fit_successfully = n_ok, units_failed_or_unavailable = n_failed,
      failure_reasons = as.character(toJSON(failure_reasons, auto_unbox = TRUE))
    )
  }
  fit_dt <- rbindlist(fit_rows, fill = TRUE)
  strata_dt <- rbindlist(stratum_rows, fill = TRUE)
  age_dt <- rbindlist(age_rows, fill = TRUE)
  fail_dt <- rbindlist(failures, fill = TRUE)
  availability_dt <- rbindlist(availability_rows, fill = TRUE)
  if (is.null(fit_dt) || !ncol(fit_dt)) fit_dt <- data.table(country = character(), unit = character(), scenario = character(),
                                             weight_set = character(), logRRR = numeric(), variance = numeric())
  if (is.null(strata_dt) || !ncol(strata_dt)) strata_dt <- data.table(country = character(), unit = character(), scenario = character(),
                                                   disease = character(), stratum = character(), logRR = numeric())
  if (is.null(age_dt) || !ncol(age_dt)) age_dt <- data.table(country = character(), unit = character(), scenario = character(),
                                             age_group = character(), logRRR = numeric(), variance = numeric())
  if (is.null(fail_dt) || !ncol(fail_dt)) fail_dt <- data.table(country = character(), unit = character(), scenario = character(), failure = character())
  if (is.null(availability_dt) || !ncol(availability_dt)) availability_dt <- data.table(country = character(), scenario = character(), units_in_calendar_window = integer(), units_fit_successfully = integer(), units_failed_or_unavailable = integer())
  fwrite(fit_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_standardized_per_unit.csv")), na = "")
  fwrite(strata_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_strata_per_unit.csv")), na = "")
  fwrite(age_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_ageband_per_unit.csv")), na = "")
  fwrite(fail_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_failures.csv")), na = "")
  fwrite(availability_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_availability.csv")), na = "")
  pooled <- list()
  pool_groups <- if (nrow(fit_dt)) unique(fit_dt[, .(scenario, scenario_label, weight_set)]) else data.table()
  if (nrow(pool_groups)) for (i in seq_len(nrow(pool_groups))) {
    g <- pool_groups[i]
    z <- fit_dt[scenario == g$scenario & weight_set == g$weight_set]
    for (pool_method in c("DL", "REML_HK")) {
      p <- pool_meta(z$logRRR, z$variance, method = pool_method)
      pooled[[length(pooled) + 1L]] <- data.table(
        country = country_name, scenario = g$scenario, scenario_label = g$scenario_label,
        weight_set = g$weight_set, pooling = pool_method, units_fit = nrow(z), k = p$k,
        logRRR = p$estimate, se = p$se, RRR = p$RR, RR_lo = p$RR95CI[1L], RR_hi = p$RR95CI[2L],
        I2 = p$I2, tau2 = p$tau2, failure = if (is.null(p$failure)) NA_character_ else p$failure
      )
    }
  }
  pooled_dt <- rbindlist(pooled, fill = TRUE)
  if (is.null(pooled_dt) || !ncol(pooled_dt)) pooled_dt <- data.table(country = character(), scenario = character(), scenario_label = character(),
                                                                         weight_set = character(), pooling = character(),
                                                                         units_fit = integer(), k = integer(), logRRR = numeric(),
                                                                         se = numeric(), RRR = numeric(), RR_lo = numeric(), RR_hi = numeric())
  fwrite(pooled_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_pooled.csv")), na = "")
  age_pooled <- list()
  age_groups <- if (nrow(age_dt)) unique(age_dt[, .(scenario, scenario_label, age_group)]) else data.table()
  if (nrow(age_groups)) for (i in seq_len(nrow(age_groups))) {
    g <- age_groups[i]
    z <- age_dt[scenario == g$scenario & age_group == g$age_group]
    for (pool_method in c("DL", "REML_HK")) {
      p <- pool_meta(z$logRRR, z$variance, method = pool_method)
      age_pooled[[length(age_pooled) + 1L]] <- data.table(
        country = country_name, scenario = g$scenario, scenario_label = g$scenario_label,
        age_group = g$age_group, pooling = pool_method, units_fit = nrow(z), k = p$k,
        logRRR = p$estimate, se = p$se, RRR = p$RR, RR_lo = p$RR95CI[1L], RR_hi = p$RR95CI[2L],
        I2 = p$I2, tau2 = p$tau2, failure = if (is.null(p$failure)) NA_character_ else p$failure
      )
    }
  }
  age_pooled_dt <- rbindlist(age_pooled, fill = TRUE)
  if (is.null(age_pooled_dt) || !ncol(age_pooled_dt)) age_pooled_dt <- data.table(country = character(), scenario = character(),
                                                                                  age_group = character(), pooling = character(),
                                                                                  units_fit = integer(), k = integer(), logRRR = numeric(),
                                                                                  se = numeric(), RRR = numeric(), RR_lo = numeric(), RR_hi = numeric())
  fwrite(age_pooled_dt, file.path(args$outdir, paste0("clinical_", country_name, "_case_mix_ageband_pooled.csv")), na = "")
  list(per_unit = fit_dt, strata = strata_dt, ageband = age_dt, ageband_pooled = age_pooled_dt,
       eligibility = case_mix_eligibility, availability = availability_dt,
       failures = fail_dt, pooled = pooled_dt, weights = weights)
}

main <- function() {
  setDTthreads(1L)
  log_event("clinical_runner_start", country = args$country, mx_cohort = args$mx_cohort, stage = args$stage, root = args$root,
            batch_index = args$batch_index, batch_count = args$batch_count,
            brazil_aggregate = args$brazil_aggregate, mexico_aggregate = args$mexico_aggregate,
            brazil_calendar = args$calendar, mexico_calendar = args$mx_calendar,
            output = args$outdir, force = args$force)
  engine_hash <- unname(tools::md5sum(args$engine))
  analysis_plan <- jsonlite::fromJSON(args$analysis_plan, simplifyVector = FALSE)
  contexts <- list()
  if (args$country %in% c("brazil", "both")) {
    cal_loaded <- read_calendar(args$calendar)
    clinical <- read_clinical(args$brazil_aggregate)
    scenarios <- source_scenarios(include_latebilling = args$include_latebilling)
    panels <- lapply(scenarios, function(s) build_daily_panel(clinical, cal_loaded$data, s))
    contexts$brazil <- list(clinical = clinical, calendar = cal_loaded$data, calendar_audit = cal_loaded$audit,
                            panels = panels, units = sort(unique(cal_loaded$data$unit)),
                            aggregate_path = args$brazil_aggregate, calendar_path = args$calendar)
  }
  if (args$country %in% c("mexico", "both")) {
    cal_loaded <- read_calendar(args$mx_calendar)
    clinical <- read_mexico_clinical(args$mexico_aggregate, cohort = args$mx_cohort)
    scenario <- mexico_primary_scenario(cohort = args$mx_cohort)
    panels <- list(build_daily_panel(clinical, cal_loaded$data, scenario))
    contexts$mexico <- list(clinical = clinical, calendar = cal_loaded$data, calendar_audit = cal_loaded$audit,
                            panels = panels, units = sort(unique(cal_loaded$data$unit)),
                            aggregate_path = args$mexico_aggregate, calendar_path = args$mx_calendar)
  }

  weight_infos <- list()
  if (args$stage %in% c("standardized", "all")) {
    for (country_name in names(contexts)) {
      ctx <- contexts[[country_name]]
      weight_infos[[country_name]] <- case_mix_inputs(ctx$clinical, ctx$panels, ctx$calendar)
    }
  }
  selected_units <- list()
  for (country_name in names(contexts)) {
    ctx <- contexts[[country_name]]
    selected_units[[country_name]] <- select_unit_batch(ctx$units, args$batch_index, args$batch_count, args$max_units)
    contexts[[country_name]]$selected_units <- selected_units[[country_name]]
  }
  weight_signatures <- lapply(weight_infos, function(x) list(md5 = x$weight_signature_md5, weights = as.list(x$weights)))
  cache_manifest <- initialize_cache_manifest(analysis_plan, weight_signatures)
  log_event("clinical_cache_fingerprint", signature_md5 = cache_manifest$signature_md5,
            config_md5 = cache_manifest$config_md5, source_hashes = cache_manifest$source_hashes,
            weight_signatures = weight_signatures, batch_index = args$batch_index, batch_count = args$batch_count)

  results <- list()
  for (country_name in names(contexts)) {
    ctx <- contexts[[country_name]]
    units <- selected_units[[country_name]]
    for (obj in ctx$panels) {
      derived_name <- paste0(country_name, "_clinical_daily_", obj$audit$scenario, ".csv.gz")
      fwrite(obj$panel, file.path(args$outdir, "derived", derived_name), na = "")
      car_summary <- obj$selected[, .(admission_records = sum(n)), by = .(source_year, dx3, car_int_raw, car_int_group)]
      fwrite(car_summary, file.path(args$outdir, paste0(country_name, "_car_int_distribution_", obj$audit$scenario, ".csv")), na = "")
      legacy <- obj$panel[, .(
        matched_days = .N,
        legacy_K80_records = sum(legacy_y_K80, na.rm = TRUE),
        reconstructed_K80_records = sum(y_K80_all, na.rm = TRUE),
        legacy_urinary_records = sum(legacy_y_urinary, na.rm = TRUE),
        reconstructed_N20_N23_records = sum(y_N20_N23, na.rm = TRUE),
        source_2021_admission_2020_rows_included = obj$audit$source_2021_admission_2020_rows_included,
        source_2023_rows_included_for_original_admission_dates = obj$audit$source_2023_rows_included_for_original_admission_dates
      )]
      legacy[, `:=`(country = country_name, scenario = obj$audit$scenario, scenario_label = obj$audit$scenario_label)]
      fwrite(legacy, file.path(args$outdir, paste0(country_name, "_combined_original_outcomes_", obj$audit$scenario, ".csv")), na = "")
      log_event("clinical_panel_ready", country = country_name, scenario = obj$audit$scenario, audit = obj$audit)
    }
    endpoints_result <- NULL; interactions_result <- NULL; standardized_result <- NULL
    if (args$stage %in% c("endpoint", "all")) endpoints_result <- run_endpoint_stage(ctx$panels, units, force = args$force)
    if (args$stage %in% c("interaction", "all")) interactions_result <- run_interaction_stage(ctx$panels, units, force = args$force)
    if (args$stage %in% c("standardized", "all")) {
      standardized_result <- run_standardized_stage(ctx$clinical, ctx$panels, ctx$calendar, units,
                                                    force = args$force, weight_info = weight_infos[[country_name]],
                                                    cache_signature = cache_manifest$signature_md5)
    }
    contexts[[country_name]]$results <- list(endpoint = endpoints_result, interaction = interactions_result,
                                             standardized = standardized_result)
  }

  analysis_audits <- lapply(contexts, function(ctx) {
    list(
      source_aggregate = ctx$aggregate_path,
      source_calendar = ctx$calendar_path,
      source_coverage = lapply(ctx$panels, `[[`, "audit"),
      calendar_exposure_audit = ctx$calendar_audit,
      units_in_calendar = uniqueN(ctx$calendar$unit),
      units_considered = length(ctx$selected_units),
      selected_unit_ids = ctx$selected_units
    )
  })
  audit <- list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    country = args$country,
    mexico_cohort = args$mx_cohort,
    stage = args$stage,
    unit_batch = list(index = args$batch_index, count = args$batch_count,
                      method = "deterministic round-robin over sorted unit IDs after full-cohort weight calculation"),
    repair_engine = args$engine,
    repair_engine_md5 = engine_hash,
    source_file_hashes = cache_manifest$source_hashes,
    cache_signature_md5 = cache_manifest$signature_md5,
    configuration_md5 = cache_manifest$config_md5,
    reference_md5 = cache_manifest$reference_md5,
    standardization_weight_signatures = weight_signatures,
    analysis_plan = args$analysis_plan,
    analysis_plan_primary_target = analysis_plan$primary_target,
    analysis_plan_primary_design = analysis_plan$primary_design,
    analysis_plan_reference = analysis_plan$reference,
    analysis_plan_selection_policy = analysis_plan$selection_policy,
    countries = analysis_audits,
    exposure = "Complete actual-calendar tmean mean for admission day and two preceding days, rebuilt from each full calendar before outcome-period filtering.",
    diagnosis_definitions = list(
      normalization = "Uppercase diagnosis, remove punctuation, retain up to four alphanumeric characters; no fourth character is imputed when absent.",
      strict_primary = "N20-only versus all K80.",
      confirmed_calculi = "N20+N21 versus all K80.",
      composition = "N23 versus all K80.",
      inflammatory_calculous_biliary = "N20-only versus K80.0 cholecystitis plus K80.3 cholangitis; K80.3 is not labeled as acute.",
      combined_original = "N20-N23 versus K80 on the same dates and source-year/role cohort as strict endpoints."
    ),
    country_specific_cohorts = list(
      brazil = "Matched-source replication uses source years 2014-2020 and DT_INTER years 2014-2020. Late-billing reconstruction uses source years 2014-2021 for the same admission-date window; no completeness is assumed.",
      mexico = if (args$mx_cohort == "baseline") {
        paste0("Baseline uses source_role=primary_panel_year only, source years 2015-2018 and 2020-2023, and admission years 2015-2018 and 2020-2022. ",
               "Recovered 2019 rows and 2023 admission-year extension are excluded; primary source-year 2023 billings for 2021-2022 remain included.")
      } else {
        "Complete Mexico cohort uses all source roles and source years 2015-2022 for admission years 2015-2022."
      }
    ),
    admission_type = list(
      brazil_elective = "CAR_INT 01 only.", brazil_urgent = "CAR_INT 02 only.",
      brazil_other_unknown = "Other and unknown CAR_INT values remain separate in descriptive CSVs and are excluded from elective/urgent contrasts.",
      mexico = "CAR_INT/TIPSERV admission interaction is unavailable and not inferred."
    ),
    sex = "Brazil uses SEXO code 1=male and 3=female; Mexico uses the aggregate worker's male/female harmonization. Unknown sex values are excluded from the male/female contrast.",
    model = list(
      family = "quasi-Poisson with log link.",
      offset = "Brazil uses log regional population; Mexico has no population offset.",
      temperature_basis = "Unit-specific natural spline with three degrees of freedom; p95 and median shared across outcomes on each country-specific eligible calendar sample.",
      separate_calendar = "Disease-specific exposure, day-of-week, month, and 12-df time spline.",
      shared_calendar = "Shared day-of-week, month, and 12-df time spline sensitivity; disease-specific exposure response.",
      covariance = "Paired actual-calendar Bartlett HAC lag 14 primary; lag 28, paired-day HC0, and quasi-model covariance sensitivities where available.",
      pooling = "Per-unit log-RRR pooled with DerSimonian-Laird and REML/Hartung-Knapp; between-unit independence limitation remains."
    ),
    interactions = list(
      brazil_admission_type = "Four endpoint fit: K80/N20-only by CAR_INT 01/02; report urgent-minus-elective log-RRR and joint 3-coefficient temperature-basis Wald test.",
      sex = "Four endpoint fit: K80/N20-only by male/female; report male-minus-female log-RRR and joint 3-coefficient temperature-basis Wald test."
    ),
    case_mix_standardization = list(
      age_strata = c("18-44", "45-64", "65plus"), sex_strata = c("male", "female"),
      minimum_cases = 30, minimum_days = 400,
      weights = "Country-specific common pooled observed K80+N20-N23 adult case-mix weights, frozen before stratum effect fitting; not population-standardized rates.",
      covariance = "Twelve per-stratum daily influence scalars aligned by actual date and jointly combined with Bartlett HAC lag 14.",
      sensitivities = c("uniform one-sixth weights", "age-band-specific contrasts"),
      failure_reporting = "Units below any case/day threshold or with model failures are written to country-specific clinical_case_mix_failures.csv; unavailable unit coverage is not imputed."
    ),
    software = list(R = R.version.string, data_table = as.character(packageVersion("data.table")))
  )
  audit_path <- file.path(args$outdir, paste0("clinical_", args$country, "_analysis_audit.json"))
  jsonlite::write_json(audit, audit_path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  endpoint_rows <- sum(vapply(contexts, function(ctx) if (is.null(ctx$results$endpoint)) 0L else nrow(ctx$results$endpoint$per_unit), integer(1)))
  interaction_rows <- sum(vapply(contexts, function(ctx) if (is.null(ctx$results$interaction)) 0L else nrow(ctx$results$interaction$per_unit), integer(1)))
  standardized_rows <- sum(vapply(contexts, function(ctx) if (is.null(ctx$results$standardized)) 0L else nrow(ctx$results$standardized$per_unit), integer(1)))
  log_event("clinical_runner_complete", country = args$country, stage = args$stage,
            endpoint_rows = endpoint_rows, interaction_rows = interaction_rows,
            standardized_rows = standardized_rows, audit = audit_path)
}

main()
