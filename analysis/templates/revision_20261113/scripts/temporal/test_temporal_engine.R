#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
engine_path <- if (length(args)) args[[1L]] else
  file.path("work", "revision_20261113", "temporal", "scripts", "repair_engine_temporal.R")
source(engine_path)

set.seed(241113L)
n <- 3L * 365L
dates <- seq(as.Date("2017-01-01"), by = "day", length.out = n)
t <- seq_len(n)
exposure <- 20 + 6 * sin(2 * pi * t / 365.25) + 0.7 * cos(2 * pi * t / 29)
holiday <- as.integer(as.integer(format(dates, "%j")) %% 29L == 0L)
basis <- new_temperature_basis(exposure, df = 3L)
p95 <- unname(quantile(exposure, 0.95, type = 7))
median <- stats::median(exposure)
data <- data.frame(
  date = dates, E3 = exposure, national_holiday = holiday,
  y_K80 = stats::rpois(n, lambda = exp(1.4 + 0.03 * exposure / 10)),
  y_urinary = stats::rpois(n, lambda = exp(1.3 + 0.05 * exposure / 10))
)

flexible <- fit_stacked_contrast(
  data = data, outcomes = c("y_K80", "y_urinary"), exposure = "E3", date = "date",
  basis = basis, p95 = p95, median = median, nuisance = "flexible",
  hac_lags = c(14L, 28L), full_covariance = TRUE
)
stopifnot(flexible$flexible_time_df == 18L)
stopifnot(flexible$national_holiday_days == sum(holiday == 1L))
stopifnot(identical(flexible$p95, p95), identical(flexible$median, median))
stopifnot(grepl("national_holiday", flexible$audit$formula, fixed = TRUE))
stopifnot(grepl("ns(time, 18)", flexible$audit$formula, fixed = TRUE))
stopifnot("y_urinary_vs_y_K80" %in% names(flexible$contrasts))
stopifnot(all(c("hac14", "hac28", "quasi", "day_hc0") %in%
  names(flexible$contrasts[["y_urinary_vs_y_K80"]]$coefficient_contrast$by_covariance)))

yearmonth <- fit_stacked_contrast(
  data = data, outcomes = c("y_K80", "y_urinary"), exposure = "E3", date = "date",
  basis = basis, p95 = p95, median = median, nuisance = "yearmonth",
  hac_lags = c(14L, 28L), full_covariance = TRUE
)
stopifnot(grepl("ym", yearmonth$audit$formula, fixed = TRUE))
stopifnot(identical(yearmonth$p95, p95), identical(yearmonth$median, median))

missing_holiday <- tryCatch({
  fit_stacked_contrast(
    data = data[setdiff(names(data), "national_holiday")],
    outcomes = c("y_K80", "y_urinary"), exposure = "E3", date = "date",
    basis = basis, p95 = p95, median = median, nuisance = "flexible"
  )
  NULL
}, error = identity)
stopifnot(inherits(missing_holiday, "error"), grepl("holiday indicator", conditionMessage(missing_holiday)))
cat("temporal engine smoke passed; flexible df=18, target overrides retained, HAC14/HAC28 emitted, year-month retained\n")
