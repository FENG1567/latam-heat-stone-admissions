options(warn = 1)
warning_log <- Sys.getenv("CURVE_WARNING_LOG")
if (!nzchar(warning_log)) stop("CURVE_WARNING_LOG is required")
warning_n <- 0L
withCallingHandlers({
  sys.source("/data/admin/fengshuo/latam_stone_repair_20261112/revision_20261113/scripts/marginals_curves/audit_curve_pooling.R", envir = .GlobalEnv)
}, warning = function(w) {
  warning_n <<- warning_n + 1L
  call_text <- if (is.null(conditionCall(w))) "<no call>" else paste(deparse(conditionCall(w)), collapse = " ")
  line <- sprintf("[%03d] %s | %s | %s", warning_n, paste(class(w), collapse = ","), call_text, conditionMessage(w))
  write(line, file = warning_log, append = TRUE)
  invokeRestart("muffleWarning")
})
cat(sprintf("Captured %d warnings into %s\n", warning_n, warning_log))
