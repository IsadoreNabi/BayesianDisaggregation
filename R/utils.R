# ----- Internal package state -------------------------------------------------
# Private environment for package-level state (no writes to package namespace)
.bd_env <- new.env(parent = emptyenv())
.bd_env$log_level <- "INFO"

# Log levels (ordering)
.level_order <- c("TRACE" = 1, "DEBUG" = 2, "INFO" = 3, "WARN" = 4, "ERROR" = 5)

#' Enable logging at a specific level
#'
#' Sets the package-wide logging verbosity.
#'
#' @param level Character scalar. One of "TRACE", "DEBUG", "INFO", "WARN", "ERROR".
#' @export
log_enable <- function(level = "INFO") {
  if (!is.character(level) || length(level) != 1L || !level %in% names(.level_order)) {
    stop('Invalid log level. Choose one of: "TRACE","DEBUG","INFO","WARN","ERROR".')
  }
  .bd_env$log_level <- level
  invisible(level)
}

#' Log message with timestamp
#'
#' Internal helper that prints a timestamped message when the current log
#' level is at least \code{level}.
#'
#' @param level Character level: "TRACE","DEBUG","INFO","WARN","ERROR".
#' @param ...   Message components (will be concatenated with spaces).
#' @keywords internal
log_msg <- function(level = "INFO", ...) {
  cur <- if (!is.null(.bd_env$log_level)) .bd_env$log_level else "INFO"
  if (.level_order[[level]] >= .level_order[[cur]]) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s][%s] %s", ts, level, paste(..., collapse = " ")))
  }
}

# ----- Mathematical utilities -------------------------------------------------

safe_div <- function(p, q) {
  p_safe <- pmax(p, .Machine$double.eps)
  q_safe <- pmax(q, .Machine$double.eps)
  p_safe / q_safe
}

kl_divergence <- function(P_row, Q_row) {
  sum(P_row * log(safe_div(P_row, Q_row)))
}

total_variation <- function(P_row, Q_row) {
  0.5 * sum(abs(P_row - Q_row))
}

robust_cor <- function(x, y) {
  c_p <- suppressWarnings(cor(x, y, method = "pearson"))
  c_s <- suppressWarnings(cor(x, y, method = "spearman"))
  if (is.na(c_p) && is.na(c_s)) return(0)
  if (is.na(c_p)) return(c_s)
  if (is.na(c_s)) return(c_p)
  if (abs(c_s) > abs(c_p)) c_s else c_p
}

#' Row-wise L1 normalization (safe)
#' Ensures each row sums to 1; protects against NaN/Inf and near-zero sums.
#' @keywords internal
#' @noRd
row_norm1 <- function(X, eps = .Machine$double.eps) {
  rs <- rowSums(X)
  rs[!is.finite(rs) | rs < eps] <- 1
  X / rs
}

to_num_commas <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- stringr::str_replace_all(as.character(x), "\\.", "")
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

# If NOTES about globals appear in checks, keep this; otherwise you can remove it.
utils::globalVariables(c(
  ".", "Year", "CPI", "Industry", "Weight_raw", "Weight",
  "i", "composite", "value", "name", "cfg_id", "desc"
))