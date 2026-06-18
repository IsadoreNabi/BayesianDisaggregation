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
#' @return (Invisibly) the level set.
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

#' Row-wise L1 normalization (safe)
#' Ensures each row sums to 1; protects against NaN/Inf and near-zero sums.
#' @keywords internal
#' @noRd
row_norm1 <- function(X, eps = .Machine$double.eps) {
  rs <- rowSums(X)
  rs[!is.finite(rs) | rs < eps] <- 1
  X / rs
}

#' Parse localized numerics (thousands as ".", decimal as ",")
#' @keywords internal
#' @noRd
to_num_commas <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- stringr::str_replace_all(as.character(x), "\\.", "")
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

# NSE column names used by the tidyverse data readers (read_cpi / read_weights_matrix).
utils::globalVariables(c(".", "Year", "CPI", "Industry", "Weight_raw", "Weight"))
