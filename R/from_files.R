# =============================================================================
# Convenience reader: build the state-space disaggregation directly from the
# CPI and VAB-weights Excel files, reusing the existing robust readers and
# aligning on the common years.
# =============================================================================

#' Align CPI and VAB weights on their common years
#'
#' Reads the CPI and the weights matrix, intersects their years, and returns the
#' aligned aggregate vector and weight matrix ready for the disaggregation
#' engines. Both engines (state-space and conjugate) consume this output.
#'
#' @param path_cpi Path to the CPI Excel file (see \code{\link{read_cpi}}).
#' @param path_weights Path to the VAB-weights Excel file
#'   (see \code{\link{read_weights_matrix}}).
#' @return A list with \code{cpi} (numeric, length \eqn{T}), \code{W}
#'   (\eqn{T \times K}, rows sum to 1), \code{years}, \code{industries}.
#' @seealso \code{\link{disaggregate_from_files}}
#' @export
align_disagg_inputs <- function(path_cpi, path_weights) {
  cpi_df <- read_cpi(path_cpi)
  Wraw   <- read_weights_matrix(path_weights)
  Pfull  <- Wraw$P

  years_common <- sort(unique(intersect(as.integer(rownames(Pfull)), cpi_df$Year)))
  if (length(years_common) < 2L) {
    stop("Insufficient common years between CPI and weights data.", call. = FALSE)
  }
  W   <- Pfull[as.character(years_common), , drop = FALSE]
  cpi <- cpi_df$CPI[match(years_common, cpi_df$Year)]

  list(cpi = as.numeric(cpi), W = W, years = years_common,
       industries = colnames(W))
}

#' Evidence-based disaggregation directly from Excel files
#'
#' Thin convenience wrapper: reads and aligns the CPI and VAB-weight files
#' (\code{\link{align_disagg_inputs}}) and runs the canonical state-space engine
#' (\code{\link{disaggregate_statespace}}).
#'
#' @param path_cpi Path to the CPI Excel file (index levels, re-indexed to the
#'   same base as the production prices; see the package vignette and the data
#'   note on \code{convert_to_index}).
#' @param path_weights Path to the VAB-weights Excel file.
#' @param ... Passed to \code{\link{disaggregate_statespace}} (sampler controls,
#'   \code{priors}, \code{student_obs}, ...).
#' @return A \code{"disagg_statespace"} object.
#' @examples
#' \dontrun{
#' cpi_file <- system.file("extdata", "CPI.xlsx", package = "BayesianDisaggregation")
#' w_file   <- system.file("extdata", "WEIGHTS.xlsx", package = "BayesianDisaggregation")
#' fit <- disaggregate_from_files(cpi_file, w_file, chains = 2, iter = 800)
#' }
#' @seealso \code{\link{disaggregate_statespace}}, \code{\link{align_disagg_inputs}}
#' @export
disaggregate_from_files <- function(path_cpi, path_weights, ...) {
  al <- align_disagg_inputs(path_cpi, path_weights)
  disaggregate_statespace(cpi = al$cpi, W = al$W, years = al$years,
                          industries = al$industries, ...)
}
