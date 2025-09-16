#' Read CPI data from an Excel file
#'
#' Loads and normalizes a CPI time series from an Excel worksheet. The function
#' detects the date/year column and the CPI/value column by pattern-matching on
#' lower-cased header names, parses localized numerics (via \code{to_num_commas()}),
#' collapses duplicate years by averaging, and returns a clean, sorted data frame.
#'
#' @param path_cpi Character path to the CPI Excel file.
#'
#' @details
#' \strong{Column detection.} Headers are lower-cased and matched with:
#' \itemize{
#'   \item Date/year: patterns \code{"date|fecha|year|anio|ano"}.
#'   \item CPI/value: patterns \code{"cpi|indice|price"}.
#' }
#' If either column cannot be identified, the function errors.
#'
#' \strong{Cleaning.}
#' \itemize{
#'   \item Year is extracted as the first 4 digits of the date-like column.
#'   \item CPI is parsed with \code{to_num_commas()} (handles commas/thousands).
#'   \item \code{NA} rows are dropped; duplicates in \code{Year} are averaged.
#'   \item Output is sorted by \code{Year} ascending.
#' }
#'
#' @return A \code{data.frame} with two columns:
#' \itemize{
#'   \item \code{Year} (\code{integer})
#'   \item \code{CPI}  (\code{numeric})
#' }
#'
#' @examples
#' \donttest{
#' # Suppose "CPI.xlsx" has columns like "Fecha" and "Indice":
#' # df <- read_cpi("CPI.xlsx")
#' # head(df)
#' }
#'
#' @seealso \code{\link{read_weights_matrix}}, \code{\link{bayesian_disaggregate}}
#' @export
read_cpi <- function(path_cpi) {
  log_msg("INFO", "Reading CPI:", path_cpi)
  df <- readxl::read_excel(path_cpi)
  cn <- tolower(names(df))
  col_date <- which.max(stringr::str_detect(cn, "date|fecha|year|anio|ano"))
  col_cpi  <- which.max(stringr::str_detect(cn, "cpi|indice|price"))

  if (length(col_date) != 1 || length(col_cpi) != 1) {
    stop("Could not identify date and CPI columns in the CPI file")
  }

  df <- df %>%
    dplyr::transmute(
      Year = as.integer(substr(as.character(.[[col_date]]), 1, 4)),
      CPI  = to_num_commas(.[[col_cpi]])
    ) %>%
    dplyr::filter(!is.na(Year), !is.na(CPI)) %>%
    dplyr::arrange(Year)

  if (anyDuplicated(df$Year)) {
    log_msg("WARN", "CPI with duplicate years; aggregating by average")
    df <- df %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(CPI = mean(CPI), .groups = "drop")
  }

  df
}

#' Read a weights matrix from an Excel file
#'
#' Loads a sector-by-year weight table, normalizes weights to the simplex per year,
#' and returns a list with the \eqn{T \times K} prior matrix \code{P}, the sector
#' names, and the year vector. The first column is assumed to contain sector names
#' (renamed to \code{Industry}); all other columns are treated as years.
#'
#' @param path_weights Character path to the weights Excel file.
#'
#' @details
#' \strong{Expected layout.} One sheet with:
#' \itemize{
#'   \item First column: sector names (any header; renamed to \code{Industry}).
#'   \item Remaining columns: years; the function extracts a 4-digit year from each
#'         header using \code{stringr::str_extract(Year, "\\\\d{4}")}.
#' }
#' Values are parsed with \code{to_num_commas()}, missing rows are dropped, and
#' weights are normalized within each year to sum to 1. Any absent (sector, year)
#' entry becomes 0 when pivoting wide. Finally, rows are re-normalized with
#' \code{row_norm1()} for numerical safety.
#'
#' \strong{Safeguards.}
#' \itemize{
#'   \item Rows with all-missing/zero after parsing are dropped by the filters.
#'   \item If no valid year columns are found, the function errors.
#' }
#'
#' @return A \code{list} with:
#' \describe{
#'   \item{\code{P}}{\eqn{T \times K} numeric matrix of prior weights (rows sum to 1).}
#'   \item{\code{industries}}{Character vector of sector names (length \eqn{K}).}
#'   \item{\code{years}}{Integer vector of years (length \eqn{T}).}
#' }
#'
#' @examples
#' \donttest{
#' # Suppose "WEIGHTS.xlsx" has first column "Industry" and year columns "2019","2020",...
#' # w <- read_weights_matrix("WEIGHTS.xlsx")
#' # str(w)
#' }
#'
#' @seealso \code{\link{read_cpi}}, \code{\link{bayesian_disaggregate}}
#' @export
read_weights_matrix <- function(path_weights) {
  log_msg("INFO", "Reading weights matrix:", path_weights)
  raw <- readxl::read_excel(path_weights)
  names(raw)[1] <- "Industry"

  long <- raw %>%
    tidyr::pivot_longer(-Industry, names_to = "Year", values_to = "Weight_raw") %>%
    dplyr::mutate(
      Year   = as.integer(stringr::str_extract(Year, "\\d{4}")),
      Weight = to_num_commas(Weight_raw)
    ) %>%
    dplyr::filter(!is.na(Year), !is.na(Weight))

  long <- long %>%
    dplyr::group_by(Year) %>%
    dplyr::mutate(Weight = Weight / sum(Weight, na.rm = TRUE)) %>%
    dplyr::ungroup()

  wide <- long %>%
    dplyr::select(Industry, Year, Weight) %>%
    tidyr::pivot_wider(names_from = Industry, values_from = Weight, values_fill = 0) %>%
    dplyr::arrange(Year)

  P <- as.matrix(wide[, -1, drop = FALSE])
  rownames(P) <- wide$Year
  P <- row_norm1(P)

  list(P = P, industries = colnames(P), years = as.integer(rownames(P)))
}