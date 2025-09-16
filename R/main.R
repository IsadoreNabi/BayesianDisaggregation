#' Run Bayesian disaggregation
#'
#' Performs Bayesian disaggregation of an aggregated time series (e.g., CPI)
#' into \eqn{K} components using one of four deterministic update rules:
#' \emph{weighted}, \emph{multiplicative}, \emph{dirichlet}, \emph{adaptive}.
#'
#' @param path_cpi Path to the CPI Excel file. Must contain at least the
#'   columns \code{date} and \code{value} (case/locale tolerant in \code{read_cpi()}).
#' @param path_weights Path to the Excel file with the baseline weight matrix
#'   (prior): either \eqn{T \times K} (years in rows, sectors in columns) or a
#'   length-\eqn{K} vector (constant across time). Rows are renormalized to the simplex.
#' @param method Disaggregation method: \code{"weighted"}, \code{"multiplicative"},
#'   \code{"dirichlet"}, or \code{"adaptive"}.
#' @param lambda Weight for the \code{"weighted"} method in \eqn{[0,1]}. Ignored otherwise.
#' @param gamma Uncertainty factor for the \code{"dirichlet"} method (\eqn{> 0}).
#' @param coh_mult Multiplier for the coherence increment \eqn{\Delta\rho}.
#' @param coh_const Constant offset for coherence, truncated to \eqn{[0,1]}.
#' @param stab_a Sensitivity for row-sum deviation penalty \eqn{| \sum w - 1 |}.
#' @param stab_b Sensitivity for negative-values penalty (count of negatives).
#' @param stab_kappa Sensitivity for temporal variation (average \eqn{|\Delta|}).
#' @param likelihood_pattern Temporal spreading pattern for the likelihood:
#'   \code{"constant"}, \code{"recent"}, \code{"linear"}, or \code{"bell"}.
#'
#' @details
#' Assumptions: (i) prior/posterior rows lie on the simplex; (ii) no MCMC is used,
#' updates are analytic/deterministic; (iii) \code{read_*} helpers coerce benign
#' formatting issues and error on malformed inputs.
#'
#' @return A \code{list} with:
#' \describe{
#'   \item{\code{years}}{Integer vector of years used.}
#'   \item{\code{industries}}{Character vector of sector/column names.}
#'   \item{\code{prior}}{Tibble prior \eqn{T \times (1+K)} with \code{Year} then sectors.}
#'   \item{\code{likelihood_t}}{Tibble likelihood over time (same shape as \code{prior}).}
#'   \item{\code{likelihood}}{Tibble \code{Sector}, \code{L} (length \eqn{K}).}
#'   \item{\code{posterior}}{Tibble posterior \eqn{T \times (1+K)} (rows sum to 1).}
#'   \item{\code{metrics}}{Tibble with hyperparameters + \code{coherence}, \code{stability},
#'         \code{interpretability}, \code{efficiency}, \code{composite}, \code{T}, \code{K}.}
#' }
#'
#' @examples
#' \donttest{
#' # Minimal synthetic run (no files):
#' T <- 6; K <- 4
#' P <- matrix(rep(1/K, T*K), nrow = T)
#' L <- runif(K); L <- L/sum(L)
#' LT <- spread_likelihood(L, T, "recent")
#' W  <- posterior_weighted(P, LT, lambda = 0.7)
#' }
#'
#' @seealso \code{\link{read_cpi}}, \code{\link{read_weights_matrix}},
#'   \code{\link{compute_L_from_P}}, \code{\link{spread_likelihood}},
#'   \code{\link{posterior_weighted}}, \code{\link{posterior_multiplicative}},
#'   \code{\link{posterior_dirichlet}}, \code{\link{posterior_adaptive}},
#'   \code{\link{coherence_score}}, \code{\link{stability_composite}},
#'   \code{\link{interpretability_score}}
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr mutate relocate
#' @export
bayesian_disaggregate <- function(path_cpi, path_weights,
                                  method = c("weighted", "multiplicative", "dirichlet", "adaptive"),
                                  lambda = 0.7, gamma = 0.1,
                                  coh_mult = 3.0, coh_const = 0.5,
                                  stab_a = 1000, stab_b = 10, stab_kappa = 50,
                                  likelihood_pattern = "recent") {
  method <- match.arg(method)
  if (!is.numeric(lambda) || lambda < 0 || lambda > 1) stop("lambda must be in [0,1].")
  if (!is.numeric(gamma)  || gamma <= 0)                stop("gamma must be > 0.")
  likelihood_pattern <- match.arg(likelihood_pattern, c("constant","recent","linear","bell"))
  coh_const <- max(0, min(1, coh_const))

  # Read inputs
  cpi  <- read_cpi(path_cpi)
  Wraw <- read_weights_matrix(path_weights)
  Pfull <- Wraw$P

  # Intersect years
  years_common <- intersect(as.integer(rownames(Pfull)), cpi$Year)
  years_common <- sort(unique(years_common))
  if (length(years_common) < 2) {
    stop("Insufficient common years between CPI and weights data.")
  }

  P <- Pfull[as.character(years_common), , drop = FALSE]
  Tn <- nrow(P); K <- ncol(P)

  # Likelihood (vector) and temporal spreading
  L  <- compute_L_from_P(P)
  LT <- spread_likelihood(L, Tn, pattern = likelihood_pattern)

  # Posterior by method
  W <- switch(method,
    weighted       = posterior_weighted(P, LT, lambda),
    multiplicative = posterior_multiplicative(P, LT),
    dirichlet      = posterior_dirichlet(P, LT, gamma),
    adaptive       = posterior_adaptive(P, LT)
  )

  # Metrics
  coh   <- coherence_score(P, W, L, mult = coh_mult, const = coh_const)
  stab  <- stability_composite(W, a = stab_a, b = stab_b, kappa = stab_kappa)
  interp<- interpretability_score(P, W)
  eff   <- c(weighted = 0.90, multiplicative = 0.80, dirichlet = 0.75, adaptive = 0.65)[[method]]
  composite <- 0.30*coh + 0.25*stab + 0.25*interp + 0.20*eff

  years <- as.integer(rownames(P))
  prior_tbl <- tibble::as_tibble(P) %>% dplyr::mutate(Year = years) %>% dplyr::relocate(Year)
  post_tbl  <- tibble::as_tibble(W) %>% dplyr::mutate(Year = years) %>% dplyr::relocate(Year)
  like_tbl  <- tibble::as_tibble(LT) %>% dplyr::mutate(Year = years) %>% dplyr::relocate(Year)
  L_tbl     <- tibble::tibble(Sector = colnames(P), L = as.numeric(L))

  metrics <- tibble::tibble(
    method = method, lambda = lambda, gamma = gamma,
    coh_mult = coh_mult, coh_const = coh_const,
    stab_a = stab_a, stab_b = stab_b, stab_kappa = stab_kappa,
    likelihood_pattern = likelihood_pattern,
    coherence = round(coh, 4), stability = round(stab, 4),
    interpretability = round(interp, 4), efficiency = eff,
    composite = round(composite, 4), T = Tn, K = K
  )

  list(
    years = years, industries = colnames(P),
    prior = prior_tbl, likelihood_t = like_tbl, likelihood = L_tbl,
    posterior = post_tbl, metrics = metrics
  )
}

#' Save disaggregation results to disk
#'
#' Writes CSV extracts and a single Excel workbook with the key outputs from
#' \code{\link{bayesian_disaggregate}}.
#'
#' @param res A result object returned by \code{bayesian_disaggregate()}.
#' @param out_dir Directory where files will be written. Created if missing.
#'
#' @return (Invisibly) the path to the Excel file written.
#' @examples
#' \dontrun{
#' res <- bayesian_disaggregate("CPI.xlsx","WEIGHTS.xlsx")
#' save_results(res, "out")
#' }
#' @importFrom utils write.csv
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @export
save_results <- function(res, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # CSVs
  utils::write.csv(res$prior,        file.path(out_dir, "prior_P.csv"),                 row.names = FALSE)
  utils::write.csv(res$likelihood_t, file.path(out_dir, "likelihood_time_Lt.csv"),      row.names = FALSE)
  utils::write.csv(res$likelihood,   file.path(out_dir, "likelihood_vector_L.csv"),     row.names = FALSE)
  utils::write.csv(res$posterior,    file.path(out_dir, "posterior_W.csv"),             row.names = FALSE)
  utils::write.csv(res$metrics,      file.path(out_dir, "metrics_summary.csv"),         row.names = FALSE)

  # Excel workbook
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Metrics");       openxlsx::writeData(wb, "Metrics",      res$metrics)
  openxlsx::addWorksheet(wb, "Prior_P");       openxlsx::writeData(wb, "Prior_P",      res$prior)
  openxlsx::addWorksheet(wb, "Likelihood_t");  openxlsx::writeData(wb, "Likelihood_t", res$likelihood_t)
  openxlsx::addWorksheet(wb, "Likelihood_L");  openxlsx::writeData(wb, "Likelihood_L", res$likelihood)
  openxlsx::addWorksheet(wb, "Posterior_W");   openxlsx::writeData(wb, "Posterior_W",  res$posterior)

  out_xlsx <- file.path(out_dir, "Resumen_Disagg.xlsx")
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

  log_msg("INFO", "Results saved to:", normalizePath(out_xlsx))
  invisible(out_xlsx)
}