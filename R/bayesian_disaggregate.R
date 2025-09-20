#' Run Bayesian disaggregation
#'
#' Performs Bayesian disaggregation of an aggregated time series (e.g., CPI)
#' into K components using one of four deterministic updating rules:
#' \emph{weighted}, \emph{multiplicative}, \emph{dirichlet}, \emph{adaptive}.
#'
#' @param path_cpi Path to an Excel file containing the aggregate index (CPI).
#'   Must include columns \code{date} and \code{value}. Values are normalized
#'   to the unit simplex if necessary.
#' @param path_weights Path to an Excel file with the baseline weight matrix
#'   (prior), either \eqn{T \times K} or length-\eqn{K} (constant across time).
#'   Columns represent sectors; each row should sum to 1 (small deviations
#'   are corrected).
#' @param method Disaggregation method: \code{"weighted"}, \code{"multiplicative"},
#'   \code{"dirichlet"}, or \code{"adaptive"}.
#' @param lambda Weight parameter for the \code{"weighted"} method (0–1). Ignored otherwise.
#' @param gamma Uncertainty factor for the \code{"dirichlet"} method (>0; smaller = sharper).
#' @param coh_mult Multiplier for the coherence increment \eqn{\Delta \rho}.
#' @param coh_const Constant offset for coherence in \eqn{[0,1]}.
#' @param stab_a Sensitivity parameter for the sum-deviation penalty \eqn{\lvert \sum w - 1 \rvert}.
#' @param stab_b Sensitivity parameter for the negative-values penalty.
#' @param stab_kappa Sensitivity parameter for temporal variation (mean \eqn{\lvert \Delta \rvert}).
#' @param likelihood_pattern Temporal spreading pattern for the likelihood:
#'   \code{"constant"}, \code{"recent"}, \code{"linear"}, or \code{"bell"}.
#'
#' @details
#' \strong{Assumptions:} (i) prior and posterior weights lie on the simplex;
#' (ii) no MCMC; all updates are analytic and deterministic; (iii) missing values
#' are either conservatively imputed or rejected by the \code{read_*} helpers.
#'
#' \strong{Methods:}
#' \itemize{
#'   \item \emph{weighted}: \eqn{W = \lambda P + (1-\lambda)L}, then row-normalize.
#'   \item \emph{multiplicative}: \eqn{W \propto P \odot L}{W \propto P * L}, then row-normalize.
#'   \item \emph{dirichlet}: posterior mean with \eqn{\alpha = (P+L)/\gamma}.
#'   \item \emph{adaptive}: sector-wise mixing from prior volatility (see \code{posterior_adaptive}).
#' }
#'
#' \strong{Scores:} Coherence scales the correlation increment \eqn{\Delta \rho}
#' with \code{coh_mult} and \code{coh_const}. Stability combines an exponential
#' penalty on numerical errors with temporal smoothness controlled by \code{stab_kappa}.
#'
#' @return A \code{list} with:
#' \describe{
#'   \item{\code{prior}}{\eqn{T \times K} prior matrix (rows sum to 1).}
#'   \item{\code{likelihood}}{\eqn{T \times K} likelihood matrix (temporal pattern applied).}
#'   \item{\code{posterior}}{\eqn{T \times K} posterior matrix (rows sum to 1).}
#'   \item{\code{metrics}}{\code{list} with \code{coherence}, \code{stability},
#'         \code{interpretability}, \code{efficiency}, and \code{composite}.}
#'   \item{\code{diagnostics}}{\code{list} with numerical deviations, negatives count,
#'         temporal variation, and logging messages.}
#' }
#'
#' @section File formats:
#' \subsection{CPI (\code{path_cpi})}{
#'   First sheet by default. Required columns: \code{date}, \code{value}.
#'   Dates may be \code{Date}/\code{POSIXct} or ISO text. \code{value} is rescaled if not in [0,1].
#' }
#' \subsection{Weights (\code{path_weights})}{
#'   Either (i) \eqn{T \times K} matrix with sector columns; or (ii) length-\eqn{K}
#'   vector (constant across time). Small deviations from 1 in row sums are corrected.
#' }
#'
#' @examples
#' \donttest{
#' # Synthetic example without files:
#' T <- 6; K <- 4
#' P <- matrix(rep(1/K, T*K), nrow=T)  # flat prior
#' L <- runif(K); L <- L/sum(L)
#' LT <- spread_likelihood(L, T, "recent")
#' W  <- posterior_weighted(P, LT, lambda = 0.7)
#'
#' # Metrics
#' dims <- evaluate_all_dims(P, W, L, mult = 3, const = 0.5, a=1000, b=10, kappa=50)
#' str(dims)
#' }
#'
#' @seealso \code{\link{posterior_weighted}}, \code{\link{posterior_multiplicative}},
#'   \code{\link{posterior_dirichlet}}, \code{\link{posterior_adaptive}},
#'   \code{\link{coherence_score}}, \code{\link{stability_composite}},
#'   \code{\link{interpretability_score}}, \code{\link{spread_likelihood}}
#' @export