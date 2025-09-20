#' Coherence score (prior → posterior alignment improvement)
#'
#' Measures how much the posterior \eqn{W} improves alignment with the
#' sectoral signal \eqn{L} relative to the prior \eqn{P}. We compute the
#' correlation increment \eqn{\Delta\rho = \max(0,\rho(W,L)-\rho(P,L))} using
#' \code{robust_cor} (chooses Pearson/Spearman by larger absolute value), then
#' map it to \eqn{[0,1]} with \code{mult} and \code{const}:
#' \deqn{\mathrm{score}=\min\{1,\ \mathrm{const}+\mathrm{mult}\cdot\Delta\rho\}.}
#'
#' @param P Prior matrix (\eqn{T \times K}); rows should sum to 1 (approximately).
#' @param W Posterior matrix (\eqn{T \times K}); rows should sum to 1 (approximately).
#' @param L Likelihood vector (length \eqn{K}), non-negative and summing to 1.
#' @param mult Non-negative multiplier applied to the correlation increment (default \code{3.0}).
#' @param const Constant offset in \eqn{[0,1]} (default \code{0.5}).
#'
#' @return Scalar coherence score in \eqn{[0,1]}.
#'
#' @examples
#' T <- 6; K <- 4
#' P <- matrix(runif(T*K), T); P <- P/rowSums(P)
#' W <- matrix(runif(T*K), T); W <- W/rowSums(W)
#' L <- runif(K); L <- L/sum(L)
#' coherence_score(P, W, L)
#'
#' @export
coherence_score <- function(P, W, L, mult = 3.0, const = 0.5) {
  if (!is.matrix(P) || !is.matrix(W)) stop("P and W must be matrices.")
  if (!all(dim(P) == dim(W))) stop("P and W must have the same dimensions.")
  if (!is.numeric(L) || length(L) != ncol(P)) stop("L must be numeric of length K = ncol(P).")
  if (!is.numeric(mult) || mult < 0) stop("mult must be a non-negative number.")
  if (!is.numeric(const) || const < 0 || const > 1) stop("const must be in [0,1].")
  if (any(!is.finite(P)) || any(!is.finite(W)) || any(!is.finite(L))) {
    stop("P, W, and L must contain only finite values.")
  }

  pm <- colMeans(P)
  wm <- colMeans(W)
  rho_prior <- robust_cor(pm, L)
  rho_post  <- robust_cor(wm, L)
  delta <- max(0, rho_post - rho_prior)
  pmin(1.0, mult * delta + const)
}

#' Numerical stability (exponential penalty)
#'
#' Penalizes numerical issues in \eqn{W}: (i) deviation of each row sum from 1,
#' measured by the mean absolute deviation, and (ii) count of negative entries.
#' The score is \eqn{\exp(-(a\,\mathrm{msd} + b\,\#\mathrm{neg}))}.
#'
#' @param W Posterior matrix (\eqn{T \times K}).
#' @param a Sensitivity for row-sum deviation (default \code{1000}).
#' @param b Sensitivity for negative counts (default \code{10}).
#'
#' @return Scalar numerical stability score in \eqn{[0,1]}.
#'
#' @examples
#' W <- matrix(runif(20), 5); W <- W/rowSums(W)
#' numerical_stability_exp(W)
#'
#' @export
numerical_stability_exp <- function(W, a = 1000, b = 10) {
  if (!is.matrix(W)) stop("W must be a matrix.")
  if (any(!is.finite(W))) stop("W must contain only finite values.")
  if (!is.numeric(a) || a < 0 || !is.numeric(b) || b < 0) {
    stop("a and b must be non-negative numbers.")
  }

  msd  <- mean(abs(rowSums(W) - 1))
  nneg <- sum(W < 0, na.rm = TRUE)
  as.numeric(exp(-(a * msd + b * nneg)))
}

#' Temporal stability (smoothness over time)
#'
#' Rewards smooth posteriors by penalizing average absolute period-to-period
#' changes per sector. Maps to \eqn{[0,1]} via \eqn{1/(1 + \kappa \cdot \mathrm{mv})}.
#'
#' @param W Posterior matrix (\eqn{T \times K}).
#' @param kappa Sensitivity to average absolute change (default \code{50}).
#'
#' @return Scalar temporal stability score in \eqn{[0,1]}. If \eqn{T<2}, returns \code{0.8}.
#'
#' @examples
#' W <- matrix(runif(30), 6); W <- W/rowSums(W)
#' temporal_stability(W, kappa = 40)
#'
#' @export
temporal_stability <- function(W, kappa = 50) {
  if (!is.matrix(W)) stop("W must be a matrix.")
  if (any(!is.finite(W))) stop("W must contain only finite values.")
  if (!is.numeric(kappa) || kappa < 0) stop("kappa must be a non-negative number.")
  if (nrow(W) < 2) return(0.8)

  # mean absolute delta per sector, then average across sectors
  mv <- mean(apply(W, 2, function(col) mean(abs(diff(col)), na.rm = TRUE)), na.rm = TRUE)
  1 / (1 + kappa * mv)
}

#' Composite stability score (numerical and temporal)
#'
#' Convex combination of numerical and temporal stability:
#' \deqn{0.6 \cdot \mathrm{numerical\_stability\_exp}(W) + 0.4 \cdot \mathrm{temporal\_stability}(W).}
#'
#' @param W Posterior matrix (\eqn{T \times K}).
#' @param a Sensitivity for row-sum deviation in numerical part (default \code{1000}).
#' @param b Sensitivity for negatives in numerical part (default \code{10}).
#' @param kappa Sensitivity for temporal smoothness (default \code{50}).
#'
#' @return Scalar composite stability score in \eqn{[0,1]}.
#'
#' @examples
#' W <- matrix(runif(20), 5); W <- W/rowSums(W)
#' stability_composite(W)
#'
#' @seealso \code{\link{numerical_stability_exp}}, \code{\link{temporal_stability}}
#' @export
stability_composite <- function(W, a = 1000, b = 10, kappa = 50) {
  0.6 * numerical_stability_exp(W, a, b) +
    0.4 * temporal_stability(W, kappa)
}

#' Interpretability score (structure preservation + plausibility)
#'
#' Balances two ideas: (i) \emph{preservation} of the average sectoral structure
#' (correlation between \code{colMeans(P)} y \code{colMeans(W)}; truncated at 0),
#' and (ii) \emph{plausibility} of relative changes
#' \eqn{|W\_bar - P\_bar| / (P\_bar + \varepsilon)}, summarized by the 90th
#' percentile (or the maximum). The score is
#' \deqn{0.6\cdot \mathrm{preservation} + 0.4\cdot \frac{1}{1+2\cdot \mathrm{change}}.}
#'
#' @param P Prior matrix (\eqn{T \times K}).
#' @param W Posterior matrix (\eqn{T \times K}).
#' @param use_q90 If \code{TRUE} (default), use 90th percentile of relative changes;
#'   if \code{FALSE}, use the maximum.
#'
#' @return Scalar interpretability score in \eqn{[0,1]}.
#'
#' @examples
#' T <- 6; K <- 5
#' P <- matrix(runif(T*K), T); P <- P/rowSums(P)
#' W <- matrix(runif(T*K), T); W <- W/rowSums(W)
#' interpretability_score(P, W)
#'
#' @importFrom stats quantile
#' @export
interpretability_score <- function(P, W, use_q90 = TRUE) {
  if (!is.matrix(P) || !is.matrix(W)) stop("P and W must be matrices.")
  if (!all(dim(P) == dim(W))) stop("P and W must have the same dimensions.")
  if (any(!is.finite(P)) || any(!is.finite(W))) stop("P and W must contain only finite values.")

  pm <- colMeans(P); wm <- colMeans(W)
  preservation <- max(0, robust_cor(pm, wm))

  eps <- 1e-10
  relchg <- abs(wm - pm) / (pm + eps)
  change <- if (use_q90) stats::quantile(relchg, 0.9, na.rm = TRUE) else max(relchg, na.rm = TRUE)
  plausibility <- 1 / (1 + 2 * change)

  0.6 * preservation + 0.4 * plausibility
}