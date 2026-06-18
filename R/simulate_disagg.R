# =============================================================================
# Synthetic data generator for the state-space disaggregation model.
# Used by the recovery tests, the goldens and the documentation examples. The
# DGP is the SAME as the Stan model (random-walk-with-drift in log phi, partial
# pooling on drift and innovation scale, VAB-weighted aggregate observation), so
# parameter/path recovery is a well-posed question.
# =============================================================================

#' Simulate from the state-space disaggregation DGP
#'
#' Generates a synthetic aggregate index \code{cpi}, the (known) VAB weights
#' \code{W}, and the latent sectoral price-index paths \code{phi_true} from the
#' same data-generating process as \code{\link{disaggregate_statespace}}. The
#' innovation scale is kept modest so the log random walk stays in a numerically
#' stable region over the simulated horizon (the same care taken in the sibling
#' OU simulator).
#'
#' @param T Integer; number of periods.
#' @param K Integer; number of sectors.
#' @param phi1_center Numeric; central level of the initial cross-section.
#' @param omega_struct Numeric; cross-sectional log-level dispersion at t = 1.
#' @param delta_mu,delta_sigma Common drift and its cross-sector dispersion.
#' @param tau_mu,tau_sigma Geometric mean innovation scale and log-dispersion
#'   (so \eqn{\tau_k = \tau_{mu}\exp(\tau_{sigma} z_k)}).
#' @param sigma_cpi Observation noise scale on the aggregate.
#' @param nu Student-t degrees of freedom of the observation (\code{Inf} = Gaussian).
#' @param seed Integer RNG seed.
#'
#' @return A list with \code{cpi} (length \eqn{T}), \code{W} (\eqn{T \times K},
#'   rows sum to 1), \code{phi_true} (\eqn{T \times K}), \code{agg_true}
#'   (length \eqn{T}), and \code{params} (the true scalar/vector parameters).
#' @examples
#' sim <- simulate_disagg(T = 20, K = 3, seed = 42)
#' str(sim$params)
#' @export
simulate_disagg <- function(T = 40L, K = 5L, phi1_center = 100,
                            omega_struct = 0.3,
                            delta_mu = 0.02, delta_sigma = 0.01,
                            tau_mu = 0.04, tau_sigma = 0.3,
                            sigma_cpi = 1.0, nu = Inf, seed = 1234L) {
  T <- as.integer(T); K <- as.integer(K)
  if (T < 2L || K < 1L) stop("Need T >= 2 and K >= 1.", call. = FALSE)
  set.seed(seed)

  delta <- delta_mu + delta_sigma * stats::rnorm(K)
  tau   <- tau_mu * exp(tau_sigma * stats::rnorm(K))

  log_phi <- matrix(NA_real_, T, K)
  log_phi[1, ] <- log(phi1_center) + omega_struct * stats::rnorm(K)
  for (t in 2:T) {
    log_phi[t, ] <- log_phi[t - 1, ] + delta + tau * stats::rnorm(K)
  }
  phi_true <- exp(log_phi)

  # VAB weights: a smooth, persistent simplex per period (Dirichlet-ish via
  # a slowly drifting log-share), so they are realistic but known.
  base_share <- stats::runif(K, 0.5, 1.5)
  log_share <- matrix(log(base_share), T, K, byrow = TRUE) +
    apply(matrix(stats::rnorm(T * K, sd = 0.02), T, K), 2, cumsum)
  W <- exp(log_share)
  W <- W / rowSums(W)

  agg_true <- rowSums(W * phi_true)
  noise <- if (is.finite(nu)) sigma_cpi * stats::rt(T, df = nu) else stats::rnorm(T, sd = sigma_cpi)
  cpi <- agg_true + noise
  cpi <- pmax(cpi, .Machine$double.eps)             # keep strictly positive

  list(
    cpi = as.numeric(cpi), W = W, phi_true = phi_true, agg_true = agg_true,
    params = list(delta = delta, tau = tau, delta_mu = delta_mu,
                  delta_sigma = delta_sigma, tau_mu = tau_mu, tau_sigma = tau_sigma,
                  omega_struct = omega_struct, sigma_cpi = sigma_cpi, nu = nu,
                  phi1_center = phi1_center, T = T, K = K)
  )
}
