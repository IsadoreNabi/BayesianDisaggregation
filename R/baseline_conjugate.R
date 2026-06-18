# =============================================================================
# Closed-form conjugate baseline: linear-Gaussian state-space disaggregation.
#
# The canonical engine (disaggregate_statespace) works on log phi with Student-t
# observations and a partial-pooling hierarchy, and is sampled by HMC. This
# baseline is its closest CLOSED-FORM analogue: a linear-Gaussian random walk in
# LEVELS with a scalar VAB-weighted observation, whose exact posterior is the
# Kalman/RTS smoother. The contrast is deliberate and documented:
#   * shared: random-walk transition + aggregate observation = genuine evidence;
#   * traded away for closed form: positivity (levels may go slightly negative),
#     heavy-tail robustness, and the cross-sector/innovation-scale hierarchy.
# Retained for comparison only; the canonical engine is the reusable one.
# =============================================================================

#' Kalman filter + RTS smoother for a diagonal random walk with scalar observation
#'
#' State \eqn{\varphi_t\in\mathbb{R}^K} follows \eqn{\varphi_t=\varphi_{t-1}+
#' \mathcal{N}(0,Q)}; observation \eqn{y_t=w_t^\top\varphi_t+\mathcal{N}(0,r)}.
#' Returns the smoothed means and marginal variances.
#' @keywords internal
#' @noRd
.kalman_rts <- function(y, W, m0, P0, Q, r) {
  Tn <- length(y); K <- ncol(W)
  a_pred <- vector("list", Tn); P_pred <- vector("list", Tn)
  a_filt <- vector("list", Tn); P_filt <- vector("list", Tn)

  a_pred[[1]] <- m0; P_pred[[1]] <- P0
  for (t in seq_len(Tn)) {
    w <- W[t, ]
    Pw <- P_pred[[t]] %*% w
    F  <- as.numeric(crossprod(w, Pw)) + r
    v  <- y[t] - sum(w * a_pred[[t]])
    Kg <- Pw / F
    a_filt[[t]] <- a_pred[[t]] + as.numeric(Kg) * v
    P_filt[[t]] <- P_pred[[t]] - (Pw %*% t(Pw)) / F
    P_filt[[t]] <- (P_filt[[t]] + t(P_filt[[t]])) / 2          # symmetrize
    if (t < Tn) {
      a_pred[[t + 1]] <- a_filt[[t]]
      P_pred[[t + 1]] <- P_filt[[t]] + Q
    }
  }

  a_sm <- vector("list", Tn); P_sm <- vector("list", Tn)
  a_sm[[Tn]] <- a_filt[[Tn]]; P_sm[[Tn]] <- P_filt[[Tn]]
  for (t in (Tn - 1):1) {
    if (t < 1) break
    J <- P_filt[[t]] %*% solve(P_pred[[t + 1]])
    a_sm[[t]] <- a_filt[[t]] + J %*% (a_sm[[t + 1]] - a_pred[[t + 1]])
    P_sm[[t]] <- P_filt[[t]] + J %*% (P_sm[[t + 1]] - P_pred[[t + 1]]) %*% t(J)
    P_sm[[t]] <- (P_sm[[t]] + t(P_sm[[t]])) / 2
  }

  mean_mat <- do.call(rbind, lapply(a_sm, as.numeric))
  var_mat  <- do.call(rbind, lapply(P_sm, function(P) pmax(diag(P), 0)))
  list(mean = mean_mat, var = var_mat, P_sm = P_sm)
}

#' Conjugate (closed-form) disaggregation baseline
#'
#' Exact linear-Gaussian state-space posterior (Kalman/RTS smoother) for the
#' sectoral price-index levels \eqn{\varphi_{t,k}} given the aggregate index and
#' the VAB weights. Optionally returns joint posterior draws via the
#' Durbin-Koopman simulation smoother (so the draws can also feed
#' \code{bayesianOU::fit_ou_nested_mi}), and a pointwise Gaussian log-likelihood.
#'
#' @param cpi Numeric vector (length \eqn{T}); observed aggregate index (levels).
#' @param W Numeric matrix (\eqn{T \times K}); VAB weights (rows sum to 1).
#' @param years,industries Optional period and sector labels.
#' @param q_frac Random-walk innovation sd as a fraction of \code{sd(cpi)}
#'   (state noise). Default \code{0.10}.
#' @param r_frac Observation sd as a fraction of \code{sd(cpi)}. Default \code{0.05}.
#' @param p0_frac Initial cross-sectional sd as a fraction of \code{cpi[1]}.
#'   Default \code{0.30}.
#' @param n_draws Integer; number of joint posterior draws (simulation smoother).
#'   \code{0} (default) returns only the smoothed mean and bands.
#' @param seed Integer RNG seed (used only when \code{n_draws > 0}).
#'
#' @return An object of class \code{"disagg_conjugate"}: a list with
#'   \code{phi_summary} (median = smoothed mean, \code{q2.5}/\code{q97.5} from the
#'   marginal Gaussian), \code{agg_summary}, \code{loglik} (total Gaussian
#'   log-likelihood of \code{cpi}), \code{phi_draws} (\code{[T, K, n_draws]} or
#'   \code{NULL}), \code{years}, \code{industries}, \code{config}.
#'
#' @seealso \code{\link{disaggregate_statespace}} (canonical engine).
#' @examples
#' sim <- simulate_disagg(T = 25, K = 4, seed = 7)
#' bl  <- disaggregate_conjugate(sim$cpi, sim$W)
#' dim(bl$phi_summary$median)
#' @export
disaggregate_conjugate <- function(cpi, W, years = NULL, industries = NULL,
                                   q_frac = 0.10, r_frac = 0.05, p0_frac = 0.30,
                                   n_draws = 0L, seed = 1234L) {
  cpi <- as.numeric(cpi)
  if (length(cpi) < 2L || any(!is.finite(cpi)) || any(cpi <= 0)) {
    stop("`cpi` must be a finite, strictly positive numeric vector of length >= 2.",
         call. = FALSE)
  }
  W <- row_norm1(as.matrix(W))
  Tn <- length(cpi); K <- ncol(W)
  if (nrow(W) != Tn) stop("`W` rows must match length(cpi).", call. = FALSE)

  industries <- industries %||% colnames(W) %||% paste0("sector_", seq_len(K))
  years <- years %||% seq_len(Tn)

  sd_cpi <- stats::sd(cpi)
  r  <- max(.Machine$double.eps, (r_frac * sd_cpi)^2)
  Q  <- diag((q_frac * sd_cpi)^2, K)
  m0 <- rep(cpi[1], K)
  P0 <- diag((p0_frac * cpi[1])^2, K)

  sm <- .kalman_rts(cpi, W, m0, P0, Q, r)
  med <- sm$mean; sdv <- sqrt(sm$var)
  lo <- med - 1.959964 * sdv; hi <- med + 1.959964 * sdv
  dimnames(med) <- dimnames(lo) <- dimnames(hi) <- list(years, industries)

  agg_mean <- rowSums(W * med)
  agg_var  <- vapply(seq_len(Tn), function(t) {
    as.numeric(crossprod(W[t, ], sm$P_sm[[t]] %*% W[t, ]))
  }, numeric(1))
  agg_summary <- cbind(q2.5 = agg_mean - 1.959964 * sqrt(agg_var),
                       median = agg_mean,
                       q97.5 = agg_mean + 1.959964 * sqrt(agg_var))
  rownames(agg_summary) <- years

  # Gaussian total log-likelihood of the aggregate at the smoothed fit.
  loglik <- sum(stats::dnorm(cpi, agg_mean, sqrt(agg_var + r), log = TRUE))

  phi_draws <- NULL
  n_draws <- as.integer(n_draws)
  if (n_draws > 0L) {
    phi_draws <- .simulation_smoother(cpi, W, m0, P0, Q, r, sm$mean, n_draws, seed)
    dimnames(phi_draws) <- list(years, industries, NULL)
  }

  out <- list(
    phi_summary = list(median = med, q2.5 = lo, q97.5 = hi),
    agg_summary = agg_summary,
    loglik = loglik,
    phi_draws = phi_draws,
    cpi = cpi, W = W, years = years, industries = industries,
    config = list(T = Tn, K = K, q_frac = q_frac, r_frac = r_frac,
                  p0_frac = p0_frac, n_draws = n_draws)
  )
  class(out) <- c("disagg_conjugate", "list")
  out
}

#' Durbin-Koopman simulation smoother for the linear-Gaussian RW state space
#'
#' Returns joint posterior draws \code{[T, K, n_draws]} of \eqn{\varphi} with the
#' correct cross-time/cross-sector covariance (not marginal approximations).
#' @keywords internal
#' @noRd
.simulation_smoother <- function(y, W, m0, P0, Q, r, a_sm_real, n_draws, seed) {
  set.seed(seed)
  Tn <- length(y); K <- ncol(W)
  cP0 <- chol(P0); cQ <- chol(Q)
  out <- array(NA_real_, dim = c(Tn, K, n_draws))
  for (d in seq_len(n_draws)) {
    phi_plus <- matrix(NA_real_, Tn, K)
    phi_plus[1, ] <- m0 + as.numeric(crossprod(cP0, stats::rnorm(K)))
    for (t in 2:Tn) phi_plus[t, ] <- phi_plus[t - 1, ] + as.numeric(crossprod(cQ, stats::rnorm(K)))
    y_plus <- vapply(seq_len(Tn), function(t) {
      sum(W[t, ] * phi_plus[t, ]) + stats::rnorm(1, 0, sqrt(r))
    }, numeric(1))
    sm_plus <- .kalman_rts(y_plus, W, m0, P0, Q, r)
    out[, , d] <- phi_plus - sm_plus$mean + a_sm_real
  }
  out
}

#' @export
print.disagg_conjugate <- function(x, ...) {
  cat("<disagg_conjugate>  closed-form linear-Gaussian baseline (Kalman/RTS)\n")
  cat(sprintf("  periods T = %d, sectors K = %d, joint draws = %d\n",
              x$config$T, x$config$K, x$config$n_draws))
  cat(sprintf("  aggregate Gaussian log-likelihood = %.2f\n", x$loglik))
  invisible(x)
}
