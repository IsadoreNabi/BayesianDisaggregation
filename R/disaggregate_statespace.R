# =============================================================================
# Canonical evidence-based disaggregation engine (state-space, Model A).
# This is the single, reusable implementation: convergenceDFM imports it instead
# of carrying its own deterministic blend.
# =============================================================================

#' Default prior scales for the state-space disaggregation
#'
#' Weakly-informative scales on the raw index level, derived from the observed
#' aggregate so the model is location/scale aware. Override any field by passing
#' a (partial) named list to \code{disaggregate_statespace(priors = ...)}.
#'
#' @param cpi Numeric vector; the observed aggregate index (level).
#' @return Named list of prior scales (see \code{\link{disaggregate_statespace}}).
#' @export
disagg_default_priors <- function(cpi) {
  if (!is.numeric(cpi) || length(cpi) < 2L || any(!is.finite(cpi)) || any(cpi <= 0)) {
    stop("`cpi` must be a finite, strictly positive numeric vector of length >= 2.",
         call. = FALSE)
  }
  sd_cpi <- stats::sd(cpi)
  list(
    phi1_center     = cpi[1],            # anchor the t = 1 cross-section at the aggregate
    s_omega_struct  = 0.5,               # cross-sectional log-level dispersion (estimable)
    s_delta_mu      = 0.1,               # common annual log-drift sd
    s_delta_sigma   = 0.1,               # drift dispersion across sectors
    log_tau_loc     = log(0.05),         # innovation scale ~ 5% per period (log scale)
    s_log_tau_mu    = 1.0,
    s_log_tau_sigma = 0.5,
    s_sigma_cpi     = max(.Machine$double.eps, 0.1 * sd_cpi)  # observation noise scale
  )
}

#' Evidence-based Bayesian disaggregation (state-space; canonical engine)
#'
#' Disaggregates an observed aggregate index (CPI) into \eqn{K} latent sectoral
#' price indices \eqn{\varphi_{t,k}} with a Bayesian state-space model in which
#' the aggregate enters as a genuine observation density (not a renormalization
#' identity). The model couples a random-walk-with-drift transition in
#' \eqn{\log\varphi} (partial pooling on the drift and the innovation scale), an
#' estimable cross-sectional concentration, and a Student-t (or Gaussian)
#' observation \eqn{cpi_t \mid \varphi \sim \mathrm{Student\text{-}t}(\nu,
#' \sum_k W_{t,k}\varphi_{t,k}, \sigma)}. See \code{vignette("evidence-based-disaggregation")}.
#'
#' The returned posterior draws of \eqn{\varphi} (a \code{[T, K, draws]} array)
#' are exactly the multiple-imputation input consumed by
#' \code{bayesianOU::fit_ou_nested_mi()}, propagating the disaggregation
#' uncertainty into the downstream nested-OU analysis (Rubin's rule).
#'
#' @param cpi Numeric vector (length \eqn{T}); the observed aggregate index in
#'   levels (e.g. CPI re-indexed to a common base). Strictly positive.
#' @param W Numeric matrix (\eqn{T \times K}); the (known) VAB aggregation
#'   weights, rows summing to 1 (small deviations are renormalized). Sector
#'   columns must align with the desired output ordering.
#' @param years Optional integer vector (length \eqn{T}); period labels.
#' @param industries Optional character vector (length \eqn{K}); sector labels.
#'   Defaults to \code{colnames(W)} when present.
#' @param student_obs Logical; if \code{TRUE} (default) the observation is
#'   Student-t (robust to aggregate outliers), otherwise Gaussian.
#' @param priors Optional named list overriding \code{\link{disagg_default_priors}}.
#' @param chains,iter,warmup,thin Sampler controls (HMC/NUTS). Defaults
#'   \code{4 / 2000 / 1000 / 1}.
#' @param cores Integer; parallel chains. Default \code{min(chains, detectCores())}.
#' @param adapt_delta,max_treedepth NUTS tuning. Defaults \code{0.95 / 12}.
#' @param seed Integer RNG seed. Default \code{1234}.
#' @param init Sampler init; a numeric scalar is an init radius (cmdstanr) or is
#'   translated to \code{init_r} (rstan). Default \code{0.5}.
#' @param keep_fit Logical; keep the raw Stan fit object in the result. Default
#'   \code{TRUE} (needed to draw further quantities or run LOO).
#' @param verbose Logical; print progress. Default \code{FALSE}.
#'
#' @return An object of class \code{"disagg_statespace"}: a list with
#'   \describe{
#'     \item{phi_draws}{\code{[T, K, draws]} numeric array of posterior draws of
#'       \eqn{\varphi} (the multiple-imputation input for the nested OU).}
#'     \item{phi_summary}{List of \eqn{T \times K} matrices \code{median},
#'       \code{q2.5}, \code{q97.5} (credible bands per sector and period).}
#'     \item{agg_summary}{\eqn{T \times 3} matrix: posterior median and 95\% band
#'       of the fitted aggregate \eqn{\sum_k W\varphi} (against which \code{cpi}
#'       is the evidence).}
#'     \item{years, industries}{Period and sector labels.}
#'     \item{diagnostics}{\code{rhat_max}, \code{divergences}.}
#'     \item{stan_fit}{The Stan fit (if \code{keep_fit}).}
#'     \item{config}{Sampler/prior configuration and \code{T}, \code{K}.}
#'   }
#'
#' @seealso \code{\link{disaggregate_conjugate}} (closed-form Bayesian baseline),
#'   \code{\link{disaggregate_from_files}}, \code{\link{simulate_disagg}}.
#' @examples
#' \dontrun{
#' set.seed(1)
#' sim <- simulate_disagg(T = 30, K = 4)
#' fit <- disaggregate_statespace(sim$cpi, sim$W, chains = 2, iter = 800)
#' dim(fit$phi_draws)            # T x K x draws
#' }
#' @export
disaggregate_statespace <- function(cpi, W, years = NULL, industries = NULL,
                                    student_obs = TRUE, priors = NULL,
                                    chains = 4L, iter = 2000L, warmup = 1000L,
                                    thin = 1L,
                                    cores = NULL, adapt_delta = 0.95,
                                    max_treedepth = 12L, seed = 1234L,
                                    init = 0.5, keep_fit = TRUE, verbose = FALSE) {
  # ---- validate inputs --------------------------------------------------------
  cpi <- as.numeric(cpi)
  if (length(cpi) < 2L || any(!is.finite(cpi)) || any(cpi <= 0)) {
    stop("`cpi` must be a finite, strictly positive numeric vector of length >= 2.",
         call. = FALSE)
  }
  W <- as.matrix(W)
  if (!is.numeric(W) || any(!is.finite(W)) || any(W < 0)) {
    stop("`W` must be a finite, non-negative numeric matrix.", call. = FALSE)
  }
  Tn <- length(cpi); K <- ncol(W)
  if (nrow(W) != Tn) {
    stop(sprintf("`W` has %d rows but `cpi` has length %d; they must match.",
                 nrow(W), Tn), call. = FALSE)
  }
  W <- row_norm1(W)                                   # enforce the simplex per period

  industries <- industries %||% colnames(W) %||% paste0("sector_", seq_len(K))
  if (length(industries) != K) {
    stop("`industries` must have length K = ncol(W).", call. = FALSE)
  }
  years <- years %||% seq_len(Tn)
  if (length(years) != Tn) stop("`years` must have length T = length(cpi).", call. = FALSE)

  cores <- as.integer(cores %||% min(as.integer(chains), parallel::detectCores()))

  # ---- assemble priors and Stan data -----------------------------------------
  pr <- disagg_default_priors(cpi)
  if (!is.null(priors)) {
    if (!is.list(priors)) stop("`priors` must be a named list.", call. = FALSE)
    unknown <- setdiff(names(priors), names(pr))
    if (length(unknown)) {
      stop("Unknown prior fields: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    pr[names(priors)] <- priors
  }

  stan_dat <- .build_stan_dat(cpi, W, pr, student_obs)

  # ---- sample -----------------------------------------------------------------
  fit <- .run_stan_disagg(stan_dat, chains = chains, iter = iter, warmup = warmup,
                          thin = thin, cores = cores, adapt_delta = adapt_delta,
                          max_treedepth = max_treedepth, seed = seed, init = init,
                          verbose = verbose)

  # ---- extract ----------------------------------------------------------------
  phi_arr <- .matrix_param_to_array(fit, "phi", Tn, K)
  if (is.null(phi_arr)) stop("Could not extract `phi` draws from the Stan fit.", call. = FALSE)
  dimnames(phi_arr) <- list(years, industries, NULL)

  phi_summary <- .summarize_phi(phi_arr, years, industries)
  agg_summary <- .summarize_agg(fit, Tn, years)

  out <- list(
    phi_draws   = phi_arr,
    phi_summary = phi_summary,
    agg_summary = agg_summary,
    cpi         = cpi,
    W           = W,
    years       = years,
    industries  = industries,
    diagnostics = .stan_diagnostics(fit),
    stan_fit    = if (keep_fit) fit else NULL,
    config      = list(T = Tn, K = K, student_obs = isTRUE(student_obs),
                       priors = pr, chains = chains, iter = iter, warmup = warmup,
                       thin = thin, seed = seed, n_draws = dim(phi_arr)[3])
  )
  class(out) <- c("disagg_statespace", "list")
  out
}

#' Assemble the Stan data list from inputs and resolved priors
#'
#' Single source of truth for the \code{stan_dat} layout, shared by the fit
#' driver and the golden generator (so the frozen fixture uses exactly the data
#' the model is fitted with).
#' @keywords internal
#' @noRd
.build_stan_dat <- function(cpi, W, pr, student_obs) {
  list(
    T = length(cpi), K = ncol(W), cpi = cpi, W = W,
    phi1_center = pr$phi1_center, s_omega_struct = pr$s_omega_struct,
    s_delta_mu = pr$s_delta_mu, s_delta_sigma = pr$s_delta_sigma,
    log_tau_loc = pr$log_tau_loc, s_log_tau_mu = pr$s_log_tau_mu,
    s_log_tau_sigma = pr$s_log_tau_sigma, s_sigma_cpi = pr$s_sigma_cpi,
    student_obs = as.integer(isTRUE(student_obs))
  )
}

#' Summarize phi draws into median / 95% credible bands (T x K matrices)
#' @keywords internal
#' @noRd
.summarize_phi <- function(phi_arr, years, industries) {
  d <- dim(phi_arr); Tn <- d[1]; K <- d[2]
  med <- matrix(NA_real_, Tn, K, dimnames = list(years, industries))
  lo <- med; hi <- med
  for (k in seq_len(K)) {
    q <- t(apply(phi_arr[, k, , drop = TRUE], 1, stats::quantile,
                 probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    lo[, k] <- q[, 1]; med[, k] <- q[, 2]; hi[, k] <- q[, 3]
  }
  list(median = med, q2.5 = lo, q97.5 = hi)
}

#' Summarize the fitted aggregate sum_k W phi (median + 95% band)
#' @keywords internal
#' @noRd
.summarize_agg <- function(fit, Tn, years) {
  M <- .draws_matrix(fit, "agg")
  if (is.null(M)) return(NULL)
  q <- t(apply(M, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  rownames(q) <- years; colnames(q) <- c("q2.5", "median", "q97.5")
  q
}

# ---- S3 methods -------------------------------------------------------------

#' @export
print.disagg_statespace <- function(x, ...) {
  cat("<disagg_statespace>  evidence-based state-space disaggregation\n")
  cat(sprintf("  periods T = %d, sectors K = %d, posterior draws = %d\n",
              x$config$T, x$config$K, x$config$n_draws))
  cat(sprintf("  observation: %s\n",
              if (x$config$student_obs) "Student-t (robust)" else "Gaussian"))
  cat(sprintf("  diagnostics: rhat_max = %.3f, divergences = %s\n",
              x$diagnostics$rhat_max %||% NA_real_,
              format(x$diagnostics$divergences %||% NA_real_)))
  invisible(x)
}

#' @export
summary.disagg_statespace <- function(object, ...) {
  med <- object$phi_summary$median
  structure(list(
    T = object$config$T, K = object$config$K,
    industries = object$industries,
    phi_median_range = range(med, na.rm = TRUE),
    diagnostics = object$diagnostics,
    config = object$config
  ), class = "summary.disagg_statespace")
}

#' @export
print.summary.disagg_statespace <- function(x, ...) {
  cat("Evidence-based state-space disaggregation (summary)\n")
  cat(sprintf("  T = %d periods, K = %d sectors\n", x$T, x$K))
  cat(sprintf("  phi (median) range: [%.3f, %.3f]\n",
              x$phi_median_range[1], x$phi_median_range[2]))
  cat(sprintf("  rhat_max = %.3f, divergences = %s\n",
              x$diagnostics$rhat_max %||% NA_real_,
              format(x$diagnostics$divergences %||% NA_real_)))
  invisible(x)
}
