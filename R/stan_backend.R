# =============================================================================
# Stan backend plumbing for the canonical (state-space) disaggregation engine.
# Mirrors the conventions used by the sibling package bayesianOU so that both
# behave identically with respect to backend selection, model-file resolution
# and backend-agnostic draw extraction.
# =============================================================================

#' Null-coalescing helper
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Locate the canonical state-space Stan model file
#'
#' Resolves the path to \code{inst/stan/disagg_statespace.stan} (the single
#' source of truth for the evidence-based engine), both for the installed
#' package and during development (source-tree fallbacks).
#'
#' @return Character path to the \code{.stan} file.
#' @keywords internal
#' @noRd
.disagg_stan_file_path <- function() {
  p <- system.file("stan", "disagg_statespace.stan", package = "BayesianDisaggregation")
  if (nzchar(p) && file.exists(p)) return(p)

  candidates <- c(
    file.path("inst", "stan", "disagg_statespace.stan"),
    file.path("..", "inst", "stan", "disagg_statespace.stan"),
    file.path("..", "..", "inst", "stan", "disagg_statespace.stan")
  )
  for (cand in candidates) if (file.exists(cand)) return(normalizePath(cand))

  stop("Could not locate 'disagg_statespace.stan'. Reinstall the package or run ",
       "from the package root.", call. = FALSE)
}

#' Stan code for the canonical disaggregation model
#'
#' Returns the complete Stan code of the evidence-based state-space model read
#' from the canonical file (single source of truth; no embedded duplicate).
#'
#' @return Character string with the Stan model code.
#' @examples
#' code <- disagg_stan_code()
#' cat(substr(code, 1, 200))
#' @export
disagg_stan_code <- function() {
  paste(readLines(.disagg_stan_file_path(), warn = FALSE), collapse = "\n")
}

#' Detect an available Stan backend (cmdstanr preferred, then rstan)
#' @keywords internal
#' @noRd
check_stan_backend <- function(verbose = FALSE) {
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    ok <- tryCatch(!is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
                   error = function(e) FALSE)
    if (ok) {
      if (verbose) message("Using cmdstanr backend")
      return("cmdstanr")
    }
  }
  if (requireNamespace("rstan", quietly = TRUE)) {
    if (verbose) message("Using rstan backend")
    return("rstan")
  }
  "none"
}

#' Compile and sample the canonical state-space model (cmdstanr or rstan)
#' @keywords internal
#' @noRd
.run_stan_disagg <- function(stan_dat, chains, iter, warmup, thin, cores,
                             adapt_delta, max_treedepth, seed, init, verbose) {
  backend <- check_stan_backend(verbose)
  if (backend == "none") {
    stop("A Stan backend is required (install 'cmdstanr' or 'rstan').", call. = FALSE)
  }
  par_chains <- max(1L, min(as.integer(chains), as.integer(cores)))

  if (backend == "cmdstanr") {
    if (verbose) message("Compiling state-space model with cmdstanr")
    cache_dir <- tryCatch({
      d <- tools::R_user_dir("BayesianDisaggregation", "cache")
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(d)) d else tempdir()
    }, error = function(e) tempdir())
    mod <- cmdstanr::cmdstan_model(.disagg_stan_file_path(), dir = cache_dir,
                                   pedantic = FALSE)
    if (verbose) message("Running MCMC sampling")
    return(mod$sample(
      data = stan_dat, chains = chains, parallel_chains = par_chains,
      iter_warmup = warmup, iter_sampling = iter - warmup, thin = thin,
      seed = seed, refresh = if (verbose) 200 else 0,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth, init = init
    ))
  }

  if (verbose) message("Compiling state-space model with rstan")
  sm <- rstan::stan_model(model_code = disagg_stan_code())
  rstan_init <- "random"; rstan_init_r <- 2
  if (!is.null(init)) {
    if (is.numeric(init) && length(init) == 1L) rstan_init_r <- init else rstan_init <- init
  }
  rstan::sampling(
    sm, data = stan_dat, chains = chains, iter = iter, warmup = warmup,
    thin = thin, seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = if (verbose) 200 else 0, init = rstan_init, init_r = rstan_init_r
  )
}

#' Draws as a [draws x columns] matrix for a parameter, backend-agnostic
#' @keywords internal
#' @noRd
.draws_matrix <- function(fit, p) {
  if (inherits(fit, "CmdStanMCMC")) {
    df <- tryCatch(fit$draws(p, format = "df"), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    keep <- grep(sprintf("^%s(\\[|$)", p), names(df), value = TRUE)
    keep <- setdiff(keep, c(".chain", ".iteration", ".draw"))
    if (length(keep) == 0) return(NULL)
    return(as.matrix(df[, keep, drop = FALSE]))
  }
  out <- tryCatch(rstan::extract(fit, pars = p, permuted = TRUE)[[1]],
                  error = function(e) NULL)
  if (is.null(out)) return(NULL)
  if (is.null(dim(out))) return(matrix(out, ncol = 1))
  if (length(dim(out)) == 2L) return(out)
  matrix(out, nrow = dim(out)[1])
}

#' Reshape draws of a T x K matrix parameter into a [T, K, D] array
#'
#' Columns of \code{.draws_matrix(fit, name)} are named \code{name[t,k]} (t
#' fastest, matching cmdstanr/rstan order). Returns a \code{[T, K, D]} array with
#' the imputations on the third margin, the contract consumed by
#' \code{bayesianOU::fit_ou_nested_mi}.
#' @keywords internal
#' @noRd
.matrix_param_to_array <- function(fit, name, Tn, K) {
  M <- .draws_matrix(fit, name)
  if (is.null(M)) return(NULL)
  cn <- colnames(M)
  tt <- as.integer(sub(sprintf("^%s\\[(\\d+),(\\d+)\\]$", name), "\\1", cn))
  kk <- as.integer(sub(sprintf("^%s\\[(\\d+),(\\d+)\\]$", name), "\\2", cn))
  D <- nrow(M)
  arr <- array(NA_real_, dim = c(Tn, K, D))
  for (j in seq_along(cn)) arr[tt[j], kk[j], ] <- M[, j]
  arr
}

#' Backend-agnostic convergence diagnostics (rhat_max, divergences)
#' @keywords internal
#' @noRd
.stan_diagnostics <- function(fit) {
  if (inherits(fit, "CmdStanMCMC")) {
    rhat_max <- tryCatch({
      s <- fit$summary()
      max(s$rhat, na.rm = TRUE)
    }, error = function(e) NA_real_)
    div <- tryCatch({
      dg <- fit$diagnostic_summary(quiet = TRUE)
      sum(dg$num_divergent)
    }, error = function(e) NA_real_)
    return(list(rhat_max = rhat_max, divergences = div))
  }
  rhat_max <- tryCatch(max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE),
                       error = function(e) NA_real_)
  div <- tryCatch({
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  }, error = function(e) NA_real_)
  list(rhat_max = rhat_max, divergences = div)
}
