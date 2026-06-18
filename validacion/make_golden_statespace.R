# =============================================================================
# Golden generator for the state-space engine (pattern borrowed from bayesianOU).
#
# Freezes, for a small fitted model: (i) a handful of posterior PARAMETER draws,
# (ii) the exact stan_dat, and (iii) the reference pointwise log-likelihood
# recomputed by the model. The companion test recomputes log_lik via
# generate_quantities on the frozen draws+data and demands a BIT-FOR-BIT match,
# so any change to the Stan computation path is caught while serialization
# rounding is isolated (the comparison is gq-output vs gq-output).
#
# Run:  NOT_CRAN=true Rscript validacion/make_golden_statespace.R
# =============================================================================
suppressMessages({
  library(cmdstanr); cmdstanr::set_cmdstan_path(path.expand("~/.cmdstan/cmdstan-2.38.0"))
  library(pkgload); library(posterior)
})
# Run from the package root (Rscript validacion/make_golden_statespace.R).
pkg <- if (dir.exists("R")) "." else stop("Run from the package root.")
pkgload::load_all(pkg, quiet = TRUE)

set.seed(2024)
sim <- simulate_disagg(T = 18, K = 3, seed = 2024)
fit <- disaggregate_statespace(sim$cpi, sim$W, chains = 2, iter = 600, warmup = 300,
                               seed = 2024, student_obs = TRUE)
sf <- fit$stan_fit

param_vars <- c("z_phi1", "omega_struct", "delta_mu", "delta_sigma", "z_delta",
                "log_tau_mu", "log_tau_sigma", "z_tau", "z_eta", "sigma_cpi", "nu_tilde")
draws_par <- sf$draws(variables = param_vars, format = "draws_array")
# keep a small, evenly-spaced set of draws
nd <- posterior::ndraws(draws_par)
keep <- unique(round(seq(1, nd, length.out = 12)))
draws_small <- posterior::subset_draws(draws_par, draw = keep)

pr <- fit$config$priors
stan_dat <- BayesianDisaggregation:::.build_stan_dat(sim$cpi, sim$W, pr, TRUE)

# Reference log_lik via generate_quantities on the FROZEN draws+data, i.e. the
# SAME path the test uses. This isolates code equivalence from CSV serialization
# rounding (both sides pass through the identical write/read), per bayesianOU
# D-IMPL-1. Reading log_lik off the original fit (already CSV-rounded) would
# spuriously fail an exact comparison against a full-precision recompute.
mod <- cmdstanr::cmdstan_model(BayesianDisaggregation:::.disagg_stan_file_path(),
                               dir = tempdir(), pedantic = FALSE)
gq_ref <- mod$generate_quantities(fitted_params = draws_small, data = stan_dat)
ll_ref <- gq_ref$draws(variables = "log_lik", format = "draws_matrix")

golden <- list(
  draws = draws_small, stan_dat = stan_dat, log_lik_ref = ll_ref,
  meta = list(T = stan_dat$T, K = stan_dat$K, n_draws = length(keep),
              created = as.character(Sys.Date()))
)
out <- file.path(pkg, "tests", "testthat", "fixtures", "golden_statespace.rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(golden, out)
cat(sprintf("Wrote %s  (T=%d K=%d draws=%d)\n", out, stan_dat$T, stan_dat$K, length(keep)))
