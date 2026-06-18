# Golden regression for the state-space computation path.
# Recomputes log_lik via generate_quantities on FROZEN parameter draws + data and
# demands a bit-for-bit match against the frozen reference. Catches any change to
# the Stan model's computed quantities. gq is deterministic (no RNG in the
# generated quantities block), so an exact match is the correct expectation.
# Gated: needs cmdstanr (gq) and the fixture; auto-skips on CRAN.

test_that("state-space generated quantities reproduce the golden log_lik bit-for-bit", {
  skip_on_cran()
  fx <- testthat::test_path("fixtures", "golden_statespace.rds")
  skip_if(!file.exists(fx), "golden fixture absent (run validacion/make_golden_statespace.R)")
  skip_if_not_installed("cmdstanr")
  if (check_stan_backend() != "cmdstanr") skip("golden recompute requires cmdstanr")

  golden <- readRDS(fx)
  mod <- cmdstanr::cmdstan_model(BayesianDisaggregation:::.disagg_stan_file_path(),
                                 dir = tempdir(), pedantic = FALSE)
  gq <- mod$generate_quantities(fitted_params = golden$draws, data = golden$stan_dat)
  ll_new <- gq$draws(variables = "log_lik", format = "draws_matrix")

  expect_equal(dim(ll_new), dim(golden$log_lik_ref))
  expect_identical(as.numeric(ll_new), as.numeric(golden$log_lik_ref))
})
