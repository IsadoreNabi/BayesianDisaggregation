# Recovery of the state-space engine on data simulated from its own DGP.
# Heavy (compiles + samples): gated behind NOT_CRAN and an explicit flag.
#   NOT_CRAN=true BAYESDISAGG_RUN_RECOVERY=1 Rscript -e 'devtools::test(...)'
#
# Rigor by layers: the AGGREGATE is strongly identified (must be recovered
# tightly); the SECTORAL split is weakly identified by construction (the bands
# are honestly wide, so coverage is high and the point correlation is only
# moderate) — we assert the identified quantity tightly and the unidentified one
# conservatively, never claiming sectoral precision that the data cannot give.

test_that("state-space recovers the aggregate and yields the phi_draws contract", {
  skip_on_cran()
  if (Sys.getenv("BAYESDISAGG_RUN_RECOVERY") != "1") skip("heavy recovery gated")
  if (check_stan_backend() == "none") skip("no Stan backend available")

  sim <- simulate_disagg(T = 35, K = 4, seed = 101)
  fit <- disaggregate_statespace(sim$cpi, sim$W, chains = 2, iter = 900, warmup = 450,
                                 seed = 101, adapt_delta = 0.97, max_treedepth = 12)

  # contract consumed by bayesianOU::fit_ou_nested_mi
  expect_equal(dim(fit$phi_draws)[1:2], c(35L, 4L))
  expect_true(dim(fit$phi_draws)[3] >= 100L)
  expect_true(all(is.finite(fit$phi_draws)))

  # convergence
  expect_lt(fit$diagnostics$rhat_max, 1.15)

  # aggregate: strongly identified
  expect_gt(cor(fit$agg_summary[, "median"], sim$cpi), 0.95)
  expect_gt(cor(fit$agg_summary[, "median"], sim$agg_true), 0.95)

  # sectoral: weakly identified -> wide bands -> conservative coverage
  cov95 <- mean(sim$phi_true >= fit$phi_summary$q2.5 & sim$phi_true <= fit$phi_summary$q97.5)
  expect_gt(cov95, 0.70)
})
