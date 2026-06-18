# Closed-form conjugate baseline (Kalman/RTS) — pure R, always runs.

test_that("simulate_disagg is deterministic and well-formed", {
  s1 <- simulate_disagg(T = 20, K = 3, seed = 5)
  s2 <- simulate_disagg(T = 20, K = 3, seed = 5)
  expect_identical(s1$cpi, s2$cpi)
  expect_identical(s1$phi_true, s2$phi_true)
  expect_equal(dim(s1$phi_true), c(20L, 3L))
  expect_true(all(s1$cpi > 0))
  expect_true(all(abs(rowSums(s1$W) - 1) < 1e-10))   # weights on the simplex
  expect_equal(s1$agg_true, rowSums(s1$W * s1$phi_true), tolerance = 1e-12)
})

test_that("conjugate disaggregation has correct structure and tracks the aggregate", {
  sim <- simulate_disagg(T = 25, K = 4, seed = 7)
  bl  <- disaggregate_conjugate(sim$cpi, sim$W)
  expect_s3_class(bl, "disagg_conjugate")
  med <- bl$phi_summary$median
  expect_equal(dim(med), c(25L, 4L))
  # credible bands ordered
  expect_true(all(bl$phi_summary$q2.5 <= med))
  expect_true(all(med <= bl$phi_summary$q97.5))
  # the smoothed aggregate reproduces sum_k W phi and tracks the CPI tightly
  expect_equal(bl$agg_summary[, "median"], rowSums(sim$W * med), tolerance = 1e-8)
  expect_gt(cor(bl$agg_summary[, "median"], sim$cpi), 0.95)
  expect_true(is.finite(bl$loglik))
})

test_that("K = 1 with unit weight: smoother essentially observes the level", {
  set.seed(1)
  Tn <- 30
  cpi <- 100 + cumsum(rnorm(Tn, 0, 1))
  W <- matrix(1, Tn, 1)
  bl <- disaggregate_conjugate(cpi, W, r_frac = 0.02)
  # with a direct (unit-weight) observation and small obs noise the smoothed
  # state must hug the data
  expect_gt(cor(bl$phi_summary$median[, 1], cpi), 0.98)
})

test_that("simulation smoother returns reproducible joint draws centred on the smooth", {
  sim <- simulate_disagg(T = 20, K = 3, seed = 9)
  b1 <- disaggregate_conjugate(sim$cpi, sim$W, n_draws = 200, seed = 123)
  b2 <- disaggregate_conjugate(sim$cpi, sim$W, n_draws = 200, seed = 123)
  expect_equal(dim(b1$phi_draws), c(20L, 3L, 200L))
  expect_identical(b1$phi_draws, b2$phi_draws)               # seeded -> reproducible
  # Durbin-Koopman draws are unbiased for the smoothed mean
  draw_mean <- apply(b1$phi_draws, c(1, 2), mean)
  expect_equal(draw_mean, b1$phi_summary$median, tolerance = 0.5)
})

test_that("conjugate draws satisfy the bayesianOU phi_draws array contract", {
  sim <- simulate_disagg(T = 18, K = 2, seed = 3)
  bl <- disaggregate_conjugate(sim$cpi, sim$W, n_draws = 10, seed = 1)
  expect_true(is.array(bl$phi_draws))
  expect_length(dim(bl$phi_draws), 3L)
  expect_equal(dim(bl$phi_draws)[1:2], c(18L, 2L))
})
