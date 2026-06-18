# Input validation and small utilities — pure R, always runs (no Stan needed:
# the validation happens before any sampling).

test_that("disaggregate_statespace rejects malformed inputs before sampling", {
  W <- matrix(1/3, 5, 3)
  expect_error(disaggregate_statespace(c(1), W), "length >= 2")
  expect_error(disaggregate_statespace(c(-1, 2, 3, 4, 5), W), "strictly positive")
  expect_error(disaggregate_statespace(rep(100, 4), W), "must match")          # T mismatch
  expect_error(disaggregate_statespace(rep(100, 5), W, industries = c("a", "b")),
               "length K")
  expect_error(disaggregate_statespace(rep(100, 5), W, years = 1:4), "length T")
  expect_error(disaggregate_statespace(rep(100, 5), W, priors = list(nope = 1)),
               "Unknown prior")
})

test_that("disagg_default_priors validates the aggregate and scales with it", {
  expect_error(disagg_default_priors(c(1)), "length >= 2")
  expect_error(disagg_default_priors(c(0, 1, 2)), "strictly positive")
  pr <- disagg_default_priors(c(100, 110, 120, 130))
  expect_equal(pr$phi1_center, 100)
  expect_true(pr$s_sigma_cpi > 0)
  expect_named(pr, c("phi1_center", "s_omega_struct", "s_delta_mu", "s_delta_sigma",
                     "log_tau_loc", "s_log_tau_mu", "s_log_tau_sigma", "s_sigma_cpi"))
})

test_that("disaggregate_conjugate rejects malformed inputs", {
  W <- matrix(1/3, 5, 3)
  expect_error(disaggregate_conjugate(c(1), W), "length >= 2")
  expect_error(disaggregate_conjugate(c(-1, 1, 1, 1, 1), W), "strictly positive")
  expect_error(disaggregate_conjugate(rep(100, 4), W), "match length")
})

test_that("row_norm1 normalizes rows and is robust to degenerate rows", {
  rn <- BayesianDisaggregation:::row_norm1
  M <- matrix(c(1, 1, 2, 2), 2, byrow = TRUE)
  expect_equal(rowSums(rn(M)), c(1, 1))
  M0 <- matrix(c(0, 0, 1, 1), 2, byrow = TRUE)          # first row all-zero
  out <- rn(M0)
  expect_true(all(is.finite(out)))
})

test_that("to_num_commas parses localized numerics", {
  f <- BayesianDisaggregation:::to_num_commas
  expect_equal(f("1.234,5"), 1234.5)
  expect_equal(f(c(10, 20)), c(10, 20))
})
