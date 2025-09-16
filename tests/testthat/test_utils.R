context("Utility functions")

test_that("safe_div prevents division by zero", {
  expect_equal(safe_div(c(1, 2, 3), c(1, 1, 1)), c(1, 2, 3))
  expect_equal(safe_div(c(1, 2, 0), c(1, 1, 0)), c(1, 2, .Machine$double.eps/.Machine$double.eps))
})

test_that("kl_divergence computes correctly", {
  p <- c(0.5, 0.5)
  q <- c(0.5, 0.5)
  expect_equal(kl_divergence(p, q), 0)
  
  p <- c(1, 0)
  q <- c(0.5, 0.5)
  expect_true(kl_divergence(p, q) > 0)
})

test_that("row_norm1 normalizes rows to sum to 1", {
  M <- matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  normalized <- row_norm1(M)
  expect_equal(rowSums(normalized), c(1, 1))
})