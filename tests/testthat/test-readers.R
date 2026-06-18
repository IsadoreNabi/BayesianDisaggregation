# Data readers and alignment against the bundled extdata — always runs.

cpi_file <- function() system.file("extdata", "CPI.xlsx", package = "BayesianDisaggregation")
w_file   <- function() system.file("extdata", "WEIGHTS.xlsx", package = "BayesianDisaggregation")

test_that("read_cpi returns a clean, sorted Year/CPI frame", {
  skip_if(cpi_file() == "")
  df <- read_cpi(cpi_file())
  expect_true(all(c("Year", "CPI") %in% names(df)))
  expect_false(is.unsorted(df$Year))
  expect_true(all(is.finite(df$CPI)))
  expect_true(all(c(1982, 1983, 1984) %in% df$Year))    # base period present
})

test_that("read_weights_matrix returns a simplex-per-year prior matrix", {
  skip_if(w_file() == "")
  w <- read_weights_matrix(w_file())
  expect_true(is.matrix(w$P))
  expect_true(all(abs(rowSums(w$P) - 1) < 1e-8))
  expect_equal(length(w$industries), ncol(w$P))
  expect_equal(length(w$years), nrow(w$P))
})

test_that("align_disagg_inputs intersects on common years", {
  skip_if(cpi_file() == "" || w_file() == "")
  al <- align_disagg_inputs(cpi_file(), w_file())
  expect_equal(length(al$cpi), nrow(al$W))
  expect_true(all(abs(rowSums(al$W) - 1) < 1e-8))
  expect_true(length(al$years) >= 2)
  expect_equal(al$industries, colnames(al$W))
})
