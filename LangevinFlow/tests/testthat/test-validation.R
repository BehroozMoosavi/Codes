test_that("ula() rejects bad inputs", {
  g <- function(x) x

  expect_error(ula(numeric(0), g, 0.1, 100), "non-empty")
  expect_error(ula(c(NA, 1), g, 0.1, 100), "finite")
  expect_error(ula(c(1, 1), g, -0.1, 100), "positive")
  expect_error(ula(c(1, 1), g, 0.1, -10), "positive integer")
  expect_error(ula(c(1, 1), g, 0.1, 1.5), "positive integer")
  expect_error(ula(c(1, 1), g, 0.1, 100, beta = 0), "positive")
  expect_error(ula(c(1, 1), g, 0.1, 100, burn_in = -1), "non-negative")
  expect_error(ula(c(1, 1), g, 0.1, 100, burn_in = 100),
               "strictly less than")
  expect_error(ula(c(1, 1), "not a function", 0.1, 100), "must be a function")

  bad_g <- function(x) c(x, 0)  # wrong length
  expect_error(ula(c(1, 1), bad_g, 0.1, 100), "length 2")
})

test_that("mala() rejects bad inputs", {
  U <- function(x) 0.5 * sum(x^2)
  g <- function(x) x

  expect_error(mala(c(1, 1), U = "no", grad_u = g, 0.1, 100),
               "must be a function")
  bad_U <- function(x) c(1, 2)  # not scalar
  expect_error(mala(c(1, 1), U = bad_U, grad_u = g, 0.1, 100),
               "single finite numeric")
})
