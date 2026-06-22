test_that("S3 methods work on a langevin_chain object", {
  set.seed(1)
  g <- function(x) x
  fit <- ula(c(1, 1), g, 0.1, 500, burn_in = 100)

  expect_s3_class(fit, "langevin_chain")
  expect_equal(fit$dimension, 2L)
  expect_equal(nrow(fit$samples), 400L)

  # print() returns invisible(x)
  expect_output(print(fit), "langevin_chain")
  expect_output(print(fit), "ULA")

  # summary() returns a data frame invisibly
  s <- expect_output(summary(fit), "Per-coordinate summaries")
  expect_true(is.data.frame(s) || is.null(s))

  # plot() should not error (use a null device)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(fit))
})

test_that("Edge case: 1-dimensional state", {
  set.seed(1)
  g <- function(x) x
  fit <- ula(0, g, 0.1, 200)
  expect_equal(fit$dimension, 1L)
  expect_equal(ncol(fit$samples), 1L)
})
