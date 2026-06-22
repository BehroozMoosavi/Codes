test_that("ULA on standard Gaussian recovers mean and variance", {
  # Target: N(0, I_2). U(x) = 0.5 * ||x||^2, grad_U(x) = x.
  set.seed(42)
  g <- function(x) x
  fit <- ula(init_x = c(0, 0), grad_u = g, step_size = 0.05,
             n_iter = 6000, burn_in = 1000)

  mu  <- colMeans(fit$samples)
  v   <- apply(fit$samples, 2, var)

  # Tolerances are loose to keep the test fast and stable.
  expect_lt(max(abs(mu)), 0.20)
  expect_lt(max(abs(v - 1)), 0.30)
  expect_equal(fit$algorithm, "ULA")
  expect_true(is.na(fit$acceptance_rate))
})

test_that("MALA on standard Gaussian recovers mean, variance, and has plausible acceptance rate", {
  set.seed(42)
  U <- function(x) 0.5 * sum(x^2)
  g <- function(x) x
  fit <- mala(init_x = c(0, 0), U = U, grad_u = g,
              step_size = 0.4, n_iter = 4000, burn_in = 1000)

  mu <- colMeans(fit$samples)
  v  <- apply(fit$samples, 2, var)

  expect_lt(max(abs(mu)), 0.20)
  expect_lt(max(abs(v - 1)), 0.30)
  expect_equal(fit$algorithm, "MALA")
  expect_true(fit$acceptance_rate > 0.3 && fit$acceptance_rate < 0.95)
})

test_that("MALA on non-isotropic Gaussian recovers covariance", {
  # Target: N(mu, Sigma) with Sigma diagonal = (1, 4).
  # U(x) = 0.5 * (x - mu)^T Sigma^{-1} (x - mu).
  set.seed(7)
  mu_true <- c(1, -1)
  prec    <- diag(c(1, 1 / 4))
  U <- function(x) {
    d <- x - mu_true
    0.5 * sum(d * (prec %*% d))
  }
  g <- function(x) as.numeric(prec %*% (x - mu_true))

  fit <- mala(init_x = c(0, 0), U = U, grad_u = g,
              step_size = 0.5, n_iter = 5000, burn_in = 1000)

  mu_hat <- colMeans(fit$samples)
  v_hat  <- apply(fit$samples, 2, var)

  expect_lt(abs(mu_hat[1] - mu_true[1]), 0.15)
  expect_lt(abs(mu_hat[2] - mu_true[2]), 0.30)
  expect_lt(abs(v_hat[1]  - 1), 0.35)
  expect_lt(abs(v_hat[2]  - 4), 1.2)
})
