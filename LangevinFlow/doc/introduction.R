## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(LangevinFlow)

# Target: N(0, I_2). U(x) = 0.5 * ||x||^2.
U      <- function(x) 0.5 * sum(x^2)
grad_u <- function(x) x

set.seed(1)
fit_ula <- ula(init_x = c(3, -3), grad_u = grad_u,
               step_size = 0.05, n_iter = 5000, burn_in = 1000)
summary(fit_ula)

## ----fig.alt="ULA trace plots"------------------------------------------------
plot(fit_ula)

## -----------------------------------------------------------------------------
set.seed(1)
fit_mala <- mala(init_x = c(3, -3), U = U, grad_u = grad_u,
                 step_size = 0.4, n_iter = 5000, burn_in = 1000)
fit_mala$acceptance_rate
summary(fit_mala)

## -----------------------------------------------------------------------------
Sigma     <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
Sigma_inv <- solve(Sigma)
mu_true   <- c(2, -1)

U      <- function(x) 0.5 * as.numeric(t(x - mu_true) %*% Sigma_inv %*% (x - mu_true))
grad_u <- function(x) as.numeric(Sigma_inv %*% (x - mu_true))

set.seed(2)
fit <- mala(c(0, 0), U, grad_u, step_size = 0.4,
            n_iter = 8000, burn_in = 2000)
colMeans(fit$samples)
cov(fit$samples)

