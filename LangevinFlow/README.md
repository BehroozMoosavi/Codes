# LangevinFlow

Langevin diffusion Markov chain Monte Carlo samplers for R, with a C++
backend via `Rcpp` and `RcppArmadillo`.

## Installation

Development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("BehroozMoosavi/LangevinFlow")
```

## Usage

The package exposes two samplers, both targeting
$\pi(x) \propto \exp(-\beta U(x))$ where the user supplies the potential
$U(x) = -\log \pi(x)$ (and, for MALA, its gradient).

### Unadjusted Langevin Algorithm

```r
library(LangevinFlow)

# Standard 2D Gaussian: U(x) = 0.5 * ||x||^2
grad_u <- function(x) x

set.seed(1)
fit <- ula(init_x = c(3, -3), grad_u = grad_u,
           step_size = 0.05, n_iter = 5000, burn_in = 1000)
summary(fit)
plot(fit)
```

### Metropolis-Adjusted Langevin Algorithm

```r
U      <- function(x) 0.5 * sum(x^2)
grad_u <- function(x) x

set.seed(1)
fit <- mala(init_x = c(3, -3), U = U, grad_u = grad_u,
            step_size = 0.4, n_iter = 5000, burn_in = 1000)
fit$acceptance_rate
```

## When to use which

| Sampler | Asymptotically exact? | Tuning sensitivity | Best for |
|---|---|---|---|
| `ula()` | No — biased by `step_size` | Less | Quick exploration; stochastic optimization |
| `mala()` | Yes | More (target ~0.574 acceptance) | Accurate inference |

## References

- Roberts, G. O., & Tweedie, R. L. (1996). Exponential convergence of
  Langevin distributions and their discrete approximations. *Bernoulli*,
  2(4), 341–363.
- Roberts, G. O., & Rosenthal, J. S. (1998). Optimal scaling of discrete
  approximations to Langevin diffusions. *JRSS-B*, 60(1), 255–268.
