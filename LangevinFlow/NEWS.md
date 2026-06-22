# LangevinFlow 0.1.0

* First CRAN submission.
* Implements the Unadjusted Langevin Algorithm (`ula()`) and the
  Metropolis-Adjusted Langevin Algorithm (`mala()`) with a C++ backend
  using `Rcpp` and `RcppArmadillo`.
* Provides an S3 class `"langevin_chain"` with `print()`, `summary()`,
  and `plot()` methods.
* All random number generation is routed through R's RNG state, so
  `set.seed()` controls reproducibility of the chains.
