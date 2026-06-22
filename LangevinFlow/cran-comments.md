## Submission notes

This is a first submission of the LangevinFlow package.

## Test environments

* local: macOS Sequoia 15.5 (aarch64-apple-darwin20), R 4.4.3
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Downstream dependencies

There are currently no downstream dependencies.

## Notes for reviewers

* All random number generation uses R's RNG (via `Rcpp::rnorm` /
  `Rcpp::runif`), so `set.seed()` controls reproducibility from the R
  side, as exercised by `tests/testthat/test-reproducibility.R`.
* All examples run in well under 5 seconds on a modest laptop.
* The package writes nothing to disk and changes no global options.
* `Rcpp::checkUserInterrupt()` is called inside the main C++ sampling
  loops to ensure long runs are interruptible.
