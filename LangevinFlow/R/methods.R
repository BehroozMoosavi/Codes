#' @export
print.langevin_chain <- function(x, ...) {
  cat("<langevin_chain>\n")
  cat(sprintf("  Algorithm:       %s\n", x$algorithm))
  cat(sprintf("  Dimension:       %d\n", x$dimension))
  cat(sprintf("  Iterations:      %d (burn-in: %d)\n",
              x$n_iter, x$burn_in))
  cat(sprintf("  Step size:       %g\n", x$step_size))
  cat(sprintf("  Inverse temp.:   %g\n", x$beta))
  if (!is.na(x$acceptance_rate)) {
    cat(sprintf("  Acceptance rate: %.3f\n", x$acceptance_rate))
  }
  if (!is.na(x$elapsed_secs)) {
    cat(sprintf("  Elapsed:         %.3f sec\n", x$elapsed_secs))
  }
  invisible(x)
}

#' Summarize a Langevin Chain
#'
#' Computes per-coordinate posterior summaries (mean, standard deviation,
#' and selected quantiles) from a fitted Langevin chain.
#'
#' @param object An object of class \code{"langevin_chain"}.
#' @param probs Numeric vector of probabilities for quantile computation.
#' @param ... Currently unused.
#' @return Invisibly, a data frame of summaries. Printed by default.
#' @export
summary.langevin_chain <- function(object,
                                   probs = c(0.025, 0.5, 0.975),
                                   ...) {
  S <- object$samples
  d <- ncol(S)

  means <- colMeans(S)
  sds   <- apply(S, 2, stats::sd)
  qs    <- t(apply(S, 2, stats::quantile, probs = probs, na.rm = TRUE))

  df <- data.frame(
    mean = means,
    sd   = sds,
    qs,
    check.names = FALSE,
    row.names   = if (is.null(colnames(S))) {
      paste0("x", seq_len(d))
    } else {
      colnames(S)
    }
  )

  print(object)
  cat("\nPer-coordinate summaries:\n")
  print(signif(df, 4))
  invisible(df)
}

#' Plot a Langevin Chain
#'
#' Trace plot for each dimension, plus a marginal histogram. For
#' high-dimensional chains, only the first \code{max_dim} coordinates are
#' shown.
#'
#' @param x An object of class \code{"langevin_chain"}.
#' @param max_dim Maximum number of coordinates to display. Defaults to 4.
#' @param ... Passed to \code{\link[graphics]{plot}}.
#' @return Called for side effect. Returns \code{x} invisibly.
#' @export
plot.langevin_chain <- function(x, max_dim = 4L, ...) {
  S <- x$samples
  d <- ncol(S)
  k <- min(d, max_dim)

  old_par <- graphics::par(mfrow = c(k, 2), mar = c(3, 3, 2, 1),
                           mgp = c(1.8, 0.6, 0))
  on.exit(graphics::par(old_par), add = TRUE)

  iter <- seq_len(nrow(S))
  for (j in seq_len(k)) {
    graphics::plot(iter, S[, j], type = "l",
                   xlab = "iteration",
                   ylab = sprintf("x[%d]", j),
                   main = sprintf("Trace: dim %d", j), ...)
    graphics::abline(h = mean(S[, j]), lty = 2)
    graphics::hist(S[, j], breaks = 40,
                   xlab = sprintf("x[%d]", j),
                   main = sprintf("Marginal: dim %d", j))
  }

  invisible(x)
}
