#' Metropolis-Adjusted Langevin Algorithm
#'
#' Asymptotically exact sampler for \eqn{\pi(x) \propto \exp(-\beta U(x))}.
#' At each step a Langevin proposal is generated as in \code{\link{ula}}
#' and then accepted or rejected by the Metropolis-Hastings rule so that
#' \eqn{\pi} is exactly invariant. Optimal acceptance rates are typically
#' near 0.574 in high dimension (Roberts & Rosenthal, 1998).
#'
#' @param init_x Numeric vector. Starting state of the chain.
#' @param U Function. Takes a numeric vector of length
#'   \code{length(init_x)} and returns the scalar potential
#'   \eqn{U(x) = -\log \pi(x)} (additive constants do not matter).
#' @param grad_u Function. Returns the gradient of \code{U}; same shape
#'   convention as for \code{\link{ula}}.
#' @param step_size Positive numeric scalar. Discretization step
#'   \eqn{\gamma}.
#' @param n_iter Positive integer. Number of iterations.
#' @param beta Positive numeric scalar. Inverse temperature; defaults to 1.
#' @param burn_in Non-negative integer, strictly less than \code{n_iter}.
#'   Defaults to 0.
#'
#' @return An object of class \code{"langevin_chain"} (see \code{\link{ula}}
#'   for structure). The \code{acceptance_rate} component is populated and
#'   \code{accepted} is a logical vector indicating per-iteration outcomes
#'   (post-burn-in).
#'
#' @references
#' Roberts, G. O., & Rosenthal, J. S. (1998). Optimal scaling of discrete
#' approximations to Langevin diffusions. \emph{Journal of the Royal
#' Statistical Society, Series B}, 60(1), 255-268.
#'
#' @seealso \code{\link{ula}}
#'
#' @examples
#' # Standard 2D Gaussian: U(x) = 0.5 * ||x||^2.
#' set.seed(1)
#' U      <- function(x) 0.5 * sum(x^2)
#' grad_u <- function(x) x
#' fit <- mala(init_x = c(3, -3), U = U, grad_u = grad_u,
#'             step_size = 0.3, n_iter = 2000, burn_in = 500)
#' summary(fit)
#'
#' @export
mala <- function(init_x, U, grad_u, step_size, n_iter, beta = 1,
                 burn_in = 0L) {
  init_x    <- .check_init(init_x)
  step_size <- .check_pos_scalar(step_size, "step_size")
  beta      <- .check_pos_scalar(beta, "beta")
  n_iter    <- .check_pos_int(n_iter, "n_iter")

  if (!is.numeric(burn_in) || length(burn_in) != 1L || burn_in < 0 ||
      burn_in != round(burn_in)) {
    stop("'burn_in' must be a non-negative integer.", call. = FALSE)
  }
  burn_in <- as.integer(burn_in)
  if (burn_in >= n_iter) {
    stop("'burn_in' must be strictly less than 'n_iter'.", call. = FALSE)
  }

  U      <- .check_fun(U,      "U",      init_x, expect_scalar = TRUE)
  grad_u <- .check_fun(grad_u, "grad_u", init_x, expect_scalar = FALSE)

  t0 <- proc.time()[["elapsed"]]
  res <- mala_sampler_cpp(
    init_x    = init_x,
    U         = U,
    grad_u    = grad_u,
    step_size = step_size,
    n_iter    = n_iter,
    beta      = beta
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  .new_langevin_chain(
    samples         = res$samples,
    algorithm       = "MALA",
    step_size       = step_size,
    beta            = beta,
    n_iter          = n_iter,
    burn_in         = burn_in,
    acceptance_rate = res$acceptance_rate,
    accepted        = as.logical(res$accepted),
    elapsed         = elapsed
  )
}
