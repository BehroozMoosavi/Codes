#' Unadjusted Langevin Algorithm
#'
#' Draws a Markov chain whose stationary distribution approximates
#' \eqn{\pi(x) \propto \exp(-\beta U(x))} using the Euler-Maruyama
#' discretization of the overdamped Langevin diffusion. The discretization
#' introduces a bias that vanishes as \code{step_size} tends to zero; for
#' an asymptotically exact sampler, use \code{\link{mala}}.
#'
#' @param init_x Numeric vector. Starting state of the chain. Its length
#'   defines the dimension.
#' @param grad_u Function. Takes a numeric vector of length
#'   \code{length(init_x)} and returns the gradient of
#'   \eqn{U(x) = -\log \pi(x)}, as a numeric vector of the same length.
#'   See the sign convention in \code{?LangevinFlow}.
#' @param step_size Positive numeric scalar. Discretization step
#'   \eqn{\gamma}.
#' @param n_iter Positive integer. Number of iterations to run.
#' @param beta Positive numeric scalar. Inverse temperature; defaults to 1.
#' @param burn_in Non-negative integer. Number of initial samples to
#'   discard. Must be strictly less than \code{n_iter}. Defaults to 0.
#'
#' @return An object of class \code{"langevin_chain"}, which is a list with
#'   components:
#' \describe{
#'   \item{samples}{Numeric matrix of post-burn-in samples; rows index
#'     iterations and columns index dimensions.}
#'   \item{algorithm}{Character string \code{"ULA"}.}
#'   \item{step_size, beta, n_iter, burn_in, dimension}{Echoed inputs.}
#'   \item{acceptance_rate}{\code{NA} for ULA (no accept/reject step).}
#'   \item{elapsed_secs}{Wall-clock runtime of the sampler in seconds.}
#' }
#'
#' @references
#' Roberts, G. O., & Tweedie, R. L. (1996). Exponential convergence of
#' Langevin distributions and their discrete approximations.
#' \emph{Bernoulli}, 2(4), 341-363.
#'
#' @seealso \code{\link{mala}}
#'
#' @examples
#' # Standard 2D Gaussian target: U(x) = 0.5 * ||x||^2, grad_U(x) = x.
#' set.seed(1)
#' grad_u <- function(x) x
#' fit <- ula(init_x = c(3, -3), grad_u = grad_u,
#'            step_size = 0.05, n_iter = 2000, burn_in = 500)
#' summary(fit)
#'
#' @export
ula <- function(init_x, grad_u, step_size, n_iter, beta = 1, burn_in = 0L) {
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

  grad_u <- .check_fun(grad_u, "grad_u", init_x, expect_scalar = FALSE)

  t0 <- proc.time()[["elapsed"]]
  samples <- ula_sampler_cpp(
    init_x    = init_x,
    grad_u    = grad_u,
    step_size = step_size,
    n_iter    = n_iter,
    beta      = beta
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  .new_langevin_chain(
    samples         = samples,
    algorithm       = "ULA",
    step_size       = step_size,
    beta            = beta,
    n_iter          = n_iter,
    burn_in         = burn_in,
    acceptance_rate = NA_real_,
    accepted        = NULL,
    elapsed         = elapsed
  )
}
