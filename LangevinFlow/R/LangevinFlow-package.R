#' LangevinFlow: Langevin Diffusion Samplers with a C++ Backend
#'
#' Provides Markov chain Monte Carlo samplers based on overdamped Langevin
#' diffusion. The unadjusted variant (\code{\link{ula}}) implements the
#' Euler-Maruyama discretization directly; the Metropolis-adjusted variant
#' (\code{\link{mala}}) corrects the discretization bias with an
#' accept/reject step.
#'
#' @section Sign convention:
#' All samplers expect the user to supply derivatives of the
#' \emph{potential} \eqn{U(x) = -\log \pi(x)}, not the log-density itself.
#' That is, if your target density is \eqn{\pi(x) \propto \exp(-\beta U(x))},
#' you pass \code{U} and \eqn{\nabla U}. This matches the physics
#' convention in the Langevin SDE
#' \deqn{dX_t = -\nabla U(X_t)\,dt + \sqrt{2/\beta}\,dW_t.}
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{ula}}: Unadjusted Langevin Algorithm.
#'   \item \code{\link{mala}}: Metropolis-Adjusted Langevin Algorithm.
#' }
#'
#' @references
#' Roberts, G. O., & Tweedie, R. L. (1996). Exponential convergence of
#' Langevin distributions and their discrete approximations.
#' \emph{Bernoulli}, 2(4), 341-363.
#'
#' Roberts, G. O., & Rosenthal, J. S. (1998). Optimal scaling of
#' discrete approximations to Langevin diffusions. \emph{Journal of the
#' Royal Statistical Society, Series B}, 60(1), 255-268.
#'
#' @keywords internal
#' @aliases LangevinFlow-package
"_PACKAGE"

## usethis namespace: start
#' @useDynLib LangevinFlow, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
