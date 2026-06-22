// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat ula_sampler_cpp(arma::vec init_x,
                          Function grad_u,
                          double step_size,
                          int n_iter,
                          double beta) {
   int d = init_x.n_elem;
   arma::mat samples(n_iter, d);
   arma::vec current_x = init_x;

   // Diffusion scale: sqrt(2 * gamma / beta)
   const double noise_scale = std::sqrt(2.0 * step_size / beta);

   for (int i = 0; i < n_iter; ++i) {
     // CRAN: allow user to break long C++ loops
     if (i % 1000 == 0) Rcpp::checkUserInterrupt();

     // Evaluate gradient of U at current state
     NumericVector grad_res = grad_u(current_x);
     arma::vec grad = as<arma::vec>(grad_res);

     if (grad.n_elem != (unsigned int)d) {
       Rcpp::stop("Gradient function returned a vector of wrong length.");
     }

     // R-linked RNG (respects set.seed() on the R side)
     arma::vec noise = as<arma::vec>(rnorm(d));

     // Euler-Maruyama step
     current_x = current_x - step_size * grad + noise_scale * noise;

     // Divergence guard

     if (!current_x.is_finite()) {
       Rcpp::stop("Sampler diverged at iteration %d. "
                    "Try a smaller step_size.", i + 1);
     }

     samples.row(i) = current_x.t();
   }

   return samples;
 }
