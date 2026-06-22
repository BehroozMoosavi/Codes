// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// Log-density of multivariate normal N(mean, sigma^2 I), up to constants
// that cancel in the acceptance ratio.
static inline double log_q_proposal(const arma::vec& y,
                                    const arma::vec& mean,
                                    double sigma_sq) {
  arma::vec diff = y - mean;
  return -0.5 * arma::dot(diff, diff) / sigma_sq;
}

// [[Rcpp::export]]
List mala_sampler_cpp(arma::vec init_x,
                      Function U,
                      Function grad_u,
                       double step_size,
                       int n_iter,
                       double beta) {
   int d = init_x.n_elem;
   arma::mat samples(n_iter, d);
   LogicalVector accepted(n_iter);

   arma::vec current_x = init_x;
   double current_U = as<double>(U(current_x));
   arma::vec current_grad = as<arma::vec>(grad_u(current_x));

   if (current_grad.n_elem != (unsigned int)d) {
     Rcpp::stop("Gradient function returned a vector of wrong length.");
   }

   const double sigma_sq = 2.0 * step_size / beta;
   const double sigma = std::sqrt(sigma_sq);

   int n_accept = 0;

   for (int i = 0; i < n_iter; ++i) {
     if (i % 1000 == 0) Rcpp::checkUserInterrupt();

     // Forward proposal mean: x - gamma * grad_U(x)
     arma::vec fwd_mean = current_x - step_size * current_grad;
     arma::vec noise = as<arma::vec>(rnorm(d));
     arma::vec proposal = fwd_mean + sigma * noise;

     if (!proposal.is_finite()) {
       // Reject silently, record current state
       samples.row(i) = current_x.t();
       accepted[i] = false;
       continue;
     }

     double prop_U = as<double>(U(proposal));
     arma::vec prop_grad = as<arma::vec>(grad_u(proposal));

     // Reverse proposal mean: y - gamma * grad_U(y)
     arma::vec rev_mean = proposal - step_size * prop_grad;

     // log alpha = -beta * (U(y) - U(x)) + log q(x|y) - log q(y|x)
     double log_alpha = -beta * (prop_U - current_U)
       + log_q_proposal(current_x, rev_mean, sigma_sq)
       - log_q_proposal(proposal, fwd_mean, sigma_sq);

       // Accept/reject using R's RNG
       double u = as<double>(runif(1));
       if (std::log(u) < log_alpha) {
         current_x = proposal;
         current_U = prop_U;
         current_grad = prop_grad;
         accepted[i] = true;
         ++n_accept;
       } else {
         accepted[i] = false;
       }

       samples.row(i) = current_x.t();
   }

   double accept_rate = (double)n_accept / (double)n_iter;

   return List::create(
     Named("samples") = samples,
     Named("accepted") = accepted,
     Named("acceptance_rate") = accept_rate
   );
 }
