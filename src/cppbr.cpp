#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cppbr(int iter, arma::vec y, arma::mat X, double prior_a, 
                    double prior_b, arma::vec prior_mu, arma::vec prior_prec) {
  
  // First get the number of cases n
  int n = y.size();
  // And the number of variables k
  int k = X.n_cols;

  // Then calculate the OLS estimates for beta
  arma::mat v = arma::inv(X.t()*X);
  arma::vec bhat = v*X.t()*y;
  
  // Use this estimate for the posterior of sigma-squared
  
  arma::vec posterior_b = prior_b + ((y-X*bhat).t()*(y-X*bhat))/2;
  double posterior_a = prior_a + n/2;
  
  // For the posteriors of the betas we need the prior precision matrix first
  arma::mat prior_tau = arma::diagmat(prior_prec); // diagonal precision matrix

  // Calculate the posterior mean vector
  arma::vec posterior_mu = arma::inv(X.t()*X + prior_tau)*(prior_tau*prior_mu + X.t()*y);
  
  
  // Before we start sampling, create the output matrix!
  arma::mat out((iter),(k+2)); // columns for variance, betas and r-squared
  
  
  // Sample from marginal distributions
  for (int i=0; i<iter; i++) {
    
      
    // Posterior for var
    double varn = 1/(R::rgamma(posterior_a, 1/posterior_b[0]));
    
    
    // Posteriors for coefficients
    // First calculate posterior sigma (matrix)
    arma::mat posterior_sigma = varn * arma::inv(X.t()*X + prior_tau);
      
    // Sample from multivariate normal
    arma::mat Y = arma::randn(1, k);
    arma::vec betan = (posterior_mu.t() + Y * arma::chol(posterior_sigma)).t();
    
    // Calculate R-squared
    arma::vec idvec = arma::vec(n); // initialise identity vector
    std::fill(idvec.begin(), idvec.end(), 1); // fill with 1s
   
    arma::vec res = y-X*betan; // Save Residuals
    arma::vec resmean = idvec.t()*res*arma::inv(idvec.t()*idvec); // calculate residual mean
    arma::vec resvar = (res-resmean[0]).t()*(res-resmean[0])/(n-1); // calculate var
    
    arma::vec ymean = idvec.t()*y*arma::inv(idvec.t()*idvec); //calculate outcome mean
    arma::vec yvar = (y-ymean[0]).t()*(y-ymean[0])/(n-1); // calculate var
   
    double rsq = 1 - resvar[0]/yvar[0]; // R-squared
    
    
    // Store variables in output matrix
    out(i,0) = varn;
    for (int o = 1; o<k+1; o++){
      out(i,o) = betan[(o-1)];
    }
    out(i,k+1) = rsq;
  }
  
  
  // Return to R
  return out;
}