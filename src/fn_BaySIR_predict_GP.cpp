// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace Rcpp;
using namespace arma;


vec mvrnorm(vec &mu, mat &Sigma) {
  
  vec X; X.randn(mu.size());
  mat R = chol(Sigma);
  return(trans(R) * X + mu);

}



// [[Rcpp::export]]
vec BaySIR_predict_GP (vec &beta, vec &mu, double sigma_beta, double rho,
  mat &X, mat &X_pred, int T, int T_pred) {
  

  mat Sigma_beta(T + T_pred + 1, T + T_pred + 1);
  Sigma_beta.fill(0.0);
  
  int i, j;

  for (i = 0; i < T + T_pred + 1; i++) {
    for (j = 0; j < T + T_pred + 1; j++) {
      Sigma_beta(i, j) = pow(rho, (double) abs(i - j));
    }
  }
  
  Sigma_beta = pow(sigma_beta, 2) * Sigma_beta;


  mat Sigma_beta_11 = Sigma_beta.submat(0, 0, T, T);
  mat Sigma_beta_21 = Sigma_beta.submat(T + 1, 0, T + T_pred, T);
  mat Sigma_beta_22 = Sigma_beta.submat(T + 1, T + 1, T + T_pred, T + T_pred);
  mat inv_Sigma_beta_11 = inv(Sigma_beta_11);
  
  vec beta_pred_mean = X_pred * mu + 
            Sigma_beta_21 * inv_Sigma_beta_11 * (log(beta) - X * mu);

  mat beta_pred_cov = Sigma_beta_22 - Sigma_beta_21 * inv_Sigma_beta_11 * trans(Sigma_beta_21);
  
  // Rcout << beta_pred_cov << std::endl; 

  vec beta_pred = exp(mvrnorm(beta_pred_mean, beta_pred_cov));

  return(beta_pred);

}


