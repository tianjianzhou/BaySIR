// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "Rcpp.h"

using namespace Rcpp;
using namespace arma;




double logit(double x) {

  return(log(x / (1.0 - x)));

}


double probit(double x) {

  return(R::qnorm(x, 0.0, 1.0, 1, 0));

}




double cloglog(double x) {

  return(log(-log(1.0 - x)));
   
}


double gamma_link(double gamma, int link) {

  if (link == 2) {
    return(probit(gamma));
  } 
  else if (link == 3) {
    return(cloglog(gamma));
  }
  else {
    return(logit(gamma));
  }

}



vec mvrnorm(vec &mu, mat &Sigma) {
  
  vec X; X.randn(mu.size());
  mat R = chol(Sigma);
  return(trans(R) * X + mu);

}






void update_I_U_0 (vec &S, vec &I_U, vec &I_D, vec &R_U, vec &R_D, 
  vec &B, double N,
  vec &beta, 
  vec &gamma, mat &Y, vec &eta, double sigma_gamma,
  double kappa, double alpha, 
  int T, double nu_1, double nu_2, int link,
  double Tmp){
  
  double h, A_t;
  int t, stop_trigger = 0;

  vec S_pro(T + 1); 
  vec I_U_pro(T + 1);
  vec R_U_pro(T + 1);
  
  vec gamma_pro(T + 1); 
  vec gamma_tilde_pro(T + 1);

  vec Y_eta = Y * eta;
  
  h = R::rnorm(0.0, 0.05);
  
  I_U_pro(0) = I_U(0) * exp(h);
  S_pro(0) = S(0) - (I_U_pro(0) - I_U(0));
  R_U_pro(0) = 0;

  gamma_pro(0) = B(0) / ( (1.0 - alpha) * I_U_pro(0) );
  gamma_tilde_pro(0) = gamma_link(gamma_pro(0), link);

  for (t = 0; t < T; t++) {
    
    if ( (1.0 - alpha) * I_U_pro(t) <= B(t) ) {
      stop_trigger = 1;
      break;
    }

    A_t = beta(t) * S_pro(t) * (I_U_pro(t) + kappa * I_D(t)) / N;

    S_pro(t + 1) = S_pro(t) - A_t;
    I_U_pro(t + 1) = (1.0 - alpha) * I_U_pro(t) + A_t - B(t);
    R_U_pro(t + 1) = R_U_pro(t) + alpha * I_U_pro(t);
    
    gamma_pro(t + 1) = B(t + 1) / ( (1.0 - alpha) * I_U_pro(t + 1) );
    gamma_tilde_pro(t + 1) = gamma_link(gamma_pro(t + 1), link);
    
  }

  if ( (1.0 - alpha) * I_U_pro(T) <= B(T) ) {
    stop_trigger = 1;
  }


  if (stop_trigger == 0) {
    
    double loglik_diff = 0, logprior_diff = 0, log_Jacobian = 0, u;

    for (t = 0; t <= T; t++) {

      loglik_diff = loglik_diff - (0.5 / pow(sigma_gamma, 2)) * 
        ( pow(gamma_tilde_pro(t) - Y_eta(t), 2) - 
          pow(gamma_link(gamma(t), link) - Y_eta(t), 2) );
    
    }

    loglik_diff = loglik_diff / Tmp;
    
    logprior_diff = logprior_diff + 
      (nu_1 - 1.0) * h - 
      (nu_2 / I_D(0)) * (I_U_pro(0) - I_U(0));
    
    log_Jacobian = h;

    u = R::runif(0.0, 1.0);

    // Accept proposal
    if(log(u) < loglik_diff + logprior_diff + log_Jacobian){

      for (t = 0; t <= T; t++) {
        
        S(t) = S_pro(t);
        I_U(t) = I_U_pro(t);
        R_U(t) = R_U_pro(t);

        gamma(t) = gamma_pro(t);

      }

    }

  }

  return;

}







void update_beta_t (vec &S, vec &I_U, vec &I_D, vec &R_U, vec &R_D, 
  vec &B, double N,
  vec &beta, mat &X, vec &mu, double sigma_beta, double rho,
  vec &gamma, mat &Y, vec &eta, double sigma_gamma,
  double kappa, double alpha,
  int t, int T, int link,
  double Tmp) {
  

  // Propose beta_t_pro from its prior
  vec X_mu = X * mu;
  double beta_t_pro;

  if (t == 0) {
    beta_t_pro = exp( R::rnorm( X_mu(0) + rho * (log(beta(1)) - X_mu(1)), 
      sigma_beta * sqrt(1.0 - pow(rho, 2)) ) );
  }
  else if (t == T) {
    beta_t_pro = exp( R::rnorm( X_mu(T) + rho * (log(beta(T - 1)) - X_mu(T - 1)), 
      sigma_beta * sqrt(1.0 - pow(rho, 2)) ) );
  }
  else {
    beta_t_pro = exp( R::rnorm( X_mu(t) + rho * (log(beta(t - 1)) - X_mu(t - 1)) / (1 + pow(rho, 2)) +
      rho * (log(beta(t + 1)) - X_mu(t + 1)) / (1 + pow(rho, 2)),
      sigma_beta * sqrt( (1.0 - pow(rho, 2)) / (1.0 + pow(rho, 2)) ) ) );
  }
  
  // Update S, I_U, R_U and gamma accordingly
  vec S_pro(T + 1);
  vec I_U_pro(T + 1);
  vec R_U_pro(T + 1);

  vec gamma_pro(T + 1); 
  vec gamma_tilde_pro(T + 1);

  vec Y_eta = Y * eta;

  double A_t;

  int tt, stop_trigger = 0;

  S_pro(0) = S(0);
  I_U_pro(0) = I_U(0);
  R_U_pro(0) = R_U(0);
  
  gamma_pro(0) = B(0) / ( (1.0 - alpha) * I_U_pro(0) );
  gamma_tilde_pro(0) = gamma_link(gamma_pro(0), link);
  

  for (tt = 0; tt < T; tt++) {

    if ( (1.0 - alpha) * I_U_pro(tt) <= B(tt) ) {
      stop_trigger = 1;
      break;
    }
    
    if (tt == t) {
      A_t = beta_t_pro * S_pro(tt) * (I_U_pro(tt) + kappa * I_D(tt)) / N;
    }
    else {
      A_t = beta(tt) * S_pro(tt) * (I_U_pro(tt) + kappa * I_D(tt)) / N;
    }
    
    S_pro(tt + 1) = S_pro(tt) - A_t;
    I_U_pro(tt + 1) = (1.0 - alpha) * I_U_pro(tt) + A_t - B(tt);
    R_U_pro(tt + 1) = R_U_pro(tt) + alpha * I_U_pro(tt);
    
    gamma_pro(tt + 1) = B(tt + 1) / ( (1.0 - alpha) * I_U_pro(tt + 1) );
    gamma_tilde_pro(tt + 1) = gamma_link(gamma_pro(tt + 1), link);

  }

  if ( (1.0 - alpha) * I_U_pro(T) <= B(T) ) {
    stop_trigger = 1;
  }

  // Evaluate proposal likelihood
  if (stop_trigger == 0) {
    
    double loglik_diff = 0, u;

    for (tt = 0; tt <= T; tt++) {

      loglik_diff = loglik_diff - (0.5 / pow(sigma_gamma, 2)) * 
        ( pow(gamma_tilde_pro(tt) - Y_eta(tt), 2) - 
          pow(gamma_link(gamma(tt), link) - Y_eta(tt), 2) );
    
    }
    
    loglik_diff = loglik_diff / Tmp;

    u = R::runif(0.0, 1.0);

    // Accept proposal
    if(log(u) < loglik_diff){
      
      beta(t) = beta_t_pro;

      for (tt = 0; tt <= T; tt++) {
        
        S(tt) = S_pro(tt);
        I_U(tt) = I_U_pro(tt);
        R_U(tt) = R_U_pro(tt);
        
        gamma(tt) = gamma_pro(tt);
        
      }

    }

  }

  return;

}





void update_mu_sigma_beta_rho (vec &beta, mat &X, vec &mu, double *sigma_beta, double *rho,
  int Q, int T,
  vec &mu_tilde, vec &sigma_mu) {

  ////////// Update mu
  int q, t;
  double sigma_beta_cur, sigma_beta_pro;
  double rho_cur, rho_pro;

  sigma_beta_cur = *sigma_beta;
  rho_cur = *rho;

  mat mu_post_cov;
  vec mu_post_mean;

  mat inv_Sigma_mu;
  inv_Sigma_mu = eye(Q, Q);

  for (q = 0; q < Q; q++) {
    inv_Sigma_mu(q, q) = 1.0 / pow(sigma_mu(q), 2);
  }

  mat inv_Sigma_beta(T + 1, T + 1);
  inv_Sigma_beta.fill(0.0);

  inv_Sigma_beta(0, 0) = 1.0;
  inv_Sigma_beta(0, 1) = -rho_cur;
  for (t = 1; t < T; t++) {
    inv_Sigma_beta(t, t - 1) = -rho_cur;
    inv_Sigma_beta(t, t) = 1.0 + pow(rho_cur, 2);
    inv_Sigma_beta(t, t + 1) = -rho_cur;
  }
  inv_Sigma_beta(T, T - 1) = -rho_cur;
  inv_Sigma_beta(T, T) = 1.0;
  
  // Check this!!!
  inv_Sigma_beta = inv_Sigma_beta / ( 1.0 - pow(rho_cur, 2) );
  
  // Rcout << inv_Sigma_mu << std::endl; 
  // Rcout << pow(1.0 / sigma_beta_cur, 2) << std::endl; 
  // Rcout << inv_Sigma_beta(0, 0) << inv_Sigma_beta(0, 1) << std::endl; 
  // Rcout << trans(X) * inv_Sigma_beta * X << std::endl; 

  mu_post_cov = inv(inv_Sigma_mu + 
    pow(1.0 / sigma_beta_cur, 2) * trans(X) * inv_Sigma_beta * X);
  
  mu_post_mean = mu_post_cov * (inv_Sigma_mu * mu_tilde + 
    pow(1.0 / sigma_beta_cur, 2) * trans(X) * inv_Sigma_beta * log(beta));

  mu = mvrnorm(mu_post_mean, mu_post_cov);

  

  ////////// Update sigma_beta

  double sigma_beta_post_shape, sigma_beta_post_rate;
  mat RSS = trans(log(beta) - X * mu) * inv_Sigma_beta * (log(beta) - X * mu);
  
  // Rcout << RSS << std::endl;

  sigma_beta_post_shape = 11.0 + 0.5 * ((double) T + 1.0);
  sigma_beta_post_rate = 1.0 + 0.5 * RSS(0, 0);
  
  sigma_beta_pro = 1.0 / sqrt(R::rgamma(sigma_beta_post_shape, 1.0 / sigma_beta_post_rate));
  
  *sigma_beta = sigma_beta_pro;

  

  ////////// Update rho

  double h;

  h = R::rnorm(0.0, 0.05);

  rho_pro = rho_cur * exp(h);
  
  if (rho_pro < 0.98) {

    mat inv_Sigma_beta_pro(T + 1, T + 1);
    inv_Sigma_beta_pro.fill(0.0);
    
    inv_Sigma_beta_pro(0, 0) = 1.0;
    inv_Sigma_beta_pro(0, 1) = -rho_pro;
    for (t = 1; t < T; t++) {
      inv_Sigma_beta_pro(t, t - 1) = -rho_pro;
      inv_Sigma_beta_pro(t, t) = 1.0 + pow(rho_pro, 2);
      inv_Sigma_beta_pro(t, t + 1) = -rho_pro;
    }
    inv_Sigma_beta_pro(T, T - 1) = -rho_pro;
    inv_Sigma_beta_pro(T, T) = 1.0;
    
    inv_Sigma_beta_pro = inv_Sigma_beta_pro / ( 1.0 - pow(rho_pro, 2) );
    
    mat RSS_pro = trans(log(beta) - X * mu) * inv_Sigma_beta_pro * (log(beta) - X * mu);
    
    double loglik_diff, logprior_diff, log_Jacobian, u;

    loglik_diff = - 0.5 * (double) T * 
                        ( log( 1.0 - pow(rho_pro, 2) ) - log( 1.0 - pow(rho_cur, 2) ) ) - 
                    0.5 * pow(1.0 / sigma_beta_pro, 2) * ( RSS_pro(0, 0) - RSS(0, 0) );

    logprior_diff = (4.0 - 1.0) * h + 
                    (1.0 - 1.0) * ( log(1.0 - rho_pro) - log(1.0 - rho_cur) );

    log_Jacobian = h;

    u = R::runif(0.0, 1.0);

    // Accept proposal
    if (log(u) < loglik_diff + logprior_diff + log_Jacobian) {
      
      *rho = rho_pro;

    }

  }

  
  return;

}











void update_alpha (vec &S, vec &I_U, vec &I_D, vec &R_U, vec &R_D, 
  vec &B, double N,
  vec &beta, 
  vec &gamma, mat &Y, vec &eta, double sigma_gamma,
  double kappa, double *alpha, int T,
  double nu_alpha_1, double nu_alpha_2, int link,
  double Tmp) {
  
  double alpha_cur = *alpha;
  double alpha_pro; 
  
  vec S_pro(T + 1);
  vec I_U_pro(T + 1);
  vec I_D_pro(T + 1);
  vec R_U_pro(T + 1);
  vec R_D_pro(T + 1);
  
  vec gamma_pro(T + 1); 
  vec gamma_tilde_pro(T + 1);

  vec Y_eta = Y * eta;

  double h, A_t;

  int t, stop_trigger = 0;
  
  h = R::rnorm(0, 0.05);
  
  alpha_pro = alpha_cur * exp(h);
  
  if (alpha_pro > 1) {
    stop_trigger = 1;
  }
  else {
    
    S_pro(0) = S(0);
    I_U_pro(0) = I_U(0);
    I_D_pro(0) = I_D(0);
    R_U_pro(0) = R_U(0);
    R_D_pro(0) = R_D(0);
    
    gamma_pro(0) = B(0) / ( (1.0 - alpha_pro) * I_U_pro(0) );
    gamma_tilde_pro(0) = gamma_link(gamma_pro(0), link);
    
    for (t = 0; t < T; t++) {
    
      if ( (1.0 - alpha_pro) * I_U_pro(t) <= B(t) ) {
        stop_trigger = 1;
        break;
      }
      
      A_t = beta(t) * S_pro(t) * (I_U_pro(t) + kappa * I_D_pro(t)) / N;
      
      S_pro(t + 1) = S_pro(t) - A_t;
      I_U_pro(t + 1) = (1 - alpha_pro) * I_U_pro(t) + A_t - B(t);
      I_D_pro(t + 1) = (1 - alpha_pro) * I_D_pro(t) + B(t);
      R_U_pro(t + 1) = R_U_pro(t) + alpha_pro * I_U_pro(t);
      R_D_pro(t + 1) = R_D_pro(t) + alpha_pro * I_D_pro(t);
    
      gamma_pro(t + 1) = B(t + 1) / ( (1.0 - alpha_pro) * I_U_pro(t + 1) );
      gamma_tilde_pro(t + 1) = gamma_link(gamma_pro(t + 1), link);

    }

    if ( (1.0 - alpha_pro) * I_U_pro(T) <= B(T) ) {
      stop_trigger = 1;
    }

  }
  

  if (stop_trigger == 0) {
    
    double loglik_diff = 0, logprior_diff = 0, log_Jacobian = 0, u;

    for (t = 0; t <= T; t++) {

      loglik_diff = loglik_diff - (0.5 / pow(sigma_gamma, 2)) * 
        ( pow(gamma_tilde_pro(t) - Y_eta(t), 2) - 
          pow(gamma_link(gamma(t), link) - Y_eta(t), 2) );
    
    }
    
    loglik_diff = loglik_diff / Tmp;
      
    logprior_diff = logprior_diff - 
      (nu_alpha_1 - 1.0) * h - 
      nu_alpha_2 * (1.0 / alpha_pro - 1.0 / alpha_cur);

    // log_Jacobian = log(alpha_pro) - log(alpha_cur);
    log_Jacobian = h;

    u = R::runif(0.0, 1.0);

    // Accept proposal
    if(log(u) < loglik_diff + logprior_diff + log_Jacobian){
      
      *alpha = alpha_pro;

      for (t = 0; t <= T; t++) {
        
        S(t) = S_pro(t);
        I_U(t) = I_U_pro(t);
        I_D(t) = I_D_pro(t);
        R_U(t) = R_U_pro(t);
        R_D(t) = R_D_pro(t);
        
        gamma(t) = gamma_pro(t);
        
      }

    }

  }
  
  return;

}








void update_eta (vec &gamma, mat &Y, vec &eta, double sigma_gamma, int K,
  vec &eta_tilde, double sigma_eta,
  int T, int link, double Tmp) {

  vec gamma_tilde(T + 1);
  
  int t;
  for (t = 0; t <= T; t++) {
    gamma_tilde(t) = gamma_link(gamma(t), link);
  }

  // vec eta_pro;
  mat eta_post_cov;
  vec eta_post_mean;

  eta_post_cov = inv(pow(1.0 / sigma_eta, 2) * eye(K, K) + 
    pow(1.0 / sigma_gamma, 2) * trans(Y) * Y / Tmp);

  eta_post_mean = eta_post_cov * ( pow(1.0 / sigma_eta, 2) * eye(K, K) * eta_tilde + 
    pow(1.0 / sigma_gamma, 2) * trans(Y) * gamma_tilde / Tmp );
  
  // Rcout << eta_tilde << "\n";
  // Rcout << eta_post_mean << "\n";

  eta = mvrnorm(eta_post_mean, eta_post_cov);
  
  // int k;
  // for (k = 0; k < K; k++) {
  //   eta[k] = eta_pro(k);
  // }

  return;

}




void update_sigma_gamma (vec &gamma, mat &Y, vec &eta, double *sigma_gamma,
  int T, int link, double Tmp) {
  
  double sigma_gamma_post_shape, sigma_gamma_post_rate;
  int t;

  vec Y_eta = Y * eta;
  
  sigma_gamma_post_shape = 1.0 + 0.5 * ((double) T + 1.0) / Tmp;
  sigma_gamma_post_rate = 1.0;

  for (t = 0; t <= T; t++) {
    sigma_gamma_post_rate = sigma_gamma_post_rate + 
      0.5 * pow(gamma_link(gamma(t), link) - Y_eta(t), 2) / Tmp;
  }
  
  *sigma_gamma = 1.0 / sqrt(R::rgamma(sigma_gamma_post_shape, 1.0 / sigma_gamma_post_rate));

  return;

}







double calc_loglik (vec &gamma, mat &Y, vec &eta, double sigma_gamma, int link, int T) {
  
  double loglik = 0;
  int t;
  vec Y_eta = Y * eta;

  for (t = 0; t <= T; t++) {
    loglik = loglik - log(sigma_gamma) - 
      (0.5 / pow(sigma_gamma, 2)) * pow(gamma_link(gamma(t), link) - Y_eta(t), 2);
  }

  return loglik;
  
}







// [[Rcpp::export]]
RcppExport SEXP BaySIR_MCMC_cpp (SEXP S, SEXP I_U, SEXP I_D, SEXP R_U, SEXP R_D,
  SEXP B, SEXP N,
  SEXP beta, SEXP X, SEXP mu, SEXP sigma_beta, SEXP rho,
  SEXP gamma, SEXP Y, SEXP eta, SEXP sigma_gamma, 
  SEXP kappa, SEXP alpha, 
  SEXP Q, SEXP K, SEXP T, 
  SEXP nu_1, SEXP nu_2, 
  SEXP mu_tilde, SEXP sigma_mu,
  SEXP eta_tilde, SEXP sigma_eta,
  SEXP nu_alpha_1, SEXP nu_alpha_2, 
  SEXP Delta, SEXP M,
  SEXP link,
  SEXP niter) {
  
  // All variables
  NumericMatrix S_(S);
  NumericMatrix I_U_(I_U);
  NumericMatrix I_D_(I_D);
  NumericMatrix R_U_(R_U);
  NumericMatrix R_D_(R_D);
  
  NumericMatrix beta_(beta);
  NumericMatrix mu_(mu);
  NumericVector sigma_beta_(sigma_beta);
  NumericVector rho_(rho);
  
  NumericMatrix gamma_(gamma);
  NumericMatrix eta_(eta);
  NumericVector sigma_gamma_(sigma_gamma);

  NumericVector alpha_(alpha);
  
  // Observables
  NumericVector B_(B);

  double N_ = as <double> (N);

  NumericMatrix X_(X);
  NumericMatrix Y_(Y);

  double kappa_ = as <double> (kappa);
  
  int Q_ = as <int> (Q);
  int K_ = as <int> (K);
  int T_ = as <int> (T);

  double nu_1_ = as <double> (nu_1);
  double nu_2_ = as <double> (nu_2);

  double nu_alpha_1_ = as <double> (nu_alpha_1);
  double nu_alpha_2_ = as <double> (nu_alpha_2);

  NumericVector mu_tilde_(mu_tilde); 
  NumericVector sigma_mu_(sigma_mu);

  NumericVector eta_tilde_(eta_tilde); 
  double sigma_eta_ = as <double> (sigma_eta);
  
  NumericVector Delta_(Delta);
  int M_ = as <int> (M);
  int link_ = as <int> (link);
  int niter_ = as <int> (niter);



  // initialize arma vecs
  // variables
  mat S_arma(S_.begin(), T_ + 1, M_, false);
  mat I_U_arma(I_U_.begin(), T_ + 1, M_, false);
  mat I_D_arma(I_D_.begin(), T_ + 1, M_, false);
  mat R_U_arma(R_U_.begin(), T_ + 1, M_, false);
  mat R_D_arma(R_D_.begin(), T_ + 1, M_, false);

  mat beta_arma(beta_.begin(), T_ + 1, M_, false);
  mat mu_arma(mu_.begin(), Q_, M_, false);
  vec sigma_beta_arma(sigma_beta_.begin(), M_, false);
  vec rho_arma(rho_.begin(), M_, false);
  
  mat gamma_arma(gamma_.begin(), T_ + 1, M_, false);
  mat eta_arma(eta_.begin(), K_, M_, false);
  vec sigma_gamma_arma(sigma_gamma_.begin(), M_, false);

  vec alpha_arma(alpha_.begin(), M_, false);

  // observables
  vec B_arma(B_.begin(), T_ + 1, false);

  mat X_arma(X_.begin(), T_ + 1, Q_, false);
  mat Y_arma(Y_.begin(), T_ + 1, K_, false);

  vec mu_tilde_arma(mu_tilde_.begin(), Q_, false); 
  vec sigma_mu_arma(sigma_mu_.begin(), Q_, false); 

  vec eta_tilde_arma(eta_tilde_.begin(), K_, false); 
  
  vec Delta_arma(Delta_.begin(), M_, false);
  
  vec loglik(M_);

  int i, t, m;
  double u, log_prob_swap, temp_swap;

  vec S_arma_m(T_ + 1);
  vec I_U_arma_m(T_ + 1);
  vec I_D_arma_m(T_ + 1);
  vec R_U_arma_m(T_ + 1);
  vec R_D_arma_m(T_ + 1);
  
  vec beta_arma_m(T_ + 1);
  vec mu_arma_m(Q_);
  
  vec gamma_arma_m(T_ + 1);
  vec eta_arma_m(K_);


  for (i = 0; i < niter_; i++){
    
    
    for (m = 0; m < M_; m++) {
      
      S_arma_m = S_arma.col(m);
      I_U_arma_m = I_U_arma.col(m);
      I_D_arma_m = I_D_arma.col(m);
      R_U_arma_m = R_U_arma.col(m);
      R_D_arma_m = R_D_arma.col(m);
      
      beta_arma_m = beta_arma.col(m);
      mu_arma_m = mu_arma.col(m);
      
      gamma_arma_m = gamma_arma.col(m);
      eta_arma_m = eta_arma.col(m);

      update_I_U_0(S_arma_m, I_U_arma_m, I_D_arma_m, R_U_arma_m, R_D_arma_m, 
        B_arma, N_,
        beta_arma_m, gamma_arma_m, Y_arma, eta_arma_m, sigma_gamma_arma(m), 
        kappa_, alpha_arma(m), T_, nu_1_, nu_2_, link_, Delta_arma(m));

      for (t = 0; t <= T_; t++) {
        
        update_beta_t(S_arma_m, I_U_arma_m, I_D_arma_m, R_U_arma_m, R_D_arma_m, 
          B_arma, N_,
          beta_arma_m, X_arma, mu_arma_m, sigma_beta_arma(m), rho_arma(m), 
          gamma_arma_m, Y_arma, eta_arma_m, sigma_gamma_arma(m),
          kappa_, alpha_arma(m), t, T_, link_, Delta_arma(m));

      }
      
      update_mu_sigma_beta_rho (beta_arma_m, X_arma, mu_arma_m, 
        &sigma_beta_arma(m), &rho_arma(m),
        Q_, T_, mu_tilde_arma, sigma_mu_arma);

      update_alpha(S_arma_m, I_U_arma_m, I_D_arma_m, R_U_arma_m, R_D_arma_m, 
        B_arma, N_,
        beta_arma_m, gamma_arma_m, Y_arma, eta_arma_m, sigma_gamma_arma(m),
        kappa_, &alpha_arma(m), T_,
        nu_alpha_1_, nu_alpha_2_, link_, Delta_arma(m));
  
      update_eta(gamma_arma_m, Y_arma, eta_arma_m, sigma_gamma_arma(m), K_, 
        eta_tilde_arma, sigma_eta_, T_, link_, Delta_arma(m));
      
      update_sigma_gamma(gamma_arma_m, Y_arma, eta_arma_m, &sigma_gamma_arma(m), T_, 
        link_, Delta_arma(m));
      
      loglik(m) = calc_loglik(gamma_arma_m, Y_arma, eta_arma_m, sigma_gamma_arma(m), 
        link_, T_);
      
      S_arma.col(m) = S_arma_m;
      I_U_arma.col(m) = I_U_arma_m;
      I_D_arma.col(m) = I_D_arma_m;
      R_U_arma.col(m) = R_U_arma_m;
      R_D_arma.col(m) = R_D_arma_m;
      
      beta_arma.col(m) = beta_arma_m;
      mu_arma.col(m) = mu_arma_m;
      
      gamma_arma.col(m) = gamma_arma_m;
      eta_arma.col(m) = eta_arma_m;

    }
    
    
    for (m = 0; m < M_ - 1; m++) {
      
      u = R::runif(0.0, 1.0);

      log_prob_swap = (1 / Delta_arma(m) - 1 / Delta_arma(m+1)) * 
        (loglik(m + 1) - loglik(m));
      
      if (log(u) < log_prob_swap) {
      // if (log(u) < log_prob_swap) {
        
        S_arma.swap_cols(m, m + 1);
        I_U_arma.swap_cols(m, m + 1);
        I_D_arma.swap_cols(m, m + 1);
        R_U_arma.swap_cols(m, m + 1);
        R_D_arma.swap_cols(m, m + 1);
        
        beta_arma.swap_cols(m, m + 1);
        mu_arma.swap_cols(m, m + 1);
        gamma_arma.swap_cols(m, m + 1);
        eta_arma.swap_cols(m, m + 1);
        
        temp_swap = sigma_beta_arma(m);
        sigma_beta_arma(m) = sigma_beta_arma(m + 1);
        sigma_beta_arma(m + 1) = temp_swap;

        temp_swap = rho_arma(m);
        rho_arma(m) = rho_arma(m + 1);
        rho_arma(m + 1) = temp_swap;

        temp_swap = sigma_gamma_arma(m);
        sigma_gamma_arma(m) = sigma_gamma_arma(m + 1);
        sigma_gamma_arma(m + 1) = temp_swap;

        temp_swap = alpha_arma(m);
        alpha_arma(m) = alpha_arma(m + 1);
        alpha_arma(m + 1) = temp_swap;

      }
      

    }

  }
  
  List BaySIR_MCMC_output;
  
  BaySIR_MCMC_output["S"] = S_arma;
  BaySIR_MCMC_output["I_U"] = I_U_arma;
  BaySIR_MCMC_output["I_D"] = I_D_arma;
  BaySIR_MCMC_output["R_U"] = R_U_arma;
  BaySIR_MCMC_output["R_D"] = R_D_arma;

  BaySIR_MCMC_output["beta"] = beta_arma;
  BaySIR_MCMC_output["mu"] = mu_arma;
  BaySIR_MCMC_output["sigma_beta"] = sigma_beta_arma;
  BaySIR_MCMC_output["rho"] = rho_arma;

  BaySIR_MCMC_output["gamma"] = gamma_arma;
  BaySIR_MCMC_output["eta"] = eta_arma;
  BaySIR_MCMC_output["sigma_gamma"] = sigma_gamma_arma;
  
  BaySIR_MCMC_output["alpha"] = alpha_arma;

  return(wrap(BaySIR_MCMC_output));

}





// [[Rcpp::export]]
RcppExport SEXP BaySIR_predict_GP (SEXP beta, SEXP mu, SEXP sigma_beta, SEXP rho,
  SEXP X, SEXP X_pred, SEXP T, SEXP T_pred, SEXP Q) {
  
  // All variables
  NumericVector beta_(beta);
  NumericVector mu_(mu);
  double sigma_beta_ = as <double> (sigma_beta);
  double rho_ = as <double> (rho);
  NumericMatrix X_(X);
  NumericMatrix X_pred_(X_pred);
  int T_ = as <int> (T);
  int T_pred_ = as <int> (T_pred);
  int Q_ = as <int> (Q);

  // initialize arma vecs
  // variables
  vec beta_arma(beta_.begin(), T_ + 1, false);
  vec mu_arma(mu_.begin(), Q_, false);
  mat X_arma(X_.begin(), T_ + 1, Q_, false);
  mat X_pred_arma(X_pred_.begin(), T_pred_, Q_, false);
  
  mat Sigma_beta(T_ + T_pred_ + 1, T_ + T_pred_ + 1);
  Sigma_beta.fill(0.0);
  

  int i, j;

  for (i = 0; i < T_ + T_pred_ + 1; i++) {
    for (j = 0; j < T_ + T_pred_ + 1; j++) {
      Sigma_beta(i, j) = pow(rho_, (double) abs(i - j));
    }
  }
  
  Sigma_beta = pow(sigma_beta_, 2) * Sigma_beta;


  mat Sigma_beta_11 = Sigma_beta.submat(0, 0, T_, T_);
  mat Sigma_beta_21 = Sigma_beta.submat(T_ + 1, 0, T_ + T_pred_, T_);
  mat Sigma_beta_22 = Sigma_beta.submat(T_ + 1, T_ + 1, T_ + T_pred_, T_ + T_pred_);
  mat inv_Sigma_beta_11 = inv(Sigma_beta_11);
  
  vec beta_pred_mean = X_pred_arma * mu_arma + 
            Sigma_beta_21 * inv_Sigma_beta_11 * (log(beta_arma) - X_arma * mu_arma);

  mat beta_pred_cov = Sigma_beta_22 - Sigma_beta_21 * inv_Sigma_beta_11 * trans(Sigma_beta_21);
  
  // Rcout << beta_pred_cov << std::endl; 

  vec beta_pred = exp(mvrnorm(beta_pred_mean, beta_pred_cov));

  return(wrap(beta_pred));

}


