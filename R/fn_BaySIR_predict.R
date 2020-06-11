

BaySIR_predict = function(T_pred = 10, MCMC_spls,
  B, I_D_0, N, 
  confirmed_cases_cum = NULL,
  X_pred = NULL, Y_pred = NULL,
  X = NULL, Y = NULL, 
  kappa = 1) {
  
  ############################################################
  ## Setup
  ############################################################
  if (missing(MCMC_spls)) {
    stop("Must supply the MCMC samples for the parameters (BaySIR_MCMC).")
  }
  
  if (missing(B) | missing(I_D_0)) {
    
    if (!is.null(confirmed_cases_cum)) {
      B = confirmed_cases_cum[-1] - head(confirmed_cases_cum, -1)
      I_D_0 = confirmed_cases_cum[1]
    } else {
      stop("Must specify (B, I_D_0) or confirmed_cases_cum.")
    }

  }

  if (any(B < 0)) stop("Error: some daily confirmed cases (B) are negative.")
  
  if (I_D_0 <= 0) stop("Error: initial number of documented infections (I_D[0]) is not positive.")

  if (missing(N)) stop("Must specify N (population size).")
  
  

  # B: length T + 1 vector, avoid gamma being 0
  B[B == 0] = 1
  T = length(B) - 1

  # may change 
  if (is.null(X_pred)) X_pred = cbind(rep(1, T_pred), T + (1:T_pred))
  if (is.null(Y_pred)) Y_pred = as.matrix(rep(1, T_pred))

  if (is.null(X)) X = cbind(rep(1, T + 1), 0:T)
  if (is.null(Y)) Y = as.matrix(rep(1, T + 1))

  if (!is.matrix(X_pred)) stop("X_pred must be a matrix.")
  if (!is.matrix(Y_pred)) stop("Y_pred must be a matrix.")
  if (!is.matrix(X)) stop("X must be a matrix.")
  if (!is.matrix(Y)) stop("Y must be a matrix.")
  
  if (dim(X_pred)[1] != T_pred) stop("Dimension of X_pred does not match with T_pred.")
  if (dim(Y_pred)[1] != T_pred) stop("Dimension of Y_pred does not match with T_pred.")
  if (dim(X)[1] != (T + 1)) stop("Dimension of X does not match with the time series.")
  if (dim(Y)[1] != (T + 1)) stop("Dimension of Y does not match with the time series.")

  if (dim(X_pred)[2] != dim(X)[2]) stop("Dimensions of X_pred and X do not match.")
  if (dim(Y_pred)[2] != dim(Y)[2]) stop("Dimensions of Y_pred and Y do not match.")

  Q = dim(X)[2]
  K = dim(Y)[2]

  
  ############################################################
  ## Retrieve MCMC samples
  ############################################################
  S_spls = MCMC_spls$S_spls
  I_U_spls = MCMC_spls$I_U_spls
  I_D_spls = MCMC_spls$I_D_spls
  R_U_spls = MCMC_spls$R_U_spls
  R_D_spls = MCMC_spls$R_D_spls

  mu_spls = MCMC_spls$mu_spls
  beta_spls = MCMC_spls$beta_spls
  sigma_beta_spls = MCMC_spls$sigma_beta_spls
  rho_spls = MCMC_spls$rho_spls

  eta_spls = MCMC_spls$eta_spls
  gamma_spls = MCMC_spls$gamma_spls
  sigma_lambda_spls = MCMC_spls$sigma_lambda_spls
  
  alpha_spls = MCMC_spls$alpha_spls

  niter = dim(S_spls)[2]
  
  ############################################################
  ## Initialize posterior predictive samples
  ############################################################
  B_pred_spls = matrix(0, T_pred, niter)

  S_pred_spls = matrix(0, T_pred + 1, niter)
  I_U_pred_spls = matrix(0, T_pred + 1, niter)
  I_D_pred_spls = matrix(0, T_pred + 1, niter)
  R_U_pred_spls = matrix(0, T_pred + 1, niter)
  R_D_pred_spls = matrix(0, T_pred + 1, niter)
  
  gamma_pred_spls = matrix(0, T_pred, niter)
  beta_pred_spls = matrix(0, T_pred, niter)
  
  ############################################################
  ## Sample from posterior predictive distribution
  ############################################################
  for (i in 1:niter) {
    
    beta = beta_spls[ , i] # length T + 1 vector
    
    mu = mu_spls[ , i]
    sigma_beta = sigma_beta_spls[i]
    rho = rho_spls[i]
    
    beta_pred_spls[ , i] = .Call("BaySIR_predict_GP", 
                                 beta = beta, mu = mu, 
                                 sigma_beta = sigma_beta, rho = rho, 
                                 X = X, X_pred = X_pred, 
                                 T = T, T_pred = T_pred)

    alpha = alpha_spls[i]

    eta = eta_spls[ , i]
    sigma_lambda = sigma_lambda_spls[i]

    
    S_last = S_spls[T + 1, i]
    I_U_last = I_U_spls[T + 1, i]
    I_D_last = I_D_spls[T + 1, i]
    R_U_last = R_U_spls[T + 1, i]
    R_D_last = R_D_spls[T + 1, i]

    beta_last = beta_spls[T + 1, i]
    B_last = B[T + 1]
    
    
    A_t = beta_last * S_last * (I_U_last + kappa * I_D_last) / N

    S_pred_spls[1, i] = S_last - A_t
    I_U_pred_spls[1, i] = (1 - alpha) * I_U_last + A_t - B_last
    I_D_pred_spls[1, i] = (1 - alpha) * I_D_last + B_last
    R_U_pred_spls[1, i] = R_U_last + alpha * I_U_last
    R_D_pred_spls[1, i] = R_D_last + alpha * I_D_last

    epsilon_pred = rnorm(T_pred, 0, sigma_lambda)
    gamma_pred_spls[ , i] = 1 - exp(-exp(c(Y_pred %*% eta + epsilon_pred)))


    for (t in 1:T_pred) {

      A_t = beta_pred_spls[t, i] * S_pred_spls[t, i] * (I_U_pred_spls[t, i] + kappa * I_D_pred_spls[t, i]) / N
      B_pred_spls[t, i] = gamma_pred_spls[t, i] * (1 - alpha) * I_U_pred_spls[t, i]

      S_pred_spls[t + 1, i] = S_pred_spls[t, i] - A_t
      I_U_pred_spls[t + 1, i] = (1 - alpha) * I_U_pred_spls[t, i] + A_t - B_pred_spls[t, i]
      I_D_pred_spls[t + 1, i] = (1 - alpha) * I_D_pred_spls[t, i] + B_pred_spls[t, i]
      R_U_pred_spls[t + 1, i] = R_U_pred_spls[t, i] + alpha * I_U_pred_spls[t, i]
      R_D_pred_spls[t + 1, i] = R_D_pred_spls[t, i] + alpha * I_D_pred_spls[t, i]

    }
    

  }
  
  pred_spls = list()
  
  pred_spls$S_pred_spls = S_pred_spls
  pred_spls$I_U_pred_spls = I_U_pred_spls
  pred_spls$I_D_pred_spls = I_D_pred_spls
  pred_spls$R_U_pred_spls = R_U_pred_spls
  pred_spls$R_D_pred_spls = R_D_pred_spls
  pred_spls$B_pred_spls = B_pred_spls
  pred_spls$beta_pred_spls = beta_pred_spls
  
  return(pred_spls)

}



