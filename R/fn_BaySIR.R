

# library(Rcpp)
# library(RcppArmadillo)
# sourceCpp("fn_BaySIR_MCMC.cpp")



BaySIR_MCMC = function(B, I_D_0, N, 
  confirmed_cases_cum = NULL,
  X = NULL, Y = NULL, 
  kappa = 1,
  niter = 1000, burnin = 10000, thin = 20, Delta = (1.5)^(0:9)) {
  
  if (missing(B) | missing(I_D_0)) {
    
    if (!is.null(confirmed_cases_cum)) {
      B = confirmed_cases_cum[-1] - head(confirmed_cases_cum, -1)
      I_D_0 = confirmed_cases_cum[1]
    } else {
      stop("Must specify (B, I_D_0) or confirmed_cases_cum.")
    }

  }

  if (any(B < 0)) {
    stop("Error: some daily confirmed cases (B) are negative.")
  }

  if (I_D_0 <= 0) {
    stop("Error: initial number of documented infections (I_D[0]) is not positive.")
  }

  if (missing(N)) {
    stop("Must specify N (population size).")
  }
  
  # B: length T + 1 vector, avoid gamma being 0
  B[B == 0] = 1
  T = length(B) - 1

  # may change 
  if (is.null(X)) {
    X = cbind(rep(1, T + 1), 0:T)
  }
  
  if (is.null(Y)) {
    Y = as.matrix(rep(1, T + 1))
  }
  # Y = cbind(rep(1, T + 1), 0:T) 
  # Y = as.matrix(log(z)) 

  if (!is.matrix(X)) {
    stop("X must be a matrix.")
  }

  if (!is.matrix(Y)) {
    stop("Y must be a matrix.")
  }

  if (dim(X)[1] != (T + 1)) {
    stop("Dimension of X does not match with the time series.")
  }
  
  if (dim(Y)[1] != (T + 1)) {
    stop("Dimension of Y does not match with the time series.")
  }

  Q = dim(X)[2]
  K = dim(Y)[2]
  
  M = length(Delta)
  
  #################################################################
  # Setting hyperparameters
  #################################################################
  nu_1 = 5
  nu_2 = 1

  nu_alpha_1 = 175
  nu_alpha_2 = 25

  mu_tilde = c(log(2.5 * nu_alpha_2 / nu_alpha_1),rep(0, Q-1))
  sigma_mu = c(0.3, rep(1, Q-1))

  eta_tilde = c(log(-log(1 - 0.5)), rep(0, K-1))
  sigma_eta = 1

  #################################################################
  # Initialize Markov chains
  #################################################################
  S_spls = matrix(0, T + 1, niter)
  I_U_spls = matrix(0, T + 1, niter)
  I_D_spls = matrix(0, T + 1, niter)
  R_U_spls = matrix(0, T + 1, niter)
  R_D_spls = matrix(0, T + 1, niter)

  mu_spls = matrix(0, Q, niter)
  beta_spls = matrix(0, T + 1, niter)
  sigma_beta_spls = rep(0, niter)
  rho_spls = rep(0, niter)
  
  eta_spls = matrix(0, K, niter)
  gamma_spls = matrix(0, T + 1, niter)
  sigma_lambda_spls = rep(0, niter)

  alpha_spls = rep(0, niter)
  
  #################################################################
  # Set initial values
  #################################################################
  alpha_spls[1] = nu_alpha_2 / nu_alpha_1 

  I_D_spls[1, 1] = I_D_0
  R_D_spls[1, 1] = 0

  for (t in 1:T) {
    
    I_D_spls[t + 1, 1] = (1 - alpha_spls[1]) * I_D_spls[t, 1] + B[t]
    R_D_spls[t + 1, 1] = R_D_spls[t, 1] + alpha_spls[1] * I_D_spls[t, 1]

  }

  I_U_spls[1, 1] = 5 * I_D_spls[1, 1]
  S_spls[1, 1] = N - I_U_spls[1, 1] - I_D_spls[1, 1]
  R_U_spls[1, 1] = 0
  
  rho_spls[1] = 0.9
  sigma_beta_spls[1] = 0.1
  
  # basic 
  BRN_start = 2.5
  BRN_end = 2

  while (any(I_U_spls[ , 1] <= 0) | 
         any(gamma_spls[ , 1] <= 0) |
         any(gamma_spls[ , 1] >= 1)) {

    mu_spls[ , 1] = c(log(BRN_start * alpha_spls[1]), (log(BRN_end * alpha_spls[1]) - log(BRN_start * alpha_spls[1])) / T)
    beta_spls[ , 1] = exp(c(X %*% mu_spls[ , 1] + rnorm(T + 1, 0, sigma_beta_spls[1])))

    for (t in 1:T) {
    
      A_t = beta_spls[t, 1] * S_spls[t, 1] * (I_U_spls[t, 1] + kappa * I_D_spls[t, 1]) / N
      S_spls[t + 1, 1] = S_spls[t, 1] - A_t 
      I_U_spls[t + 1, 1] = (1 - alpha_spls[1]) * I_U_spls[t, 1] + A_t - B[t]
      R_U_spls[t + 1, 1] = R_U_spls[t, 1] + alpha_spls[1] * I_U_spls[t, 1]
  
    }
    
    gamma_spls[ , 1] = B / ((1 - alpha_spls[1]) * I_U_spls[ , 1])

    BRN_start = BRN_start + 0.05
    BRN_end = BRN_end + 0.05
    # print(sprintf("Try again, BRN_start = %.2f, BRN_end = %.2f", BRN_start, BRN_end))

  }
  
  
  eta_spls[ , 1] = c(log(-log(1 - 0.2)), rep(0, K-1))
  sigma_lambda_spls[1] = 0.01
  

  

  
  #################################################################
  # Set parallel tempering chains
  #################################################################
  S_PT = matrix(0, T + 1, M)
  I_U_PT = matrix(0, T + 1, M)
  I_D_PT = matrix(0, T + 1, M)
  R_U_PT = matrix(0, T + 1, M)
  R_D_PT = matrix(0, T + 1, M)
  
  mu_PT = matrix(0, Q, M)
  beta_PT = matrix(0, T + 1, M)
  sigma_beta_PT = rep(0, M)
  rho_PT = rep(0, M)

  eta_PT = matrix(0, K, M)
  gamma_PT = matrix(0, T + 1, M)
  sigma_lambda_PT = rep(0, M)

  alpha_PT = rep(0, M)
  
  S_PT[] = S_spls[ , 1]
  I_U_PT[] = I_U_spls[ , 1]
  I_D_PT[] = I_D_spls[ , 1]
  R_U_PT[] = R_U_spls[ , 1]
  R_D_PT[] = R_D_spls[ , 1]
  
  mu_PT[] = mu_spls[ , 1]
  beta_PT[] = beta_spls[ , 1]
  sigma_beta_PT[] = sigma_beta_spls[1]
  rho_PT[] = rho_spls[1]
  
  eta_PT[] = eta_spls[ , 1]
  gamma_PT[] = gamma_spls[ , 1]
  sigma_lambda_PT[] = sigma_lambda_spls[1]
  
  alpha_PT[] = alpha_spls[1]
  

  print(sprintf("BaySIR: Posterior simulation started. Date: %s.", date()))
  
  #################################################################
  # MCMC burn-in
  #################################################################
  print(sprintf("BaySIR: Burn-in started (%d iterations). Date: %s.", burnin, date()))
  
  output_C = .Call("BaySIR_MCMC_cpp", S = S_PT, 
                      I_U = I_U_PT,
                      I_D = I_D_PT,
                      R_U = R_U_PT,
                      R_D = R_D_PT,
                      B = B,
                      N = N, 
                      X = X, 
                      mu = mu_PT, 
                      beta = beta_PT, 
                      sigma_beta = sigma_beta_PT,
                      rho = rho_PT,
                      Y = Y,
                      eta = eta_PT, 
                      gamma = gamma_PT, 
                      sigma_lambda = sigma_lambda_PT,
                      kappa = as.double(kappa),
                      alpha = alpha_PT,
                      Q = Q, 
                      K = K, 
                      T = T, 
                      nu_1 = as.double(nu_1), 
                      nu_2 = as.double(nu_2),
                      nu_alpha_1 = nu_alpha_1,
                      nu_alpha_2 = nu_alpha_2,
                      mu_tilde = mu_tilde, 
                      sigma_mu = as.double(sigma_mu),
                      eta_tilde = eta_tilde,
                      sigma_eta = as.double(sigma_eta),
                      Delta = Delta, M = M,
                      niter = burnin)
  
  S_PT = output_C$S
  I_U_PT = output_C$I_U
  I_D_PT = output_C$I_D
  R_U_PT = output_C$R_U
  R_D_PT = output_C$R_D
  
  mu_PT = output_C$mu
  beta_PT = output_C$beta
  sigma_beta_PT = output_C$sigma_beta
  rho_PT = output_C$rho
  
  eta_PT = output_C$eta
  gamma_PT = output_C$gamma
  sigma_lambda_PT = output_C$sigma_lambda

  alpha_PT = output_C$alpha
  
  S_spls[ , 1] = S_PT[ , 1]
  I_U_spls[ , 1] = I_U_PT[ , 1]
  I_D_spls[ , 1] = I_D_PT[ , 1]
  R_U_spls[ , 1] = R_U_PT[ , 1]
  R_D_spls[ , 1] = R_D_PT[ , 1]
  
  mu_spls[ , 1] = mu_PT[ , 1]
  beta_spls[ , 1] = beta_PT[ , 1]
  sigma_beta_spls[1] = sigma_beta_PT[1]
  rho_spls[1] = rho_PT[1]
  
  eta_spls[ , 1] = eta_PT[ , 1]
  gamma_spls[ , 1] = gamma_PT[ , 1]
  sigma_lambda_spls[1] = sigma_lambda_PT[1]

  alpha_spls[1] = alpha_PT[1]

  print(sprintf("BaySIR: Burn-in finished. Date: %s.", date()))
  
  #################################################################
  # MCMC iterations
  #################################################################
  print(sprintf("BaySIR: MCMC started (%d iterations, %d thin). Date: %s.", 
        niter, thin, date()))

  for (i in 2:niter) {
    
    output_C = .Call("BaySIR_MCMC_cpp", S = S_PT, 
                        I_U = I_U_PT,
                        I_D = I_D_PT,
                        R_U = R_U_PT,
                        R_D = R_D_PT,
                        B = B,
                        N = N, 
                        X = X, 
                        mu = mu_PT, 
                        beta = beta_PT, 
                        sigma_beta = sigma_beta_PT,
                        rho = rho_PT,
                        Y = Y,
                        eta = eta_PT, 
                        gamma = gamma_PT, 
                        sigma_lambda = sigma_lambda_PT,
                        kappa = kappa,
                        alpha = alpha_PT,
                        Q = Q, 
                        K = K, 
                        T = T, 
                        nu_1 = nu_1, 
                        nu_2 = nu_2,
                        nu_alpha_1 = nu_alpha_1, 
                        nu_alpha_2 = nu_alpha_2,
                        mu_tilde = mu_tilde, 
                        sigma_mu = sigma_mu,
                        eta_tilde = eta_tilde,
                        sigma_eta = sigma_eta,
                        Delta = Delta, M = M,
                        niter = thin)
  
    S_PT = output_C$S
    I_U_PT = output_C$I_U
    I_D_PT = output_C$I_D
    R_U_PT = output_C$R_U
    R_D_PT = output_C$R_D
    
    mu_PT = output_C$mu
    beta_PT = output_C$beta
    sigma_beta_PT = output_C$sigma_beta
    rho_PT = output_C$rho
    
    eta_PT = output_C$eta
    gamma_PT = output_C$gamma
    sigma_lambda_PT = output_C$sigma_lambda

    alpha_PT = output_C$alpha
    
    S_spls[ , i] = S_PT[ , 1]
    I_U_spls[ , i] = I_U_PT[ , 1]
    I_D_spls[ , i] = I_D_PT[ , 1]
    R_U_spls[ , i] = R_U_PT[ , 1]
    R_D_spls[ , i] = R_D_PT[ , 1]
    
    mu_spls[ , i] = mu_PT[ , 1]
    beta_spls[ , i] = beta_PT[ , 1]
    sigma_beta_spls[i] = sigma_beta_PT[1]
    rho_spls[i] = rho_PT[1]
    
    eta_spls[ , i] = eta_PT[ , 1]
    gamma_spls[ , i] = gamma_PT[ , 1]
    sigma_lambda_spls[i] = sigma_lambda_PT[1]

    alpha_spls[i] = alpha_PT[1]
    
    if ((i %% 200) == 0) {
      print(sprintf("BaySIR: MCMC %.2f %% finished. Date: %s.", 
            round(i / niter * 100, 2), date()))
    }
    
  }

  print(sprintf("BaySIR: MCMC finished. Date: %s.", date()))

  MCMC_spls = list()

  MCMC_spls$S_spls = S_spls
  MCMC_spls$I_U_spls = I_U_spls
  MCMC_spls$I_D_spls = I_D_spls
  MCMC_spls$R_U_spls = R_U_spls
  MCMC_spls$R_D_spls = R_D_spls

  MCMC_spls$mu_spls = mu_spls
  MCMC_spls$beta_spls = beta_spls
  MCMC_spls$sigma_beta_spls = sigma_beta_spls
  MCMC_spls$rho_spls = rho_spls

  MCMC_spls$eta_spls = eta_spls
  MCMC_spls$gamma_spls = gamma_spls
  MCMC_spls$sigma_lambda_spls = sigma_lambda_spls
  
  MCMC_spls$alpha_spls = alpha_spls

  return(MCMC_spls)

}

