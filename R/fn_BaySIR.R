#' BaySIR
#'
#' @description
#' The BaySIR package includes two major functions, which can be used for 
#' estimating the time-varying effective reproduction number of COVID-19 in a
#' specific country/region, and predicting future case counts. 
#' Data of daily confirmed cases are needed. 
#' The two functions are:
#' \describe{
#' \item{\code{\link{BaySIR_MCMC}}}{}
#' \item{\code{\link{BaySIR_predict}}}{}
#' }
#' Use \code{?Function_Name} or \code{help(Function_Name)} to retrieve the 
#' documentation for a specific function. E.g., \code{?BaySIR_MCMC}.
#'
#' @docType package
#' @name BaySIR
#' @section Author:
#' Tianjian Zhou, \email{tjzhou95@gmail.com}
#' @seealso \code{\link{BaySIR_MCMC}} for MCMC sampling, 
#'   \code{\link{BaySIR_predict}} for posterior prediction.
NULL



###########################################################################
## 1. BaySIR: Function for MCMC sampling
###########################################################################
#' BaySIR MCMC sampling 
#' 
#' @description
#' Implementing the parallel-tempering Markov chain Monte Carlo 
#' (PTMCMC) sampling from the posterior distribution of the parameters.
#'
#' @param B A length \code{T + 1} vector of daily new confirmed cases 
#'          \code{B[0], ..., B[T]}. \code{B[t]} represents the increment in
#'          confirmed cases between day \code{t} and day \code{t + 1}.
#' @param I_D_0 The total number of confirmed cases on day 0.
#' @param N The population size.
#' @param confirmed_cases_cum Optional. A length \code{T + 2} vector of cumulative 
#'          confirmed case counts. If \code{B} and \code{I_D_0} have already been 
#'          specified, then \code{confirmed_cases_cum} will be ignored. If not both 
#'          \code{B} and \code{I_D_0} are specified and \code{confirmed_cases_cum}
#'          is supplied, then \code{confirmed_cases_cum}
#'          will be used to calculate \code{B} and \code{I_D_0}.
#' @param X Optional. A \code{(T + 1) * Q} matrix, covariates related to the disease 
#'          transmission rate. Default is an intercept term plus a time trend, 
#'          \code{X[t, ] = (1, t)}.
#' @param Y Optional. A \code{(T + 1) * K} matrix, covariates related to the diagnosis 
#'          rate. Default is only an intercept term, \code{Y[t, ] = 1}.
#' @param niter Optional. The number of MCMC samples to return. Default is \code{1000}. 
#'              The total number of MCMC iterations to be run is 
#'              \code{niter * thin + burnin}.
#' @param burnin Optional. The number of MCMC iterations that will be discarded as 
#'               initial burn-in. Default is \code{10000}.
#' @param thin Optional, meaning keep 1 draw every \code{thin} MCMC iterations. 
#'             Default is \code{20}.
#' @param link Optional. The link function to be used for the diagnosis rate (the rate 
#'             is between 0 and 1, and the link function transforms it to the real line).
#'             \code{1}, \code{2} and \code{3} represent the logit, probit and 
#'             complementary log-log link, respectively. Default is \code{1}.
#' @param Delta Optional. A monotonically increasing vector (each element is larger 
#'              than the previous) defining the temperatures of the parallel Markov 
#'              chains (parallel tempering). The first element must be \code{1},  
#'              corresponding to the original posterior. Default is \code{1.5^(0:9)}.
#' @param nu_alpha_1 Optional. Prior shape parameter for the length of the infectious
#'        period (\code{1 / alpha}). Default is \code{325.5}.
#' @param nu_alpha_2 Optional. Prior rate parameter for the length of the infectious
#'        period (\code{1 / alpha}). Default is \code{35}, meaning the prior mean for
#'        the infectious period is 9.3 days. \code{nu_alpha_1} and \code{nu_alpha_2}
#'        may be changed to reflect different modeling assumptions. For example, if 
#'        one defines the infectious period as the time from infection to recovery/death,
#'        then a larger prior mean (~ 25 days) for \code{1 / alpha} should be chosen.
#'
#' @return A list of the following:
#' \describe{
#' \item{\code{MCMC_spls}}{Again, a list of MCMC samples for the parameters.
#' \itemize{
#'   \item \code{R_eff_spls} A \code{(T + 1) * niter} matrix, MCMC samples of the 
#'                       effective reproduction number. The \code{l}-th column,
#'                       \code{R_eff_spls[ , l]}, represents the posterior samples of 
#'                       the effective reproduction numbers for \code{T + 1} days at
#'                       the \code{l}-th MCMC iteration.
#'   \item \code{THETA_spls} MCMC samples of a specific parameter \code{THETA}. 
#'                       Here, \code{THETA} can be \code{S}, \code{I_U}, \code{I_D},
#'                       \code{R_U}, \code{R_D}, 
#'                       \code{beta}, \code{mu}, \code{sigma_beta}, \code{rho},
#'                       \code{gamma}, \code{eta}, \code{sigma_gamma}, and 
#'                       \code{alpha}.
#' }}
#' \item{\code{MCMC_summary}}{Posterior summaries for the parameters.
#' \itemize{
#'   \item \code{R_eff} A \code{(T + 1) * 3} matrix, posterior summaries of the 
#'                      effective reproduction number. Columns 1, 2 and 3 correspond to
#'                      the posterior medians, 2.5% quantiles and 97.5% quantiles of 
#'                      the effective reproduction numbers for \code{T + 1} days.
#'                      For example, to access the posterior medians of R_eff, use
#'                      \code{MCMC_summary$R_eff[ , 1]}.
#'   \item \code{THETA} Posterior summaries of a specific parameter \code{THETA}. 
#'                       Here, \code{THETA} can be \code{S}, \code{I_U}, \code{I_D},
#'                       \code{R_U}, \code{R_D}, 
#'                       \code{beta}, \code{mu}, \code{sigma_beta}, \code{rho},
#'                       \code{gamma}, \code{eta}, \code{sigma_gamma}, and 
#'                       \code{alpha}.
#' }}
#' }
#' @examples
#' library(BaySIR)
#'   
#' # read data
#' data(data_sim_1)
#' B = data_sim_1$B
#' I_D_0 = data_sim_1$I_D[1]
#' N = data_sim_1$N
#' 
#' # run MCMC (may take a few minutes, depending on computer. ~ 2 mins on Macbook Pro)
#' result_list = BaySIR_MCMC(B = B, I_D_0 = I_D_0, N = N)
#' 
#' # for testing purpose, use smaller number of MCMC burn-in/iterations
#' # result_list = BaySIR_MCMC(B = B, I_D_0 = I_D_0, N = N, burnin = 1000, thin = 2, niter = 500)
#' 
#' # retrieve posterior summaries for the effective reproduction number
#' result_list$MCMC_summary$R_eff
#'
#' # End(Not run)

BaySIR_MCMC = function(B, I_D_0, N, 
  confirmed_cases_cum = NULL,
  X = NULL, Y = NULL, 
  kappa = 1,
  niter = 1000, burnin = 20000, thin = 30, 
  link = 1, Delta = (1.5)^(0:9),
  nu_alpha_1 = 325.5, nu_alpha_2 = 35) {
  
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
  
  if (T < 4) {
    stop("Must have COVID data for at least 5 days.")
  }

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

  mu_tilde = c(log(2.5 * nu_alpha_2 / nu_alpha_1),rep(0, Q-1))
  sigma_mu = c(0.3, rep(1, Q-1))
  
  eta_tilde = rep(0, K)
  # eta_tilde = c(log(-log(1 - 0.5)), rep(0, K-1))
  sigma_eta = 1

  #################################################################
  # Initialize Markov chains
  #################################################################
  S_spls = matrix(0, T + 1, niter)
  I_U_spls = matrix(0, T + 1, niter)
  I_D_spls = matrix(0, T + 1, niter)
  R_U_spls = matrix(0, T + 1, niter)
  R_D_spls = matrix(0, T + 1, niter)
  
  beta_spls = matrix(0, T + 1, niter)
  mu_spls = matrix(0, Q, niter)
  sigma_beta_spls = rep(0, niter)
  rho_spls = rep(0, niter)
  
  gamma_spls = matrix(0, T + 1, niter)
  eta_spls = matrix(0, K, niter)
  sigma_gamma_spls = rep(0, niter)

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
  sigma_beta_spls[1] = 0.02
  
  # Set initial values for the daily transmission rate, beta
  # Use a piecewise linear form  
  n_T_fragment = ceiling((T + 1) / 50)
  n_T_breakpoint = n_T_fragment + 1
  
  T_breakpoint = 
    ceiling((0:(n_T_breakpoint - 1) * T / (n_T_breakpoint - 1))) + 1

  BRN_all = list()
  for(ntb in 1:n_T_breakpoint) {
    if(ntb == 1) {
      BRN_all[[ntb]] = seq(1, 9, 0.5)
    } else {
      BRN_all[[ntb]] = seq(1, 7, 0.5)
    }
  }

  BRN_all = as.matrix(expand.grid(BRN_all))


  for (kk in 1:nrow(BRN_all)) {
    
    BRN = BRN_all[kk, ]

    for(tbb in 1:(n_T_breakpoint - 1)) {
      beta_spls[T_breakpoint[tbb]:T_breakpoint[tbb + 1], 1] = 
        seq(BRN[tbb], BRN[tbb + 1], 
            length.out = T_breakpoint[tbb + 1] - T_breakpoint[tbb] + 1) * 
          alpha_spls[1]
    }
    
    for (t in 1:T) {
    
      A_t = beta_spls[t, 1] * S_spls[t, 1] * (I_U_spls[t, 1] + kappa * I_D_spls[t, 1]) / N
      S_spls[t + 1, 1] = S_spls[t, 1] - A_t 
      I_U_spls[t + 1, 1] = (1 - alpha_spls[1]) * I_U_spls[t, 1] + A_t - B[t]
      R_U_spls[t + 1, 1] = R_U_spls[t, 1] + alpha_spls[1] * I_U_spls[t, 1]
  
    }
    
    gamma_spls[ , 1] = B / ((1 - alpha_spls[1]) * I_U_spls[ , 1])
    
    if(all(I_U_spls[ , 1] > 0) & 
       all(gamma_spls[ , 1] > 0) &
       all(gamma_spls[ , 1] < 1)) {
      break
    }
    
  }
  
  if (any(I_U_spls[ , 1] <= 0) | 
      any(gamma_spls[ , 1] <= 0) |
      any(gamma_spls[ , 1] >= 1)) {
    stop("Parameter initialization fails...")
  }
  
  times = 1:(T+1)
  lmfit_logbeta = lm(log(beta_spls[ , 1]) ~ times)
  mu_spls[ , 1] = lmfit_logbeta$coefficients

  eta_spls[ , 1] = c(-1.4, rep(0, K-1))
  sigma_gamma_spls[1] = 0.01

  print("Parameter initialization successful.")

  

  
  #################################################################
  # Set parallel tempering chains
  #################################################################
  S_PT = matrix(0, T + 1, M)
  I_U_PT = matrix(0, T + 1, M)
  I_D_PT = matrix(0, T + 1, M)
  R_U_PT = matrix(0, T + 1, M)
  R_D_PT = matrix(0, T + 1, M)
  
  beta_PT = matrix(0, T + 1, M)
  mu_PT = matrix(0, Q, M)
  sigma_beta_PT = rep(0, M)
  rho_PT = rep(0, M)
  
  gamma_PT = matrix(0, T + 1, M)
  eta_PT = matrix(0, K, M)
  sigma_gamma_PT = rep(0, M)

  alpha_PT = rep(0, M)
  
  S_PT[] = S_spls[ , 1]
  I_U_PT[] = I_U_spls[ , 1]
  I_D_PT[] = I_D_spls[ , 1]
  R_U_PT[] = R_U_spls[ , 1]
  R_D_PT[] = R_D_spls[ , 1]
  
  beta_PT[] = beta_spls[ , 1]
  mu_PT[] = mu_spls[ , 1]
  sigma_beta_PT[] = sigma_beta_spls[1]
  rho_PT[] = rho_spls[1]
  
  gamma_PT[] = gamma_spls[ , 1]
  eta_PT[] = eta_spls[ , 1]
  sigma_gamma_PT[] = sigma_gamma_spls[1]
  
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
                      beta = beta_PT,
                      X = X, 
                      mu = mu_PT, 
                      sigma_beta = sigma_beta_PT,
                      rho = rho_PT,
                      gamma = gamma_PT, 
                      Y = Y,
                      eta = eta_PT, 
                      sigma_gamma = sigma_gamma_PT,
                      kappa = as.double(kappa),
                      alpha = alpha_PT,
                      Q = Q, 
                      K = K, 
                      T = T, 
                      nu_1 = as.double(nu_1), 
                      nu_2 = as.double(nu_2),
                      mu_tilde = mu_tilde, 
                      sigma_mu = as.double(sigma_mu),
                      eta_tilde = eta_tilde,
                      sigma_eta = as.double(sigma_eta),
                      nu_alpha_1 = nu_alpha_1,
                      nu_alpha_2 = nu_alpha_2,
                      Delta = Delta, M = M,
                      link = link,
                      niter = burnin)
  
  S_PT = output_C$S
  I_U_PT = output_C$I_U
  I_D_PT = output_C$I_D
  R_U_PT = output_C$R_U
  R_D_PT = output_C$R_D
  
  beta_PT = output_C$beta
  mu_PT = output_C$mu
  sigma_beta_PT = output_C$sigma_beta
  rho_PT = output_C$rho
  
  gamma_PT = output_C$gamma
  eta_PT = output_C$eta
  sigma_gamma_PT = output_C$sigma_gamma

  alpha_PT = output_C$alpha
  
  S_spls[ , 1] = S_PT[ , 1]
  I_U_spls[ , 1] = I_U_PT[ , 1]
  I_D_spls[ , 1] = I_D_PT[ , 1]
  R_U_spls[ , 1] = R_U_PT[ , 1]
  R_D_spls[ , 1] = R_D_PT[ , 1]
  
  beta_spls[ , 1] = beta_PT[ , 1]
  mu_spls[ , 1] = mu_PT[ , 1]
  sigma_beta_spls[1] = sigma_beta_PT[1]
  rho_spls[1] = rho_PT[1]
  
  gamma_spls[ , 1] = gamma_PT[ , 1]
  eta_spls[ , 1] = eta_PT[ , 1]
  sigma_gamma_spls[1] = sigma_gamma_PT[1]

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
                        beta = beta_PT, 
                        X = X, 
                        mu = mu_PT, 
                        sigma_beta = sigma_beta_PT,
                        rho = rho_PT,
                        gamma = gamma_PT, 
                        Y = Y,
                        eta = eta_PT, 
                        sigma_gamma = sigma_gamma_PT,
                        kappa = kappa,
                        alpha = alpha_PT,
                        Q = Q, 
                        K = K, 
                        T = T, 
                        nu_1 = nu_1, 
                        nu_2 = nu_2,
                        mu_tilde = mu_tilde, 
                        sigma_mu = sigma_mu,
                        eta_tilde = eta_tilde,
                        sigma_eta = sigma_eta,
                        nu_alpha_1 = nu_alpha_1, 
                        nu_alpha_2 = nu_alpha_2,
                        Delta = Delta, M = M,
                        link = link,
                        niter = thin)
  
    S_PT = output_C$S
    I_U_PT = output_C$I_U
    I_D_PT = output_C$I_D
    R_U_PT = output_C$R_U
    R_D_PT = output_C$R_D
    
    beta_PT = output_C$beta
    mu_PT = output_C$mu
    sigma_beta_PT = output_C$sigma_beta
    rho_PT = output_C$rho
    
    gamma_PT = output_C$gamma
    eta_PT = output_C$eta
    sigma_gamma_PT = output_C$sigma_gamma

    alpha_PT = output_C$alpha
    
    S_spls[ , i] = S_PT[ , 1]
    I_U_spls[ , i] = I_U_PT[ , 1]
    I_D_spls[ , i] = I_D_PT[ , 1]
    R_U_spls[ , i] = R_U_PT[ , 1]
    R_D_spls[ , i] = R_D_PT[ , 1]
    
    beta_spls[ , i] = beta_PT[ , 1]
    mu_spls[ , i] = mu_PT[ , 1]
    sigma_beta_spls[i] = sigma_beta_PT[1]
    rho_spls[i] = rho_PT[1]
    
    gamma_spls[ , i] = gamma_PT[ , 1]
    eta_spls[ , i] = eta_PT[ , 1]
    sigma_gamma_spls[i] = sigma_gamma_PT[1]

    alpha_spls[i] = alpha_PT[1]
    
    if ((i %% 200) == 0) {
      print(sprintf("BaySIR: MCMC %.2f %% finished. Date: %s.", 
            round(i / niter * 100, 2), date()))
    }
    
  }

  print(sprintf("BaySIR: MCMC finished. Date: %s.", date()))
  
  R_eff_spls = sweep(beta_spls * S_spls / N, 2, alpha_spls, "/")

  ## Posterior samples
  MCMC_spls = list()

  MCMC_spls$S_spls = S_spls
  MCMC_spls$I_U_spls = I_U_spls
  MCMC_spls$I_D_spls = I_D_spls
  MCMC_spls$R_U_spls = R_U_spls
  MCMC_spls$R_D_spls = R_D_spls
  
  MCMC_spls$beta_spls = beta_spls
  MCMC_spls$mu_spls = mu_spls
  MCMC_spls$sigma_beta_spls = sigma_beta_spls
  MCMC_spls$rho_spls = rho_spls
  
  MCMC_spls$gamma_spls = gamma_spls
  MCMC_spls$eta_spls = eta_spls
  MCMC_spls$sigma_gamma_spls = sigma_gamma_spls
  
  MCMC_spls$alpha_spls = alpha_spls

  MCMC_spls$R_eff_spls = R_eff_spls

  ## Posterior summaries
  MCMC_summary = list()

  MCMC_summary$S = t(apply(S_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$I_U = t(apply(I_U_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$I_D = t(apply(I_D_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$R_U = t(apply(R_U_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$R_D = t(apply(R_D_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  
  MCMC_summary$beta = t(apply(beta_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$mu = t(apply(mu_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$sigma_beta = quantile(sigma_beta_spls, probs = c(0.5, 0.025, 0.975))
  MCMC_summary$rho = quantile(rho_spls, probs = c(0.5, 0.025, 0.975))
  
  MCMC_summary$gamma = t(apply(gamma_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$eta = t(apply(eta_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  MCMC_summary$sigma_gamma = quantile(sigma_gamma_spls, probs = c(0.5, 0.025, 0.975))
  
  MCMC_summary$alpha = quantile(alpha_spls, probs = c(0.5, 0.025, 0.975))

  MCMC_summary$R_eff = t(apply(R_eff_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  

  result_list = list()
  result_list$MCMC_spls = MCMC_spls
  result_list$MCMC_summary = MCMC_summary

  return(result_list)

}











###########################################################################
## 2. BaySIR: Function for posterior prediction
###########################################################################
#' BaySIR posterior prediction 
#' 
#' @description
#' Sampling from the posterior predictive distribution of future observations,
#' in particular, future case counts.
#'
#' @param T_pred The number of future days that you would like to predict. 
#'               Will predict B[T + 1], ..., B[T + T_pred]. Default is 10 days.
#' @param MCMC_spls A list of the MCMC samples for the parameters obtained 
#'                  from \code{\link{BaySIR_MCMC}}.
#' @param B A length \code{T + 1} vector of daily new confirmed cases 
#'          \code{B[0], ..., B[T]}. \code{B[t]} represents the increment in
#'          confirmed cases between day \code{t} and day \code{t + 1}.
#' @param I_D_0 The total number of confirmed cases on day 0.
#' @param N The population size.
#' @param confirmed_cases_cum Optional. A length \code{T + 2} vector of cumulative 
#'          confirmed case counts. If \code{B} and \code{I_D_0} have already been 
#'          specified, then \code{confirmed_cases_cum} will be ignored. If not both 
#'          \code{B} and \code{I_D_0} are specified and \code{confirmed_cases_cum}
#'          is supplied, then \code{confirmed_cases_cum}
#'          will be used to calculate \code{B} and \code{I_D_0}.
#' @param X_pred Optional. A \code{T_pred * Q} matrix, covariates related to the disease 
#'               transmission rate for future days \code{T + 1, ..., T + T_pred}. 
#'               Default is an intercept term plus a time trend, 
#'               \code{X_pred[t, ] = (1, T + t)}.
#' @param Y_pred Optional. A \code{T_pred * K} matrix, covariates related to the diagnosis 
#'               ratefor future days \code{T + 1, ..., T + T_pred}. Default 
#'               contains only an intercept term, \code{Y_pred[t, ] = 1}.
#' @param X Optional. A \code{(T + 1) * Q} matrix, covariates related to the disease 
#'          transmission rate. Default is an intercept term plus a time trend, 
#'          \code{X[t, ] = (1, t)}.
#' @param Y Optional. A \code{(T + 1) * K} matrix, covariates related to the diagnosis 
#'          rate. Default is only an intercept term, \code{Y[t, ] = 1}.
#' @param link Optional. The link function to be used for the diagnosis rate (the rate 
#'             is between 0 and 1, and the link function transforms it to the real line).
#'             \code{1}, \code{2} and \code{3} represent the logit, probit and 
#'             complementary log-log link, respectively. Default is \code{1}.
#'
#' @return A list of the following:
#' \describe{
#' \item{\code{pred_spls}}{Again, a list of samples for the parameters from their 
#'                         posterior predictive distributions.
#' \itemize{
#'   \item \code{B_pred_spls} A \code{T_pred * niter} matrix, samples of the 
#'                       B (case counts) for future days 
#'                       \code{T + 1, ..., T + T_pred} from its posterior predictive
#'                       distribution. Each column corresponds to a sample, and each
#'                       row corresponds to a day.
#'   \item \code{R_eff_pred_spls} A \code{T_pred * niter} matrix, samples of the 
#'                       effective reproduction number for future days 
#'                       \code{T + 1, ..., T + T_pred} from its posterior predictive
#'                       distribution. Each column corresponds to a sample, and each
#'                       row corresponds to a day.
#'   \item \code{THETA_pred_spls} Samples of a specific time-dependent parameter 
#'                       \code{THETA} from its posterior predictive distribution
#'                       for future days. 
#'                       Here, \code{THETA} can be \code{S}, \code{I_U},
#'                       \code{beta}, etc.
#' }}
#' \item{\code{pred_summary}}{Summaries for the posterior predictive distributions
#'                            of the parameters.
#' \itemize{
#'   \item \code{B} A \code{T_pred * 3} matrix, summaries for the posterior 
#'                  predictive distribution of B (case counts)
#'                  for future days 
#'                  \code{T + 1, ..., T + T_pred}. Columns 1, 2 and 3 correspond to
#'                  the posterior medians, 2.5% quantiles and 97.5% quantiles.
#'                  For example, to access the posterior median of B on day
#'                  \code{T + j}, use
#'                  \code{pred_summary$B[j, 1]}.
#'   \item \code{R_eff} A \code{T_pred * 3} matrix, summaries for the posterior 
#'                      predictive distribution of the 
#'                      effective reproduction number for future days 
#'                       \code{T + 1, ..., T + T_pred}. Columns 1, 2 and 3 correspond to
#'                      the posterior medians, 2.5% quantiles and 97.5% quantiles.
#'                      For example, to access the posterior median of R_eff on day
#'                      \code{T + j}, use
#'                      \code{pred_summary$R_eff[j, 1]}.
#'   \item \code{THETA} Summaries of a specific time-dependent parameter 
#'                      \code{THETA} from its posterior predictive distribution
#'                      for future days. 
#'                      Here, \code{THETA} can be \code{S}, \code{I_U},
#'                      \code{beta}, etc.
#' }}
#' }
#' @examples
#' library(BaySIR)
#'   
#' # read data
#' data(data_sim_1)
#' B = data_sim_1$B
#' I_D_0 = data_sim_1$I_D[1]
#' N = data_sim_1$N
#' 
#' # run MCMC (may take a few minutes, depending on computer. ~ 2 mins on Macbook Pro)
#' result_list = BaySIR_MCMC(B = B, I_D_0 = I_D_0, N = N)
#' 
#' # sample from posterior predictive distribution (for future 30 days)
#' predict_list = BaySIR_predict(T_pred = 30, MCMC_spls = result_list$MCMC_spls, B = B, I_D_0 = I_D_0, N = N)
#' 
#' # posterior median of future B's
#' predict_list$pred_summary$B[ , 1]
#' 
#' # posterior summary for the future effective reproduction numbers
#' predict_list$pred_summary$R_eff
#'
#' # End(Not run)

BaySIR_predict = function(T_pred = 10, MCMC_spls,
  B, I_D_0, N, 
  confirmed_cases_cum = NULL,
  X_pred = NULL, Y_pred = NULL,
  X = NULL, Y = NULL, 
  kappa = 1, link = 1) {
  
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

  beta_spls = MCMC_spls$beta_spls
  mu_spls = MCMC_spls$mu_spls
  sigma_beta_spls = MCMC_spls$sigma_beta_spls
  rho_spls = MCMC_spls$rho_spls

  gamma_spls = MCMC_spls$gamma_spls
  eta_spls = MCMC_spls$eta_spls
  sigma_gamma_spls = MCMC_spls$sigma_gamma_spls
  
  alpha_spls = MCMC_spls$alpha_spls

  niter = dim(S_spls)[2]

  if (dim(mu_spls)[1] != dim(X)[2]) {
    stop("Dimension of the posterior samples of mu does not match with X.")
  }

  if (dim(mu_spls)[1] != dim(X_pred)[2]) {
    stop("Dimension of the posterior samples of mu does not match with X_pred.")
  }
  
  if (dim(eta_spls)[1] != dim(Y)[2]) {
    stop("Dimension of the posterior samples of eta does not match with Y.")
  }

  if (dim(eta_spls)[1] != dim(Y_pred)[2]) {
    stop("Dimension of the posterior samples of eta does not match with Y_pred.")
  }
  
  ############################################################
  ## Initialize posterior predictive samples
  ############################################################
  B_pred_spls = matrix(0, T_pred, niter)

  S_pred_spls = matrix(0, T_pred + 1, niter)
  I_U_pred_spls = matrix(0, T_pred + 1, niter)
  I_D_pred_spls = matrix(0, T_pred + 1, niter)
  R_U_pred_spls = matrix(0, T_pred + 1, niter)
  R_D_pred_spls = matrix(0, T_pred + 1, niter)
  
  beta_pred_spls = matrix(0, T_pred, niter)
  gamma_pred_spls = matrix(0, T_pred, niter)
  
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
                                 T = T, T_pred = T_pred, Q = Q)

    alpha = alpha_spls[i]

    eta = eta_spls[ , i]
    sigma_gamma = sigma_gamma_spls[i]

    
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

    epsilon_pred = rnorm(T_pred, 0, sigma_gamma)
    
    if (link == 2) {
      # probit
      gamma_pred_spls[ , i] = pnorm(c(Y_pred %*% eta + epsilon_pred))
    } else if (link == 3) {
      # cloglog
      gamma_pred_spls[ , i] = 1 - exp(-exp(c(Y_pred %*% eta + epsilon_pred)))
    } else {
      # logit
      gamma_pred_spls[ , i] = 1 / (1 + exp(-c(Y_pred %*% eta + epsilon_pred)))
    }
    


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

  R_eff_pred_spls = sweep(beta_pred_spls * S_pred_spls[1:T_pred, ] / N, 2, alpha_spls, "/")
  
  pred_spls = list()
  
  pred_spls$S_pred_spls = S_pred_spls
  pred_spls$I_U_pred_spls = I_U_pred_spls
  pred_spls$I_D_pred_spls = I_D_pred_spls
  pred_spls$R_U_pred_spls = R_U_pred_spls
  pred_spls$R_D_pred_spls = R_D_pred_spls
  pred_spls$B_pred_spls = B_pred_spls
  pred_spls$beta_pred_spls = beta_pred_spls
  pred_spls$R_eff_pred_spls = R_eff_pred_spls

  pred_summary = list()
  pred_summary$S = t(apply(S_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$I_U = t(apply(I_U_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$I_D = t(apply(I_D_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$R_U = t(apply(R_U_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$R_D = t(apply(R_D_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$B = t(apply(B_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$beta = t(apply(beta_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  pred_summary$R_eff = t(apply(R_eff_pred_spls, 1, function(x) quantile(x, probs = c(0.5, 0.025, 0.975))))
  
  result_list = list()
  result_list$pred_spls = pred_spls
  result_list$pred_summary = pred_summary

  return(result_list)

}



