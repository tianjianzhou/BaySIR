# BaySIR
**Semiparametric Bayesian Inference for the Transmission Dynamics of COVID-19 with a State-Space Model**

I am actively updating the package. Please visit again for the latest version!

Link to manuscript: http://arxiv.org/abs/2006.05581

## Updates
**[6/29]** Changed the prior mean for the infectious period to 9.3 days (based on [He et al., Nature Medicine](https://www.nature.com/articles/s41591-020-0869-5)). Previously it was 7 days. Also, added different specifications for the link function of the diagnosis rate (logit, probit, cloglog). Previously cloglog was used as default. Please install again for the updated version.

**[6/23]** There has been a mistake so that the earlier version of this package does not implement the parallel tempering procedure (now fixed of course). Please install again for the updated version.


## Installation

First of all, you need [R](https://www.r-project.org/). For Windows users, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is needed. For Mac users, [Xcode](https://apps.apple.com/us/app/xcode/id497799835) may be needed.

The `BaySIR` package has two dependencies: `Rcpp` and `RcppArmadillo`.
BLAS/LAPACK libraries are needed by `RcppArmadillo`. Windows 7 users may have difficulties installing these.
The package was tested on MacOS (High Sierra 10.13.6) and Windows 10.

The `BaySIR` package can be easily installed with the `devtools` package in R. After `devtools` has been installed, run the following commands in R to install the `BaySIR` package.
```
library(devtools)
install_github("tianjianzhou/BaySIR")
```

## Examples
Please refer to [Documentation](https://github.com/tianjianzhou/BaySIR/blob/master/README.md#documentation) for details about the BaySIR functions.

### Example 1: Simulated Data

```
library(BaySIR)
  
# read data
data(data_sim_1)
B = data_sim_1$B
I_D_0 = data_sim_1$I_D[1]
N = data_sim_1$N

# run MCMC (may take a few minutes, depending on computer. ~ 2 mins on Macbook Pro)
result_list = BaySIR_MCMC(B = B, I_D_0 = I_D_0, N = N)

# for testing purpose, use smaller number of MCMC burn-in/iterations
# result_list = BaySIR_MCMC(B = B, I_D_0 = I_D_0, N = N, burnin = 1000, thin = 2, niter = 500)

# posterior summary for the effective reproduction number
result_list$MCMC_summary$R_eff

# sample from posterior predictive distribution (for future 30 days)
predict_list = BaySIR_predict(T_pred = 30, MCMC_spls = result_list$MCMC_spls, B = B, I_D_0 = I_D_0, N = N)

# posterior median of future B's
predict_list$pred_summary$B[ , 1]

# posterior summary for the future effective reproduction numbers
predict_list$pred_summary$R_eff
```

### Example 2: Real Data (Illinois) 
Data from [JHU CSSE](https://github.com/CSSEGISandData/COVID-19)
```
library(BaySIR)

# read data
data = read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv", header = TRUE)

confirmed_cases_cum = colSums(as.matrix(data[data$Province_State == "Illinois", -(1:11)]))
confirmed_cases_cum = unname(confirmed_cases_cum[confirmed_cases_cum > 100])

# population in IL
N = 12671821

# run MCMC
result_list = BaySIR_MCMC(confirmed_cases_cum = confirmed_cases_cum, N = N)

# for testing purpose, use smaller number of MCMC burn-in/iterations
# result_list = BaySIR_MCMC(confirmed_cases_cum = confirmed_cases_cum, N = N, burnin = 1000, thin = 2, niter = 500)

# posterior summary for the effective reproduction number
result_list$MCMC_summary$R_eff

# sample from posterior predictive distribution (for future 30 days)
predict_list = BaySIR_predict(T_pred = 30, MCMC_spls = result_list$MCMC_spls, confirmed_cases_cum = confirmed_cases_cum, N = N)

# posterior median of future B's
predict_list$pred_summary$B[ , 1]

# posterior summary for the future effective reproduction numbers
predict_list$pred_summary$R_eff
```

## Documentation

### 1. BaySIR_MCMC
- Usage: `BaySIR_MCMC(B, I_D_0, N, ...)`

- Arguments
  - Required
    - `B`: a length `T + 1` vector of daily new confirmed cases `B[0], ..., B[T]`.
    - `I_D_0`: the total number of confirmed cases on day 0.
    - `N`: the population size.
  - Optional
    - `confirmed_cases_cum`: an optional length `T + 2` vector of cumulative confirmed case counts. If `B` and `I_D_0` have already been specified, then `confirmed_cases_cum` will be ignored. If not both `B` and `I_D_0` are specified and `confirmed_cases_cum` is supplied, then `confirmed_cases_cum` will be used to calculate `B` and `I_D_0`.
    - `X`: a `(T + 1) * Q` matrix, covariates related to the disease transmission rate. Default is an intercept term plus a time trend, `X[t, ] = (1, t)`. It is possible to include other covariates. For example, `X[t, ] = (1, t, 1(stay-at-home order on day t), t * 1(stay-at-home order on day t))` and thus `Q = 4`.
    - `Y`: a `(T + 1) * K` matrix, covariates related to the diagnosis rate. Default contains only an intercept term, `Y[t, ] = 1`. It is possible to include other covariates. For example, `Y[t, ] = (log number of test on day t)` or `Y[t, ] = (1, t)`.
    - `niter`: number of MCMC samples to return. Default is `1000`. The total number of MCMC iterations to be run is `niter * thin + burnin`.
    - `burnin`: number of MCMC iterations that will be discarded as initial burn-in. Default is `10000`.
    - `thin`: meaning keep 1 draw every `thin` MCMC iterations. Default is `20`.
    - `Delta`: a monotonically increasing vector (each element is larger than the previous) defining the temperatures of the parallel Markov chains (parallel tempering). The first element must be 1, corresponding to the original posterior. Default is `1.5^(0:9)`.

- Output
  - A list of the following:
    - `MCMC_spls`: a list of the MCMC samples for the parameters.
    - `MCMC_summary`: a list of the posterior summaries for the parameters. For each parameter, its posterior median, 2.5% quantile and 97.5% quantile are reported.


### 2. BaySIR_predict
- Usage: `BaySIR_predict(T_pred = 10, MCMC_spls, B, I_D_0, N, ...)`

- Arguments
  - Required
    - `MCMC_spls`: a list of the MCMC samples for the parameters obtained from `BaySIR_MCMC(...)`.
    - `B`: a length `T + 1` vector of daily new confirmed cases `B[0], ..., B[T]`.
    - `I_D_0`: the total number of confirmed cases on day 0.
    - `N`: the population size.
  - Optional
    - `T_pred`: the number of future days that you would like to predict. Will predict `B[T + 1], ..., B[T + T_pred]`. Default is 10 days.
    - `confirmed_cases_cum`: an optional length `T + 2` vector of cumulative confirmed case counts. If `B` and `I_D_0` have already been specified, then `confirmed_cases_cum` will be ignored. If not both `B` and `I_D_0` are specified and `confirmed_cases_cum` is supplied, then `confirmed_cases_cum` will be used to calculate `B` and `I_D_0`.
    - `X_pred`: a `T_pred * Q` matrix, covariates related to the disease transmission rate for future days `T + 1, ..., T + T_pred`. Default is an intercept term plus a time trend, `X_pred[t, ] = (1, T + t)`. Note that if policy indicator is used as a covariate, we may not know what the policy will look like in the future and have to impute its future value. 
    - `Y_pred`: a `T_pred * K` matrix, covariates related to the diagnosis rate for future days `T + 1, ..., T + T_pred`. Default contains only an intercept term, `Y_pred[t, ] = 1`. Note that if number of tests is used as a covariate, we do not know what the number of tests will be in the future and have to impute its future value. 
    - `X`: a `(T + 1) * Q` matrix, covariates related to the disease transmission rate. Default is an intercept term plus a time trend, `X[t, ] = (1, t)`.
    - `Y`: a `(T + 1) * K` matrix, covariates related to the diagnosis rate. Default contains only an intercept term, `Y[t, ] = 1`.
    
- Output
  - A list of the following:
    - `pred_spls`: a list of the MC samples for the parameters from their posterior predictive distributions.
    - `pred_summary`: a list of the summaries for the posterior predictive distributions of the parameters. For each parameter, its posterior (predictive) median, 2.5% quantile and 97.5% quantile are reported.
