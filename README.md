# BaySIR
Semiparametric Bayesian Inference for the Transmission Dynamics of COVID-19 with a State-Space Model

## Installation
The `BaySIR` package has two dependencies: `Rcpp` and `RcppArmadillo`.

The `BaySIR` package can be easily installed with the `devtools` package in R. After `devtools` has been installed, run the following commands in R to install the `BaySIR` package.
```
library(devtools)
install_github("tianjianzhou/BaySIR")
```


## Functions

### 1. BaySIR_MCMC
- Usage: `BaySIR_MCMC(B, I_D_0, N, ...)`

- Arguments
  - Required
    - `B`: a length `T + 1` vector of daily new confirmed cases `B[0], ..., B[T]`.
    - `I_D_0`: the total number of confirmed cases on day 0.
    - `N`: the population size.
  - Optional
    - `confirmed_cases_cum`: an optional length `T + 2` vector of cumulative confirmed case counts. If `B` and `I_D_0` have already been specified, then `confirmed_cases_cum` will be ignored. If not both `B` and `I_D_0` are specified and `confirmed_cases_cum` is supplied, then `confirmed_cases_cum` will be used to calculate `B` and `I_D_0`.
    - `X`: a `(T + 1) * Q` matrix, covariates related to the disease transmission rate. Default is an intercept term plus a time trend, `X[t, ] = (1, t)`.
    - `Y`: a `(T + 1) * K` matrix, covariates related to the diagnosis rate. Default contains only an intercept term, `Y[t, ] = 1`.
    - `niter`: number of MCMC samples to return. Default is `1000`. The total number of MCMC iterations to be run is `niter * thin + burnin`.
    - `burnin`: number of MCMC iterations that will be discarded as initial burn-in. Default is `10000`.
    - `thin`: meaning keep 1 draw every `thin` MCMC iterations. Default is `20`.
    - `Delta`: a monotonically increasing vector (each element is larger than the previous) defining the temperatures of the parallel Markov chains (parallel tempering). The first element must be 1, corresponding to the original posterior. Default is `1.5^(0:9)`.

- Output
  - A list of the following:
    - `MCMC_spls`: a list of the MCMC samples for the parameters.
    - `MCMC_summary`: a list of the posterior summaries for the parameters. For each parameter, its posterior median, 2.5% quantile and 97.5% quantile are reported.


- Examples
  ```
  library(BaySIR)
  
  data(data_sim_1)
  B = data_sim_1$B
  I_D_0 = data_sim_1$I_D[1]
  N = data_sim_1$N
  
  result_list = BaySIR_MCMC(B = B, I_D_0 = I_D_0, N = N)
  result_list$MCMC_summary$R_eff
  ```
