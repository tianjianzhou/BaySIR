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

### BaySIR_MCMC
- Usage: `BaySIR_MCMC(B, I_D_0, N, ...)`

- Arguments
  - Required
    - `B`: a length `T + 1` vector of daily new confirmed cases `B[0], ..., B[T]`
    - `I_D_0`: the total number of confirmed cases on day 0
    - `N`: the population size
  - Optional
    - `confirmed_cases_cum`: a length `T + 2` vector. If `confirmed_cases_cum` is specified and not both `B` and `I_D_0`
    - `X`: a `(T + 1) * Q` matrix, covariates related to the disease transmission rate. Default is an intercept term plus a time trend
    - `Y`: a `(T + 1) * K` matrix, covariates related to the diagnosis rate. Default contains only an intercept term
    - `niter`
    - `burnin`
    - `thin`
    - `Delta`

- Output
  - A list of posterior samples
