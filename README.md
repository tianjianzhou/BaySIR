# BaySIR
Semiparametric Bayesian Inference for the Transmission Dynamics of COVID-19 with a State-Space Model

## Installation
The `BaySIR` package has two dependencies: `Rcpp` and `RcppArmadillo`.

The `BaySIR` package can be easily installed with the `devtools` package in R. After `devtools` has been installed, run the following commands in R to install the `RNDClone` package.
```
library(devtools)
install_github("tianjianzhou/BaySIR")
```


## Usage

#### BaySIR_MCMC
```
library(BaySIR)
MCMC_spls = BaySIR_MCMC(B, I_D_0, N)
```
- `B` is a length `T + 1` vector of daily new confirmed cases `B[0], ..., B[T]`
- `I_D_0` is the total number of confirmed cases on day 0
- `N` is the population size.
