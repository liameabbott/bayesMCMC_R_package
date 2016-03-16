\name{bayesMCMC-package}
\alias{bayesMCMC-package}
\alias{bayesMCMC}
\docType{package}
\title{
 Fitting a Bayesian Linear Regression Model
}
\description{
This package implements the Metropolis-Hastings and Gibbs sampler   Markov Chain Monte Carlo (MCMC) algorithms to estimate coefficients for a Bayesian multiple linear regression model.
}
\details{
The gibbsSampler function takes the number of MCMC iterations (nMC), the response variable values, covariate values, and values for prior parameters a, b, and kappa as input. All response and covariate values must be numeric.

The function returns a matrix with nMC rows and p+1 columns (where there are p columns in the covariate matrix). Each row is a vector of parameter estimates produced by a single iteration of the algorithm. The first p columns contain estimates of the beta coefficients, and the p+1 column contains an estimate of the sigma-squared parameter.

The RWMfun function takes the number of MCMC iterations (nMC), the response variable values, covariate values, and the magnitude of the size of the steps to take through the parameter space as input.

The function returns a matrix with nMC rows and p+1 columns (where there are p columns in the covariate matrix). Each row is a vector of parameter estimates produced by a single iteration of the algorithm. The first p columns contain estimates of the beta coefficients, and the p+1 column contains an estimate of the sigma-squared parameter.

Also included in the package are versions of each function written in C++ (gibbsSamplerC and RWMfunC). Those versions should be faster than the versions written natively in R.
}
\author{
Liam Abbott, David Reynolds

Maintainer: Liam Abbott <test@test.com>
}

