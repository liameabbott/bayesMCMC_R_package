# bayesMCMC

This is an R package co-authered by David Reynolds (https://github.com/DBomber60) and myself as a joint final project for BIOSTAT 615, 'Statistical Computing in C++ and R' in the Fall 2015 term at the University of Michigan.

The package implements the Metropolis-Hastings and Gibbs sampler Markov Chain Monte Carlo (MCMC) algorithms to estimate the coefficients for a Bayesian multiple linear regression model.

Functions within the package are implemented using both R and C++ (using the Rcpp library) code.

To install the package:
* Download the installation package in this repository (bayesMCMC_1.0.tar.gz)
* In the terminal, navigate to the directory containing the installation package and type 'R CMD INSTALL bayesMCMC_1.0.tar.gz'
* Open R and type 'library(bayesMCMC)' to make sure the package was installed correctly
* After loading the package in an R session, type '? bayesMCMC' to view an overview of the functionality of the package, including short descriptions of the functions included in the package.
