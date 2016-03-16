\name{gibbsSampler}
\alias{gibbsSampler}
\title{
Gibbs Sampler - Bayes Parameter Estimation 
}
\description{
This function uses the Gibbs sampler algorithm to obtain sample parameter draws from the posterior density, which can be used to estimate the coefficients of a linear regression model. The function is coded in R.
}
\usage{
gibbsSampler(nMC, y, X, a=1.1, b=1.1, kappa=0.1)
}
\arguments{
	\item{nMC}{
		An integer that tells the function how many MCMC iterations to run.
	}
	\item{y}{
		A numeric vector containing the response variable value for each observation.
	}
	\item{X}{
		A numeric matrix containing covariate values for each observation. 
	}
	\item{a}{
		A numeric value that is an estimate of a prior distribution parameter. The default value is 1.1.
	}
	\item{b}{
		A numeric value that is an estimate of a prior distribution parameter. The default value is 1.1.
	}
	\item{kappa}{
		A numeric value is an estimate of a prior distribution parameter. The default value is 0.1.
	}
}
\details{
The gibbsSampler function takes the number of MCMC iterations (nMC), the response variable values, covariate values, and values for prior parameters a, b, and kappa as input. All response and covariate values must be numeric.

The function returns a matrix with nMC rows and p+1 columns (where there are p columns in the covariate matrix). Each row is a vector of parameter estimates produced by a single iteration of the algorithm. The first p columns contain estimates of the beta coefficients, and the p+1 column contains an estimate of the sigma-squared parameter.
}
\author{
Liam Abbott, David Reynolds
}
\examples{
 # simulate data
set.seed(2)
n = 10000
beta_star = c(1, -12, 50, 22)
p = length(beta_star)
sigmasq_true = 0.1
X = cbind(rep(1, n), matrix(rnorm((p-1)*n), ncol=(p-1)))
y = X \%*\% beta_star + sqrt(sigmasq_true)*rnorm(n, 1)

# set length of Markov chain
nMC = 100000

# define parameters of the prior distribution
a = 1.1
b = 1.1
kappa = 0.1

# call function
results = gibbsSampler(nMC=nMC, y=y, X=X, a=a, b=b, kappa=kappa)
}
