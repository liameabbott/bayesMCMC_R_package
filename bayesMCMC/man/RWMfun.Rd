\name{RWMfun}
\alias{RWMfun}
\title{
Random Walk Metropolis - Bayes Parameter Estimation 
}
\description{
This function implements a simple Metropolis-Hastings algorithm, which we call the 'Random Walk Metropolis', to obtain sample parameter draws from the posterior density. These sample parameter draws can be used to estimate the coefficients of a linear regression model. The function is coded in R.
}
\usage{
RWMfun(nMC, y, X, stepsize=exp(-5))
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
	\item{stepsize}{
		A numeric value that indicates the size of the steps to take  when seeking a new proposal value in the Markov chain. This defaults to exp(-5).
	}
}
\details{
The RWMfun function takes the number of MCMC iterations (nMC), the response variable values, covariate values, and the magnitude of the size of the steps to take through the parameter space as input.

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

# call function
results = RWMfun(nMC=nMC, y=y, X=X, stepsize=exp(-5))
}
