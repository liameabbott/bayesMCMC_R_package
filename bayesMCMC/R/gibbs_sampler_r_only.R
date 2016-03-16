# Gibbs sampler implemented in R

gibbsSampler = function(nMC, y, X, a=1.1, b=1.1, kappa=0.1) {
  
  # define variables
  n = dim(X)[1]
  p = dim(X)[2]
  
  # create matrix to store posterior density draws
  Res = matrix(NA, ncol=p+1, nrow=nMC)
  
  # initialize the sampler
  beta = as.matrix(rep(0, p))
  M = kappa*diag(p) + t(X) %*% X
  invM = solve(M)
  R = chol(M)
  
  # run the sampler
  for (jj in 1:nMC)
  {
    
    # update sigma
    A = a + 0.5*(n+p)
    B = b + 0.5*kappa*crossprod(beta, beta) + 0.5*crossprod(y - X %*% beta, y - X %*% beta)
    sigmasq = 1/rgamma(1, shape=A, rate=B)
    
    # update betas
    beta = invM %*% (t(X) %*% y) + sqrt(sigmasq)*solve(R, rnorm(p))
    
    # record results
    Res[jj,] = c(as.vector(beta), sigmasq)
  }
  
  # return results
  return (Res)
  
}


