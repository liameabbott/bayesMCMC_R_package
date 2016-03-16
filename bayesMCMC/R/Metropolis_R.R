

RWMfun = function(nMC,y,X,stepsize=exp(-5)){
  
  ## log-density function
  lpi = function(beta,sigmasq,y,X) {
      # prior
      tau_0=1e-3; tau_1=1e-3; sigmasq_0=3;
      val=-(1+.5*(length(y)+1))*log(sigmasq);
      val = val-0.5*tau_0*(crossprod(beta))-0.5*(sigmasq_0/sigmasq);
      val = val-(.5/sigmasq)*crossprod(y-X%*%beta);
      return(val);
  }
  
  #### Metropolis algorithm with proposal N(0,exp(lsig))
  p = ncol(X)
  prop_mat = diag(p)*stepsize;
  beta=numeric(p)
  Res=matrix(NA,nrow=nMC,ncol=p+1)
  beta = numeric(p); sigmasq = 1
  log_pi=lpi(beta,sigmasq,y,X)
  for (jj in 1:nMC){
    beta_prop =beta + prop_mat%*%rnorm(p,0,1)
    sigmasq_prop = sigmasq + stepsize*rnorm(1);
    if (sigmasq_prop <=0){
      Acc=0 }else{
        log_pi_prop=lpi(beta_prop, sigmasq_prop, y,X);
        Acc=min(1,exp(log_pi_prop-log_pi))
        if (runif(1)<=Acc){
          beta=beta_prop;
          sigmasq = sigmasq_prop
          log_pi=log_pi_prop;
        }
      }
    Res[jj,]=c(beta,sigmasq)
  }
  return(Res)
}