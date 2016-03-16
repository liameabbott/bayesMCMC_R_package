#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double lpi(NumericVector y, NumericVector beta, double sigmasq, NumericMatrix X) {
  double val = -(1+0.5*(y.size()+1))*log(sigmasq);
  // cross product of beta
  double crossprod_beta = 0;
  for (int i=0; i<beta.size(); i++) {
    crossprod_beta += beta[i]*beta[i];}
  // cross product of y-XB
  NumericVector mu(y.size());
  for (int i=0; i<y.size(); i++) {
    mu[i] = 0;
    for (int j=0; j<beta.size(); j++) {
      mu[i] += X(i,j)*beta[j];}}
  NumericVector diff(y.size());
  for (int i=0; i<y.size(); i++) {
    diff[i] = y[i]-mu[i];}
  double crossprod_d = 0;
  for (int i=0; i<mu.size(); i++) {
    crossprod_d += diff[i]*diff[i];}
  val += -0.5*0.001*crossprod_beta-0.5*(3/sigmasq);
  val -= (0.5/sigmasq)*crossprod_d;
  return val; }

// [[Rcpp::export]]
NumericVector upd(int p, NumericVector oldbeta, NumericMatrix cholesky, NumericVector rnorms) {
  NumericVector result(p);
  for (int r=0; r<p; r++) {
    result[r] += oldbeta[r];
    for (int c=0; c<p; c++) {
      result[r] += cholesky(r,c)*rnorms[c];
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix RWMfunC(int nMC, NumericVector y, NumericMatrix X, double stepsize = exp(-5)) {
  // beta_prop
  int p = X.ncol();
  NumericVector beta(p);
  NumericMatrix C(p,p);
  for (int r=0; r<p; r++) {
    for (int c=0; c<p; c++) {
      if (r==c) C(r,c)=stepsize;
    }
  }
  double sigmasq = 1;
  NumericMatrix Res2(nMC,p+1);
  double log_lpi = lpi(y, beta, sigmasq, X);
  
  //NumericVector betaprop=beta;
  for (int k=0; k<nMC; k++) {
    double sigmasq_prop = 0;
    sigmasq_prop += sigmasq + stepsize*rnorm(1,0,1)[0];
    if (sigmasq_prop>0) {
      NumericVector betaprop(p);
      betaprop = upd(p, beta, C, rnorm(p,0,1));
      double new_lpi = lpi(y, betaprop, sigmasq_prop, X);
      double compare = exp(new_lpi-log_lpi);
      if (compare >1) compare = 1;
      double u = runif(1,0,1)[0];
      if (u < compare) {
        beta = betaprop;
        log_lpi = new_lpi;
        sigmasq = sigmasq_prop;
      }
    }
    Res2(k,p) = sigmasq;
    for (int c=0; c<p; c++) {
      Res2(k,c) = beta[c];
    }
  }
  return Res2;
  }


/*** R

*/
