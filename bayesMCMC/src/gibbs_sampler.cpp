// Liam Abbott, leabbott
// BIOSTAT 615

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd gibbsSamplerC(int nMC,
                       VectorXd y,
                       MatrixXd X,
                       double a = 1.1,
                       double b = 1.1,
                       double kappa = 0.1) {

    int n = X.rows();
    int p = X.cols();
    
    VectorXd beta = VectorXd::Zero(p);
    MatrixXd Res(nMC, p+1);
    
    // take Cholesky decomposition of M matrix
    MatrixXd M = ( kappa * MatrixXd::Identity(p, p) ) + ( X.transpose() * X );
    LLT<MatrixXd> lltOfM(M);
    MatrixXd L = lltOfM.matrixL();
    
    // to be used in draws of sigmasq estimates
    double A = a + 0.5*(n+p);
    
    // to be used in draws of beta estimates
    VectorXd beta_means = M.inverse() * (X.transpose() * y);
    
    // run Markov Chain
    for ( int i=0; i<nMC; ++i ) {
        
        // UPDATE SIGMASQ ESTIMATE
        
        // update scale parameter
        VectorXd mat = y - (X * beta);
        double beta_val = (beta.transpose() * beta).sum();
        double mat_val = (mat.transpose() * mat).sum();
        double B = ( b + 0.5 * kappa * beta_val ) + ( 0.5 * mat_val );
        
        // new sigmasq inverse gamma draw
        double gamma_draw = rgamma(1, A, 1/B)[0];
        double new_sigmasq = 1.0/gamma_draw;
        
        
        // UPDATE BETA ESTIMATES
        
        for ( int c=0; c<p; ++c ) {
            beta(c) = beta_means(c) + sqrt(new_sigmasq) * L.inverse().sum() * rnorm(1, 0, 1)[0];
        }
        
        
        // STORE RESULTS
        
        Res(i, p) = new_sigmasq;
        Res.row(i).leftCols(p) = beta;
        
    }
    
    return Res;
    
}
