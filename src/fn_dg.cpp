#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec fn_dg(double t, double gamma_all, int len_par) {
  
  arma::vec res(len_par, fill::zeros);
  
  res[6] = pow(t, gamma_all - 1) * (1 + gamma_all * log(t));
    
  return res;
  
} 
