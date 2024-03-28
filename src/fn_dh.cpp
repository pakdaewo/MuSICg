#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec fn_dh(double t, double gamma_all, int len_par) {
  
  arma::vec res(len_par, fill::zeros);
  
  if (t == 0) {
    res[6] = 0;
  } else {
    res[6] = pow(t, gamma_all) * log(t);
  }
  
  return res;
  
} 
