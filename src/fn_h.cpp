#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double fn_h(double t, double gamma_all) {
  return pow(t, gamma_all);
} 