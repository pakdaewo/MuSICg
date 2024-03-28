#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double fn_g(double t, double gamma_all) {
  return gamma_all * pow(t, gamma_all - 1.0);
}