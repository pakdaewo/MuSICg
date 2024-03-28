#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double fn_Sto(double t, arma::vec q, double gamma_all, arma::vec px, arma::vec pw) {
  
  int npx = px.n_elem;
  
  double num = 0;
  double den = 0;
  arma::vec expQ(4, fill::zeros);
  
  arma::vec s = t * (px + 1)/2;
  
  for (int i = 0; i < npx; ++i) {
    expQ = fn_expQ(0, 0, s[i], q, gamma_all);
    num += pw[i] * expQ[1] * q[3] * fn_g(s[i], gamma_all) * t/2;
  }
  
  arma::vec r = (1 - px)/(1 + px);
  for (int j = 0; j < npx; ++j) {
    expQ = fn_expQ(0, 0, r[j], q, gamma_all);
    den += pw[j] * expQ[1] * q[3] * fn_g(r[j], gamma_all) * 2 / pow(1 + px[j], 2.0);
  }
  
  return num/den;
  
}