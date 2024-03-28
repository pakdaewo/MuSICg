#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec fn_fto(arma::vec t, arma::vec q, double gamma_all, arma::vec px, arma::vec pw) {
  
  int nt = t.n_elem;
  int npx = px.n_elem;
  
  arma::vec num(nt);
  double den = 0;
  arma::vec expQ(4, fill::zeros);
  
  for (int i = 0; i < nt; ++i) {
    expQ = fn_expQ(0, 0, t[i], q, gamma_all);
    num[i] = expQ[1] * q[3] * fn_g(t[i], gamma_all);
  }
  
  arma::vec r = (1 - px)/(1 + px);
  for (int j = 0; j < npx; ++j) {
    expQ = fn_expQ(0, 0, r[j], q, gamma_all);
    den += pw[j] * expQ[1] * q[3] * fn_g(r[j], gamma_all) * 2 / pow(1 + px[j], 2.0);
  }
  
  return num/den;
  
}