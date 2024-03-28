#include <RcppArmadillo.h>
#include "utilities.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec fn_expQ(int cstat, double t1, double t2, arma::vec q, double gamma_all) {

  arma::vec valq = fn_egval(q);
  arma::mat vecq = fn_egvec(q);
  double u = fn_h(t2, gamma_all) - fn_h(t1, gamma_all);
  arma::mat exp_D = diagmat(exp(valq * u));
  arma::mat expQ = vecq * exp_D * inv(vecq);

  arma::vec res(4);
  
  arma::vec test(4, fill::zeros);
  
  res[0] = expQ(cstat, 0);
  res[1] = expQ(cstat, 1);
  res[2] = expQ(cstat, 2);
  res[3] = expQ(cstat, 3);
  
  // for (int jj = 0; jj < 4; ++jj) {
  //   for (int kk = 0; kk < 4; ++kk) {
  //     test[jj] += expQ(jj, kk);
  //   }
  //   if (test[jj] > 1.000001) res[jj] = NA_REAL;
  // }

  
  return res;

}
