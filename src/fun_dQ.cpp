#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fun_dQ(arma::vec q, arma::vec xi, double gamma_all, int len_par) {
  
  int len_q = q.n_elem;
  int nx = xi.n_elem;
  int nn = len_q + 1 + nx * len_q;
  
  arma::vec diag_e = ones(len_q);
  arma::mat dQe = join_cols(diagmat(diag_e), kron(xi, diagmat(diag_e)));
  int k = 0;
  arma::mat res(len_par, 3, fill::zeros);
  for (int i = 0; i < len_par; ++i) {
    
    if ((i != 6) & (k < nn)) {
      res(i, 0) = dQe(k, 1) * q[1];
      res(i, 1) = dQe(k, 4) * q[4];
      res(i, 2) = dQe(k, 5) * q[5];
      k += 1;
    }
    
    if (i == 6) {
      res(i, 0) = 0;
      res(i, 2) = 0;
      res(i, 1) = 0;
    }
    
  }

  return res;

} 
