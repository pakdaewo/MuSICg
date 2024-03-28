#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fn_egvec(arma::vec q) {
  
  arma::mat egvec(4, 4);
  
  egvec(0, 0) = 1;
  egvec(1, 0) = 1;
  egvec(2, 0) = 1;
  egvec(3, 0) = 1;
  
  egvec(0, 1) = (-(q[0] + q[1] - q[2] - q[3] - q[4] + sqrt(pow(q[0], 2.0) + 2*q[0]*q[1] + pow(q[1], 2.0) + 2*q[0]*q[2] - 2*q[1]*q[2] + pow(q[2], 2.0) - 2*q[0]*q[3] - 2*q[1]*q[3] + 2*q[2]*q[3] + pow(q[3], 2.0) - 2*q[0]*q[4] - 2*q[1]*q[4] + 2*q[2]*q[4] + 2*q[3]*q[4] + pow(q[4], 2.0)))/(2*q[2]));
  egvec(1, 1) = 1;
  egvec(2, 1) = 0;
  egvec(3, 1) = 0;
  
  egvec(0, 2) = (-(q[0] + q[1] - q[2] - q[3] - q[4] - sqrt(pow(q[0], 2.0) + 2*q[0]*q[1] + pow(q[1], 2.0) + 2*q[0]*q[2] - 2*q[1]*q[2] + pow(q[2], 2.0) - 2*q[0]*q[3] - 2*q[1]*q[3] + 2*q[2]*q[3] + pow(q[3], 2.0) - 2*q[0]*q[4] - 2*q[1]*q[4] + 2*q[2]*q[4] + 2*q[3]*q[4] + pow(q[4], 2.0)))/(2*q[2]));
  egvec(1, 2) = 1;
  egvec(2, 2) = 0;
  egvec(3, 2) = 0;
  
  egvec(0, 3) = (- q[0]*q[3])/(-q[1]*q[2] - q[0]*q[3] - q[1]*q[3] - q[0]*q[4] - q[1]*q[4] + q[0]*q[5] + q[1]*q[5] + q[2]*q[5] + q[3]*q[5] + q[4]*q[5] - pow(q[5], 2.0));
  egvec(1, 3) = (-(-q[0]*q[3] - q[1]*q[3] + q[3]*q[5]))/(q[1]*q[2] + q[0]*q[3] + q[1]*q[3] + q[0]*q[4] + q[1]*q[4] - q[0]*q[5] - q[1]*q[5] - q[2]*q[5] - q[3]*q[5] - q[4]*q[5] + pow(q[5], 2.0));
  egvec(2, 3) = 1;
  egvec(3, 3) = 0;
  
  return egvec;
}