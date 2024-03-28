#include <RcppArmadillo.h>
#include "utilities.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat fn_dexpQ(int cstat, double t1, double t2, arma::vec q, arma::vec Xi, double gamma_all, int len_par) {
  
  int len_q = q.n_elem;
  arma::vec diag_e = ones(len_q);
  arma::vec valq = fn_egval(q);
  arma::mat vecq = fn_egvec(q);
  double u = fn_h(t2, gamma_all) - fn_h(t1, gamma_all);
  arma::mat exp_D = diagmat(exp(valq * u));
  arma::mat expQ = vecq * exp_D * inv(vecq);
  arma::mat dQe = join_cols(diagmat(diag_e), kron(Xi, diagmat(diag_e)));

  int nn = dQe.n_rows;
  arma::mat Q(4,4, fill::zeros);
  arma::mat dQi(4,4, fill::zeros);
  arma::mat Q_expQ(4,4, fill::zeros);
  arma::mat G(4,4, fill::zeros);
  arma::mat V(4,4, fill::zeros);
  arma::mat matt(4,4, fill::zeros);
  arma::mat res(len_par,4, fill::zeros);
  arma::vec expdiff(4);
  double dh1 = 0;
  double dh2 = 0;

  int r = 0;
  for (int i = 0; i < len_par; ++i) {

    if ((i != 6) & (r < nn)) {
      dQi.zeros();
      G.zeros();
      V.zeros();
      
      dQi(0, 0) = (- dQe(r, 0) * q[0] - dQe(r, 1) * q[1]);
      dQi(0, 1) = dQe(r, 0) * q[0];
      dQi(0, 3) = dQe(r, 1) * q[1];

      dQi(1, 0) = dQe(r, 2) * q[2];
      dQi(1, 1) = (- dQe(r, 2) * q[2] - dQe(r, 3) * q[3] - dQe(r, 4) * q[4]);
      dQi(1, 2) = dQe(r, 3) * q[3];
      dQi(1, 3) = dQe(r, 4) * q[4];

      dQi(2, 2) = (- dQe(r, 5) * q[5]);
      dQi(2, 3) = dQe(r, 5) * q[5];

      G = inv(vecq) * dQi * vecq;
      expdiff = exp(valq * u);

      for (int j = 0; j < 4; ++j) {
        for (int k = j; k < 4; ++k) {
          if (j == k) {
            V(j, j) = u * expdiff[j];
          } else {
            V(j, k) = V(k, j) = (expdiff[j] - expdiff[k]) / (valq[j] - valq[k]);
          }
        }
      }

      matt = vecq * (G % V) * inv(vecq);

      res(i, 0) = matt(cstat, 0);
      res(i, 1) = matt(cstat, 1);
      res(i, 2) = matt(cstat, 2);
      res(i, 3) = matt(cstat, 3);
      // res.row(i) = matt.rows(cstat);

      r += 1;
    }

    if (i == 6) {
      Q(0, 0) = (- q[0] - q[1]);
      Q(0, 1) = q[0];
      Q(0, 3) = q[1];

      Q(1, 0) = q[2];
      Q(1, 1) = (- q[2] - q[3] - q[4]);
      Q(1, 2) = q[3];
      Q(1, 3) = q[4];

      Q(2, 2) = (-q[5]);
      Q(2, 3) = q[5];

      Q_expQ = Q * expQ;
      if (t1 == 0) dh1 = 0; else dh1 = pow(t1, gamma_all) * log(t1);
      if (t2 == 0) dh2 = 0; else dh2 = pow(t2, gamma_all) * log(t2);

      res(i, 0) = (dh2 - dh1) * Q_expQ(cstat, 0);
      res(i, 1) = (dh2 - dh1) * Q_expQ(cstat, 1);
      res(i, 2) = (dh2 - dh1) * Q_expQ(cstat, 2);
      res(i, 3) = (dh2 - dh1) * Q_expQ(cstat, 3);
      
      // res.row(i) = (dh2 - dh1) * Q_expQ.rows(cstat);
      // 
    }
  }

  return res;

}
