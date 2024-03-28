#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"
#include "utilities3.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec fun_lik_incident_grad_i(arma::vec ids, arma::vec ID, arma::vec visits, arma::vec states, arma::mat Qe, arma::mat X_ID, arma::vec cen, double gamma_all, int len_par) {
  
  arma::vec visit_i = visits.elem(find(ID == ids[0]));
  arma::vec state_i = states.elem(find(ID == ids[0]));
  arma::vec cens_i = cen.elem(find(ID == ids[0]));
  
  arma::vec Q_i(6);
  for (int k1 = 0; k1 < 6; ++k1) {
    Q_i[k1] = Qe(ids[1]-1, k1);
  }
  
  int nx = X_ID.n_cols;
  arma::vec X_i(nx);
  for (int k2 = 0; k2 < nx; ++k2) {
    X_i[k2] = X_ID(ids[1]-1, k2);
  }
  
  double cen_i = cens_i[0];
  double nv = visit_i.n_elem;
  
  arma::mat times_i(nv-1, 2, fill::zeros);
  arma::mat trans_i(nv-1, 2, fill::zeros);
  for (int i = 0; i < (nv - 1); ++i) {
    times_i(i, 0) = visit_i[i];
    times_i(i, 1) = visit_i[i+1];
    trans_i(i, 0) = state_i[i];
    trans_i(i, 1) = state_i[i+1];
  }
  
  arma::vec sum_log_tranp_i(len_par, fill::zeros);
  double before_time;
  double after_time;
  int before_state;
  int after_state;
  arma::vec expQ;
  arma::mat dexpQ(len_par, 4);
  double gt_die;
  arma::vec gt_grad_die(len_par);
  arma::mat dQ(len_par, 3);
  
  for (int j = 0; j < (nv-1); ++j) {

    before_time = times_i(j, 0);
    after_time = times_i(j, 1);
    before_state = trans_i(j, 0);
    after_state = trans_i(j, 1);

    gt_die = fn_g(after_time, gamma_all);
    gt_grad_die = fn_dg(after_time, gamma_all, len_par);
    dQ = fun_dQ(Q_i, X_i, gamma_all, len_par);
    expQ = fn_expQ(before_state, before_time, after_time, Q_i, gamma_all);
    dexpQ = fn_dexpQ(before_state, before_time, after_time, Q_i, X_i, gamma_all, len_par);
    double den;
    double num1;
    double num2;
    double num3;
    
    for (int j2 = 0; j2 < len_par; ++j2) {
    if (after_state != 3) {
      sum_log_tranp_i[j2] += dexpQ(j2, after_state)/expQ[after_state];
    } else {
      if (cen_i == 1) {
        den = gt_die * (expQ[0] * Q_i[1] + expQ[1] * Q_i[4] + expQ[2] * Q_i[5]);
        num1 = (dexpQ(j2, 0) * Q_i[1] * gt_die + expQ[0] * dQ(j2, 0) * gt_die + expQ[0] * Q_i[1] * gt_grad_die[j2]);
        num2 = (dexpQ(j2, 1) * Q_i[4] * gt_die + expQ[1] * dQ(j2, 1) * gt_die + expQ[1] * Q_i[4] * gt_grad_die[j2]);
        num3 = (dexpQ(j2, 2) * Q_i[5] * gt_die + expQ[2] * dQ(j2, 2) * gt_die + expQ[2] * Q_i[5] * gt_grad_die[j2]);
        sum_log_tranp_i[j2] += (num1 + num2 + num3)/den;
      } else {
          sum_log_tranp_i[j2] += (dexpQ(j2, 0) + dexpQ(j2, 1) + dexpQ(j2, 2))/(expQ[0] + expQ[1] + expQ[2]);
        }
      }
    }
    
  }
    
  return sum_log_tranp_i;
  
}

  