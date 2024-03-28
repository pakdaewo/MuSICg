#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"
#include "utilities3.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double fun_lik_incident(arma::vec uID, arma::vec ID, arma::vec visits, arma::vec states, arma::mat Qe, arma::mat X_ID, arma::vec cen, double gamma_all) {
  
  int nID = uID.n_elem;
  arma::vec Q_i(6);
  double cen_i;
  int nv = 0;
  
  double sum_log_tranp_i = 0;
  double before_time;
  double after_time;
  int before_state;
  int after_state;
  arma::vec expQ(4);
  double gt_die;
  
  for (int r = 0; r < nID; ++r) {
    
    arma::vec visit_i = visits.elem(find(ID == uID[r]));
    arma::vec state_i = states.elem(find(ID == uID[r]));
    arma::vec cens_i = cen.elem(find(ID == uID[r]));
    
    for (int k1 = 0; k1 < 6; ++k1) {
      Q_i[k1] = Qe(r, k1);
    }
    
    cen_i = cens_i[0];
    nv = visit_i.n_elem;
    
    arma::mat times_i(nv-1, 2, fill::zeros);
    arma::mat trans_i(nv-1, 2, fill::zeros);
    
    for (int i = 0; i < (nv - 1); ++i) {
      times_i(i, 0) = visit_i[i];
      times_i(i, 1) = visit_i[i+1];
      trans_i(i, 0) = state_i[i];
      trans_i(i, 1) = state_i[i+1];
    }
    
    for (int j = 0; j < (nv-1); ++j) {
      before_time = times_i(j, 0);
      after_time = times_i(j, 1);
      before_state = trans_i(j, 0);
      after_state = trans_i(j, 1);
      
      expQ = fn_expQ(before_state, before_time, after_time, Q_i, gamma_all);
      
      if (after_state != 3) {
        sum_log_tranp_i += log(expQ[after_state]);
      } else {
        if (cen_i == 1) {
          gt_die = fn_g(after_time, gamma_all);
          sum_log_tranp_i += log(gt_die * (expQ[0] * Q_i[1] + expQ[1] * Q_i[4] + expQ[2] * Q_i[5]));
        } else {
          sum_log_tranp_i += log(expQ[0] + expQ[1] + expQ[2]);
        }
      }
      
    }
    
  }
  
  return sum_log_tranp_i;
  
}

