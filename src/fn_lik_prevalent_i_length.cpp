#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"
#include "utilities3.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double fn_lik_prevalent_i_length(int id, arma::vec u_ID, arma::vec v_ID, arma::mat Qe, arma::vec cen, double gamma_all, 
                                 double tau, arma::vec px, arma::vec pw, arma::vec px21, arma::vec px22, arma::vec pw2) {
  
  double ui = u_ID[id-1];
  double vi = v_ID[id-1];
  
  arma::vec q(6);
  for (int k1 = 0; k1 < 6; ++k1) {
    q[k1] = Qe(id-1, k1);
  }
  
  double lam23 = q[5];
  int ceni = cen[id-1];
  arma::vec s = ui * (px + 1)/2;
  int ns = s.n_elem;
  
  arma::vec fto = fn_fto(s, q, gamma_all, px, pw);
  
  double St0to = 0;
  double ft0to = 0;
  double num = 0;
  double fs = 0;
  
  for (int i = 0; i < ns; ++i) {
    St0to = exp(-lam23 * pow(ui + vi, gamma_all) + lam23 * pow(s[i], gamma_all));
    ft0to = lam23 * gamma_all * pow(ui + vi, gamma_all - 1) * St0to;
    
    if (ceni == 1) fs = ft0to;
    if (ceni == 0) fs = St0to;
    
    num += pw[i] * fs * 1/tau * fto[i] * ui/2;
  }
  
  arma::vec as1 = tau/2 * (1 - px21);
  arma::vec as2 = (1 - px22)/(1 + px22);
  
  int nas = as1.n_elem;
  double St0to_den = 0;
  double ga_den = 0;
  double den = 0;
  
  arma::vec fto_den = fn_fto(as2, q, gamma_all, px, pw);
  
  for (int j = 0; j < nas; ++j) {
    
    St0to_den = exp(-lam23 * pow(as1[j] + as2[j], gamma_all) + lam23 * pow(as2[j], gamma_all));
    den += pw2[j] * St0to_den * 1/tau * fto_den[j] * (tau/2) * 2 / pow(px22[j] + 1, 2);
  }
  
  if(den > 1.0) den = 1;
  
  return log(num) - log(den);
}
