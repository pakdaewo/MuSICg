#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"
#include "utilities3.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double fn_lik_prevalent_exp(arma::vec u_ID, arma::vec v_ID, arma::mat Qe, arma::vec cen, double gamma_all, 
                          double eta, arma::vec px, arma::vec pw, arma::vec px21, arma::vec px22, arma::vec pw2) {
  
  int nID = u_ID.n_elem;
  double ui;
  double vi;
  arma::vec q(6);
  double lam23;
  int ceni;
  int ns = px.n_elem;
  int nas = px21.n_elem;
  
  arma::vec s(ns);
  arma::vec as1(nas);
  arma::vec as2(nas);
  arma::vec fto(ns);
  arma::vec fto_den(nas);
    
  double St0to = 0;
  double ft0to = 0;
  double ga = 0;
  double num = 0;
  double fs = 0;
  double lik = 0;
  double St0to_den = 0;
  double ga_den = 0;
  double den = 0;
  
  for (int r = 0; r < nID; ++r) {
    
    num = 0;
    den = 0;
    
    ui = u_ID[r];
    vi = v_ID[r];
    
    for (int k1 = 0; k1 < 6; ++k1) {
      q[k1] = Qe(r, k1);
    }
    
    lam23 = q[5];
    
    ceni = cen[r];
    s = ui * (px + 1)/2;

    fto = fn_fto(s, q, gamma_all, px, pw);
    
    for (int i = 0; i < ns; ++i) {
      St0to = exp(-lam23 * pow(ui + vi, gamma_all) + lam23 * pow(s[i], gamma_all));
      ft0to = lam23 * gamma_all * pow(ui + vi, gamma_all - 1) * St0to;
      ga = eta *  exp(-eta * (ui - s[i]));
      
      if (ceni == 1) fs = ft0to;
      if (ceni == 0) fs = St0to;
      
      num += pw[i] * fs * ga * fto[i] * ui/2;
    }
    
    as1 = (1 - px21)/(1 + px21);
    as2 = (1 - px22)/(1 + px22);
    
    fto_den = fn_fto(as2, q, gamma_all, px, pw);
    
    for (int j = 0; j < nas; ++j) {
      
      St0to_den = exp(-lam23 * pow(as1[j] + as2[j], gamma_all) + lam23 * pow(as2[j], gamma_all));
      ga_den = eta * exp(-eta * as1[j]);
      den += pw2[j] * St0to_den * ga_den * fto_den[j] * 2 / pow(px21[j] + 1, 2) * 2 / pow(px22[j] + 1, 2);
    }
    
    if(den > 1.0) den = 1;
    
    lik += log(num) - log(den);
    
  }
  
  return lik;
  
}
