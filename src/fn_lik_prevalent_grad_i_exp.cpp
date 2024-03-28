#include <RcppArmadillo.h>
#include "utilities.hpp"
#include "utilities2.hpp"
#include "utilities3.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec fn_lik_prevalent_grad_i_exp(int id, arma::vec u_ID, arma::vec v_ID, arma::mat X_ID, arma::mat Qe, arma::vec cen, double gamma_all, 
                           double eta, arma::vec px, arma::vec pw, arma::vec px21, arma::vec px22, arma::vec pw2, int len_par) {
  
  double ui = u_ID[id-1];
  double vi = v_ID[id-1];
  
  arma::vec q(6);
  for (int k1 = 0; k1 < 6; ++k1) {
    q[k1] = Qe(id-1, k1);
  }

  int nx = X_ID.n_cols;
  arma::vec xi(nx);
  for (int k2 = 0; k2 < nx; ++k2) {
    xi[k2] = X_ID(id-1, k2);
  }
  
  double ceni = cen[id-1];
  
  arma::vec s = ui * (px + 1)/2;
  int ns = s.n_elem;
  arma::vec r = (1 - px)/(1 + px);
  int nr = r.n_elem;
  
  // denominator of fto
  double fto_den = 0;
  double h12 = 0;
  double P01 = 0;
  arma::vec fto_dden(len_par, fill::zeros);
  arma::vec expQ(4, fill::zeros);
  arma::mat dexpQ(len_par, 4, fill::zeros); // here?
  arma::vec dP01(len_par, fill::zeros);
  arma::vec dh12(len_par, fill::zeros);
  
  for (int i = 0; i < nr; ++i) {
    
    expQ = fn_expQ(0, 0, r[i], q, gamma_all);
    P01 = expQ[1];
    dexpQ = fn_dexpQ(0, 0, r[i], q, xi, gamma_all, len_par);
    for (int i2 = 0; i2 < len_par; ++i2) {
      dP01[i2] = dexpQ(i2, 1);
    }
    
    h12 = q[3] * fn_g(r[i], gamma_all);
    // h12 = q[3] * gamma_all * pow(r[i], gamma_all - 1);
    dh12[3] = h12;
    dh12[6] = h12 * (1/gamma_all + log(r[i]));
    
    for (int i3 = 0; i3 < nx; ++i3) {
      dh12[7 + 6 * i3 + 3] = h12 * xi[i3];
    }
    
    fto_den += pw[i] * P01 * h12 * 2 / pow(1 + px[i], 2);
    fto_dden += pw[i] * (dP01 * h12 + P01 * dh12) * 2 / pow(1 + px[i], 2);
    
  }
  
  // num_part
  double si = 0;
  double num_part = 0;
  arma::vec dnum_part(len_par, fill::zeros);
  
  double h23 = 0;
  double Sc23 = 0;
  double fto_num = 0;
  double fto = 0;
  double ft0_to = 0;
  double St0_to = 0;
  double ga_h = 0;
  double ga_S = 0;
  double ga = 0;
  double num = 0;
  arma::vec dh23(len_par, fill::zeros);
  arma::vec dSc23(len_par, fill::zeros);
  arma::vec fto_dnum(len_par, fill::zeros);
  arma::vec dfto(len_par, fill::zeros);
  arma::vec dft0_to(len_par, fill::zeros);
  arma::vec dSt0_to(len_par, fill::zeros);
  arma::vec dga_h(len_par, fill::zeros);
  arma::vec dga_S(len_par, fill::zeros);
  arma::vec dga(len_par, fill::zeros);
  arma::vec dnum(len_par, fill::zeros);
  
  for (int j = 0; j < ns; ++j) {
    
    si = s[j];
    expQ = fn_expQ(0, 0, si, q, gamma_all);
    P01 = expQ[1];
    dexpQ = fn_dexpQ(0, 0, si, q, xi, gamma_all, len_par);
    for (int i2 = 0; i2 < len_par; ++i2) {
      dP01[i2] = dexpQ(i2, 1);
    }
    
    h12 = q[3] * fn_g(si, gamma_all);
    dh12[3] = h12;
    dh12[6] = h12 * (1/gamma_all + log(si));
    for (int i3 = 0; i3 < nx; ++i3) {
      dh12[7 + 6 * i3 + 3] = h12 * xi[i3];
    }
    
    fto_num = P01 * h12;
    fto_dnum = dP01 * h12 + P01 * dh12;
    
    fto = fto_num/fto_den;
    dfto = fto_dnum/fto_den - fto_num * fto_dden / pow(fto_den, 2);
    
    if(ceni == 1) {
      
      h23 = q[5] * fn_g(ui + vi, gamma_all);
      dh23[5] = h23;
      dh23[6] = h23 * (1/gamma_all + log(ui + vi));
      for (int i4 = 0; i4 < nx; ++i4) {
        dh23[7 + 6 * i4 + 5] = h23 * xi[i4];
      }
      
      Sc23 = exp(-q[5] * pow(ui + vi, gamma_all) + q[5] * pow(si, gamma_all));
      dSc23[5] = Sc23 * (-q[5] * pow(ui + vi, gamma_all) + q[5] * pow(si, gamma_all));
      dSc23[6] = Sc23 * (- q[5] * pow(ui + vi, gamma_all) * log(ui + vi) + q[5] * pow(si, gamma_all) * log(si));
      for (int i4 = 0; i4 < nx; ++i4) {
        dSc23[7 + 6 * i4 + 5] = Sc23 * (-q[5] * pow(ui + vi, gamma_all) + q[5] * pow(si, gamma_all)) * xi[i4];
      }
      
      ft0_to = h23 * Sc23;
      dft0_to = dh23 * Sc23 + h23 * dSc23;
      
      ga_h = eta;
      dga_h[len_par - 1] = 1;
      
      ga_S = exp(- eta * (ui - si));
      dga_S[len_par - 1] = -ga_S * (ui - si);
      
      ga = ga_h * ga_S;
      dga = dga_h * ga_S + ga_h * dga_S;
        
      num = ft0_to * ga;
      dnum = dft0_to * ga + ft0_to * dga;
      
    } else {
      Sc23 = exp(-q[5] * pow(ui + vi, gamma_all) + q[5] * pow(si, gamma_all));
      dSc23[5] = Sc23 * (-q[5] * pow(ui + vi, gamma_all) + q[5] * pow(si, gamma_all));
      dSc23[6] = Sc23 * (- q[5] * pow(ui + vi, gamma_all) * log(ui + vi) + q[5] * pow(si, gamma_all) * log(si));
      for (int i4 = 0; i4 < nx; ++i4) {
        dSc23[7 + 6 * i4 + 5] = Sc23 * (-q[5] * pow(ui + vi, gamma_all) + q[5] * pow(si, gamma_all)) * xi[i4];
      }
      
      St0_to = Sc23;
      dSt0_to = dSc23;
      
      ga_h = eta;
      dga_h[len_par - 1] = 1;
      
      ga_S = exp(- eta * (ui - si));
      dga_S[len_par - 1] = -ga_S * (ui - si);
      
      ga = ga_h * ga_S;
      dga = dga_h * ga_S + ga_h * dga_S;
      
      num = St0_to * ga;
      dnum = dSt0_to * ga + St0_to * dga;
      
    }
    
    num_part += pw[j] * num * fto * ui/2;
    
    for (int i5 = 0; i5 < len_par; ++i5) {
      dnum_part[i5] += pw[j] * (dnum[i5] * fto + num * dfto[i5]) * ui/2;
    }

  }

  // den_part
  arma::vec as1 = (1 - px21)/(1 + px21);
  arma::vec as2 = (1 - px22)/(1 + px22);
  int nas = as1.n_elem;

  double den_part = 0;
  arma::vec dden_part(len_par, fill::zeros);

  for (int jj = 0; jj < nas; ++jj) {

    expQ = fn_expQ(0, 0, as2[jj], q, gamma_all);
    P01 = expQ[1];
    dexpQ = fn_dexpQ(0, 0, as2[jj], q, xi, gamma_all, len_par);
    for (int i2 = 0; i2 < len_par; ++i2) {
      dP01[i2] = dexpQ(i2, 1);
    }

    h12 = q[3] * fn_g(as2[jj], gamma_all);
    dh12[3] = h12;
    dh12[6] = h12 * (1/gamma_all + log(as2[jj]));
    for (int i3 = 0; i3 < nx; ++i3) {
      dh12[7 + 6 * i3 + 3] = h12 * xi[i3];
    }

    fto_num = P01 * h12;
    fto_dnum = dP01 * h12 + P01 * dh12;

    fto = fto_num/fto_den;
    dfto = fto_dnum/fto_den - fto_num * fto_dden / pow(fto_den, 2);

    Sc23 = exp(-q[5] * pow(as1[jj] + as2[jj], gamma_all) + q[5] * pow(as2[jj], gamma_all));
    dSc23[5] = Sc23 * (- q[5] * pow(as1[jj] + as2[jj], gamma_all) + q[5] * pow(as2[jj], gamma_all));
    dSc23[6] = Sc23 * (- q[5] * pow(as1[jj] + as2[jj], gamma_all) * log(as1[jj] + as2[jj]) + q[5] * pow(as2[jj], gamma_all) * log(as2[jj]));
    for (int i4 = 0; i4 < nx; ++i4) {
      dSc23[7 + 6 * i4 + 5] = Sc23 * (-q[5] * pow(as1[jj] + as2[jj], gamma_all) + q[5] * pow(as2[jj], gamma_all)) * xi[i4];
    }

    St0_to = Sc23;
    dSt0_to = dSc23;

    ga_h = eta;
    dga_h[len_par - 1] = 1;

    ga_S = exp(- eta * as1[jj]);
    dga_S[len_par - 1] = -ga_S * as1[jj];

    ga = ga_h * ga_S;
    dga = dga_h * ga_S + ga_h * dga_S;

    num = St0_to * ga;
    dnum = dSt0_to * ga + St0_to * dga;

    den_part += pw2[jj] * num * fto * 2/pow(px21[jj] + 1, 2) * 2/pow(px22[jj] + 1, 2);

    for (int i5 = 0; i5 < len_par; ++i5) {
      dden_part[i5] += pw2[jj] * (dnum[i5] * fto + num * dfto[i5]) * 2/pow(px21[jj] + 1, 2) * 2/pow(px22[jj] + 1, 2);
    }

  }

  return dnum_part/num_part - dden_part/den_part;

}

  

