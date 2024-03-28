#ifndef __UTILITIES__
#define __UTILITIES__

arma::vec fn_egvec(arma::vec q);
arma::mat fn_egval(arma::vec q);
double fn_h(double t, double gamma_all);
double fn_g(double t, double gamma_all);
arma::vec fn_dh(double t, double gamma_all, int len_par);
arma::vec fn_dg(double t, double gamma_all, int len_par);
arma::mat fun_dQ(arma::vec q, arma::vec xi, double gamma_all, int len_par);

#endif //__UTILITIES__