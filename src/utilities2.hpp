#ifndef __UTILITIES2__
#define __UTILITIES2__

arma::vec fn_expQ(int cstat, double t1, double t2, arma::vec q, double gamma_all);
arma::mat fn_dexpQ(int cstat, double t1, double t2, arma::vec q, arma::vec Xi, double gamma_all, int len_par);

#endif //__UTILITIES2__