#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector fn_sumxtx(NumericMatrix x) {

  int nc = x.ncol();
  int nr = x.nrow();

  NumericVector res(nc * nc);

  for (int i = 0; i < nr; ++i) {

    for (int j = 0; j < nc; ++j) {

      for (int k = 0; k < nc; ++k) {

        res(nc * j + k) += x(i, j) * x(i, k);
      }

    }

  }

  return res;
}
