#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec fn_egval(arma::vec q) {
  
  
  arma::vec egval(4);


  egval[0] = 0;
  egval[1] = 0.5 * (- sqrt(pow(q[0] + q[1] + q[2] + q[3] + q[4], 2.0) - 4 *(q[0]*q[3] + q[0]*q[4] + q[1]*q[2] + q[1]*q[3] + q[1]*q[4])) - (q[0] + q[1] + q[2] + q[3] + q[4]));
  egval[2] = 0.5 * (sqrt(pow(q[0] + q[1] + q[2] + q[3] + q[4], 2.0) - 4 *(q[0]*q[3] + q[0]*q[4] + q[1]*q[2] + q[1]*q[3] + q[1]*q[4])) - (q[0] + q[1] + q[2] + q[3] + q[4]));
  egval[3] = -q[5];
  
  return egval;
}

// 
// double a = q[0];
// double b = q[1];
// double c = q[2];
// double d = q[3];
// double e = q[4];
// double f = q[5];
// 
// a = q[1];
// b = q[2];
// c = q[3];
// d = q[4];
// e = q[5];
// f = q[6];
// 
// egval1 = 0;
// egval2 = 1/2 * (- sqrt((a + b + c + d + e)^2 - 4 *(a*d + a*e + b*c + b*d + b*e)) - (a + b + c + d + e));
// egval3 = 1/2 * (sqrt((a + b + c + d + e)^2 - 4 *(a*d + a*e + b*c + b*d + b*e)) - (a + b + c + d + e));
// egval4 = -f;
// egval = c(egval1, egval2, egval3, egval4)
//   
//   egvec1 = c(1, 1, 1, 1);
// egvec2 = c(-(a + b - c - d - e + sqrt(a^2 + 2*a*b + b^2 + 2*a*c - 2*b*c + c^2 - 2*a*d - 2*b*d + 2*c*d + d^2 - 2*a*e - 2*b*e + 2*c*e + 2*d*e + e^2))/(2*c), 1, 0, 0)
//   egvec3 = c(-(a + b - c - d - e - sqrt(a^2 + 2*a*b + b^2 + 2*a*c - 2*b*c + c^2 - 2*a*d - 2*b*d + 2*c*d + d^2 - 2*a*e - 2*b*e + 2*c*e + 2*d*e + e^2))/(2*c), 1, 0, 0)
//   egvec4 = c(- a*d/(-b*c - a*d - b*d - a*e - b*e + a*f + b*f + c*f + d*f + e*f -f^2),
//              -(-a*d - b*d + d*f)/(b*c + a*d + b*d + a*e + b*e - a*f - b*f - c*f - d*f - e*f +f^2), 1, 0)
//   egvec = cbind(egvec1, egvec2, egvec3, egvec4)
//   
//   ei <- eigen(Q)
//   (Q - diag(ei$values))
//   
//   ei$vectors %*% diag(ei$values) %*% solve(ei$vectors)
//   
//   egvec %*% diag(egval) %*% solve(egvec)
//   