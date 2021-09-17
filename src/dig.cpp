#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dig_vec(NumericVector u, NumericVector v, NumericVector theta,
                      NumericVector alpha)
{ int i;
  int n = u.size();
  double dig_single (double, double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = dig_single(u[i],v[i],theta[i],alpha[i]); }
  return(out);
}

// 0<u<1, 0<v<1, theta>0 alpha>0
double dig_single(double u, double v, double theta, double alpha)
{ double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double interp_kappa_D1_single(double, double, double);
  double interp_gen_D1_single(double, double, double);
  double eps,bd,y,pdf;
  int mxiter;
  mxiter=20; eps=1.e-12; bd=5.;
  y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd, 0);
  pdf = interp_kappa_D1_single(y, (1.-u)*theta, alpha) /
    interp_gen_D1_single(y, theta, alpha);
  return(pdf);
}
