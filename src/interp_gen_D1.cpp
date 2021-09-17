#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector interp_gen_D1_vec(NumericVector x, NumericVector eta, NumericVector alpha)
{ int n = x.size();
  int i;
  double interp_gen_D1_single (double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = interp_gen_D1_single(x[i], eta[i], alpha[i]); }
  return(out);
}

double interp_gen_D1_single (double x, double eta, double alpha)
{ double igl_gen_single(double, double);
  double igl_gen_D_single(double, double);
  double dgamma(double, double, double, int);
  double res,prod;
  prod = x * eta * alpha;
  if (ISNAN(prod)) return(prod);
  if (x == 0.)
  { res = - (1. + eta / 2. * R::dgamma(0., alpha, 1., 0)); }
  else
  { res = -exp(-x) *
    (igl_gen_single(eta*x, alpha) - eta * igl_gen_D_single(eta*x, alpha)); }
  return(res);
}
