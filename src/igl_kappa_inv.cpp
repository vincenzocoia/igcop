#include <Rcpp.h>
using namespace Rcpp;

//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector igl_kappa_inv_vec(NumericVector p, NumericVector alpha)
{ int n = p.size();
  int i;
  double igl_kappa_inv_single (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_kappa_inv_single(p[i], alpha[i]); }
  return(out);
}

// The `igl_kappa_inv()` function, with scalar inputs and output.
double igl_kappa_inv_single(double p, double alpha)
{ double qgamma(double, double, double, int, int);
  double res,prod;
  prod = p * alpha;
  if (ISNAN(prod)) return(prod);
  res = R::qgamma(1. - p, alpha, 1., 1, 0);
  return(res);
}
