#include <Rcpp.h>
using namespace Rcpp;

//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector igl_kappa_D_vec(NumericVector x, NumericVector alpha)
{ int n = x.size();
  int i;
  double igl_kappa_D_single (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_kappa_D_single(x[i], alpha[i]); }
  return(out);
}

// The `igl_kappa_D()` function, with scalar inputs and output.
double igl_kappa_D_single(double x, double alpha)
{ double dgamma(double, double, double, int);
  double res,prod;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  res = -R::dgamma(x, alpha, 1., 0);
  return(res);
}
