#include <Rcpp.h>
using namespace Rcpp;

//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector igl_kappa_vec(NumericVector x, NumericVector alpha)
{ int n = x.size();
  int i;
  double igl_kappa_single (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_kappa_single(x[i], alpha[i]); }
  return(out);
}

// The `igl_kappa()` function, with scalar inputs and output.
double igl_kappa_single(double x, double alpha)
{ double pgamma(double, double, double, int, int);
  double prod;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  return(R::pgamma(x, alpha, 1., 0, 0));
}

