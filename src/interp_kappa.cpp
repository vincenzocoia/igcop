#include <Rcpp.h>
using namespace Rcpp;

//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector interp_kappa_vec(NumericVector x, NumericVector eta,
                               NumericVector alpha)
{ int n = x.size();
  int i;
  double interp_kappa_single (double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = interp_kappa_single(x[i], eta[i], alpha[i]); }
  return(out);
}

// The `interp_kappa()` function, with scalar inputs and output.
double interp_kappa_single(double x, double eta, double alpha)
{ double igl_kappa_single(double, double);
  double res;
  res = exp(-x) * igl_kappa_single(eta*x, alpha);
  return(res);
}
