#include <Rcpp.h>
using namespace Rcpp;

//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector igl_gen_D_vec(NumericVector x, NumericVector alpha)
{ int i;
  int n = x.size();
  double igl_gen_D_single (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_gen_D_single(x[i],alpha[i]); }
  return(out);
}

// The `igl_gen_D()` function, with scalar inputs and output.
double igl_gen_D_single(double x, double alpha)
{ double pgamma(double, double, double, int, int);
  double dgamma(double, double, double, int);
  double prod;
  prod = x * alpha;
  if (ISNAN(prod)) { return(prod); }
  if (x == 0.) {
    return(-R::dgamma(0., alpha, 1., 0) / 2.);
  }
  return(
    -(alpha / (x*x)) * R::pgamma(x, alpha + 1., 1., 1, 0)
  );
}
