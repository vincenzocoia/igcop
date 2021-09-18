#include <Rcpp.h>
using namespace Rcpp;

//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector interp_gen_vec(NumericVector x, NumericVector eta,
                             NumericVector alpha)
{ int i;
  int n = x.size();
  double interp_gen_single (double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = interp_gen_single(x[i],eta[i],alpha[i]); }
  return(out);
}

// The `interp_gen()` function, with scalar inputs and output.
double interp_gen_single (double x, double eta, double alpha)
{ double igl_gen_single(double, double);
  double res;
  res = exp(-x) * igl_gen_single(eta*x, alpha);
  return(res);
}


