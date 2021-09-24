#include <Rcpp.h>
using namespace Rcpp;

//' IG/IGL Generators and Related Functions: matching inputs
//'
//' These are the psi, H, and kappa functions
//' of the IG and IGL copula families, but with inputs
//' needing to be vectors of the same length.
//' These functions are called by the R functions of the
//' same name, without the `_vec` suffix.
//' @param x Function argument. Vector of non-negative values.
//' @param p Function inverse argument. Vector of values between 0 and 1.
//' @param eta,alpha Function parameters. Vector of positive values.
//' @note If calling this function manually, make sure each input
//' are vectors of a common length. 
//' @seealso `igl_gen()` and family; 
//' `dig_vec()`, `pcondig12_vec()`, and `qcondig12_vec()`. 
//' @rdname generators_vec
// [[Rcpp::export]]
NumericVector igl_gen_vec(NumericVector x, NumericVector alpha)
{ int n = x.size();
  int i;
  double igl_gen_single (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_gen_single(x[i],alpha[i]); }
  return(out);
}

// The `igl_gen()` function, with scalar inputs and output.
double igl_gen_single(double x, double alpha)
{ double pgamma(double, double, double, int, int);
  double prod,res,term1,term2;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  if (x == 0.) return(1.);
  term1 = R::pgamma(x, alpha, 1., 0, 0);
  term2 = alpha * R::pgamma(x, alpha + 1., 1., 1, 0) / x;
  res = term1+term2;
  return(res);
}

