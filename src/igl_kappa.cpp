#include <Rcpp.h>
using namespace Rcpp;

//' Functions generating the IG and IGL copulas
//'
//' Kappa function and its relatives have prefix `igl_kappa`;
//' Psi function and its relatives have prefix `igl_gen`;
//' Interpolating function H with either kappa or psi has
//' `igl` prefix replaced with `interp`. Relatives of these functions:
//' suffix `inv` indicates inverse; suffix `D` represents function
//' derivative, and `D1` derivative with respect to the first argument.
//'. Suffix `_vec` indicates that the entries must be vectors of
//' the same length; `_single` means entries must be
//' scalars.
//' @param x Function argument
//' @param p Inverse function argument
//' @param eta,alpha Function parameters
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

double igl_kappa_single(double x, double alpha)
{ double pgamma(double, double, double, int, int);
  double prod;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  return(R::pgamma(x, alpha, 1., 0, 0));
}

