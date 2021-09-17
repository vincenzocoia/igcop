#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector igl_gen_vec(NumericVector x, NumericVector alpha)
{ int n = x.size();
  int i;
  double igl_gen (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_gen_single(x[i],alpha[i]); }
  return(out);
}

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

