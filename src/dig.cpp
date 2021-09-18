#include <Rcpp.h>
using namespace Rcpp;

//' Select IG copula quantities: matching inputs
//'
//' The density function, 1|2 conditional cdf,
//' and 1|2 conditional quantile function of the IG
//' copula family. Inputs
//' need to be vectors of the same length.
//' These functions are called by the R functions of the
//' same name, without the `_vec` suffix.
//' 
//' @param u,v Copula arguments. Vector of values between 0 and 1.
//' @param p Function inverse argument. Vector of values between 0 and 1.
//' @param theta,alpha IG copula parameters. Vector of positive values.
//' @note If calling these functions manually, make sure each input
//' are vectors of a common length. 
//' @details The `qcondig12()` function needs its own Newton
//' Raphson algorithm. It also needs access to some version 
//' of `pcondig12()` and `dig()`. So, these three functions
//' are coded up in C++, each with a scalar and vector pair
//' of functions. 
//' @seealso `dig()`, `pcondig12()`, and `qcondig12()`; 
//' and `igl_gen_vec()` and family.
//' @rdname ig_cpp_vec
// [[Rcpp::export]]
NumericVector dig_vec(NumericVector u, NumericVector v, NumericVector theta,
                      NumericVector alpha)
{ int i;
  int n = u.size();
  double dig_single (double, double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = dig_single(u[i],v[i],theta[i],alpha[i]); }
  return(out);
}

// The `dig()` function, with scalar inputs and output.
double dig_single(double u, double v, double theta, double alpha)
{ double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double interp_kappa_D1_single(double, double, double);
  double interp_gen_D1_single(double, double, double);
  double eps,bd,y,pdf;
  int mxiter;
  mxiter=20; eps=1.e-12; bd=5.;
  y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd, 0);
  pdf = interp_kappa_D1_single(y, (1.-u)*theta, alpha) /
    interp_gen_D1_single(y, theta, alpha);
  return(pdf);
}
