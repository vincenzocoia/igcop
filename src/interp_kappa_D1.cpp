#include <Rcpp.h>
using namespace Rcpp;

// The `interp_kappa_D1()` function, with scalar inputs and output.
// Note that no _vec version of this function is needed, as
// it is only ever called from other _single or _algo
// C++ functions (and not from R).
double interp_kappa_D1_single(double x, double eta, double alpha)
{ double igl_kappa_single(double, double);
  double igl_kappa_D_single(double, double);
  double res;
  res = -exp(-x) * (igl_kappa_single(eta * x, alpha) -
    eta* igl_kappa_D_single(eta * x, alpha));
  return(res);
}