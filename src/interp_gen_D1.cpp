#include <Rcpp.h>
using namespace Rcpp;

// The `interp_gen_D1()` function, with scalar inputs and output.
// Note that no _vec version of this function is needed, as
// it is only ever called from other _single or _algo
// C++ functions (and not from R).
double interp_gen_D1_single (double x, double eta, double alpha)
{ double igl_gen_single(double, double);
  double igl_gen_D_single(double, double);
  double dgamma(double, double, double, int);
  double res,prod;
  prod = x * eta * alpha;
  if (ISNAN(prod)) return(prod);
  if (x == 0.)
  { res = - (1. + eta / 2. * R::dgamma(0., alpha, 1., 0)); }
  else
  { res = -exp(-x) *
    (igl_gen_single(eta*x, alpha) - eta * igl_gen_D_single(eta*x, alpha)); }
  return(res);
}
