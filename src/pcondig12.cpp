#include <Rcpp.h>
using namespace Rcpp;

//' @rdname ig_cpp_vec
// [[Rcpp::export]]
NumericVector pcondig12_vec(NumericVector u, NumericVector v,
                            NumericVector theta, NumericVector alpha)
{ int n = u.size();
  int i;
  double pcondig12_single (double, double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = pcondig12_single(u[i],v[i],theta[i],alpha[i]); }
  return(out);
}

// The `pcondig12()` function, with scalar inputs and output.
double pcondig12_single(double u, double v, double theta, double alpha)
{ double interp_gen_inv_algo (double, double, double, int, double, double);
  double interp_gen_D1_single(double, double, double);
  double eps,bd,y,pcond;
  int mxiter;
  mxiter=25; eps=1.e-13; bd=5.;
  y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd);
  pcond = 1. - (1.-u) *
    interp_gen_D1_single(y, theta*(1.-u), alpha) /
      interp_gen_D1_single(y, theta, alpha);
  return(pcond);
}


