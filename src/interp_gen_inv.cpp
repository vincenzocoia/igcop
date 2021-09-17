#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector interp_gen_inv_vec(NumericVector p, NumericVector eta,
                                 NumericVector alpha, int mxiter, double eps, double bd)
{ int i;
  int n = p.size();
  double interp_gen_inv_algo (double, double, double, int, double, double, int);
  NumericVector inv(n);
  for(i=0;i<n;i++) {
    inv[i] = interp_gen_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0);
    R_CheckUserInterrupt();
  }
  return(inv);
}

double interp_gen_inv_algo(double p, double eta, double alpha, int mxiter,
                           double eps, double bd, int iprint)
{ double x1,x2,p1,p2,diff1,diff2;
  double x,diff,g,gp,prod;
  double interp_gen_single(double, double, double);
  double interp_gen_D1_single(double, double, double);
  double igl_gen_inv_algo (double, double, int, double, double, int);
  int iter;
  prod = alpha * eta * p;
  if (ISNAN(prod)) return(prod);
  if (p <= 0.) return(DBL_MAX); // Inf
  if (p >= 1.) return(0.);
  x1 = -log(p);
  x2 = igl_gen_inv_algo(p, alpha, mxiter, eps, bd, 0) / eta;
  p1 = interp_gen_single(x1, eta, alpha);
  p2 = interp_gen_single(x2, eta, alpha);
  diff1 = fabs(p1-p); diff2 = fabs(p2-p);
  x=x1;
  if(diff2<diff1) { x=x2; }
  //x = log(x);
  iter = 0; diff = 1.;
  while(iter < mxiter && fabs(diff) > eps)
  { //ex = exp(x);
    //g = interp_gen(ex, eta, alpha) - p;
    //gp = interp_gen_D1(ex, eta, alpha) * ex;
    g = interp_gen_single(x, eta, alpha) - p;
    gp = interp_gen_D1_single(x, eta, alpha);
    diff = g / gp;
    if (diff > bd) diff = bd;
    if (diff < -bd) diff = -bd;
    if (x - diff < 0.) diff = x / 2.;
    x -= diff;
    //while (fabs(diff) > bd)
    //{ diff/=2.; x += diff; R_CheckUserInterrupt(); }
    iter++;
    if(iprint) Rprintf("%d %f %f\n", iter, x, diff);
    R_CheckUserInterrupt();
  }
  //return(exp(x));
  return(x);
}
