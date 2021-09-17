#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector igl_gen_inv_vec(NumericVector p, NumericVector alpha)
{ int n = p.size();
  int i;
  double igl_gen_inv_algo (double, double, int, double, double, int);
  NumericVector inv(n);
  double eps = 1.e-12, bd = 5.;
  int mxiter = 20;
  for(i=0;i<n;i++) {
    inv[i] = igl_gen_inv_algo(p[i],alpha[i],mxiter,eps,bd, 0);
  }
  return(inv);
}

double igl_gen_inv_algo (double p, double alpha, int mxiter, double eps,
                         double bd, int iprint)
{ double x1,x2,x3,p1,p2,p3,diff1,diff2,diff3;
  double x,best,diff,g,gp,prod;
  double qgamma(double, double, double, int, int);
  double igl_gen_single(double, double);
  double igl_gen_D_single(double, double);
  int iter;
  prod = alpha * p;
  if (ISNAN(prod)) return(prod);
  if (p == 0.) return(DBL_MAX);  // Inf
  if (p == 1.) return(0.);
  x1 = 1. / (pow(1.-p, -1./ alpha) - 1.);
  x2 = alpha / p;
  x3 = R::qgamma(p, alpha + 1., 1., 1, 0);
  p1 = igl_gen_single(x1, alpha);
  p2 = igl_gen_single(x2, alpha);
  p3 = igl_gen_single(x3, alpha);
  diff1 = fabs(p1 - p);
  diff2 = fabs(p2 - p);
  diff3 = fabs(p3 - p);
  x = x1;
  best = diff1;
  if(diff2<best) { x=x2; best=diff2; }
  if(diff3<best) { x=x3; best=diff3; }
  if (x == 0.) x = eps;
  //ex = x;
  //x = log(x);
  iter = 0; diff = 1.;
  while (iter < mxiter && fabs(diff) > eps)
  { //g = igl_gen(ex, alpha) - p;
    //gp = igl_gen_D(ex, alpha) * ex;
    g = igl_gen_single(x, alpha) - p;
    gp = igl_gen_D_single(x, alpha);
    diff = g / gp;
    if (diff > bd) diff = bd;
    if (diff < -bd) diff = -bd;
    if (x - diff < 0.) diff = x / 2.;
    x -= diff;
    //while (fabs(diff) > bd)
    //{ diff/=2.; x += diff; R_CheckUserInterrupt(); }
    //ex = exp(x);
    iter++;
    if(iprint) Rprintf("%d %f %f\n", iter, x, diff);
    R_CheckUserInterrupt();
  }
  //return(exp(x));
  return(x);
}
