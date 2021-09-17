#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector interp_kappa_inv_vec(NumericVector p, NumericVector eta,
                                   NumericVector alpha, int mxiter, double eps, double bd)
{ int n = p.size();
  int i;
  double interp_kappa_inv_algo (double, double, double, int, double, double, int);
  NumericVector inv(n);
  for(i=0;i<n;i++) {
    inv[i] = interp_kappa_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0);
  }
  return(inv);
}

double interp_kappa_inv_algo(double p, double eta, double alpha, int mxiter,
                             double eps, double bd, int iprint)
{ double igl_kappa_single(double p, double alpha);
  double igl_kappa_D_single(double p, double alpha);
  double igl_kappa_inv_single(double p, double alpha);
  double interp_kappa_single(double x, double eta, double alpha);
  double interp_kappa_D1_single(double x, double eta, double alpha);
  double x1,x2, p1,p2, diff1,diff2, x, ik, diff,g,gp,prod, logx;
  int iter;
  prod = alpha * eta * p;
  if (ISNAN(prod)) return(prod);
  if (p <= 0.) return(DBL_MAX);
  if (p >= 1.) return(0.);
  x1 = -log(p);
  x2 = igl_kappa_inv_single(p, alpha) / eta;
  p1 = interp_kappa_single(x1, eta, alpha);
  p2 = interp_kappa_single(x2, eta, alpha);
  diff1 = fabs(p1-p); diff2 = fabs(p2-p);
  x=x1;
  if(diff2<diff1) { x=x2; }
  iter = 0; diff = 1.;
  while(iter < mxiter && fabs(diff) > eps)
  { //pex = p * exp(x);
    //g = igl_kappa(eta*x, alpha) - pex;
    //gp = eta * igl_kappa_D(eta*x, alpha) - pex;
    //enx = exp(-x);
    //ik = enx * igl_kappa(eta*x, alpha);
    //g = ik - p;
    //gp = enx * eta * igl_kappa_D(eta*x, alpha) - ik;
    logx = log(x);
    ik = interp_kappa_single(x, eta, alpha);
    g = log(ik) - log(p);
    gp = interp_kappa_D1_single(x, eta, alpha) / ik * x;
    diff = g / gp;
    if (diff > bd) diff = bd;
    if (diff < -bd) diff = -bd;
    //if (x - diff < 0.) diff = x / 2.;
    logx -= diff;
    x = exp(logx);
    //while(fabs(diff) > bd)
    //{ diff /= 2.;  x += diff; R_CheckUserInterrupt(); }
    iter++;
    if(iprint) Rprintf("%d %f %f\n", iter, x, diff);
    R_CheckUserInterrupt();
  }
  return(x);
}