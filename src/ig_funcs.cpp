#include <Rcpp.h>
using namespace Rcpp;


double igl_gen_single(double x, double alpha)
{ double pgamma(double,double,double);
  double prod,res,term1,term2;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  if (x == 0.) return(1.);
  term1 = 1.-pgamma(x, alpha, 1.);
  term2 = alpha * pgamma(x, alpha+1., 1.) / x;
  res = term1+term2;
  return(res);
}

// p and alpha are vectors of length n
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

double igl_gen_D_single(double x, double alpha)
{ double pgamma(double,double,double);
  double dgamma(double,double,double);
  double prod;
  prod = x * alpha;
  if (ISNAN(prod)) { return(prod); }
  if (x == 0.) { return(-dgamma(0., alpha, 1.) / 2.); }
  return( -(alpha / (x*x)) * pgamma(x, alpha+1., 1.) );
}

// x and alpha are vectors of length n
// [[Rcpp::export]]
NumericVector igl_gen_D_vec(NumericVector x, NumericVector alpha)
{ int i;
  int n = x.size();
  double igl_gen_D_single (double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = igl_gen_D_single(x[i],alpha[i]); }
  return(out);
}

double igl_gen_inv_algo (double p, double alpha, int mxiter, double eps,
                         double bd, int iprint)
{ double x1,x2,x3,p1,p2,p3,diff1,diff2,diff3;
  double x,best,diff,g,gp,prod;
  double qgamma(double,double,double);
  double igl_gen_single(double, double);
  double igl_gen_D_single(double, double);
  int iter;
  prod = alpha * p;
  if (ISNAN(prod)) return(prod);
  if (p == 0.) return(DBL_MAX);  // Inf
  if (p == 1.) return(0.);
  x1 = 1. / (pow(1.-p, -1./ alpha) - 1.);
  x2 = alpha / p;
  x3 = qgamma(p, alpha+1., 1.);
  p1 = igl_gen_single(x1,alpha);
  p2 = igl_gen_single(x2,alpha);
  p3 = igl_gen_single(x3,alpha);
  diff1 = fabs(p1-p); diff2 = fabs(p2-p); diff3 = fabs(p3-p);
  x=x1; best=diff1;
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

// p and alpha are vectors of length n
// [[Rcpp::export]]
NumericVector igl_gen_inv_vec(NumericVector p, NumericVector alpha,
  int mxiter, double eps, double bd)
{ int n = p.size();
  int i,mxiter;
  double igl_gen_inv_algo (double, double, int, double, double, int);
  double eps,bd;
  NumericVector inv(n);
  for(i=0;i<n;i++) {
    inv[i] = igl_gen_inv_algo(p[i],alpha[i],mxiter,eps,bd, 0);
  }
  return(inv);
}

double interp_gen_single (double x, double eta, double alpha)
{ double igl_gen_single(double, double);
  double res;
  res = exp(-x) * igl_gen_single(eta*x, alpha);
  return(res);
}

// x, eta, alpha are vectors of length n
// [[Rcpp::export]]
NumericVector interp_gen_vec(NumericVector x, NumericVector eta,
                             NumericVector alpha)
{ int i;
  int n = x.size();
  double interp_gen_single (double, double, double);
  NumericVector out(n);
  for(i=0;i<n;i++)
  { out[i] = interp_gen(x[i],eta[i],alpha[i]); }
  return(out);
}

double interp_gen_D1_single (double x, double eta, double alpha)
{ double igl_gen_single(double, double);
  double igl_gen_D_single(double, double);
  double dgamma(double, double, double);
  double res,prod;
  prod = x * eta * alpha;
  if (ISNAN(prod)) return(prod);
  if (x == 0.)
  { res = - (1. + eta / 2. * dgamma(0., alpha, 1.)); }
  else
  { res = -exp(-x) *
    (igl_gen_single(eta*x, alpha) - eta * igl_gen_D_single(eta*x, alpha)); }
  return(res);
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

// p, eta and alpha are vectors of length n
// [[Rcpp::export]]
NumericVector interp_gen_inv_vec(NumericVector p, NumericVector eta,
  NumericVector alpha, int mxiter, double eps, double bd)
{ int i,mxiter;
  int n = p.size();
  double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double eps,bd;
  NumericVector inv(n);
  for(i=0;i<n;i++) {
    inv[i] = interp_gen_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0);
    R_CheckUserInterrupt();
  }
  return(inv);
}


// As per Rcpp GitHub Issue #636
// https://github.com/RcppCore/Rcpp/issues/636#issuecomment-280985661
// Attempt to avoid the message upon a package Check:
//-----
//   Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’
//
// It is good practice to use registered native symbols and to disable
//   symbol search.
//-----
void R_init_igcop(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
