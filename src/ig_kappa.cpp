#include <Rcpp.h>
using namespace Rcpp;


double  dgamma(double x, double shape, double scale)
{ double xx,tem,pdf,prod;
  prod = x * shape * scale;
  if (ISNAN(prod)) return(prod);
  if(x==0. && shape==1.) return(1./scale);
  xx = x/scale;
  tem = ((shape-1.)*log(xx) - lgamma(shape) - xx);
  pdf = exp(tem)/scale;
  return(pdf);
}

// x, shape, alpha are vectors of length n
void dgamma_void(int *n0, double *x, double *shape, 
  double *scale, double *out)
{ int i,n;
  double dgamma (double, double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = dgamma(x[i],shape[i],scale[i]); R_CheckUserInterrupt();}
}

double igl_kappa(double x, double alpha)
{ double pgamma(double,double,double);
  double res,prod;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  res = 1.-pgamma(x, alpha, 1.);
  return(res);
}

// x and alpha are vectors of length n
void igl_kappa_vec(int *n0, double *x, double *alpha, double *out)
{ int i,n;
  double igl_kappa (double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = igl_kappa(x[i], alpha[i]); R_CheckUserInterrupt();}
}

double igl_kappa_D(double x, double alpha)
{ double dgamma(double,double,double);
  double res,prod;
  prod = x * alpha;
  if (ISNAN(prod)) return(prod);
  res = -dgamma(x, alpha, 1.);
  return(res);
}

// x and alpha are vectors of length n
void igl_kappa_D_vec(int *n0, double *x, double *alpha, double *out)
{ int i,n;
  double igl_kappa_D (double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = igl_kappa_D(x[i], alpha[i]); R_CheckUserInterrupt();}
}

double igl_kappa_inv(double p, double alpha)
{ double qgamma(double,double,double);
  double res,prod;
  prod = p * alpha;
  if (ISNAN(prod)) return(prod);
  res = qgamma(1.-p, alpha, 1.);
  return(res);
}

// x and alpha are vectors of length n
void igl_kappa_inv_vec(int *n0, double *p, double *alpha, double *out)
{ int i,n;
  double igl_kappa_inv (double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = igl_kappa_inv(p[i], alpha[i]); R_CheckUserInterrupt();}
}


double interp_kappa(double x, double eta, double alpha)
{ double igl_kappa(double, double);
  double res;
  res = exp(-x) * igl_kappa(eta*x, alpha);
  return(res);
}

// x, eta, alpha are vectors of length n
void interp_kappa_vec(int *n0, double *x, double *eta, 
  double *alpha, double *out)
{ int i,n;
  double interp_kappa (double, double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = interp_kappa(x[i], eta[i], alpha[i]); R_CheckUserInterrupt();}
}

double interp_kappa_D1(double x, double eta, double alpha)
{ double igl_kappa(double, double);
  double igl_kappa_D(double, double);
  double res;
  res = -exp(-x) * (igl_kappa(eta*x, alpha) - eta* igl_kappa_D(eta*x, alpha));
  return(res);
}

// x, eta, alpha are vectors of length n
void interp_kappa_D1_vec(int *n0, double *x, double *eta, double *alpha, double *out)
{ int i,n;
  double interp_kappa_D1 (double, double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = interp_kappa_D1(x[i], eta[i], alpha[i]); R_CheckUserInterrupt();}
}


double interp_kappa_inv_algo(double p, double eta, double alpha, int mxiter,
   double eps, double bd, int iprint)
{ double igl_kappa(double p, double alpha);
  double igl_kappa_D(double p, double alpha);
  double igl_kappa_inv(double p, double alpha);
  double interp_kappa(double x, double eta, double alpha);
  double x1,x2, p1,p2, diff1,diff2, x, ik, diff,g,gp,prod, logx;
  int iter;
  prod = alpha * eta * p;
  if (ISNAN(prod)) return(prod);
  if (p <= 0.) return(DBL_MAX);
  if (p >= 1.) return(0.);
  x1 = -log(p);
  x2 = igl_kappa_inv(p, alpha) / eta;
  p1 = interp_kappa(x1, eta, alpha);
  p2 = interp_kappa(x2, eta, alpha);
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
    ik = interp_kappa(x, eta, alpha);
    g = log(ik) - log(p);
    gp = interp_kappa_D1(x, eta, alpha) / ik * x;
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

// p, eta and alpha are vectors of length n
void interp_kappa_inv(int *n0, double *p, double *eta, double *alpha,
  int *mxiter0, double *eps0, double *bd0, double *inv)
{ int i,n,mxiter;
  double interp_kappa_inv_algo (double, double, double, int, double, double, int);
  double eps,bd;
  n=*n0; mxiter=*mxiter0; eps=*eps0; bd=*bd0;
  for(i=0;i<n;i++) { 
    inv[i] = interp_kappa_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0); 
    R_CheckUserInterrupt();
  }
}

