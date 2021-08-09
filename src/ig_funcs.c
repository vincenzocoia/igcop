#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <values.h>

#ifdef MAIN
#include <malloc.h>
// sample main program for checking functions
int main(int argc, char *argv[])
{ double igl_gen(double x, double alpha);
  double interp_gen (double x, double eta, double alpha);
  double x, eta, alpha, iglf, interpf;
  int n,i,mxiter,iprint,j;
  double eps,bd,p,iglinv,interpinv;
  double igl_gen_inv_algo (double, double, int, double, double, int);
  double interp_gen_inv_algo (double, double, double, int, double, double, int);
  void interp_gen_inv (int *, double *, double *, double *, int *, double *, double *, double *);
  double *pvec,*avec,*etavec,*inv;

  n=10;
  pvec=(double *) malloc(n * sizeof(double));
  avec=(double *) malloc(n * sizeof(double));
  etavec=(double *) malloc(n * sizeof(double));
  inv=(double *) malloc(n * sizeof(double));
  for(i=1;i<=n;i++)
  { x=i*0.5; alpha=n-i+0.4;
    iglf = igl_gen(x,alpha);
    Rprintf("igl %f %f %f\n", x,alpha,iglf);
  }
  for(i=1;i<=n;i++)
  { x=i*0.1+0.91; alpha=n-i+4.4; eta=i*0.1;
    interpf = interp_gen(x,eta,alpha);
    Rprintf("interp %f %f %f %f\n", x,eta,alpha,interpf);
  }

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  for(i=1;i<=n;i++)
  { p=i/(n+1.); alpha=n-i+0.4;
    if(i>=3) iprint=0;
    iglinv = igl_gen_inv_algo(p,alpha,mxiter,eps,bd, iprint);
    iglf = igl_gen(iglinv,alpha);
    Rprintf("igl_gen_inv %f %f %f %f\n", p,alpha,iglinv,iglf);
  }

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  for(i=1;i<=n;i++)
  { p=i/(n+1.); alpha=n-i+4.4; eta=i*0.1;
    if(i>=5) iprint=0;
    interpinv = interp_gen_inv_algo(p,eta,alpha,mxiter,eps,bd, iprint);
    interpf = interp_gen(interpinv,eta,alpha);
    Rprintf("interp_gen_inv %f %f %f %f %f\n", p,eta,alpha,interpinv,interpf);
  }

  // second group of inputs
  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  Rprintf("\nAnother group of inputs\n");
  for(i=1;i<=5;i++)
  { p=1.-i*.1; alpha=1.1; eta=2.5;
    interpinv = interp_gen_inv_algo(p,eta,alpha,mxiter,eps,bd, iprint);
    interpf = interp_gen(interpinv,eta,alpha);
    Rprintf("interp_gen_inv %f %f %f %f %f\n", p,eta,alpha,interpinv,interpf);
  }

  for(j=0;j<n;j++)
  { i=j+1;; pvec[j] = i/(n+1.); avec[j]=n-i+4.4; etavec[j]=i*0.1; }
  interp_gen_inv(&n,pvec,etavec,avec,&mxiter,&eps,&bd,inv);
  for(j=0;j<n;j++) Rprintf("%f ", inv[j]); Rprintf("\n");
  free(pvec); free(avec); free(etavec); free(inv);
  return(0);
}
#endif


double igl_gen(double x, double alpha)
{ double pgamma(double,double,double);
  double prod,res,term1,term2;
  prod = x * alpha;
  if (isnan(prod)) return(prod);
  if (x == 0.) return(1.);
  term1 = 1.-pgamma(x, alpha, 1.);
  term2 = alpha * pgamma(x, alpha+1., 1.) / x;
  res = term1+term2;
  return(res);
}

// p and alpha are vectors of length n
void igl_gen_vec(int *n0, double *x, double *alpha, double *out)
{ int i,n;
  double igl_gen (double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = igl_gen(x[i],alpha[i]); R_CheckUserInterrupt();}
}

double igl_gen_D(double x, double alpha)
{ double pgamma(double,double,double);
  double dgamma(double,double,double);
  double prod;
  prod = x * alpha;
  if (isnan(prod)) { return(prod); }
  if (x == 0.) { return(-dgamma(0., alpha, 1.) / 2.); }
  return( -(alpha / (x*x)) * pgamma(x, alpha+1., 1.) );
}

// x and alpha are vectors of length n
void igl_gen_D_vec(int *n0, double *x, double *alpha, double *out)
{ int i,n;
  double igl_gen_D (double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = igl_gen_D(x[i],alpha[i]); R_CheckUserInterrupt();}
}

double igl_gen_inv_algo (double p, double alpha, int mxiter, double eps,
                         double bd, int iprint)
{ double x1,x2,x3,p1,p2,p3,diff1,diff2,diff3;
  double x,best,diff,g,gp,prod;
  double qgamma(double,double,double);
  double igl_gen(double, double);
  double igl_gen_D(double, double);
  int iter;
  prod = alpha * p;
  if (isnan(prod)) return(prod);
  if (p == 0.) return(DBL_MAX);  // Inf
  if (p == 1.) return(0.);
  x1 = 1. / (pow(1.-p, -1./ alpha) - 1.);
  x2 = alpha / p;
  x3 = qgamma(p, alpha+1., 1.);
  p1 = igl_gen(x1,alpha); p2 = igl_gen(x2,alpha); p3 = igl_gen(x3,alpha);
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
    g = igl_gen(x, alpha) - p;
    gp = igl_gen_D(x, alpha);
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
void igl_gen_inv(int *n0, double *p, double *alpha,
  int *mxiter0, double *eps0, double *bd0, double *inv)
{ int i,n,mxiter;
  double igl_gen_inv_algo (double, double, int, double, double, int);
  double eps,bd;
  n=*n0; mxiter=*mxiter0; eps=*eps0; bd=*bd0;
  for(i=0;i<n;i++) {
    inv[i] = igl_gen_inv_algo(p[i],alpha[i],mxiter,eps,bd, 0);
    R_CheckUserInterrupt();
  }
}

double interp_gen (double x, double eta, double alpha)
{ double igl_gen(double, double);
  double res;
  res = exp(-x) * igl_gen(eta*x, alpha);
  return(res);
}

// x, eta, alpha are vectors of length n
void interp_gen_vec(int *n0, double *x, double *eta,
  double *alpha, double *out)
{ int i,n;
  double interp_gen (double, double, double);
  n=*n0;
  for(i=0;i<n;i++)
  { out[i] = interp_gen(x[i],eta[i],alpha[i]); R_CheckUserInterrupt();}
}

double interp_gen_D1 (double x, double eta, double alpha)
{ double igl_gen(double, double);
  double igl_gen_D(double, double);
  double dgamma(double, double, double);
  double res,prod;
  prod = x * eta * alpha;
  if (isnan(prod)) return(prod);
  if (x == 0.)
  { res = - (1. + eta / 2. * dgamma(0., alpha, 1.)); }
  else
  { res = -exp(-x) * (igl_gen(eta*x, alpha) - eta * igl_gen_D(eta*x, alpha)); }
  return(res);
}

double interp_gen_inv_algo(double p, double eta, double alpha, int mxiter,
  double eps, double bd, int iprint)
{ double x1,x2,p1,p2,diff1,diff2;
  double x,diff,g,gp,prod;
  double interp_gen(double, double, double);
  double interp_gen_D1(double, double, double);
  double igl_gen_inv_algo (double, double, int, double, double, int);
  int iter;
  prod = alpha * eta * p;
  if (isnan(prod)) return(prod);
  if (p <= 0.) return(DBL_MAX); // Inf
  if (p >= 1.) return(0.);
  x1 = -log(p);
  x2 = igl_gen_inv_algo(p, alpha, mxiter, eps, bd, 0) / eta;
  p1 = interp_gen(x1, eta, alpha);
  p2 = interp_gen(x2, eta, alpha);
  diff1 = fabs(p1-p); diff2 = fabs(p2-p);
  x=x1;
  if(diff2<diff1) { x=x2; }
  //x = log(x);
  iter = 0; diff = 1.;
  while(iter < mxiter && fabs(diff) > eps)
  { //ex = exp(x);
    //g = interp_gen(ex, eta, alpha) - p;
    //gp = interp_gen_D1(ex, eta, alpha) * ex;
    g = interp_gen(x, eta, alpha) - p;
    gp = interp_gen_D1(x, eta, alpha);
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
void interp_gen_inv(int *n0, double *p, double *eta, double *alpha,
  int *mxiter0, double *eps0, double *bd0, double *inv)
{ int i,n,mxiter;
  double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double eps,bd;
  n=*n0; mxiter=*mxiter0; eps=*eps0; bd=*bd0;
  for(i=0;i<n;i++) {
    inv[i] = interp_gen_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0);
    R_CheckUserInterrupt();
  }
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
