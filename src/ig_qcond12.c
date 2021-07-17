#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

#ifdef MAIN2
#include <malloc.h>
// sample main program for checking qcondig12
int main(int argc, char *argv[])
{ double pcondig12(double, double, double, double);
  double qcondig12_algo(double p, double v, double theta, double alpha,
    int mxiter, double eps, double bd, int iprint);
  void qcondig12(int *, double *, double *, double *, double *, int *,
    double *, double *, double *);
  double p, v, theta, alpha, qu, pp;
  int n,i,mxiter,iprint,j;
  double eps,bd;
  double *pvec,*vvec,*thvec,*avec,*quvec;

  n=5;
  pvec=(double *) malloc(n * sizeof(double));
  vvec=(double *) malloc(n * sizeof(double));
  thvec=(double *) malloc(n * sizeof(double));
  avec=(double *) malloc(n * sizeof(double));
  quvec=(double *) malloc(n * sizeof(double));

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  for(i=1;i<=n;i++)
  { p=i/(n+1.);  v=0.3; theta=i*0.8; alpha=n-i+0.4;
    //if(i>=3) iprint=0;
    qu = qcondig12_algo(p,v,theta,alpha, mxiter,eps,bd, iprint);
    pp = pcondig12(qu,v,theta,alpha);
    printf("qcondig12_algo %f %f %f %f %f %f\n", p,v,theta,alpha,qu,pp);
  }

  for(j=0;j<n;j++)
  { i=j+1;; pvec[j]=i/(n+1.); vvec[j]=0.3; thvec[j]=i*0.8; avec[j]=n-i+0.4; }
  qcondig12(&n,pvec,vvec,thvec,avec,&mxiter,&eps,&bd,quvec);
  for(j=0;j<n;j++) printf("%f ", quvec[j]); printf("\n");
  free(pvec); free(vvec); free(thvec); free(avec); free(quvec);
  return(0);
}
#endif

// routine for linking to R
void qcondig12(int *n0, double *p, double *v, double *theta, double *alpha,
  int *mxiter0, double *eps0, double *bd0, double *qu)
{ int i,n,mxiter;
  double qcondig12_algo(double, double, double, double, int, double, double, int);
  double eps,bd;
  n=*n0; mxiter=*mxiter0; eps=*eps0; bd=*bd0;
  for(i=0;i<n;i++)
  { qu[i] = qcondig12_algo(p[i],v[i],theta[i],alpha[i],mxiter,eps,bd, 0); }
}


// 0<p<1, 0<v<1, theta>0 alpha>0
double qcondig12_algo(double p, double v, double theta, double alpha,
  int mxiter, double eps, double bd, int iprint)
{ double x0, p0, diff0;
  double mindiff, x, diff, ex, g, gp;
  int i,iter;
  double pcondig12(double, double, double, double);
  double dig(double, double, double, double);
  //double interp_gen_inv_algo (double, double, double, int, double, double, int);
  //double igl_gen(double, double);
  //double igl_gen_D(double, double);
  //double prod, y, denom;
  if (p <= 0) return(0.);
  if (p >= 1) return(1.);
  //prod = alpha * theta * v * p;
  //y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd, 0);
  // prod, y and denom are unused??
  //denom = igl_gen(theta*y, alpha) - theta * igl_gen_D(theta*y, alpha);
  // x0 = c(p, 1:99/100)
  // p0 <- pcondig12(x0, v, theta = theta, alpha = alpha)
  // i0 = which(diff0 == min(diff0))[1]
  p0 = pcondig12(p, v, theta, alpha);
  mindiff = fabs(p-p0);
  x = p;
  for(i=1.; i<100; i++)
  { x0=i/100.;
    p0 = pcondig12(x0, v, theta, alpha);
    diff0 = fabs(p - p0);
    if(diff0<mindiff) { mindiff=diff0; x=x0; }
  }
  x = -log(x);
  iter = 0; diff = 1;
  while (iter < mxiter && fabs(diff) > eps)
  { ex = exp(-x);
    g = pcondig12(ex, v, theta, alpha) - p;
    gp = -dig(ex, v, theta, alpha) * ex;
    diff = g / gp;
    if (x - diff < 0.) diff = x / 2.;
    x -= diff;
    while(fabs(diff) > bd) { diff /= 2.;  x += diff; }
    iter++;
    if(iprint) printf("%d %f %f\n", iter, x, diff);
  }
  return(exp(-x));
}

// 0<u<1, 0<v<1, theta>0 alpha>0
double pcondig12(double u, double v, double theta, double alpha)
{ double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double interp_gen_D1(double, double, double);
  double eps,bd,y,pcond;
  int mxiter;
  mxiter=20; eps=1.e-12; bd=5.;
  y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd, 0);
  pcond = 1. - (1.-u) *
        interp_gen_D1(y, theta*(1.-u), alpha) /
        interp_gen_D1(y, theta, alpha);
  return(pcond);
}

// 0<u<1, 0<v<1, theta>0 alpha>0
double dig(double u, double v, double theta, double alpha)
{ double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double interp_kappa_D1(double, double, double);
  double interp_gen_D1(double, double, double);
  double eps,bd,y,pdf;
  int mxiter;
  mxiter=20; eps=1.e-12; bd=5.;
  y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd, 0);
  pdf = interp_kappa_D1(y, (1.-u)*theta, alpha) / interp_gen_D1(y, theta, alpha);
  return(pdf);
}

