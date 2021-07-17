#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <values.h>

#ifdef MAIN
#include <malloc.h>
// sample main program for checking functions involving kappa
int main(int argc, char *argv[])
{
  double dgamma(double,double,double);
  double igl_kappa(double,double);
  double igl_kappa_D(double,double);
  double igl_kappa_inv(double,double);
  double interp_kappa(double,double,double);
  double interp_kappa_D1(double,double,double);
  double interp_kappa_inv_algo (double, double, double, int, double, double, int);
  int n,i,mxiter,iprint,j;
  double x,tem,x2, eps,bd;
  double alpha, iglk, iglk2, iglkder, iglkinv;
  double eta, interpk, interpk2, interpkder, interpkinv;
  void interp_kappa_inv (int *, double *, double *, double *, int *, double *, double *, double *);
  double *pvec,*avec,*etavec,*inv;

  n=10;
  pvec=(double *) malloc(n * sizeof(double));
  avec=(double *) malloc(n * sizeof(double));
  etavec=(double *) malloc(n * sizeof(double));
  inv=(double *) malloc(n * sizeof(double));

  for(i=1;i<=n;i++)
  { x=i/5.;
    tem=dgamma(x,x,2.);
    printf("%f %f 2 %f\n", x,x,tem);
  }
  for(i=1;i<=n;i++)
  { x=i/5.; alpha=i*.5; x2=x+.0001;
    iglk=igl_kappa(x,alpha);
    iglk2=igl_kappa(x2,alpha);
    iglkder=igl_kappa_D(x,alpha);
    iglkinv=igl_kappa_inv(iglk,alpha);
    tem=(iglk2-iglk)/0.0001;
    printf("%f %f %f %f %f %f\n", x,alpha,iglk,iglkder,tem,iglkinv);
  }

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  for(i=1;i<=n;i++)
  { x=i/5.; alpha=i*.5; x2=x+.0001; eta=i*.1;
    interpk=interp_kappa(x,eta,alpha);
    interpk2=interp_kappa(x2,eta,alpha);
    interpkder=interp_kappa_D1(x,eta,alpha);
    interpkinv=interp_kappa_inv_algo(interpk,eta,alpha,mxiter,eps,bd,iprint);
    tem=(interpk2-interpk)/0.0001;
    printf("%f %f %f %f %f %f %f\n", x,eta,alpha,interpk,interpkder,tem,interpkinv);
  }

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=0;
  for(j=0;j<n;j++)
  { i=j+1;; pvec[j] = i/(n+1.); avec[j]=n-i+4.4; etavec[j]=i*0.1;
    interpkinv=interp_kappa_inv_algo(pvec[j],etavec[j],avec[j],mxiter,eps,bd,iprint);
    printf("%f %f %f %f\n", pvec[j],etavec[j],avec[j],interpkinv);
  }
  interp_kappa_inv(&n,pvec,etavec,avec,&mxiter,&eps,&bd,inv);
  for(j=0;j<n;j++) printf("%f ", inv[j]); printf("\n");
  free(pvec); free(avec); free(etavec); free(inv);
  return(0);
}
#endif



double  dgamma(double x, double shape, double scale)
{ double xx,tem,pdf;
  xx = x/scale;
  tem = ((shape-1.)*log(xx) - lgamma(shape) - xx);
  pdf = exp(tem)/scale;
  return(pdf);
}

//' Kappa version of generating function
//'
//' The kappa version of a generating function \eqn{\psi} is
//' defined as \eqn{\kappa(x) = \psi(x) + x \psi'(x)}.
//' `igl_kappa` takes \eqn{\psi} to be the IGL generating
//' function, `igl_gen()`.
//'
//' `igl_kappa_D()` is the derivative.
//'
//' @param x Numeric argument of the kappa function. Vectorized.
//' @param p Numeric argument of the inverse function. Vectorized. Between 0 and 1.
//' @param alpha Parameter of the IGL generating function, `igl_gen()`, >0. Vectorized.
//' @rdname kappa
double igl_kappa(double x, double alpha)
{ double pgamma(double,double,double);
  double res;
  res = 1.-pgamma(x, alpha, 1.);
  return(res);
}

//' @rdname kappa
double igl_kappa_D(double x, double alpha)
{ double dgamma(double,double,double);
  double res;
  res = -dgamma(x, alpha, 1.);
  return(res);
}

//' @rdname kappa
double igl_kappa_inv(double p, double alpha)
{ double qgamma(double,double,double);
  double res;
  res = qgamma(1.-p, alpha, 1.);
  return(res);
}


double interp_kappa(double x, double eta, double alpha)
{ double igl_kappa(double, double);
  double res;
  res = exp(-x) * igl_kappa(eta*x, alpha);
  return(res);
}

double interp_kappa_D1(double x, double eta, double alpha)
{ double igl_kappa(double, double);
  double igl_kappa_D(double, double);
  double res;
  res = -exp(-x) * (igl_kappa(eta*x, alpha) - eta* igl_kappa_D(eta*x, alpha));
  return(res);
}


double interp_kappa_inv_algo(double p, double eta, double alpha, int mxiter,
   double eps, double bd, int iprint)
{ double igl_kappa(double p, double alpha);
  double igl_kappa_D(double p, double alpha);
  double igl_kappa_inv(double p, double alpha);
  double interp_kappa(double x, double eta, double alpha);
  double x1,x2, p1,p2, diff1,diff2, x, pex, diff,g,gp;
  int iter;
  //prod = alpha * eta * p;
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
  { pex = p * exp(x);
    g = igl_kappa(eta*x, alpha) - pex;
    gp = eta * igl_kappa_D(eta*x, alpha) - pex;
    diff = g / gp;
    if (x - diff < 0.) diff = x / 2.;
    x -= diff;
    while(fabs(diff) > bd) { diff /= 2;  x += diff; }
    iter++;
    if(iprint) printf("%d %f %f\n", iter, x, diff);
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
  for(i=0;i<n;i++)
  { inv[i] = interp_kappa_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0); }
}


#ifdef DONE
// convert to routine to link to R
// needed for qcondig21
interp_kappa_inv = function(p, eta, alpha, mxiter = 80, eps = 1.e-12, bd = 5) {
  l = vctrs::vec_size_common(p, eta, alpha)
  if (l == 0L) return(numeric(0L))
  args = vctrs::vec_recycle_common(p = p, eta = eta, alpha = alpha)
  with(args, {
    x = numeric(0L)
    for (i in 1:l) {
      x[i] = interp_kappa_inv_algo(
        p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
      )
    }
    x
  })
}
#endif
