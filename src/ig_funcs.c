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
    printf("igl %f %f %f\n", x,alpha,iglf);
  }
  for(i=1;i<=n;i++)
  { x=i*0.1+0.91; alpha=n-i+4.4; eta=i*0.1;
    interpf = interp_gen(x,eta,alpha);
    printf("interp %f %f %f %f\n", x,eta,alpha,interpf);
  }

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  for(i=1;i<=n;i++)
  { p=i/(n+1.); alpha=n-i+0.4;
    if(i>=3) iprint=0;
    iglinv = igl_gen_inv_algo(p,alpha,mxiter,eps,bd, iprint);
    iglf = igl_gen(iglinv,alpha);
    printf("igl_gen_inv %f %f %f %f\n", p,alpha,iglinv,iglf);
  }

  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  for(i=1;i<=n;i++)
  { p=i/(n+1.); alpha=n-i+4.4; eta=i*0.1;
    if(i>=5) iprint=0;
    interpinv = interp_gen_inv_algo(p,eta,alpha,mxiter,eps,bd, iprint);
    interpf = interp_gen(interpinv,eta,alpha);
    printf("interp_gen_inv %f %f %f %f %f\n", p,eta,alpha,interpinv,interpf);
  }

  // second group of inputs
  eps = 1.e-12; bd = 5; mxiter = 20; iprint=1;
  printf("\nAnother group of inputs\n");
  for(i=1;i<=5;i++)
  { p=1.-i*.1; alpha=1.1; eta=2.5;
    interpinv = interp_gen_inv_algo(p,eta,alpha,mxiter,eps,bd, iprint);
    interpf = interp_gen(interpinv,eta,alpha);
    printf("interp_gen_inv %f %f %f %f %f\n", p,eta,alpha,interpinv,interpf);
  }

  for(j=0;j<n;j++)
  { i=j+1;; pvec[j] = i/(n+1.); avec[j]=n-i+4.4; etavec[j]=i*0.1; }
  interp_gen_inv(&n,pvec,etavec,avec,&mxiter,&eps,&bd,inv);
  for(j=0;j<n;j++) printf("%f ", inv[j]); printf("\n");
  free(pvec); free(avec); free(etavec); free(inv);
  return(0);
}
#endif


//' Generating function for the IGL copula family
//'
//' \code{igl_gen} is the function itself, and \code{igl_gen_inv} is
//' its inverse; \code{igl_gen_D} is the
//' derivative; and \code{igl_gen_DD} is the second derivative.
//'
//' Function arguments and parameters are vectorized, except
//' for the algorithms (marked by `_algo`).
//'
//' @param x Vector of values >=0 to evaluate the function at, .
//' @param p Vector of values to evaluate the inverse function at, between
//' 0 and 1 (inclusive).
//' @param alpha Parameter of the function, alpha > 0. Vectorized.
//' @examples
//' ## Some examples of evaluating the functions.
//' arg <- c(0, 0.5, 3, Inf, NA)
//' #igl_gen(arg, alpha = 1)
//' #igl_gen_D(arg, alpha = 0.2)
//' #igl_gen_D(arg, alpha = 1)
//' #igl_gen_D(arg, alpha = 2)
//' #igl_gen_inv(c(0, 0.5, 1), alpha = 0.5)
//'
//' ## Visual
//' #foo <- function(u) igl_gen_inv(u, alpha = 0.5)
//' #curve(foo)
//' @rdname igl_gen
double igl_gen(double x, double alpha)
{ double pgamma(double,double,double);
  double res,term1,term2;
  term1 = 1.-pgamma(x, alpha, 1.);
  term2 = alpha * pgamma(x, alpha+1., 1.) / x;
  res = term1+term2;
  return(res);
}

//' @rdname igl_gen
double igl_gen_D(double x, double alpha)
{ double pgamma(double,double,double);
  return( -(alpha / (x*x)) * pgamma(x, alpha+1., 1.) );
}

//' @param mxiter Maximum number of iterations to run the Newton-Raphson
//' algorithm when computing inverse. Positive integer, default 20
//' (which is used in the calling routine).
//' @param eps The Newton-Raphson algorithm for computing an inverse will
//' stop if the step size is less than this small number. 1.e-12 used
//' in the calling routine.
//' @param bd The largest acceptable step size in the Newton-Raphson
//' algorithm. Step size is reduced if it reaches this large. `bd=5` in
//' the calling routine.
//' @rdname igl_gen
double igl_gen_inv_algo (double p, double alpha, int mxiter, double eps,
                         double bd, int iprint)
{ double x1,x2,x3,p1,p2,p3,diff1,diff2,diff3;
  double x,best,diff,ex,g,gp;
  double qgamma(double,double,double);
  double igl_gen(double, double);
  double igl_gen_D(double, double);
  int iter;
  //prod = alpha * p;
  //if (is.na(prod)) return(prod)
  if (p == 0) return(DBL_MAX);  // Inf
  if (p == 1) return(0.);
  // x_start = c(
  //1 / ((1 - p) ^ (-1 / alpha) - 1)
  x1 = 1. / (pow(1.-p, -1./ alpha) - 1.);
  x2 = alpha / p;
  x3 = qgamma(p, alpha+1., 1.);
  p1 = igl_gen(x1,alpha); p2 = igl_gen(x2,alpha); p3 = igl_gen(x3,alpha);
  // p_start = igl_gen(x_start, alpha = alpha)
  diff1 = fabs(p1-p); diff2 = fabs(p2-p); diff3 = fabs(p3-p);
  // diff_start = abs(p_start - p)
  // best_start = which(diff_start == min(diff_start))[1]
  // x = x_start[best_start]
  x=x1; best=diff1;
  if(diff2<best) { x=x2; best=diff2; }
  if(diff3<best) { x=x3; best=diff3; }
  if (x == 0.) x = eps;
  x = log(x);
  iter = 0; diff = 1.;
  while (iter < mxiter && fabs(diff) > eps)
  { ex = exp(x);
    g = igl_gen(ex, alpha) - p;
    gp = igl_gen_D(ex, alpha) * ex;
    diff = g / gp;
    // if (x - diff < 0) diff = x / 2
    x -= diff;
    while (fabs(diff) > bd)
    { diff/=2.; x += diff; }
    iter++;
    if(iprint) printf("%d %f %f\n", iter, x, diff);
  }
  return(exp(x));
}

// p and alpha are vectors of length n
void igl_gen_inv(int *n0, double *p, double *alpha,
  int *mxiter0, double *eps0, double *bd0, double *inv)
{ int i,n,mxiter;
  double igl_gen_inv_algo (double, double, int, double, double, int);
  double eps,bd;
  n=*n0; mxiter=*mxiter0; eps=*eps0; bd=*bd0;
  for(i=0;i<n;i++)
  { inv[i] = igl_gen_inv_algo(p[i],alpha[i],mxiter,eps,bd, 0); }
}

#ifdef DONE
// convert to routine to link to R
// needed for pcondigl21 pcondigl12 qcondigl12 digl pigl
igl_gen_inv = function(p, alpha, mxiter = 20, eps = 1.e-12, bd = 5){
    l = vctrs::vec_size_common(p, alpha)
    if (l == 0L) return(numeric(0L))
    args = vctrs::vec_recycle_common(p = p, alpha = alpha)
    with(args, {
        x = numeric(0L)
        for (i in 1:l) {
            x[i] = igl_gen_inv_algo(
                p[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
            )
        }
        x
    })
}
#endif


//' Interpolating Functions
//'
//' These interpolating functions, denoted \eqn{H} in the vignette,
//' depend on a generating function (of a DJ copula).
//' `interp_gen()` uses the IGL generating function \eqn{\Psi_k};
//' `interp_kappa()` uses the "kappa version" of that same function.
//'
//' Appending `_inv` to the function name indicates inverse with
//' respect to the first argument. Appending `_D1` indicates
//' derivative with respect to the first argument. Function arguments
//' and parameters are vectorized, except for the algorithms (marked by
//' `_algo`).
//'
//' @param x Vector of values >=1 to evaluate the interpolating function at.
//' @param p Vector of values between 0 and 1 to evaluate the inverse function at.
//' @param eta Vector of values >0 of second argument of the
//' interpolating function.
//' @param alpha Vector of values >0 corresponding to the \eqn{alpha} parameter
//' of the IGL generating function.
//' @rdname interpolator
double interp_gen (double x, double eta, double alpha)
{ double igl_gen(double, double);
  double res;
  res = exp(-x) * igl_gen(eta*x, alpha);
  return(res);
}

//' @rdname interpolator
double interp_gen_D1 (double x, double eta, double alpha)
{ double igl_gen(double, double);
  double igl_gen_D(double, double);
  double res;
  res = -exp(-x) * (igl_gen(eta*x, alpha) - eta* igl_gen_D(eta*x, alpha));
  return(res);
}


//' @inheritParams igl_gen_inv_algo
//' @rdname interpolator
double interp_gen_inv_algo(double p, double eta, double alpha, int mxiter,
  double eps, double bd, int iprint)
{ double x1,x2,p1,p2,diff1,diff2;
  double x,diff,ex,g,gp;
  double interp_gen(double, double, double);
  double interp_gen_D1(double, double, double);
  double igl_gen_inv_algo (double, double, int, double, double, int);
  int iter;
  //prod = alpha * eta * p;
  // if (is.na(prod)) return(prod)
  if (p <= 0.) return(DBL_MAX); // Inf
  if (p >= 1.) return(0.);
  x1 = -log(p);
  x2 = igl_gen_inv_algo(p, alpha, mxiter, eps, bd, 0) / eta;
  p1 = interp_gen(x1, eta, alpha);
  p2 = interp_gen(x2, eta, alpha);
  diff1 = fabs(p1-p); diff2 = fabs(p2-p);
  x=x1;
  if(diff2<diff1) { x=x2; }
  x = log(x);
  iter = 0; diff = 1.;
  while(iter < mxiter && fabs(diff) > eps)
  { ex = exp(x);
    g = interp_gen(ex, eta, alpha) - p;
    gp = interp_gen_D1(ex, eta, alpha) * ex;
    diff = g / gp;
    // if (x - diff < 0) diff = x / 2
    x -= diff;
    while (fabs(diff) > bd)
    { diff/=2.; x += diff; }
    iter++;
    if(iprint) printf("%d %f %f\n", iter, x, diff);
  }
  return(exp(x));
}

// p, eta and alpha are vectors of length n
void interp_gen_inv(int *n0, double *p, double *eta, double *alpha,
  int *mxiter0, double *eps0, double *bd0, double *inv)
{ int i,n,mxiter;
  double interp_gen_inv_algo (double, double, double, int, double, double, int);
  double eps,bd;
  n=*n0; mxiter=*mxiter0; eps=*eps0; bd=*bd0;
  for(i=0;i<n;i++)
  { inv[i] = interp_gen_inv_algo(p[i],eta[i],alpha[i],mxiter,eps,bd, 0); }
}

#ifdef DONE
// convert to routine to link to R
// needed for pig, pcondig21, pcondig12, dig
interp_gen_inv = function(p, eta, alpha, mxiter = 40, eps = 1.e-12, bd = 5) {
    l = vctrs::vec_size_common(p, eta, alpha)
    if (l == 0L) return(numeric(0L))
    args = vctrs::vec_recycle_common(p = p, eta = eta, alpha = alpha)
    with(args, {
        x = numeric(0L)
        for (i in 1:l) {
            x[i] = interp_gen_inv_algo(
                p[i], eta = eta[i], alpha = alpha[i], mxiter = mxiter, eps = eps, bd = bd
            )
        }
        x
    })
}
#endif
