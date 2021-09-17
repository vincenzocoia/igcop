#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector qcondig12_vec(NumericVector p, NumericVector v,
                            NumericVector theta, NumericVector alpha,
                            int mxiter, double eps, double bd)
{ int i;
  int n = p.size();
  double qcondig12_algo(double, double, double, double, int, double, double, int);
  NumericVector qu(n);
  for(i=0;i<n;i++) {
    qu[i] = qcondig12_algo(p[i],v[i],theta[i],alpha[i],mxiter,eps,bd, 0);
  }
  return(qu);
}

double qcondig12_algo(double p, double v, double theta, double alpha,
                      int mxiter, double eps, double bd, int iprint)
{ double x0, p0, diff0;
  double mindiff, x, diff, ex, g, gp, prod;
  int i,iter;
  double pcondig12_single(double, double, double, double);
  double dig_single(double, double, double, double);
  //double interp_gen_inv_algo (double, double, double, int, double, double, int);
  //double igl_gen(double, double);
  //double igl_gen_D(double, double);
  //double prod, y, denom;
  prod = alpha * theta * v * p;
  if (ISNAN(prod)) return(prod);
  if (p <= 0.) return(0.);
  if (p >= 1.) return(1.);
  //y = interp_gen_inv_algo(1.-v, theta, alpha, mxiter, eps, bd, 0);
  // prod, y and denom are unused??
  //denom = igl_gen(theta*y, alpha) - theta * igl_gen_D(theta*y, alpha);
  // x0 = c(p, 1:99/100)
  // p0 <- pcondig12(x0, v, theta = theta, alpha = alpha)
  // i0 = which(diff0 == min(diff0))[1]
  p0 = pcondig12_single(p, v, theta, alpha);
  mindiff = fabs(p-p0);
  x = p;
  for(i=1.; i<100; i++)
  { x0=i/100.;
    p0 = pcondig12_single(x0, v, theta, alpha);
    diff0 = fabs(p - p0);
    if(diff0<mindiff) { mindiff=diff0; x=x0; }
  }
  x = -log(x);
  iter = 0; diff = 1;
  while (iter < mxiter && fabs(diff) > eps)
  { ex = exp(-x);
    g = pcondig12_single(ex, v, theta, alpha) - p;
    gp = -dig_single(ex, v, theta, alpha) * ex;
    diff = g / gp;
    if (diff > bd) diff = bd;
    if (diff < -bd) diff = -bd;
    if (x - diff < 0.) diff = x / 2.;
    x -= diff;
    //while(fabs(diff) > bd)
    //{ diff /= 2.;  x += diff; R_CheckUserInterrupt(); }
    iter++;
    if(iprint) Rprintf("%d %f %f\n", iter, x, diff);
    R_CheckUserInterrupt();
  }
  return(exp(-x));
}
