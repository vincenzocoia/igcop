/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Based on Algorithm AS 91 */

//#include "Mathlib.h"
#include <math.h>

#define MAXIT 20

#define C1	0.01
#define C2	0.222222
#define C3	0.32
#define C4	0.4
#define C5	1.24
#define C6	2.2
#define C7	4.67
#define C8	6.66
#define C9	6.73

#define C10	13.32
#define C11	60.0
#define C12	70.0
#define C13	84.0
#define C14	105.0
#define C15	120.0
#define C16	127.0
#define C17	140.0
#define C18	1175.0
#define C19	210.0

#define C20	252.0
#define C21	2264.0
#define C22	294.0
#define C23	346.0
#define C24	420.0
#define C25	462.0
#define C26	606.0
#define C27	672.0
#define C28	707.0
#define C29	735.0

#define C30	889.0
#define C31	932.0
#define C32	966.0
#define C33	1141.0
#define C34	1182.0
#define C35	1278.0
#define C36	1740.0
#define C37	2520.0
#define C38	5040.0

double qgamma(double p, double alpha, double scale)
{ double pgamma(double,double,double);
  double qnorm(double,double,double);
  double lgamma(double);
  double a, b, c, ch, g, p1, v;
  double p2, q, s1, s2, s3, s4, s5, s6, t, x, xx;
  double aa = 0.6931471806;
  double e = 0.5e-6;
  double pmin = 0.000002;
  double pmax = 0.999998;
  int i;

  /* test arguments and initialise */

  //if(p < pmin || p > pmax || alpha<=0 ) DOMAIN_ERROR;

  v = 2*alpha;
  /* xx = 0.5*v; */
  xx = alpha;
  c = xx-1.0;
  g = lgamma(0.5*v);

  if(v < (-C5)*log(p)) {
      /* starting approximation for small chi-squared */
    ch = pow(p*xx*exp(g+xx*aa),1.0/xx);
    //if(ch < e) DOMAIN_ERROR;
    if(ch < e) return(-1.);
  }
  else if(v > C3) {
    /* starting approximation using Wilson and Hilferty estimate */
    x = qnorm(p, 0.0, 1.0);
    p1 = C2/v;
    ch = v*pow(x*sqrt(p1)+1.0-p1, 3.0);
      /* starting approximation for p tending to 1 */
    if( ch>C6*v+6.0 ) ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
  }
  else {
    /* starting approximation for v less than or equal to 0.32 */
    ch = C4;
    a = log(1.0-p);
    do {
      q = ch;
      p1 = 1.0+ch*(C7+ch);
      p2 = ch*(C9+ch*(C8+ch));
      t = -0.5+(C7+2.0*ch)/p1-(C9+ch*(C10+3.0*ch))/p2;
      ch = ch-(1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
    } while(fabs(q/ch-1.0) > C1);
  }

  /* algorithm as 239 and calculation of seven term taylor series */

  for( i=1 ; i<=MAXIT  ; i++ ) {
    q = ch;
    p1 = 0.5*ch;
    p2 = p-pgamma(p1, xx, 1.0);
#ifdef HAVE_ISNAN
    //if(!finite(p2)) DOMAIN_ERROR;
    if(!finite(p2)) return(-1.);
//#else
    //if(errno != 0) DOMAIN_ERROR;
#endif

    t = p2*exp(xx*aa+g+p1-c*log(ch));
    b = t/ch;
    a = 0.5*t-b*c;
    s1 = (C19+a*(C17+a*(C14+a*(C13+a*(C12+C11*a)))))/C24;
    s2 = (C24+a*(C29+a*(C32+a*(C33+C35*a))))/C37;
    s3 = (C19+a*(C25+a*(C28+C31*a)))/C37;
    s4 = (C20+a*(C27+C34*a)+c*(C22+a*(C30+C36*a)))/C38;
    s5 = (C13+C21*a+c*(C18+C26*a))/C37;
    s6 = (C15+c*(C23+C16*c))/C38;
    ch = ch+t*(1.0+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
    if(fabs(q/ch-1.0) > e) return 0.5*scale*ch;
  }
  /* possible loss of precision */
  /* errno = EDOM; */
  return 0.5*scale*ch;
}
