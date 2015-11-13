#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
#define pi 3.14159265359

using namespace std;

/* Computes the integrand we want to integrate
 * Skips the cases where the denominator
 * is sufficiently small, could cause problems.
 */
double cartesian_function(double x1, double y1, double z1, double x2, double y2, double z2) {
  double tolerance = 1.0e-6;
  double denominator = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
  if (denominator < tolerance ) {
    return 0.0;
  }

  double alpha = 2;
  double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
  double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);


  double exp_one = -2 * alpha * r1;
  double exp_two = -2 * alpha * r2;

  return exp(exp_one + exp_two) / denominator;
}

/* Computes the integrand for spherical coordinates
 */
double spherical_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
  double numerator = sin(theta1)*sin(theta2);
  double cos_beta = (cos(theta1)*cos(theta2)) + (sin(theta1)*sin(theta2)*cos(phi1-phi2));
  double discriminant = (r1*r1) + (r2*r2) - (2*r1*r2*cos_beta);
  double denominator = sqrt(discriminant);

  // Ignoring cases where our denominator is close to zero and
  // when we have complex solutions
  double tolerance = 1.0e-6;
  if (denominator < tolerance || discriminant < 0) {
    return 0.0;
  }

  return (numerator/denominator);
}
void legendre(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
     z = cos(pi * (i - 0.25)/(n + 0.5));

     /*
      ** Starting with the above approximation to the ith root
      ** we enter the mani loop of refinement bt Newtons method.
      */

     do {
       p1 =1.0;
       p2 =0.0;

       /*
        ** loop up recurrence relation to get the
        ** Legendre polynomial evaluated at x
        */

       for(j = 1; j <= n; j++) {
         p3 = p2;
         p2 = p1;
         p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
       }

       /*
        ** p1 is now the desired Legrendre polynomial. Next compute
        ** ppp its derivative by standard relation involving also p2,
        ** polynomial of one lower order.
        */

       pp = n * (z * p1 - p2)/(z * z - 1.0);
       z1 = z;
       z  = z1 - p1/pp;                   // Newton's method
     } while(fabs(z - z1) > ZERO);

     /*
      ** Scale the root to the desired interval and put in its symmetric
      ** counterpart. Compute the weight and its symmetric counterpart
      */

     *(x_low++)  = xm - xl * z;
     *(x_high--) = xm + xl * z;
     *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
     *(w_high--) = *(w_low++);
   }
}

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
                          24.01409824083091,-1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void laguerre(double x[], double w[], int n, double alf)
{
  int i, its, j;
  double ai;
  double p1, p2, p3, pp, z, z1;
  for (i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
          (1.0+3.5*ai))*(z-x[i-3])/(1.0+0.3*alf);
    }
    for (its=1;its<=MAXIT;its++) {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= EPS) break;
    }
    if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
    x[i-1]=z;
    w[i-1] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
  }
}


/* Computes the six-dimensional integral given integration points and weights.
 * Using Gaussian-Quadrature.
 */
double legendre_integral(double *x, double *w, int N, double (*func)(double, double, double, double, double, double))
{
  int i, j, k, l, m, n;

  double total = 0;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k) {
        for (l = 0; l < N; ++l) {
          for (m = 0; m < N; ++m) {
            for (n = 0; n < N; ++n) {
              total += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] * (*func)(x[i], x[j], x[k], x[l], x[m], x[n]);
            }
          }
        }
      }
    }
  }
  return total;
}

double laguerre_integral(double *xr, double *wr, double *x_theta, double *w_theta, double *x_phi, double *w_phi, int N, double (*func)(double, double, double, double, double, double)) {
  int i, j, k, l, m, n;

  double total = 0;
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k) {
        for (l = 0; l < N; ++l) {
          for (m = 0; m < N; ++m) {
            for (n = 0; n < N; ++n) {
              total += wr[i]*wr[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[n] * (*func)(xr[i], xr[j], x_theta[k], x_theta[l], x_phi[m], x_phi[n]);
            }
          }
        }
      }
    }
  }
  return total;
}
void bruteforce_MC(double lower_lim, double upper_lim, int N)
{
  auto start = std::chrono::high_resolution_clock::now();

  double total = 0.0;
  double variance = 0.0;

  double x;
  double fx;
  double *yvalues = new double[6];
  double jacobideterminant = pow(upper_lim - lower_lim, 6);
  double exact = 5 * pi * pi / (16.0*16.0);

  srand(time(NULL));
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < 6; j++) {
      yvalues[j] = lower_lim + (double) rand()/RAND_MAX*(upper_lim - lower_lim);
    }

    fx = cartesian_function(yvalues[0], yvalues[1], yvalues[2], yvalues[3], yvalues[4], yvalues[5]);
    total += fx;
    variance += fx*fx;
  }

  double mean = total / ((double) N);
  double mean_sq = (variance / (double) N);
  variance = mean_sq - mean*mean;
  double std_dev = sqrt(variance);
  auto end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  std::cout << "Brute force MC for N = " << N << " std.dev = " << jacobideterminant*(std_dev/sqrt((double) N)) << " Integral = " << mean*jacobideterminant << " absolute error = " << fabs(mean*jacobideterminant - exact) << " time " << time / 1.0e9 << std::endl;
}

void importance_sampl_MC(int N)
{

}
