#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "integration_methods.h"
#define pi 3.14159265359
// Function prototypes
int main(int argc, char *argv[])
{
  int n;
  int method_id;
  double a;
  double b;

  // Retrieving parameters
  if (argc < 5) {
    std::cout << "Usage: " << argv[0] << " N method "<< std::endl;
    std::cout << "Parameters\n" << "\t N - integration steps\n \t method - 1/2/3 { Legendre/Laguerre/MonteCarlo }\n" << "a, b - integration limits" << std::endl;
    exit(1);
  } else {
    n = atoi(argv[1]);
    method_id = atoi(argv[2]);
    a = atof(argv[3]);
    b = atof(argv[4]);
  }

  // Initializing empty arrays for weights and corresponding integration points
  // as well as integration limits
  double *x = new double[n];
  double *w = new double[n];
  double *x_theta = new double[n];
  double *w_theta = new double[n];
  double *x_phi = new double[n];
  double *w_phi = new double[n];

  switch (method_id) {
    case 1:
    {
      // Using Gauss-Legendre quadrature
      auto start = std::chrono::high_resolution_clock::now();

      legendre(a, b, x, w, n);
      double result = legendre_integral(x, w, n, &cartesian_function);
      double exact = 5 * pi * pi / (16.0*16.0);
      auto end = std::chrono::high_resolution_clock::now();
      auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      std::cout << "Gauss Legendre for " << n << " steps: " << result << " absolute error: " << fabs(result - exact) << " time: " << time/1.0e9 << " s" << std::endl;
      std::ofstream outfile;
      outfile.open("../data/gauss_legendre.dat", std::ios::app);
      outfile << n <<  " " << result << " " << time << " " << a << " " << b << "\n";
      outfile.close();
      break;
    }
    case 2:
    {
      auto start = std::chrono::high_resolution_clock::now();

      double exact = 5 * pi * pi / (16.0*16.0);
      // Using Gauss-Laguerre quadrature
      laguerre(x, w, n, 2);
      legendre(0, pi, x_theta, w_theta, n);
      legendre(0, 2*pi, x_phi, w_phi, n);
      double result = 1/1024. * laguerre_integral(x, w, x_theta, w_theta, x_phi, w_phi, n, &spherical_function);
      auto end = std::chrono::high_resolution_clock::now();
      auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      std::cout << "Gauss Laguerre for " << n << " steps: " << result << " absolute error: " << fabs(result - exact) << " time: " << time/1.0e9 << " s" << std::endl;
      break;
    }
    case 3:
      // Using MonteCarlo simulation
      bruteforce_MC(a, b, n);
      break;
    case 4:
      // Using importance sampling Monte Carlo
      importance_sampl_MC(n);
      break;
    default:
      // No valid method specified
      std::cout << "No valid method specified" << std::endl;
      std::cout << "Usage: " << argv[0] << " N method "<< std::endl;
      std::cout << "Parameters\n" << "\t N - integration steps\n \t method - 1/2/3 { Legendre/Laguerre/MonteCarlo }\n" << "a, b - integration limits" << std::endl;
      exit(1);
  }
  return 0;
}
