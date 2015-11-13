#include "library.h"

#include <chrono>
#include <cmath>
#include <iostream>

using namespace arma;

double source_term(double x_i)
{
    return 100 * exp(-10 * x_i);
}

void tdma(vec &a, vec &b, vec &c, vec &v, vec &b_tilde, int n)
{
    vec gamma;
    vec lambda;

    gamma.set_size(n);
    lambda.set_size(n);

    v[0] = 0;
    v[n-1] = 0;

    lambda[1] = b[1];
    gamma[1] = c[1] / lambda[1];
    v[1] = b_tilde[1] / lambda[1];

    /* Forward sweep */
    for(int i = 1; i < n-1; ++i)
    {
        lambda[i] = b[i] - a[i]*gamma[i-1];
        gamma[i] = c[i] / lambda[i];
        v[i] = (b_tilde[i] - a[i]*v[i-1]) / lambda[i];
    }

    /* Backward sweep */

    for(int i = n-2; i > 0; --i)
    {
        v[i] -= gamma[i] * v[i+1];
    }

    return;
}

void lu_decomp(vec &a, vec &b, vec &c, vec &v, vec &b_tilde, int n)
{

    mat A = zeros(n, n);
    mat L, U, P;
    vec y = zeros(n);

    // Filling the diagonals of A
    A.diag(-1) += -1;
    A.diag(0) += 2;
    A.diag(1) += -1;

    lu(L, U, P, A);

    // Cannot get the solve() method to work, so
    // implementing solution by hand.
    // Forward sweep
    for (int i = 1; i < n-1; ++i) {
      y(i) = b_tilde(i);
      for (int j = 1; j < i; ++j) {
        y(i) = y(i) - L(i-1, j-1) * y(j);
      }
    }

    // Backward sweep
    for (int i = n-2; i > 0; --i) {
      v(i) = y(i);
      for (int j = i+1; j < n-1; ++j) {
        v(i) = v(i) - U(i-1, j-1)*v(j);
      }
      v(i) = v(i) / U(i-1, i-1);
    }
}

void data_to_file(vec v, vec v_exact, vec grid_points, std::string file_name, int n)
{
    std::cout << "\twriting data to " << file_name << std::endl;
    std::ofstream file(file_name, ios::out);
    for(int i = 0; i < n; ++i){
        file << grid_points[i] << " " << v[i] << " " << v_exact[i] << " " << "\n";
    }
    file.close();
}

void max_relative_error(vec& v, vec& exact, int n)
{
  std::cout << "\twriting relative error to ../data/relative_error.dat" << std::endl;
  std::ofstream file("../data/relative_error.dat", ios::app);
  double eps = -100; // our final epsilon value, starting at -100 because we are comparing the logarithms of small values which tend to be negative.
  for (int i = 1; i < n-1; ++i) {
    double temp = log10(std::abs((v(i) - exact(i))/exact(i)));
    if(temp > eps)
    {
      eps = temp;
    }
  }
  file << n << " " << eps << "\n";
  file.close();
}
