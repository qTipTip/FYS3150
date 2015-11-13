#include "jacobi.h"
#include <chrono>

void Jacobi::rotate(arma::mat &A, arma::mat &R, int k, int l, int n) {

  /* Rotates around the (k, l) plane in order to zero out elements.
   * Also updates the eigenvectors */

  double c;
  double s;
  double t;
  double tau;

  if ( A(k, l) != 0.0 ) {
    tau = ( A(l, l) - A(k, k) ) / ( 2*A(k, l) );
    t = ((tau > 0) ? 1.0 / (tau + sqrt(1.0 + tau*tau)) : -1.0 / (-tau + sqrt(1.0 + tau*tau)));

    c = 1.0 / sqrt(1 + t*t);
    s = c*t;
  } else {
    c = 1.0;
    s = 0.0;
  }

  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;

  a_kk = A(k, k);
  a_ll = A(l, l);

  A(k, k) = c*c*a_kk - 2.0*c*s*A(k, l) + s*s*a_ll;
  A(l, l) = s*s*a_kk + 2.0*c*s*A(k, l) + c*c*a_ll;
  A(k, l) = 0.0;
  A(l, k) = 0.0;

  for (int i = 0; i < n; ++i) {
    if (i != k && i != l) {
      a_ik = A(i, k);
      a_il = A(i, l);
      A(i, k) = c*a_ik - s*a_il;
      A(k, i) = A(i, k);
      A(i, l) = c*a_il + s*a_ik;
      A(l, i) = A(i, l);
    }

    r_ik = R(i, k);
    r_il = R(i, l);
    R(i, k) = c*r_ik - s*r_il;
    R(i, l) = c*r_il + s*r_ik;
  }

  return;
}

double Jacobi::max_off_diag(arma::mat A, int *k, int *l, int n) {
  /*
   * returns the absolute value of the greatest element in absolute value
   * as well as updating the pointers k and l to point to the correct
   * indices. Assumes the matrix to be symmetric along the diagonal.
   */
  double max = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (fabs(A(i, j)) > max) {
        max = fabs(A(i, j));
        *l = i;
        *k = j;
      }
    }
  }
  return max;
}

void Jacobi::solve(arma::mat &A, arma::mat &R, int n) {
  /* Performs the jacobi rotation until max number of iterations has been reached,
   * or the max_off_diagonal is below some tolerance.
   */

  int k;
  int l;

  int current_iteration_count = 0;

  double epsilon = 1.0e-14;
  double max_number_of_iterations = (double) n * (double) n * (double) n;

  double maxdiag = max_off_diag(A, &k, &l, n);
 
  auto start = std::chrono::high_resolution_clock::now();
  while ((fabs(maxdiag)*fabs(maxdiag)) > epsilon && current_iteration_count <
      max_number_of_iterations){
    maxdiag = max_off_diag(A, &k, &l, n);
    rotate( A, R, k, l, n);
    current_iteration_count++;
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto total = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1.0e9;
  std::cout << "Time elapsed: " << total << " s" << std::endl;
  std::cout << "Iterations used: " << current_iteration_count << std::endl;
}
