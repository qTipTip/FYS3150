#define CATCH_CONFIG_MAIN

#include "../../test_framework/catch.hpp"
#include "../src/jacobi.h"

#include <armadillo>

TEST_CASE( "2x2 matrix rotated once zeroes out off-diagonals", "[arma::mat]" ) {

  Jacobi test_jacobi;

  arma::mat A;
  arma::mat R;

  A << 1 << 2 << arma::endr
    << 2 << 1;

  R = arma::zeros(2, 2);
  R.eye();

  test_jacobi.rotate(A, R, 0, 1, 2);

  REQUIRE ( A(0, 1) == 0 );
  REQUIRE ( A(1, 0) == 0 );
  REQUIRE ( A(0, 0) != 0 );
  REQUIRE ( A(1, 1) != 0 );

}

TEST_CASE ( "max_off_diag returns correct element and indices in a tridiagonal matrix" ) {

  Jacobi test_jacobi;

  arma::mat A;

  A << 1 << -3.7 << arma::endr
    << -3.7 << 1;

  int n = 2;
  int k;
  int l;

  double maxoff = test_jacobi.max_off_diag(A, &k, &l, n);

  REQUIRE( maxoff == 3.7 );
  REQUIRE( A(k, l) == -3.7 );

}

SCENARIO ( " Jacobi rotation preserves eigenvector orthogonality " ) {

  GIVEN ( " A 4x4 matrix with random elements " ) {

    Jacobi test_jacobi;

    arma::mat A;
    arma::mat R;

    A << 1 << 2  << 3 << 4 << arma::endr
      << 2 << 1 << 13 << 23 << arma::endr
      << 7 << 2 << 11 << 4 << arma::endr
      << 62 << 1 << 23 << 3;

    R = arma::zeros(4, 4); R.eye();

    WHEN ( "rotated at max-off-diag 100 times" ) {
      int l;
      int k;

      for(int i = 0; i < 100; ++i) {
        test_jacobi.max_off_diag(A, &k, &l, 4);
        test_jacobi.rotate(A, R, k, l, 4);
      }

      THEN ( "eigenvalues in R are orthogonal" ) {

        double tolerance = 1.0e-14;

        REQUIRE (arma::dot(R.col(0), R.col(1)) < tolerance);
        REQUIRE (arma::dot(R.col(0), R.col(2)) < tolerance);
        REQUIRE (arma::dot(R.col(0), R.col(3)) < tolerance);

        REQUIRE (arma::dot(R.col(1), R.col(2)) < tolerance);
        REQUIRE (arma::dot(R.col(1), R.col(3)) < tolerance);

        REQUIRE (arma::dot(R.col(2), R.col(3)) < tolerance);
      }
    }
  }
}

SCENARIO ( " Jacobi method for 4x4 matrix " ) {

  GIVEN ( "A 3x3 matrix with known eigenvalues ") {

    int n = 3;
    double tolerance = 1.0e-14;
    Jacobi test_jacobi;

    arma::mat A = arma::zeros(n, n);
    arma::mat R = arma::zeros(n, n);
    R.eye();

    A << 1 << 0 << 1 << arma::endr
      << 0 << 1 << 0 << arma::endr
      << 1 << 0 << 1;

    WHEN ( " eigenvalues found using the jacobi method " ) {
      test_jacobi.solve(A, R, n);

      arma::vec eigenvalues = sort(A.diag());
      THEN ( " eigenvalues correspond to exact by frobenius norm comparison" ) {

        double exact_one = 0.0;
        double exact_two = 1.0;
        double exact_three = 2.0;


        arma::vec exact_one_vec;
        arma::vec exact_two_vec;
        arma::vec exact_three_vec;
        arma::vec jaco_one;
        arma::vec jaco_two;
        arma::vec jaco_three;

        jaco_one << R(0, 0) << R(1, 0) << R(2, 0);
        jaco_two << R(0, 1) << R(1, 1) << R(2, 1);
        jaco_three << R(0, 2) << R(1, 2) << R(2, 2);


        exact_one_vec << -1 << 0 << 1;
        exact_two_vec << 0 << 1 << 0;
        exact_three_vec << 1 << 0 << 1;

        exact_one_vec = arma::normalise(exact_one_vec);
        exact_two_vec = arma::normalise(exact_two_vec);
        exact_three_vec = arma::normalise(exact_three_vec);

        REQUIRE(fabs(exact_one - eigenvalues(0)) < tolerance);
        REQUIRE(fabs(exact_two - eigenvalues(1)) < tolerance);
        REQUIRE(fabs(exact_three - eigenvalues(2)) < tolerance);


        REQUIRE(arma::norm(exact_one_vec, "fro") == arma::norm(jaco_one, "fro"));
        REQUIRE(arma::norm(exact_two_vec, "fro")== arma::norm(jaco_two, "fro"));
        REQUIRE(arma::norm(exact_three_vec, "fro") == arma::norm(jaco_three, "fro"));

      }
    }
  }
}

SCENARIO ( " Physical application of jacobi method - one electron in harmonic oscillator potential ") {

  GIVEN ( " rho_max = 5.0 , n = 100, sets up tridiag matrix A for the eigenvalue problem " ) {
    int n = 100;

    double rho_min = 0.0;
    double rho_max = 5.0;

    double step_length = (rho_max - rho_min) / (n + 1);

    double off_diag_const = -1.0 / (step_length * step_length);
    double diag_const = -2 * off_diag_const;

    arma::mat A = arma::zeros(n, n);
    arma::mat R = arma::zeros(n, n);
    arma::vec V = arma::zeros(n);
    arma::vec rho = arma::zeros(n);

    R.eye();
    A.diag(-1).fill(off_diag_const);
    A.diag(1).fill(off_diag_const);
    for (int i = 0; i < n; ++i) {
      rho(i) = rho_min + (i+1) * step_length;
      V(i) = rho(i)*rho(i);
      A.diag()(i) = V(i) + diag_const;
    }

    WHEN ( " eigenvalue problem solved using jacobi method " ){
      Jacobi solver;
      solver.solve(A, R, n);

      THEN ( " returns correct eigenvalues for the three lowest eigenstates to a given tolerance ")
      {
        arma::vec eigenvalues = sort(A.diag());
        double lambda_one = eigenvalues(0);
        double lambda_two = eigenvalues(1);
        double lambda_three = eigenvalues(2);

        double exact_one = 3.0;
        double exact_two = 7.0;
        double exact_three = 11;

        double tolerance = 1.0e-2;

        REQUIRE( fabs(lambda_one - exact_one) < tolerance);
        REQUIRE( fabs(lambda_two - exact_two) < tolerance);
        REQUIRE( fabs(lambda_three - exact_three) < tolerance);
      }
    }
  }
}
