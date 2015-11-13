# pragma once

#include <armadillo>

class Jacobi {

  public:

    void solve(arma::mat &A, arma::mat &R, int n);
    void rotate(arma::mat &A, arma::mat &R, int k, int l, int n);

    double max_off_diag(arma::mat A, int *k, int *l, int n);
};
