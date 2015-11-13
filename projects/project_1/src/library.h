# pragma once

#include <armadillo>
#include <string>


double source_term(double);
void tdma(arma::vec& a, arma::vec& b, arma::vec& c, arma::vec& v, arma::vec& b_tilde, int n);
void lu_decomp(arma::vec& a, arma::vec& b, arma::vec& c, arma::vec& v, arma::vec& b_tilde, int n);
void data_to_file(arma::vec v, arma::vec v_exact, arma::vec grid_points, std::string file_name, int n);
void max_relative_error(arma::vec& v, arma::vec& exact, int n);
