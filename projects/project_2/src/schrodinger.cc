#include "jacobi.h"
#include <armadillo>
#include <iostream>
#include <map>
#include <string>

using namespace arma;

double potential(double omega_r, double rho);
void write_to_file(int n, mat A, mat R, vec rho, double omega_r, std::string file_name);

int main(int argc, char *argv[])
{

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " n rho_max omega_r" << std::endl;
  }
  int n = atoi(argv[1]);

  double rho_min = 0.0;
  double rho_max = atof(argv[2]);
  double omega_r = atof(argv[3]);

  double step_length = (rho_max - rho_min) / (n + 1);

  double off_diag_const = -1.0 / (step_length * step_length);
  double diag_const = -2 * off_diag_const;

  mat A = zeros(n, n);
  mat R = zeros(n, n);
  vec V = zeros(n);
  vec rho = zeros(n);

  R.eye();
  A.diag(-1).fill(off_diag_const);
  A.diag(1).fill(off_diag_const);
  for (int i = 0; i < n; ++i) {
    rho(i) = rho_min + (i+1) * step_length;
    V(i) = potential(omega_r, rho(i));
    A.diag()(i) = V(i) + diag_const;
  }

  Jacobi solver;
  solver.solve(A, R, n);

  std::ostringstream stream;
  std::string file_name;

  stream << "../data/schrodinger_" << n << "_" << rho_max << "_" << omega_r << ".dat";

  file_name = stream.str();

  write_to_file(n, A, R, rho, omega_r, file_name);
  return 0;
}

inline double potential(double omega_r, double rho) {
  return omega_r*omega_r * rho*rho + (1.0 / rho);
}

void write_to_file(int n, mat A, mat R, vec rho, double omega_r, std::string file_name) {

  std::map<double, int> lambdamap;

  for (int i = 0; i < n; ++i) {
    lambdamap[A.diag(0)(i)] = i;
  }

  vec eigenvalues = sort(A.diag());
  vec eigenvector_one = R.col(lambdamap[eigenvalues(0)]);
  vec eigenvector_two = R.col(lambdamap[eigenvalues(1)]);
  vec eigenvector_three = R.col(lambdamap[eigenvalues(2)]);

  std::ofstream file(file_name, std::ios::out);
  file << "lambda " << eigenvalues(0)<< " " << eigenvalues(1)<< " " << eigenvalues(2) << "\n";
  for (int i = 0; i < n; ++i) {
    file << rho(i) << " " << eigenvector_one(i) << " " << eigenvector_two(i)  << " " <<  eigenvector_three(i) << "\n";
  }
  file.close();
  return;
}
