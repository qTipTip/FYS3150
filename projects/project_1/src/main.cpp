#include <iostream>
#include <cmath>
#include <armadillo>
#include <string>
#include <chrono>

#include "library.h"

using namespace arma;


int main(int argc, char* argv[])
{
    vec a; // Lower diagonal
    vec b; // Main diagonal
    vec c; // Upper diagonal
    vec v; // Unknown variables

    vec v_exact; // The closed form solution
    vec b_tilde; // The right hand side
    vec grid_points; // The X_i's

    int n;
    int solver;

    if(argc > 2){
        n = atoi(argv[1]);
        solver = atoi(argv[2]);
    } else {
        std::cout << "Usage: " << argv[0] << " n (time steps) solver (1: tdma/2: lu)"<< std::endl;
        exit(1);
    }

    double step_length = 1.0 / (n-1);

    std::ostringstream stream;
    std::string file_name;
    std::string timer_file_name;

    /* Initialization */
    a.set_size(n);
    b.set_size(n);
    c.set_size(n);
    v.set_size(n);

    b_tilde.set_size(n);
    grid_points.set_size(n);
    v_exact.set_size(n);

    a.fill(-1.0);
    b.fill(2.0);
    c.fill(-1.0);

    /* Populating v_exact, b_tilde and grid_points */
    double x_i;
    double exp_const = exp(-10);
    for(int i = 0; i < n; ++i)
    {
        x_i = i * step_length;
        grid_points[i] = x_i;
        v_exact[i] = 1 - (1 - exp_const) * x_i - exp(-10*x_i);
        b_tilde[i] = source_term(x_i) * step_length * step_length;
    }

    switch(solver) {
        case 1: {
            stream << "../data/tdma_" << n << ".dat";
            auto start = std::chrono::high_resolution_clock::now();
            tdma(a, b, c, v, b_tilde, n);
            auto end = std::chrono::high_resolution_clock::now();
            auto total = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

            std::ofstream file("../data/timer_tdma.dat", ios::app);
            file << n << " " << total << "\n";
            file.close();
            max_relative_error(v, v_exact, n);
            break;
                }
        case 2: {
            stream << "../data/lu_" << n << ".dat";
            auto start = std::chrono::high_resolution_clock::now();
            lu_decomp(a, b, c, v, b_tilde, n);
            auto end = std::chrono::high_resolution_clock::now();
            auto total = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            std::ofstream file("../data/timer_lu.dat", ios::app);
            file << n << " " << total << "\n";
            file.close();
            break;
                }
    }



    file_name = stream.str();
    data_to_file(v, v_exact, grid_points, file_name, n);
}


