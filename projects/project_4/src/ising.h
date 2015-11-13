#pragma once

#include <armadillo>
#include <string>

/* This class contains the methods required to simulate a magnetic system with
 * a discrete square lattice of arbitrary size for a range of temperatures.
 */
class Ising
{
  private:
    // FIELDS
    // =======================================
    // properties of the system and simulation
    int lattice_dimension;
    int number_of_spins;
    int number_of_accepted_states;
    int number_of_cycles;
    int write_intermediate_expectation;

    bool thermalized;

    double temperature;

    arma::mat spin_matrix;
    arma::vec average_values;
    arma::vec delta_energy;
    arma::vec energy_state_counter;

    std::string file_name;
    std::string file_name_expectations;

    // physical quantities of interest
    double mean_energy;
    double energy_variance;
    double mean_magnetization;
    double specific_heat;
    double susceptibility;
    double mean_absolute_magnetization;

    // auxilliary variables for computation
    double energy;
    double magnetization;

    // METHODS
    // =================
    // auxilliary methods
    void write_probabilities(int cycles_burned);
    void initialize_energy();
    void metropolis();

    int periodic(int x, int L, int offset); // returns periodic indices
    int thermalize(int max_burn_cycles);

    double get_energy(int x, int y);

  public:
    // METHODS
    // constructor
    Ising(int lattice_dimension, std::string file_name="simulation_results.dat", int write_intermediate_expectation=0);

    // initialization and simulation
    void simulate(int number_of_cycles);
    void initialize_system(double temperature);

    // results
    void print_results();
    void write_to_file();
    void write_intermediate_expectations(int current_cycle_count);
};
