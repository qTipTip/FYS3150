#include "ising.h"
#include <iomanip>
#include <time.h>

Ising::Ising(int lattice_dimension, std::string file_name, int intermediate_expectations)
{
  this->lattice_dimension = lattice_dimension;
  this->file_name = file_name;
  this->file_name_expectations = "expectations.dat";
  write_intermediate_expectation = intermediate_expectations == 1;

  thermalized = false;
  number_of_spins = lattice_dimension * lattice_dimension;

  spin_matrix = arma::mat(lattice_dimension, lattice_dimension);
  delta_energy = arma::vec(17);
  average_values = arma::vec(5);
  energy_state_counter = arma::vec(4 * number_of_spins + 1);

  energy = 0;
  magnetization = 0;

  mean_energy = 0;
  energy_variance = 0;
  mean_magnetization = 0;
  specific_heat = 0;
  susceptibility = 0;
  mean_absolute_magnetization = 0;

};

/*
 * Initializes the system with a given temperature. All average values are
 * reset to zero and the initial energy and magnetization is computed.  Based
 * on the temperature the initial configuration is either ground state, or
 * random. We also compute a lookup-table for the delta-energies, in order to
 * save significant time when simulating.
 *
 * If the system is already thermalized, this indicates that we can use the
 * previous grid configuration.
 */
void Ising::initialize_system(double temp)
{
  temperature = temp;
  number_of_accepted_states = 0;

  // energy lookup table
  for (int i = -8; i <= 8; i++) delta_energy(i+8) = 0;
  for (int i = -8; i <= 8; i += 4) delta_energy(i+8) = exp(-i/temp);
 
  if (!thermalized) {
    std::cout << "not thermalized!" << std::endl;
    // reset averages
    average_values.zeros();
    energy_state_counter.zeros();
    // initial grid_configuration and magnetization
    energy = 0;
    magnetization = 0;
    number_of_accepted_states = 0;
    for (int i = 0; i < lattice_dimension; i++)
    {
      for (int j = 0; j < lattice_dimension; ++j)
      {
        // if the temperature is sufficiently low, we might as well
        // start in a ground state. Otherwise, start in an random configuration;
        if (temperature < 1.5)
        {
          spin_matrix(i, j) = 1;
        } else {
          spin_matrix(i, j) = (rand() * 1. / RAND_MAX < 0.5) ? 1 : -1;
        }
        magnetization += spin_matrix(i, j);
      }
    }

    // initial energy
    for (int i = 0; i < lattice_dimension; ++i)
    {
      for (int j = 0; j < lattice_dimension; ++j)
      {
        energy -= spin_matrix(i, j) *(spin_matrix(periodic(i, lattice_dimension, -1), j) + spin_matrix(i, periodic(j, lattice_dimension, -1)));
      }
    }
  }
  // seed RNG
  srand(time(NULL));
};

/*
 * Starts a simulation of the Ising model for N Monte-Carlo cycles.  The method
 * first undergoes a burn-phase where a number of Monte-Carlo cycles are burned
 * in order to let the system thermalize and stabilize before we start
 * computing averages. Every 100'th Monte-Carlo cycle the expectation-values
 * are written out to a separate file in order to examine expectation values as
 * a function of the number of cycles.
 */
void Ising::simulate(int N)
{
  number_of_cycles = N;
  
  int cycles_burned = thermalize((int) 0.1 * N);
  for (int i = 0; i < N; ++i)
  {
    metropolis();
    average_values(0) += energy;
    average_values(1) += energy*energy;
    average_values(2) += magnetization;
    average_values(3) += magnetization*magnetization;
    average_values(4) += fabs(magnetization);
    energy_state_counter(energy + number_of_spins*2) += 1;

    if (write_intermediate_expectation && i != 0 && i % 100 == 0) {write_intermediate_expectations(i);}
  }

  // we divide by the total number of cycles in order to get averages
  double norm = 1 / ((double) N);
  double Eaverage = average_values[0]*norm;
  double E2average = average_values[1]*norm;
  double Maverage = average_values[2]*norm;
  double M2average = average_values[3]*norm;
  double Mabsaverage = average_values[4]*norm;

  // compute the variances
  double Evariance = (E2average - Eaverage*Eaverage) / number_of_spins;
  double Mvariance = (M2average - Maverage*Maverage) / number_of_spins;
  double M2variance = (M2average - Mabsaverage) / number_of_spins;

  // compute the physical quantitites
  mean_energy = Eaverage / number_of_spins;
  mean_magnetization = Maverage / number_of_spins;
  specific_heat = Evariance / (temperature * temperature);
  susceptibility = Mvariance / temperature;
  mean_absolute_magnetization = Mabsaverage / number_of_spins;
  energy_variance = Evariance;

  write_probabilities(N);
}

/*
 * Returns the index of the requested neighbour under the assumption that we
 * are operating with periodic boundary conditions.
 */
int Ising::periodic(int x, int L, int offset)
{
  return (x + L + offset) % L;
}

/*
 * Performs what constitutes one single monte-carlo cycle.  Picks a random
 * lattice-site: examines the change in energy if we were to flip one spin, if
 * the metropolis condition is met, it flips the spin and updates values
 * accordingly.  Otherwise, it retains the previous configuration.  Also
 * increments the variable counting the number of accepted states.
 */
void Ising::metropolis()
{
  for (int i = 0; i < number_of_spins; i++) {
    //TODO: this way of chosing random numbers is biased in the cases where
    // RAND_MAX is not a multiple of lattice_dimension - fix this;
    int rand_x = rand() % lattice_dimension;
    int rand_y = rand() % lattice_dimension;

    double delta_e = get_energy(rand_x, rand_y);
    double r = rand() * 1. / RAND_MAX;
    if(r <= delta_energy(delta_e + 8))
    {
      spin_matrix(rand_x, rand_y) *= -1;
      magnetization += 2 * spin_matrix(rand_x, rand_y);
      energy += delta_e;
      number_of_accepted_states++;
    }
  }
};

/*
 * Computes the energy at lattice site (x, y)
 */
double Ising::get_energy(int x, int y)
{
  double up = spin_matrix(x, periodic(y, lattice_dimension, 1));
  double down = spin_matrix(x, periodic(y, lattice_dimension, -1));
  double left = spin_matrix(periodic(x, lattice_dimension, -1), y);
  double right = spin_matrix(periodic(x, lattice_dimension, 1), y);
  double self = spin_matrix(x, y);
  return 2 * self * (up + down + left + right);
}

/*
 * appends results to file in the following format:
 * temperature mean_energy specific_heat susceptibility mean_absolute_magnetization accepted_states
 * File name is given in constructor
 */
void Ising::write_to_file()
{
  using namespace std;
  ofstream ofile;
  ofile.open(file_name, ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(20) << setprecision(8) << temperature;

  ofile << setw(20) << setprecision(8) << mean_energy;
  ofile << setw(20) << setprecision(8) << specific_heat;
  ofile << setw(20) << setprecision(8) << susceptibility;
  ofile << setw(20) << setprecision(8) << mean_absolute_magnetization;
  ofile << setw(20) << setprecision(8) << number_of_accepted_states / (double) number_of_cycles << endl;
  ofile.close();
}

/*
 * Computes and appends the intermediate expectation values to file.
 */
void Ising::write_intermediate_expectations(int current_cycle)
{
  double norm = 1 / ((double) current_cycle);

  double Eaverage = average_values[0]*norm;
  double E2average = average_values[1]*norm;
  double Maverage = average_values[2]*norm;
  double M2average = average_values[3]*norm;
  double Mabsaverage = average_values[4]*norm;

  double Evariance = (E2average - Eaverage*Eaverage) / number_of_spins;
  double Mvariance = (M2average - Maverage*Maverage) / number_of_spins;
  double M2variance = (M2average - Mabsaverage*Mabsaverage) / number_of_spins;

  mean_energy = Eaverage / number_of_spins;
  mean_magnetization = Maverage / number_of_spins;
  specific_heat = Evariance / (temperature * temperature);
  susceptibility = Mvariance / temperature;
  mean_absolute_magnetization = Mabsaverage / number_of_spins;

  using namespace std;
  ofstream ofile;
  ofile.open(file_name_expectations, ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(20) << setprecision(8) << temperature;
  ofile << setw(20) << setprecision(8) << mean_energy;
  ofile << setw(20) << setprecision(8) << mean_absolute_magnetization << endl;
  ofile.close();
}

/*
 * Writes the probability of encountering each energy_state in a simulation.
 */
void Ising::write_probabilities(int cycles_burned)
{
  using namespace std;
  ofstream ofile;

  ofile.open("energy_probability.dat", ios::app);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for (int i = 0; i < 2*number_of_spins + 1; ++i) {
    ofile << setw(20) << setprecision(8) << i - number_of_spins;
    ofile << setw(20) << setprecision(8) << energy_state_counter(i) / (double) number_of_cycles  << endl;
  }
  ofile.close();
}

/*
 * Computes variances for the system until the system stabilizes (thermalizes). Returns the number of cycles
 * burned in the thermalization. We zero out anything modified by the metropolis() function afterwards.
 */
int Ising::thermalize(int max_burn_cycles)
{
  //TODO: Write this properly, where it stops based on difference in variance. Currently it only burns a set number of cycles.
  double initial_energy = energy;
  double initial_magnetization = magnetization;
  
  for (int i = 0; i < max_burn_cycles; ++i)
  {
    metropolis();
  }

  // reset values
  number_of_accepted_states = 0;
  average_values.zeros();
  thermalized = true;
  return max_burn_cycles;
}
