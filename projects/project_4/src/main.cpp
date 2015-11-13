#include "ising.h"
#include <string>
#include <mpi.h>

int main(int argc, char *argv[])
{
   
  int world_size; 
  int world_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   
  // system as a whole
  int grid_size = atoi(argv[1]);
  int mc_cycles = atoi(argv[2]);
  
  double temperature = atof(argv[3]);
  double temperature_end = atof(argv[4]);
  double temperature_step = atof(argv[5]);
  std::string file_name = argv[6];
  
  double dT = (temperature_end - temperature) / (double) world_size;

  // each node
  double my_temp = temperature + world_rank * dT;
  double my_temp_end = my_temp + dT - temperature_step;
  
  if (world_rank == world_size - 1) { my_temp_end += temperature_step; } // edge case

  Ising solver = Ising(grid_size, file_name);

  while(my_temp <= my_temp_end) 
  {
    solver.initialize_system(my_temp);
    solver.simulate(mc_cycles);
    solver.write_to_file();
    my_temp += temperature_step;
  }

  MPI_Finalize();
  return 0;
}
