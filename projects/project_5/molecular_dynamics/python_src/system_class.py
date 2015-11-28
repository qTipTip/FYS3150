import os
import matplotlib.pyplot as plt
import numpy as np
class System:

    """ This class represents the setup of a system with the following paramenters:
    
        T: initial temperature
        N_each_dimension: number of unit cells in each dimension
        b: lattice constant  
    """

    def __init__(self, T=300, dt=1.0e-15, N_each_dimension=5, b=5.26, ID = 1, simulation_description='example'):
        self.T = T
        self.N_each_dimension = N_each_dimension
        self.b = b
        self.ID = ID
        self.dt = dt
        self.dir_name = '../data/%s/%d' % (simulation_description, ID)
        self.file_name = self.dir_name + '/%d_%04f_%04f.dat' % (N_each_dimension, T, b)
         
    def run_simulations(self, integration_method=2):
        build_path = '/Users/ivar/Documents/University/build-molecular-dynamics-fys3150-Desktop_Qt_5_5_1_clang_64bit-Release/'
        os.system(build_path + 'molecular-dynamics-fys3150 %d %f %f %d %.10g' % (self.N_each_dimension, self.T, self.b, integration_method, self.dt))
        
        try:
            os.stat(self.dir_name)
        except:
            os.makedirs(self.dir_name)
        
        os.system('mv movie.xyz %s' % self.dir_name)
        os.system('mv samples.dat %s' % (self.file_name))

    def read_data(self):
        self.temperature = []
        self.potential = []
        self.kinetic = []
        self.total = []
        
        with open(self.file_name, 'r') as data_file:
            for line in data_file:
                line = [float(value) for value in line.split()]
                self.temperature.append(line[0])
                self.potential.append(line[1])
                self.kinetic.append(line[2])
                self.total.append(line[2] + line[1])

if __name__ == "__main__":
    """
    TASK D
    """
    simulations = []
    initial_temp = 300 
    b = 5.26
    N = 5
    description = 'energy_fluctuations_'
    
    dt_values = 10**np.linspace(-15, -14, 40)
    for i, dt in enumerate(dt_values):
        
        euler = System(T = initial_temp, dt=dt, N_each_dimension=5, b=5.26, ID=i, simulation_description=description+'euler')
        verlet = System(T = initial_temp, dt=dt, N_each_dimension=5, b=5.26, ID=i, simulation_description=description+'verlet')
        
        simulations.append([euler, verlet])
    
    # run simulations 
    for sim in simulations:
        print sim[0].ID
        #sim[0].run_simulations(1) # euler
        #sim[1].run_simulations(2) # verlet
    
    # extract data
    euler_sd = []
    verlet_sd = []
    for sim in simulations:
        sim[0].read_data()
        sim[1].read_data()
    
        euler_sd.append(np.std(sim[0].total))
        verlet_sd.append(np.std(sim[1].total))

    # plot data
    plt.plot(dt_values, euler_sd, 'ro-', dt_values, verlet_sd, 'go-')
    plt.legend(['Euler-Cromer', 'Velocity-Verlet'], loc='best')
    plt.xlabel('$\\Delta t [s]$')
    plt.ylabel('$\\sigma_E/N_{\\mathrm{atoms}} [eV]$')
    plt.grid('on')
    plt.savefig('energy_fluctuations.pdf')
