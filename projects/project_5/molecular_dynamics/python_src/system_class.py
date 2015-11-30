import os
import matplotlib.pyplot as plt
import numpy as np
class System:

    """ This class represents the setup of a system with the following paramenters:
    
        T: initial temperature
        N_each_dimension: number of unit cells in each dimension
        b: lattice constant  
    """

    def __init__(self, T=300, dt=5.0e-15, N_each_dimension=5, b=5.26, ID = 1, simulation_description='example'):
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
        self.diffusion = []
        
        with open(self.file_name, 'r') as data_file:
            for line in data_file:
                line = [float(value) for value in line.split()]
                self.temperature.append(line[0])
                self.potential.append(line[1])
                self.kinetic.append(line[2])
                self.total.append(line[2] + line[1])
                self.diffusion.append(line[4])
        self.temperature = np.array(self.temperature)
        self.potential = np.array(self.potential)
        self.kinetic = np.array(self.kinetic)
        self.total = np.array(self.total)
        self.diffusion = np.array(self.diffusion)
