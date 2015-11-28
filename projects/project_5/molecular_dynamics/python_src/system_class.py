import os
import matplotlib.pyplot as plt
class System:

    """ This class represents the setup of a system with the following paramenters:
    
        T: initial temperature
        N_each_dimension: number of unit cells in each dimension
        b: lattice constant  
    """

    def __init__(self, T=300, N_each_dimension=5, b=5.26, ID = 1, simulation_description='example'):
        self.T = T
        self.N_each_dimension = N_each_dimension
        self.b = b
        self.ID = ID
        self.dir_name = '../data/%s/%d' % (simulation_description, ID)
        self.file_name = self.dir_name + '/%d_%04f_%04f.dat' % (N_each_dimension, T, b)
         
    def run_simulations(self):
        build_path = '/Users/ivar/Documents/University/build-molecular-dynamics-fys3150-Desktop_Qt_5_5_1_clang_64bit-Release/'
        os.system(build_path + 'molecular-dynamics-fys3150 %d %f %f' % (self.N_each_dimension, self.T, self.b))
        
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
    Simulate a set of systems for various initial temperatures (task d)
    """
    simulations = []
    for ID, T in enumerate(range(40, 200, 40)):
        simulations.append(System(T, 5, 5.26, ID, 'varying_initial_temperatures'))
    
    for sim in simulations:
        sim.read_data()
    
    plt.title('Energy conservation')
    for sim in simulations:
        plt.subplot(311)
        plt.plot(sim.kinetic, label='$T = %0.2f$' % sim.T)
        plt.xlabel='$t$'
        plt.ylabel='$E_k [eV]$'
        plt.legend(loc='best')
        plt.subplot(312)
        plt.plot(sim.potential, label='$T = %0.2f$' % sim.T)
        plt.xlabel='$t$'
        plt.ylabel='$U [eV]$'
        plt.legend(loc='best')
        plt.subplot(313)
        plt.plot(sim.total, label='$T = %0.2f$' % sim.T)
        plt.xlabel='$t$'
        plt.ylabel='$E_{\\mathrm{tot}} [eV]$'
        plt.legend(loc='best')
    plt.grid('on')
    plt.show()
