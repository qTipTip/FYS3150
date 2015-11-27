import os
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
        temperature = []
        potential = []
        kinetic = []
        total = []
        
        with open(self.file_name, 'r') as data_file:
            for line in data_file:
                line = [float(value) for value in line.split()]
                temperature.append(line[0])
                potential.append(line[1])
                kinetic.append(line[2])
                total.append(line[2] + line[1])

if __name__ == "__main__":
    """
    Simulate a set of systems for various initial temperatures (task d)
    """
    simulations = []
    for ID, T in enumerate(range(40, 200, 40)):
        simulations.append(System(T, 5, 5.26, ID, 'varying_initial_temperatures'))
    
    for sim in simulations:
        sim.read_data()
