import os
class System:

    """ This class represents the setup of a system with the following paramenters:
    
        T: initial temperature
        N_each_dimension: number of unit cells in each dimension
        b: lattice constant  
    """

    def __init__(self, T=300, N_each_dimension=5, b=5.26, ID = 1):
        self.T = T
        self.N_each_dimension = N_each_dimension
        self.b = b
        self.ID = ID
        self.dir_name = '../data/%d' % ID
        self.file_name = '%d_%04f_%04f.dat' % (N_each_dimension, T, b) + self.dir_name
         
    def run_simulations(self):
        build_path = '/Users/ivar/Documents/University/build-molecular-dynamics-fys3150-Desktop_Qt_5_5_1_clang_64bit-Release/'
        os.system(build_path + 'molecular-dynamics-fys3150 %d %f %f' % (self.N_each_dimension, self.T, self.b))
        
        try:
            os.stat(self.dir_name)
        except:
            os.mkdir(self.dir_name)
        
        os.system('mv movie.xyz %s' % self.dir_name)
        os.system('mv samples.dat %s/%s' % (self.dir_name, self.file_name))

    def read_data(self):
       pass    
if __name__ == "__main__":
    test = System(T=10, N_each_dimension=6, b=5.26)
    test.run_simulations()
