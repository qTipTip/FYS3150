import matplotlib.pyplot as plt
import numpy as np
import os

class wavefunction:
    """
    A class containing the methods needed for
    reading data and normalizing and plotting the wave-function.
    Takes file_name of solutions as well as number of steps as input.
    """

    def __init__(self, file_name, n, rho_max, omega):
        self.solution_file = file_name
        self.n = n
        self.rho_max = rho_max
        self.omega = omega
        self.h = rho_max / float(n-1)

    def get_results(self):
        """
        extracts three lowest eigenvalues, eigenvectors found in solution_file as well as rho.
        """
        self.rho = np.zeros(self.n)
        self.eigenvalues = np.zeros(self.n)
        self.eigenvectors = np.zeros((3, self.n))
        self.eigenvalues = np.zeros(self.n)

        with open(self.solution_file, 'r') as infile:
            eigenvalues = infile.readline().split()
            self.eigenvalues[0] = float(eigenvalues[1])
            self.eigenvalues[1] = float(eigenvalues[2])
            self.eigenvalues[2] = float(eigenvalues[3])

            for i, line in enumerate(infile):
                line = line.split()
                self.rho[i] = float(line[0])
                self.eigenvectors[0][i] = float(line[1])
                self.eigenvectors[1][i] = float(line[2])
                self.eigenvectors[2][i] = float(line[3])

    def normalize(self):
        """
        normalizes the wave functions
        """
        for eigenvector in self.eigenvectors:
            eigenvector = eigenvector / np.sqrt(self._trapezoidal(eigenvector**2))

    def _trapezoidal(self, vector):
        """
        Computes the area
        """
        total = 0
        for i in range(self.n-1):
            total += vector[i+1] + vector[i]
            return (self.h / 2.0) * total

    def plot_n_lowest_states(self, n):
        """
        Plots the n lowest eigenstates, in our case we only
        store the information for the three lowest.
        Saves the plot based on the string representation of the object
        """
        for i in range(n):
            plt.plot(self.rho, abs(self.eigenvectors[i])**2, label='$\\lambda_%d$' % i)
        plt.legend()
        plt.title('$\\omega_r = %g$' % self.omega, y=1.04, fontsize=20)
        plt.xlabel('$\\rho$', fontsize=20)
        plt.ylabel('$|\\psi(\\rho)|^2$', fontsize=20)
        plt.savefig(str(self))
        plt.clf()

    def __str__(self):
        return '../plots/schrodinger_%d_%d_%g.pdf' % (self.n, self.rho_max, self.omega)



if __name__ == '__main__':

    n = 500
    omega_r = [0.01, 0.5, 1, 5]
    rho_max_values = [100, 10, 7, 5]

    wavefunctions = []
    for omega, rho_max in zip(omega_r, rho_max_values):
        file_name = '../data/schrodinger_%d_%d_' % (n, rho_max) + str(omega) + ".dat"
        wavefunctions.append(wavefunction(file_name, n, rho_max, omega))
        os.system("../build/schrodinger %d %d %f" % (n, rho_max, omega))

    for psi in wavefunctions:
        psi.get_results()
        psi.plot_n_lowest_states(3)
