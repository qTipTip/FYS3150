import matplotlib.pyplot as plt
from numpy import vectorize, array, sqrt
from numpy.linalg import norm
import os

def reader(filename):
    """
    extracts rho and u(x) given datafile
    """

    with open(filename, 'r') as infile:
        rho = [0]
        u_one = [0]
        u_two = [0]
        u_three = [0]
        for line in infile:
            line = line.split()
            rho.append(float(line[0]))
            u_one.append(float(line[1]))
            u_two.append(float(line[2]))
            u_three.append(float(line[3]))

    return array(rho), array(u_one), array(u_two), array(u_three)

def trapezoidal(step_length, u):

    area = 0
    for i in range(len(u)-1):
        area += u[i+1] + u[i]
    return step_length * 0.5 * area

def normalize(step_length, u):
    """
    normalizes the vector u by dividing by dividing by the integral of u**2 computed by
    the trapezoidal rule
    """

    u /= sqrt(trapezoidal(step_length, u ** 2))

    return u


if __name__ == '__main__':

    """
    Computes the lowest ground state for various omega_r
    """
    omega_r = [0.01, 0.5, 1, 5]
    rho_max = 10
    n = 200
    h = rho_max / float(n-1)

    # Running simulations
    for omega in omega_r:
        os.system("../build/schrodinger %d %d %f" % (n, rho_max, omega))

    # Reading data and plotting lowest values for varying omega
    for omega in omega_r:
        file_name = "../data/schrodinger_%d_%d_" % (n, rho_max) + str(omega) + ".dat"
        rho, u = reader(file_name)[0:2]
        u = abs(normalize(h, u)) ** 2
        plt.plot(rho, u, label='$\omega_r = %.2f$' % omega)

    plt.legend()
    plt.xlabel('$\\rho$')
    plt.ylabel('$|\\psi(\\rho)|^2$')
    plt.savefig('../plots/schrodinger_%d_%d_' % (n, rho_max) + str(omega) + '.pdf')

    """
    Computing the three lowest eigenstates for various omega_r and appropriate rho_max
    rho_max_values = [50, 7, 5, 2]

    for omega, rho_max in zip(omega_r, rho_max_values):
        os.system("../build/schrodinger %d %d %f" % (n, rho_max, omega))

    for omega, rho_max in zip(omega_r, rho_max_values):
        file_name = "../data/schrodinger_%d_%d_" % (n, rho_max) + str(omega) + ".dat"
        rho, u_one, u_two, u_three = reader(file_name)
        u_one = abs(normalize(h, u_one)) ** 2
        u_two = abs(normalize(h, u_two)) ** 2
        u_three = abs(normalize(h, u_three)) ** 2
        plt.plot(rho, u_one, rho, u_two, rho, u_three)
        plt.legend(['$\\lambda_0$', '$\\lambda_1$', '$\\lambda_2$'])
        plt.xlabel('$\\rho$')
        plt.ylabel('$|\\psi(\\rho)|^2$')
        plt.savefig('../plots/schrodinger_%d_%d_' % (n, rho_max) + str(omega) + '.pdf')
    """

