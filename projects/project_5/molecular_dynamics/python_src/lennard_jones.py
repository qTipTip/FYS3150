import numpy as np
import matplotlib.pyplot as plt

sigma = 3.405 # angstrom
epsil = 1 # MD units
def lennard_jones(rij):


    return 4 * epsil * ( (sigma / rij) ** 12 - (sigma / rij) ** 6)

rij_values = np.linspace(3, 10, 10000)
pot_values = lennard_jones(rij_values)

plt.plot(rij_values / sigma, pot_values / epsil)
plt.grid('on')
plt.xlabel('$r_{ij} / \\sigma$')
plt.ylabel('$U_{ij} / \\epsilon$')
plt.savefig('lennard_jones.pdf')
