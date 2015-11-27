import sys
import matplotlib.pyplot as plt
import numpy as np

try:
    filename = sys.argv[1]
except:
    pass


kinetic_energy = []
potential_energy = []
temperature = []
density = []

with open(filename, 'r') as infile:
    for line in infile:
        line = [float(value) for value in line.split()]
        temperature.append(line[0])
        kinetic_energy.append(line[1])
        potential_energy.append(line[2])
        density.append(line[3])

kinetic_energy = np.array(kinetic_energy)
potential_energy = np.array(potential_energy)
temperature = np.array(temperature)
density = np.array(density)
initial_temp = 300
plt.plot(initial_temp/temperature, label='temperature')
plt.legend()
plt.show()
plt.plot(kinetic_energy, label='kinetic')
plt.legend()
plt.show()
plt.plot(potential_energy, label='potential')
plt.legend()
plt.show()
