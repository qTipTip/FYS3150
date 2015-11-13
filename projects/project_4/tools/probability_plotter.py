import matplotlib.pyplot as plt
import numpy as np


energies = range(-800, 801)
probabil = []

for line in open('energy_prob_1.dat', 'r').readlines():
    line = line.split()   
    probabil.append(float(line[-1]))

print sum(probabil)
plt.plot(energies, probabil)
plt.title('Probabilty of various energies for $T = 1$')
plt.ylabel('$P(E)$')
plt.xlabel('$E / J$')
plt.savefig('probabilites.pdf')
