import matplotlib.pyplot as plt
import numpy as np
import os

E_values = []
M_values = []

dimension = 20
temp_start = 1.0
temp_end = 1.0
dtemp = 0.1

mc_cycles = 10000
file_name = 'throwaway.dat'

os.system('rm ../data/expectations.dat')
command = '../src/a.out %d %f %f %f %d %s' % (dimension, temp_start, temp_end, dtemp, mc_cycles, file_name)
os.system(command)
for line in open('../data/expectations.dat', 'r').readlines():
    line = line.split()

    E_values.append(float(line[1]))
    M_values.append(float(line[2]))

mc_cycles = [100*i for i in range(len(E_values))]
plt.suptitle('$%d\\times%d$ Ising model with $T = %0.2f$' % (dimension,dimension, temp_start), y=1.00)
plt.subplot(211)
plt.plot(mc_cycles, E_values)
plt.xlabel('Monte-Carlo cycles')
plt.ylabel('$\\langle E \\rangle / J$')
plt.subplot(212)
plt.plot(mc_cycles, M_values)
plt.xlabel('Monte-Carlo cycles')
plt.ylabel('$\\langle |\\mathcal{M}| \\rangle / J$')
plt.tight_layout()
plt.savefig('../article/1_mc_random.pdf')
