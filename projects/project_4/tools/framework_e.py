import os
import matplotlib.pyplot as plt
temp = []
energy = []
magnet = []
heat = []
suscept = []

grid_dims = [20, 40, 60, 80]
for i, grid_dimension in enumerate(grid_dims):

    t_start = 2.0
    t_end = 2.4
    t_step = 0.0005
    cycles = 6000
    file_name = '../data/e_%d.dat' % grid_dimension
    os.system('rm ' + file_name)
    command = '../src/a.out %d %f %f %f %d %s' % (grid_dimension, t_start, t_end, t_step, cycles,  file_name)

    os.system(command)
    
    t = []
    e = []
    m = []
    c = []
    s = []

    for line in open(file_name, 'r').readlines():
        line = line.split()

        t.append(float(line[0]))
        e.append(float(line[1]))
        m.append(float(line[2]))
        c.append(float(line[3]))
        s.append(float(line[4]))

    temp.append(t)
    energy.append(e)
    magnet.append(m)
    heat.append(c)
    suscept.append(s)

for i in range(len(temp)):
    plt.plot(temp[i], energy[i], label='$L = %d$' % grid_dims[i])
plt.legend(loc='best')
plt.grid('on')
plt.xlabel('$k_BT$')
plt.ylabel('$E / J$')
plt.savefig('energy_grid_size.pdf')
plt.clf()

for i in range(len(temp)):
    plt.plot(temp[i], magnet[i], label='$L = %d$' % grid_dims[i])
plt.legend(loc='best')
plt.grid('on')
plt.xlabel('$k_BT$')
plt.ylabel('$\\langle |\\mathcal{M}| \\rangle$')
plt.savefig('magnet_grid_size.pdf')
plt.clf()

for i in range(len(temp)):
    plt.plot(temp[i], heat[i], label='$L = %d$' % grid_dims[i])
plt.legend(loc='best')
plt.grid('on')
plt.xlabel('$k_BT$')
plt.ylabel('$C_V$')
plt.savefig('heat_grid_size.pdf')
plt.clf()

for i in range(len(temp)):
    plt.plot(temp[i], suscept[i], label='$L = %d$' % grid_dims[i])
plt.legend(loc='best')
plt.grid('on')
plt.xlabel('$k_BT$')
plt.ylabel('$\chi$')
plt.savefig('suscept_grid_size.pdf')
plt.clf()
