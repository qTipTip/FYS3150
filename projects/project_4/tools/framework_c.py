import matplotlib.pyplot as plt
T = []
accepts = []
for line in open('../data/accepted_configs_as_function_of_temp.dat', 'r').readlines():
    line = line.split()
    T.append(float(line[0]) )
    accepts.append(float(line[-1]))



plt.plot(T, accepts)
plt.xlabel('$k_BT$')
plt.ylabel('Accepted configurations')
plt.title('Number of accepted configurations as a function of temperature')
plt.savefig('accepted_configs.pdf')
