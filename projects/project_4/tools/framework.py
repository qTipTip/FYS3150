from numpy import linspace, sinh, cosh, exp
from matplotlib.pyplot import plot, show, subplot, legend, xlabel, ylabel, title, suptitle, grid, savefig, tight_layout, subplots_adjust
import os

mc_cycles = 50000
T_start = 2.0
T_end = 2.4
T_step = 0.01
grid_dim = 80

temperature = []
mean_energy = []
mean_magnetization = []
specific_heat = []
susceptibility = []
absolute_magnetization = []
file_name = '../data/output.dat'
os.system('rm ../data/output.dat')
os.system('../src/a.out %d %f %f %f %d %s' % (grid_dim, T_start, T_end, T_step, mc_cycles, file_name))
for line in open('../data/output.dat', 'r').readlines():
    data = line.split()
    temperature.append(float(data[0]))
    mean_energy.append(float(data[1]))
    specific_heat.append(float(data[2]))
    susceptibility.append(float(data[3]))
    absolute_magnetization.append(float(data[4]))

suptitle('$%d\\times%d$ Ising model with $%d$ MC-cycles' % (grid_dim, grid_dim, mc_cycles), y=1.00)
subplots_adjust(top=0.5)
subplot(221)
grid('on')
plot(temperature, mean_energy)
xlabel('$k_BT$'); ylabel('$E / J$')
subplot(222)
grid('on')
plot(temperature, specific_heat)
xlabel('$k_BT$'); ylabel('$C_V / Jk_B$')
subplot(223)
grid('on')
plot(temperature, susceptibility)
xlabel('$k_BT$'); ylabel('$\\chi$')
subplot(224)
grid('on')
plot(temperature, absolute_magnetization)
xlabel('$k_BT$'); ylabel('$\\langle | M | \\rangle$')
tight_layout()
savefig('output.pdf')

n = grid_dim**2 # number of spins
partition_function = lambda temp : 4*cosh(8/temp) + 12 
exp_E     = lambda temp : -8*sinh(8/temp)/(3+cosh(8/temp)) / n
exp_M     = lambda temp : 0 / n
exp_M_abs = lambda temp : 2*(exp(8/temp) + 2)/(cosh(8/temp) + 3) / n 
C_v       = lambda temp : (64./(temp*temp))*(cosh(8/temp)/(3+cosh(8/temp)) - (sinh(8/temp)/(3+cosh(8/temp)))**2) / n
X         = lambda temp : (8./temp)*(exp(8./temp) + 1)/(3+cosh(8./temp)) / n 
X_abs     = lambda temp :(4./(temp*(cosh(8./temp)+3)))*(2*(exp(8./temp) +1 ) - ((exp(8./temp)+2)**2)/(cosh(8./temp)+3)) / n

print exp_E(1)
print exp_M(1)
print C_v(1)
print X(1)

T = 1
beta = 1.0/T;
B = 8.0*beta;

Z = 4.0*cosh(B)+12.0;
E = -(8.0*sinh(B))/(cosh(B) + 3.0);
Cv = -(((64.0)/(cosh(B) + 3.0))*(-cosh(B) + ((sinh(B)*sinh(B))/(cosh(B) + 3.0))))/(T*T);
M = 0.0;
absM = (16.0 + (8.0*exp(B)))/Z;
X = ((8.0*(exp(B) + 1.0))/(cosh(B) + 3.0))/(T);
absX = ( ((8.0*(exp(B) + 1.0))/(cosh(B) + 3.0)) - (absM*absM) )/(T);

print E/n, M/n, Cv/n, X/n
