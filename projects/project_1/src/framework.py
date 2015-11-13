""" This python program works as a framework for this project. It handles the
running of the C++ script for various grid resolutions n.
It runs both the tdma and the LU-decomposition algorithms for each time step.
It produces a plot containing the closed form solution of the Poisson equation
as well as the two numerical methods for each n. The file name will be on the form
step_size.pdf
It also extracts the largest relative error in our numerical approximation for
each n and plots these versus the step length h. File name will be on the form
relative_error.pdf
The C++ program spits out file in the format solver_timesteps.dat
Example: tdma_1000.dat
It is a three column data file with
grid_point numerical exact
"""
import os
import numpy as np
from matplotlib.pyplot import plot, savefig, legend, title, show, close, ylim, xlim, xscale, xlabel, ylabel

def plot_rel_err():
    """
    Plots the relative error in the TDMA for various n
    Generates a table of values in tex format
    """
    error_values = []
    n_values = []
    with open("../data/relative_error.dat", 'r') as errors:
        for line in errors:
            line = line.split()
            n_values.append(int(line[0]))
            error_values.append(float(line[1]))
    n_values = np.array(n_values)
    error_values = np.array(error_values)
    plot(n_values, error_values)
    xscale('log')
    legend(['relative error'])
    title('Relative error in data set')
    savefig('../plots/relative_error.pdf')
    print 'writing relative error to ../plots/relative_error.pdf'

    with open("../article/relative_error_table.tex", 'w') as out:
        output = "\\begin{tabular}{lc} \\\\\n"
        output += "\\hline\n"
        output += "$n$ & $\\varepsilon_{\\text{max}}$ \\\\\n"
        output += "\\hline\n"
        for n, eps in zip(n_values, error_values):
            output += "%d & %f" % (n, eps)
            output += "\\\\\n"
        output += "\\hline\n"
        output += "\\end{tabular}"
        out.write(output) 

def print_time(tdma, lu):
    with open('../data/timer_tdma.dat', 'r') as tdma_times, open('../data/timer_lu.dat') as lu_times:
        for line in tdma_times:
            line = line.split()
            n = int(line[0])
            time = float(line[1])
            if n in tdma:
                tdma[n].append(time)
            else:
                tdma[n] = []
                tdma[n].append(time)
        for line in lu_times:
            line = line.split()
            n = int(line[0])
            time = float(line[1])
            if n in lu:
                lu[n].append(time)
            else:
                lu[n] = []
                lu[n].append(time)
        
    for key in tdma:
        tdma[key] = sum(tdma[key]) / len(tdma[key])
    for key in lu:
        lu[key] = sum(lu[key]) / len(lu[key])
    
    with open("../article/times.tex", "w") as texfile:
        output = "\\begin{tabular}{lcc}\\\\\n"
        output += "\\hline \n"
        output += "n & tdma[s] & lu[s] \\\\"
        output += "\\hline \n"
        for x in sorted(tdma):
            if x in lu:
                output += "%d & %f & %f" % (x, tdma[x]/1.0e9, lu[x]/1.0e9)
            else:
                output += "%d & %f & N/A" % (x, tdma[x]/1.0e9)
            output += "\\\\\n"
        output += "\\hline \n"
        output += "\\end{tabular}"
        texfile.write(output) 

file_names = []
n_values = [10, 100, 1000, 10000, 100000, 1000000]
error_values = []
time_tdma = {}
time_lu = {}

try:
    os.system('rm ../data/relative_error.dat')
except:
    pass
for n in n_values:
    print 'Solving for n = %d' % n
    os.system('./a.out %d 1' % n)
    os.system('./a.out %d 2' % n)
    
    tdma_vals = []
    lu_vals = []
    exact_vals = []
    x_vals = []

    try:
        with open('../data/tdma_%d.dat' % n) as tdma:
            for line in tdma:
                line = line.split()
                x_vals.append(float(line[0])) 
                tdma_vals.append(float(line[1]))
                exact_vals.append(float(line[2]))
        plot(x_vals, tdma_vals, x_vals, exact_vals)
        ylim(0, 1)
        xlim(0, 1)
        legend(['tdma', 'analytic'])
        title('n = %d' % n)
        xlabel('$x_i$')
        ylabel('$u(x_i)$')
        savefig('../plots/tdma_%d.pdf' % n)
        close()
        print '\twriting plot to ../plots/tdma_%d.pdf' % n
    except:
        print "TDMA method failed for n = %d" % n
     
    try:
        if n < 100000 :
            with open('../data/lu_%d.dat' % n) as lu:
                for line in lu:
                    line = line.split()
                    lu_vals.append(float(line[1]))
            plot(x_vals, lu_vals, x_vals, exact_vals)
            ylim(0, 1)
            xlim(0, 1)
            legend(['lu', 'analytic'])
            title('n = %d' % n)
            xlabel('$x_i$')
            ylabel('$u(x_i)$')
            savefig('../plots/lu_%d.pdf' % n)
            close()
            print '\twriting plot to ../plots/lu_%d.pdf' % n
    except:
        print "LU method failed for n = %d" % n

print_time(time_tdma, time_lu)
plot_rel_err()
