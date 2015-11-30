from system_class import System
import os
import matplotlib.pyplot as plt
import numpy as np

def taskE(run_simulations=False):
    """
    Performs the simulations needed to compare the two
    integration methods Velocity Verlet and Euler-Cromer.
    """
    simulations = []
    initial_temp = 300 
    b = 5.26
    N = 5
    description = 'energy_fluctuations_'
    
    dt_values = 10**np.linspace(-15, -14, 40)
    for i, dt in enumerate(dt_values):
        
        euler = System(T = initial_temp, dt=dt, N_each_dimension=5, b=5.26, ID=i, simulation_description=description+'euler')
        verlet = System(T = initial_temp, dt=dt, N_each_dimension=5, b=5.26, ID=i, simulation_description=description+'verlet')
        
        simulations.append([euler, verlet])
    
    # run simulations 

    if run_simulations:
        for sim in simulations:
            sim[0].run_simulations(1) # euler
            sim[1].run_simulations(2) # verlet
    
    # extract data
    euler_sd = []
    verlet_sd = []
    for sim in simulations:
        sim[0].read_data()
        sim[1].read_data()
    
        euler_sd.append(np.std(sim[0].total))
        verlet_sd.append(np.std(sim[1].total))

    # plot data
    plt.plot(dt_values, euler_sd, 'ro-', dt_values, verlet_sd, 'go-')
    plt.legend(['Euler-Cromer', 'Velocity-Verlet'], loc='best')
    plt.xlabel('$\\Delta t [s]$')
    plt.ylabel('$\\sigma_E/N_{\\mathrm{atoms}} [eV]$')
    plt.grid('on')
    plt.savefig('energy_fluctuations.pdf')

def taskF(run_simulations=False):
    """
    This method runs simulations for a system with a fixed
    lattice constant, fixed number of unit cells, but for
    varying initial temperature. We are interested in plotting
    the ratio between the instantaneous temperature and the
    initial temperature.
    """

    simulations = []
    
    T_values = range(20, 500, 10) 
    b = 5.26
    N = 5
    description = 'temperature_ratios' 
    for i, T_i in enumerate(T_values):
        simulations.append([System(T=T_i, b=5.26, N_each_dimension=5, ID=i, simulation_description=description), T_i])

    if run_simulations:
        for sim, _ in simulations:
            sim.run_simulations()
    
    for i, sim in enumerate(simulations):
        sim[0].read_data()
        if i == int(len(simulations) / 2):
            plt.plot(sim[0].temperature / sim[1], 'black', alpha=1.0)
            continue;
        plt.plot(sim[0].temperature / sim[1], alpha=0.3)

    plt.grid('on')
    plt.xlabel('$t$')
    plt.ylabel('$T / T_i$')
    plt.savefig('temperature_ratio.pdf')

def taskG(run_simulations=False):
    """
    This method runs a set of simulations for varying initial temperature T
    and looks at the diffusion constant. This is in order to measure the melting
    point of the system.
    """

    simulations = []

    T_values = range(400, 2000, 100)
    b = 5.26
    N = 10
    description = 'diffusion_constant'
    for i, T_i in enumerate(T_values):
        simulations.append(System(T = T_i, b=5.26, N_each_dimension=N, ID=i, simulation_description=description))

    if run_simulations:
        for sim in simulations:
            sim.run_simulations()

    for sim in simulations:
        sim.read_data()
        plt.plot(sim.diffusion)
    plt.xlabel('$t$')
    plt.ylabel('$D = \\langle r^2(t) \\rangle / 6t$')
    plt.grid('on')
    plt.savefig('diffusion_constant.pdf')

taskG(True)
