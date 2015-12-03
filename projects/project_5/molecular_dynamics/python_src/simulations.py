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
    This method runs a set of simulations for varying initial temperature T and
    looks at the diffusion constant. This is in order to measure the melting
    point of the system. We pick the final D of the simulation.
    """

    simulations = []

    T_values = range(400, 800, 5)
    b = 5.26
    N = 6
    description = 'diffusion_constant'
    for i, T_i in enumerate(T_values):
        simulations.append(System(T = T_i, b=5.26, N_each_dimension=N, ID=i, simulation_description=description))

    if run_simulations:
        for sim in simulations:
            sim.run_simulations()
    
    diffusion = []
    for sim in simulations:
        sim.read_data()
        diffusion.append(sim.diffusion[-1])

    diffusion = np.array(diffusion)
    T_values = np.array(T_values)

    T_before_melt = np.nonzero(T_values < 620)
    T_after_melt = np.nonzero(T_values >= 615)
    
    print T_values[T_before_melt][-1]
    A = np.vstack([T_values[T_before_melt], np.ones(len(T_values[T_before_melt]))]).T
    m, c = np.linalg.lstsq(A, diffusion[T_before_melt])[0]
    plt.plot(T_values[T_before_melt], m*T_values[T_before_melt] + c)

    A = np.vstack([T_values[T_after_melt], np.ones(len(T_values[T_after_melt]))]).T
    m, c = np.linalg.lstsq(A, diffusion[T_after_melt])[0]
    plt.plot(T_values[T_after_melt], m*T_values[T_after_melt] + c, )

    plt.plot(T_values, diffusion, 'bo', alpha=0.5)
    plt.xlabel('$t \\quad [K]$')
    plt.ylabel('$D = \\langle r^2(t) \\rangle / 6t \\quad [m^2/s]$')
    plt.grid('on')
    plt.savefig('diffusion_constant.pdf')

taskG(False)
