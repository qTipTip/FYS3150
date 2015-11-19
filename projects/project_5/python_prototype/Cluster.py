import numpy as np
import random
import matplotlib.pyplot as plt
from  mpl_toolkits.mplot3d import Axes3D

from CelestialBody import CelestialBody
class Cluster:
    """
    Represents a system of N particles
    """
    point_of_origin = [0, 0, 0] # might need this

    def __init__(self, N=2, n_steps=6):
        self.N = N
        self.n_steps = n_steps
        self.cb_list = [CelestialBody(radius=random.randint(30, 500), M=1.0e5, n=n_steps) for i in range(N)]
         
    def visualize(self, i=0):
        """
        Shows a snapshot of the current state of the system at time t_i.
        Particles are scaled by size.
        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for body in self.cb_list:
            ax.scatter(body.r[i,0], body.r[i,1], body.r[i,2], s=body.R)
        plt.show()
        plt.clf()
    
    def animate(self):
        """
        Animates the time-evolution of the system.
        """
        for i in range(self.n_steps):
            self.visualize(i)
            for body in self.cb_list:
                body.step(1, i, self.cb_list)

    def initialize(self, radius=20):
        for body in self.cb_list:
            body.set_initial_conditions(r0=Cluster.pick_point_in_sphere(radius), v0=[0, 0, 0], a0=[0, 0, 0])

    @staticmethod
    def pick_point_in_sphere(R):
        """
        Picks a uniformly random point in a sphere of radius R. Pretends to be
        picking in a cube, then discards anything outside the sphere radius.
        """
        while True:
            x = random.random()*R
            y = random.random()*R
            z = random.random()*R
            if np.sqrt(x**2 + y**2 + z**2) < R:
                return np.array([x, y, z])
