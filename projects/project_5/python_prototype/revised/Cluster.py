from  mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil, glob
import random

from CelestialBody import CelestialBody
class Cluster:
    """
    Represents a system of N particles
    """
    point_of_origin = [0, 0, 0] # might need this

    def __init__(self, N=2, n_steps=1000):
        self.N = N
        self.n_steps = n_steps
        self.cb_list = [CelestialBody(radius=random.randint(30, 70), mass=float(random.randint(9, 11)), n=n_steps) for i in range(N)]
         
    def visualize(self, i):
        """
        Shows a snapshot of the current state of the system at time t_i.
        Particles are scaled by size.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim3d([-20, 20])
        ax.set_ylim3d([-20, 20])
        ax.set_zlim3d([-20, 20])
        for body in self.cb_list:
            ax.scatter(body.r[i,0], body.r[i,1], body.r[i,2], s=body.radius)
        plt.savefig('tmp_%04d.png' % i)
        plt.close()
    
    def animate(self):
        """
        Animates the time-evolution of the system.
        """
        subdir = 'animation'
        if os.path.isdir(subdir): # if the directory is already present, remove it
            shutil.rmtree(subdir)
        os.mkdir(subdir) # create directory
        os.chdir(subdir) # enter new directory 

        for i in range(self.n_steps):
            self.visualize(i)
            for body in self.cb_list:
                body.step(0.1, i, self.cb_list)

        cmd = 'convert -delay 4 tmp_*.png movie.gif'
        os.system(cmd)

        # Remove all the image files and change back to previous directory
        for a in glob.glob('tmp*.png'):
            os.remove(a)
        os.chdir(os.pardir)

    def initialize(self, radius=20):
        for body in self.cb_list:
            body.set_initial_conditions([Cluster.pick_point_in_sphere(radius), [0, 0, 0], [0, 0, 0]])

    @staticmethod
    def pick_point_in_sphere(R):
        """
        Picks a uniformly random point in a sphere of radius R. Pretends to be
        picking in a cube, then discards anything outside the sphere radius.
        """
        while True:

            x = -R + random.random()*2*R
            y = -R + random.random()*2*R
            z = -R + random.random()*2*R

            if np.sqrt(x**2 + y**2 + z**2) < R:
                return np.array([x, y, z])
