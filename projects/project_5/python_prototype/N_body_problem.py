"""
Solves the N-body problem using RungeKutta or Velocity-Verlet
"""
import ODESolver
import numpy as np
import random
import matplotlib.pyplot as plt
from  mpl_toolkits.mplot3d import Axes3D

def pick_point_in_sphere(R):
    """
    Picks a uniformly random point in a sphere of radius R.  Pretends to be picking in a
    cube, and discards anything outside the sphere range.
    """
    while True:
        x = random.random()*R    
        y = random.random()*R
        z = random.random()*R
        if np.sqrt(x**2 + y**2 + z**2) < R:
            return np.array([x, y, z])
    
class Cluster:
    """
    Represents a system of N particles
    """ 
    G = 6.6740831E-11
    point_of_origin = [0, 0, 0] 
    cb_list = []
    def __init__(self, N=2, n_steps=1000):
        Cluster.cb_list = [CelestialBody(n=n_steps, M=1.0e30) for i in range(N)]
        self.solvers = [ODESolver.RungeKutta4]
    
    def simulate(self, t_values):
        for i, t in enumerate(t_values[:-1]):
            dt = t_values[i+1] - t
            for body in Cluster.cb_list:
                body.step(dt, i)
    
    def visualize(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for body in Cluster.cb_list:
            print 'plotting for ', body
            ax.plot(body.r[:,0], body.r[:,1], body.r[:,2])
        plt.show() 
            
class CelestialBody:
    """
    Represents one celestial particle, either a star or a small number of
    stars. With is is associated a size (mass), a position vector, a velocity
    vector and an acceleration vector.
    """
 
    def __init__(self, radius=1, M=1.0e5, n=100):
         
        self.r = np.zeros([n+1, 3])
        self.r[0] = pick_point_in_sphere(radius)
        self.v = np.zeros([n+1, 3])
        self.a = np.zeros([n+1, 3])
        self.M = float(M)
         
    
    def __call__(self, celestial_body_list, i):
        """
        Returns the total force exerted by all the other celestial bodies on
        this body.
        """    
        total_force = 0 
        for body in celestial_body_list:
            if body == self:
                continue
            total_force += Cluster.G * self.M * body.M * (body.r[i, :] - self.r[i, :]) / float(np.linalg.norm(body.r[i, :] - self.r[i, :]) ** 3)
        return total_force 
    
    def step(self, dt, i):
        """
        Pushes the celestial body one time increment forward in time.
        """
        self.a[i] = self(Cluster.cb_list, i) / self.M
        self.v[i+1] = self.v[i] + self.a[i]*dt
        self.r[i+1] = self.r[i] + self.v[i+1]*dt

if __name__ == "__main__":

    test = Cluster(50, 100)
    test.simulate(np.linspace(0, 10, 100))
    test.visualize()
