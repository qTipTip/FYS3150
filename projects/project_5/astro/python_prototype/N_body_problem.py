"""
Solves the N-body problem using RungeKutta or Velocity-Verlet
"""
from Cluster import Cluster
if __name__ == "__main__":

    test = Cluster(N=4, n_steps=5)
    test.initialize(radius=100)
    test.animate()

