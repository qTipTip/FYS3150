"""
Solves the N-body problem using RungeKutta or Velocity-Verlet
"""
from Cluster import Cluster
if __name__ == "__main__":

    test = Cluster(N=4)
    test.initialize()
    test.visualize()

