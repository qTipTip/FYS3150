"""
Solves the N-body problem using RungeKutta or Velocity-Verlet
"""
from Cluster import Cluster
if __name__ == "__main__":
    
    test = Cluster(N=40, n_steps=100)
    test.initialize(radius=20)
    test.animate()
