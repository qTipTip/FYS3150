import numpy as np
class CelestialBody:
    """
    Represents one celestial particle. With it is associated a size, a position
    vector, a velocity vector and an acceleration vector. 
    """
    G = 6.6740831E-11
    def __init__(self, radius=20.0, M=1.0e5, n=10):
        """
        Initializes the celestial particle with vectors of given length, a mass
        and a radius.
        """
        self.r = np.zeros([n+1, 3])
        self.v = np.zeros([n+1, 3])
        self.a = np.zeros([n+1, 3])
        self.M = M
        self.R = radius

    def set_initial_conditions(self, r0, v0, a0):
        """
        Sets the initial position, velocity and acceleration of the celestial
        body.
        """
        try:
            self.r[0, :] = r0
            self.v[0, :] = v0
            self.a[0, :] = a0
        except IndexError:
            print 'Cannot assign to index - is the Celestial-Body initialized with a proper n?'
        except ValueError:
            print 'Invalid arguments - expects numpy arrays of length 3'
    
    def __call__ (self, celestial_body_list, i):
        """
        Returns the total forces exerted by all the other celestial bodies on
        this body at a timestep i.
        """
        total_force = [0, 0, 0]
        for body in celestial_body_list:
            if body == self:
                print 'Not calculating self'
                continue
            else: 
                total_force += CelestialBody.G * self.M * body.M * (self - (body, i))/float(np.linalg.norm(self - (body, i)))**3

        return total_force

    def __sub__(self, (other, i)):
        """
        Returns the distance vector between two celestial bodies at timestep i
        """
        return other.r[i,:] - self.r[i, :] 

    def step(self, dt, i, cb_list):
        """
        Pushes the celestial body one time increment - Using EulerCromer for now
        """
        self.a[i,:] = self(cb_list, i) / self.M
        self.v[i+1,:] = self.v[i,:] + self.a[i,:] * dt
        self.r[i+1,:] = self.r[i,:] + self.v[i,:] * dt

