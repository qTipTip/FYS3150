import numpy as np

class CelestialBody:
    """
    Represents one celestial particle. With it is associated a size, a position
    vector, a velocity vector and an acceleration vector. 
    """
    G = 1
    def __init__(self, radius=1, mass=10.0, n=10):
        """
        Mass is given in multiples of solar masses. Distance is given in light
        years and the gravitational constant is set to 1. radius given in km
        """
        self.r = np.zeros([n+1, 3])
        self.v = np.zeros([n+1, 3])
        self.a = np.zeros([n+1, 3])
        
        self.mass = mass
        self.radius = radius

    def set_initial_conditions(self, U0):
        """
        Initializes the position, velocity and acceleration of the celstial
        body.
        """
        r0, v0, a0 = U0 
        try:
            self.r[0, :] = r0
            self.v[0, :] = v0
            self.a[0, :] = a0
        except IndexError:
            print 'Cannot assign to index - is the Celestial-Body initialized with a proper n?'
        except ValueError:
            print 'Invalid arguments - expects numpy arrays of length 3'

    def __call__(self, celestial_body_list, i):
        """
        Returns the total forces exerted by all the other celestial bodies on
        this body at a timestep i.
        """
        total_force = [0, 0, 0]
        collision_correction = 1.0
        for body in celestial_body_list:
            if body == self:
                continue
            else: 
                total_force += CelestialBody.G * self.mass * body.mass * (self - (body, i))/float(np.linalg.norm(self - (body, i)) + collision_correction)**3
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
        self.a[i,:] = self(cb_list, i) / self.mass
        self.v[i+1,:] = self.v[i,:] + self.a[i,:] * dt
        self.r[i+1,:] = self.r[i,:] + self.v[i,:] * dt
