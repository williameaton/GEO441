import numpy as np

class Model():
    # Class holds some model parameters that are useful
    def __init__(self, x, dx, c, K, Nt, rho, dt=None):
        self.x = x                       # X-axis array
        self.dx = x[1] - x[0]            # X grid spacing - assumed to be constant
        self.c = c                       # wavespeed
        self.K = K                       # Material property
        self.rho = rho                   # Density
        self.Nt = Nt                     # Number of timesteps

        if dt != None:
            self.dt = dt
        else:
            self.dt = dx/np.max(c)           # Timestep
        self.dim = len(self.x)           # Dimension of X
