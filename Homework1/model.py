import numpy as np

class Model():

    def __init__(self, x, dx, c, K, Nt, rho):
        self.x = x                       # X-axis array
        self.dx = x[1] - x[0]            # X grid spacing - assumed to be constant
        self.c = c                       # wavespeed
        self.K = K                       # Material property
        self.rho = rho                   # Density
        self.Nt = Nt                     # Number of timesteps
        self.dt = dx/np.max(c)           # Timestep
        self.dim = len(self.x)           # Dimension of X
