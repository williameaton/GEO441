import numpy as np

class Model():

    def __init__(self, L, dx, c, K, N, rho):
        self.dx = dx                     # X grid spacing - assumed to be constant
        self.x = np.arange(0, L+dx, dx)  # X-axis array
        self.c = c                       # wavespeed
        self.K = K                       # Material property
        self.N = N                       # Number of timesteps
        self.dt = dx/np.max(c)           # Timestep
        self.rho = rho                   # Density
        self.dim = len(self.x)           # Dimension of X
