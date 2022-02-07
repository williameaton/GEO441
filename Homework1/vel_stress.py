import numpy as np
from copy import copy
from initial_conditions import set_IC


class vel_stress_form():
    def __init__(self, model, BC_left, BC_right, IC_func):
        u_init, V_init, T_init, zeros = IC_func(model.x)

        self.dt = model.dt
        self.dx = model.dx
        self.K = model.K
        self.rho = model.rho

        self.T = np.zeros((3, model.dim))
        self.T[:2,:] = T_init

        self.v = np.zeros((3, model.dim))
        self.v[:2,:] = V_init

        self.m = model
        self.IC_func = IC_func


    def march(self):

        # If left boundary = neumann
        self.T[:,0]=0
        self.T[:,-1]=0

        for i in range(1, len(self.v[0,:])-2):
            self.v[2, i] = (self.dt / (self.rho[i] * self.dx)) * (self.T[1, i + 1] - self.T[1, i - 1]) + self.v[0, i]
            j = i + 1
            self.T[2, j] = (self.K[j] * self.dt / self.dx) * (self.v[1, j + 1] - self.v[1, j - 1]) + self.T[0, j]

        # If left is dirichlet I think you want it the other way around with T first and then v

        # Update data arrays
        # Roll back arrays so that U[1,:] --> U[0,:] and U[2,:] --> U[1,:]
        self.v = np.roll(self.v, shift=-1, axis=0)
        self.T = np.roll(self.T, shift=-1, axis=0)



    def set_initial_conditions(self):
        # In case we want to reset to the initial conditions
        U, V, T, zeros = self.IC_func(self.m.x)
        self.T[:2, :] = T
        self.v[:2, :] = V