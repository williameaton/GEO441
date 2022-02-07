import numpy as np
from copy import copy
from initial_conditions import set_IC


class vel_stress_form():
    def __init__(self, model, BC_left, BC_right, IC_func, plot):

        # Boundary conditions
        self.m = model
        self.IC_func = IC_func
        self.BC_left = BC_left
        self.BC_right = BC_right

        # Model parameters
        self.dt = model.dt
        self.dx = model.dx
        self.K = model.K
        self.rho = model.rho

        # Setup initial conditions
        self.T = np.zeros((3, model.dim))
        self.v = np.zeros((3, model.dim))
        self.set_initial_conditions()

        # What to plot:
        self.plot_type = plot



    def march(self):

        if self.BC_left == "neumann":
            self.T[:,0]=0
            for i in range(0, len(self.v[0,:])-1):
                self.step_v(i, i)
                self.step_T(i, i+1)

        elif self.BC_left == "dirichlet":
            self.v[:, 0] = 0
            for i in range(0, len(self.T[0, :]) - 1):
                self.step_T(i, i)
                self.step_v(i, i+1)
        else:
            raise ValueError("BC type not supported")


        # Right Boundary condition:
        if self.BC_right == "neumann":
            self.T[:,-1] = 0
        else:
            self.v[:,-1] = 0


        # Update data arrays
        # Roll back arrays so that U[1,:] --> U[0,:] and U[2,:] --> U[1,:]
        self.v = np.roll(self.v, shift=-1, axis=0)
        self.T = np.roll(self.T, shift=-1, axis=0)

        # Update plotter - i.e. ensure plotter has correct data:
        self.plot = self.set_plot_type(self.plot_type)


    def set_initial_conditions(self):
        # In case we want to reset to the initial conditions - i.e for init() in animation
        U, V, T, zeros = self.IC_func(self.m.x)
        self.T[:2, :] = T
        self.v[:2, :] = V



    def step_v(self, k, l):
        self.v[2, l] = (self.dt / (self.rho[k] * self.dx)) * (self.T[1, k + 1] - self.T[1, k]) + self.v[0, l]

    def step_T(self, k,l):
        self.T[2, l] = (self.K[k] * self.dt / self.dx) * (self.v[1, k+1] - self.v[1, k]) + self.T[0, l]

    def set_plot_type(self, plot):
        if plot == "v":
            return self.v
        elif plot == "T":
            return self.T
        else:
            raise ValueError("Must be v or T")