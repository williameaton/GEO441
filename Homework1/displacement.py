from copy import copy
import numpy as np




class displacement_form():

    def __init__(self, m, BC_left, BC_right, IC_func):

        # Set initial conditions:
        u_init, v, T, zeros = IC_func(m.x)

        self.u  = np.zeros((3, m.dim))
        self.u[:2, :]  = u_init         # Set U_old and U to be u_initial; U_new not defined right now.

        self.m      = m
        self.IC_func = IC_func

        # Set two attributes as the corresponding functions for each boundary condition
        # This means I can call self.BCleft(0) and self.BCright(-1) and it will impliment the correct bc on the correct
        # index
        self.bcs = {"dirichlet": self.dirichlet}
        self.BCleft = self.bcs[BC_left]
        self.BCright = self.bcs[BC_right]

        # Calculate prefactor for marching eqn:
        self.prefactor = (m.c * m.dt / m.dx)**2

        self.plot = self.u


    def march(self):
        # Enforce boundary conditions:
        self.enforce_BC()

        # Calculate newest timestep
        for i in range(1, len(self.u[1,:]) - 1):
            self.u[2, i] = self.prefactor[i]*(self.u[1, i + 1] - 2 * self.u[1, i] + self.u[1, i - 1]) + 2 * self.u[1, i] - self.u[0, i]

        # Roll back arrays so that U[1,:] --> U[0,:] and U[2,:] --> U[1,:]
        self.u = np.roll(self.u, shift=-1, axis=0)

        # Can only plot u here so:
        self.plot = self.u

    def set_initial_conditions(self):
        # In case we want to reset to the initial conditions
        u_init, V, T, zeros = self.IC_func(self.m.x)
        self.u[0,:] = u_init
        self.u[1,:] = u_init
        self.u[2,:] = zeros


    def enforce_BC(self):
        # Apply correct BC functions to each index
        self.BCleft(0)
        self.BCright(-1)

    def dirichlet(self, index):
        self.u[1, index] = 0

