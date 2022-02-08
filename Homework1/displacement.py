from copy import copy
import numpy as np


class displacement_form():
    def __init__(self, m, BC_left, BC_right, label):
        # Boundary conditions:
        self.BC_left = BC_left
        self.BC_right = BC_right

        # Initialise arrays and apply initial conditions:
        self.m      = m
        self.T = np.zeros((3, m.dim))
        self.v = np.zeros((3, m.dim))
        self.u = np.zeros((3, m.dim))
        self.set_initial_conditions()

        # Calculate prefactor for marching eqn:
        #self.prefactor = (m.c * m.dt / m.dx)**2
        self.prefactor = (((self.m.dt/self.m.dx)**2)/self.m.rho)

        # Stuff for plotting
        self.plot = self.u
        self.label = label


    def march(self):
        # Function marches one timestep

        # Calculate newest timestep
        for i in range(1, len(self.u[1,:]) - 1):
            # The full inhomogenous equation is quite long so I split it here into two parts:
            homo = self.prefactor[i]*self.m.K[i]*(self.u[1, i + 1] - 2*self.u[1, i] + self.u[1, i - 1]) + 2*self.u[1, i] - self.u[0, i]

            # Inhomo = 0 unless there are gradients in kappa (only true at the one discontinuity
            inhomo = 0.25*self.prefactor[i]*(self.u[1, i + 1] - self.u[1, i - 1])*(self.m.K[i+1] - self.m.K[i-1])
            self.u[2, i] =  homo + inhomo

        # Apply Boundary conditions:
        self.apply_BC_left()
        self.apply_BC_right()

        # Roll back arrays so that U[1,:] --> U[0,:] and U[2,:] --> U[1,:]
        self.u = np.roll(self.u, shift=-1, axis=0)

        # Can only plot u here so:
        self.plot = self.u




    def set_initial_conditions(self):
        # Calculate and apply initial conditions
        U, V, T = self.IC_func()
        self.u[:2, :] = U
        self.T[:2, :] = T
        self.v[:2, :] = V



    def IC_func(self):
        # Calc initial conditions:
        u = np.exp(-0.1 * ((self.m.x - 50) ** 2))
        v = np.zeros(len(self.m.x))

        # Calculate T (apart from boundary values) by gradient of U:
        T = np.zeros((self.m.dim))
        for i in range(1, self.m.dim - 1):
            T[i] = (self.m.K[i] / (2 * self.m.dx))*(u[i + 1] - u[i - 1])

        # Now assign values at boundary - either 0 if neumann or T = adjacent element if dirichlet (this is approximation)
        if self.BC_left == "dirichlet":
            T[0] = T[1]
        else:
            T[0]= 0

        if self.BC_right == "dirichlet":
            T[-1] = T[-2]
        else:
            T[-1] = 0

        return u, v, T


    def apply_BC_left(self):
        # Apply left boundary condition
        if self.BC_left == "dirichlet":
            self.u[:,0] = 0
        else:
            self.u[:, 0] = self.u[:, 1]


    def apply_BC_right(self):
        # Apply right boundary condition
        if self.BC_right == "dirichlet":
            self.u[:, -1] = 0
        else:
            self.u[:, -1] = self.u[:, -2]