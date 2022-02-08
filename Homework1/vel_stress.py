import numpy as np
from copy import copy
#from initial_conditions import set_IC


class vel_stress_form():
    def __init__(self, model, BC_left, BC_right, plot, label):

        # Boundary conditions
        self.m = model
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
        self.u = np.zeros((3, model.dim))
        self.set_initial_conditions()

        # Stuff for plotting:
        self.plot_type = plot                           # Variable to plot:
        self.plot = self.set_plot_type(self.plot_type)  # Initialise plot variable
        self.label = label                              # Legend label


    def march(self):
        for i in range(1, self.m.dim-1):
            self.step_v(i)
            self.step_T(i)

        self.apply_boundary_conditions()

        # Update data arrays
        # Roll back arrays so that U[1,:] --> U[0,:] and U[2,:] --> U[1,:]
        self.v = np.roll(self.v, shift=-1, axis=0)
        self.T = np.roll(self.T, shift=-1, axis=0)

        # Update plotter - i.e. ensure plotter has correct data:
        self.plot = self.set_plot_type(self.plot_type)


    def set_initial_conditions(self):
        # In case we want to reset to the initial conditions - i.e for init() in animation
        U, V, T = self.IC_func()
        self.T[:2, :] = T
        self.v[:2, :] = V



    def IC_func(self):
        # Calc initial conditions:
        u = np.exp(-0.1 * ((self.m.x - 50) ** 2))
        v = np.zeros(len(self.m.x))

        # Calculate T (apart from boundary values) by gradient of U:
        T = np.zeros((self.m.dim))
        for i in range(1, self.m.dim - 1):
            T[i] = (self.K[i] / (2 * self.dx))*(u[i + 1] - u[i - 1])

        # Now assign values at boundary - either 0 if neumann or T = adjacent element if dirichlet (this is approximation)
        self.apply_boundary_conditions()

        return u, v, T



    def step_v(self, k):
        self.v[2, k] = (self.dt / (self.rho[k] * self.dx)) * (self.T[1, k + 1] - self.T[1, k-1]) + self.v[0, k]

    def step_T(self, k):
        self.T[2, k] = (self.K[k] * self.dt / self.dx) * (self.v[1, k+1] - self.v[1, k-1]) + self.T[0, k]


    def impliment_BC_left(self, zero, adjacent):
        # Either velocity or stress will be 0 at the boundary, while the other will approximately be equal to its value
        # in the adjacent grid square:
        zero[:, 0] = 0
        adjacent[:, 0] = copy(adjacent[:, 1])

    def impliment_BC_right(self, zero, adjacent):
        # Either velocity or stress will be 0 at the boundary, while the other will approximately be equal to its value
        # in the adjacent grid square:
        zero[:, -1] = 0
        adjacent[:, -1] = copy(adjacent[:, -2])


    def apply_boundary_conditions(self):
        # Left boundary condition
        if self.BC_left == "neumann":
            self.impliment_BC_left(zero=self.T, adjacent=self.v)
        elif self.BC_left == "dirichlet":
            self.impliment_BC_left(zero=self.v, adjacent=self.T)
        else:
            raise ValueError("BC type not supported")

        # Right boundary condition
        if self.BC_right == "neumann":
            self.impliment_BC_right(zero=self.T, adjacent=self.v)
        elif self.BC_right == "dirichlet":
            self.impliment_BC_right(zero=self.v, adjacent=self.T)
        else:
            raise ValueError("BC type not supported")


    def set_plot_type(self, plot):
        if plot == "v":
            return self.v
        elif plot == "T":
            return self.T
        else:
            raise ValueError("Must be v or T")
