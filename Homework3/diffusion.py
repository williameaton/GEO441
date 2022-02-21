import numpy as np
from scipy.sparse.linalg import inv
import scipy.sparse as ss

class diffusion():
    # Diffusion class requires a model object and a model for solving the diffusion equation in time
    # Method can be forward, backward of crank-nicolson and init func. essentially acts as a factory function to
    # create object using the correct method
    def __init__(self, model, method):
        self.m = model                          # Model object
        self.T   = np.zeros((2, model.dim))     # Temperature array for domain
        self.set_IC()                           # Hard coded initial conditions for this homework


        # Determine solver and create (time-independent) matrix for time-marching (self.A)
        # Point self.march to relevant marching function
        if method == "forward":
            self.A = self.create_A_matrix_forward()
        elif method == "backward":
            self.A = self.create_A_matrix_backward()
        elif method == "crank-nicolson":
            self.A = self.crank_nicholson_matrix()
        else:
            raise ValueError("Must be 'forward', 'backward' or 'crank-nicolson' ")

        self.method = method

    def set_IC(self):
        # Hard-coded initial conditions - impulse at x = 50
        self.T[:, np.where(self.m.x == 50)] = 1



    def create_A_matrix_forward(self):
        # Creates matrix for first order forward method for diffusion equation:
        # T(t+1) = (dt/(rho*cp*dx*dx))*[ K[x+1]*T[x+1] + T[x]*(- K[x+1] - K[x]) + K[x]*T[x-1] ] + T[x]
        # Approximate A[0,0] and A[m,m] as cant be calculated, as their adjacent values
        a = self.m.k[1:] * self.m.prefac[1:]
        b = (1 - (self.m.k[2:] + self.m.k[1:-1])*self.m.prefac[1:-1])
        b = np.insert(b, [0, -1], [b[0],b[-1]])
        c = (self.m.k[1:]) * self.m.prefac[:-1]

        A = ss.diags([a, b, c], offsets=[-1, 0, 1], shape=(self.m.dim, self.m.dim))
        return A



    def crank_nicholson_matrix(self):
        # See my report for definitions of each variable
        alpha = 2*self.m.rho*self.m.cp*(self.m.dx**2)/self.m.dt
        # Assuming that k is not time dependent in this homework problem:
        # \beta^n_i = \frac{k^{n}_{i+1} - k^{n}_{i-1} + k^{n-1}_{i+1} - k^{n-1}_{i-1}}{8}
        # becomes
        # \beta^n_i = \frac{k^{n}_{i+1} - k^{n}_{i-1}}{4}
        # Beta at i=0 is undefined so assume it is same as adjacent:
        beta  = (self.m.k[2:] - self.m.k[:-2])/4
        beta = np.insert(beta, [0], [beta[0]])

        # Calculate the matrix element values as defined in report:
        a = beta - self.m.k[:-1]
        b = alpha + 2*self.m.k
        c = -(beta + self.m.k[:-1])
        d = -a
        e = alpha - 2*self.m.k
        f = -c

        A = ss.diags([a,b,c], offsets=[-1, 0, 1], shape=(self.m.dim, self.m.dim))
        B = ss.diags([d,e,f], offsets=[-1, 0, 1], shape=(self.m.dim, self.m.dim))
        cn = inv(A)*B

        return cn


    def create_A_matrix_backward(self):
        # Creates matrix for first order backward method for diffusion equation:
        # Create A matrix:
        pre = self.m.dt / (self.m.rho * self.m.cp * (self.m.dx ** 2))

        # Because it uses i+1 and i-1 cant calculate end values
        C1 = (self.m.k[2:] - self.m.k[:-2]) / 4
        C1 = np.insert(C1, [0, -1], [C1[0], C1[-1]])

        a = pre[1:] * (C1[1:] - self.m.k[1:])
        b = 1 + 2 * pre[:] * self.m.k[:]
        c = -(pre[:-1] * (C1[:-1] + self.m.k[:-1]))

        A = ss.diags([a, b, c], offsets=[-1, 0, 1], shape=(self.m.dim, self.m.dim))
        A_inv = inv(A)

        return A_inv



    def march(self):
        # March one timestep using correct matrix
        self.T[1, :] = self.A*np.transpose(self.T[0, :])
        self.T = self.T[::-1, :]
        self.T[:,0] = self.T[:,-1] = 0 # Boundary condition



