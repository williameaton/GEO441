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
            self.march = self.forward_matrix_march
        elif method == "backward":
            self.A = self.create_A_matrix_backward()
            self.march = self.backward_matrix_march
        elif method == "crank-nicolson":
            self.A = self.crank_nicholson_matrix()
            self.march = self.cranknicholson_march
        else:
            raise ValueError("Must be 'forward', 'backward' or 'crank-nicolson' ")


    def update_T(self, T):
        # Function to update T if required
        self. T = T


    def set_IC(self):
        # Hard-coded initial conditions - impulse at x = 50
        self.T[:, np.where(self.m.x == 50)] = 1


    def create_A_matrix_forward(self):
        # Creates matrix for first order forward method for diffusion equation:
        # T(t+1) = (dt/(rho*cp*dx*dx))*[ K[x+1]*T[x+1] + T[x]*(- K[x+1] - K[x]) + K[x]*T[x-1] ] + T[x]

        # Insert 0 at each end of the diagonal so that the boundary condition is met
        B = ss.diags(
            [self.m.k[1:] * self.m.prefac[1:],                                                              # upper diag
            np.insert((1 - (self.m.k[2:] + self.m.k[1:-1])*self.m.prefac[1:-1]), [0, self.m.dim-2], [0,0]), # diagonal
            (self.m.k[1:]) * self.m.prefac[:-1]],                                                           # lower diag
            offsets=[-1, 0, 1], shape=(self.m.dim, self.m.dim))
        return B



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

        # implement boundary conditions:
        cn[0,:] = cn[-1,:] = 0

        return cn


    def create_A_matrix_backward(self):
        # Creates matrix for first order backward method for diffusion equation:
        # Create A matrix:
        A = np.zeros((self.m.dim, self.m.dim))

        pre = self.m.dt/(self.m.rho * self.m.cp * (self.m.dx**2))
        for i in range(1, self.m.dim - 1):

            C1 = (self.m.k[i+1] - self.m.k[i-1])/4

            A[i, i - 1]  =  pre[i] * (C1 - self.m.k[i])
            A[i, i]      =  1 + 2*pre[i]*self.m.k[i]
            A[i, i + 1]  =  -pre[i]* (C1 + self.m.k[i])



        # boundary matrix vals:
        A[0, 0] = 1 + 2 * pre[0] * self.m.k[0]
        A[0, 1] = A[1,2] # Approx?
        A[-1, -1] = 1 + 2 * pre[-1] * self.m.k[-1]
        A[-1, -2] = A[-2,-3]  # Approx?

        A_inv = np.linalg.inv(A)

        A_inv[0,0] = 0
        A_inv[-1,-1] = 0

        return A_inv



    def forward_matrix_march(self):
        self.T[0, :] = self.A*np.transpose(self.T[1, :])
        self.T = self.T[::-1, :]


    # Backward matrix march and cranknicolson are identical - could combine into one function but only a few lines of
    # code so maybe keep seperate for legibility
    def backward_matrix_march(self):
        self.T[1, :] = np.matmul(self.A, np.transpose(self.T[0, :]))
        self.T = self.T[::-1, :]

    def cranknicholson_march(self):
        self.T[1, :] = self.A*np.transpose(self.T[0, :])
        self.T = self.T[::-1, :]

