import numpy as np

class model():
    def __init__(self, dx, L, cp, rho, k, dt=None, dtc=1):
        # Model class holds parameters for model setup of diffusion simulations
        self.dx  = dx                               # Grid spacing (assumed constant)
        self.x   = np.arange(L[0], L[1]+dx, dx)     # X-dimension array
        self.cp  = self.x*0 + cp                    # Specific heat capacity
        self.rho = self.x*0 + rho                   # Density
        self.k   = self.x*0 + k                     # Conductivity
        self.dim   = len(self.x)                    # Dimension of x array (no. of grid elements)

        # Calculate dt if not specified. dtc is a constant that multiplies the timestep
        if dt!=None:
            self.dt = dt
        else:
            self.dt = dtc*(dx*dx)/( np.amax(self.k)/np.amin(self.rho)*np.amin(self.cp))

        # MOVE THIS ELSEWHER
        self.prefac = self.dt/(self.rho*self.cp*dx*dx) # 1D array

