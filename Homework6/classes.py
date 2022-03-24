import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.special as ss

class Diffusion():
    # CONSTRUCTOR
    def __init__(self, nelem, T1, q0, alpha, rho=1, cp=1, kappa=1, L=np.pi/2):
        # Domain params:
        self.n = nelem
        self.n1 = self.n+1
        self.L = L
        self.x = np.linspace(0, self.L, self.n1)
        self.h = (self.L)/(self.n)

        # Integration params:
        self.dt = 0.1
        self.alpha = alpha

        # Modal params:
        self.rho = np.zeros(self.n) + rho
        self.cp = np.zeros(self.n) + cp
        self.kappa = kappa
        self.T1 = T1
        self.q0 = q0
        self.f  = np.zeros(self.n) + 0

        # Local variables
        self.N1 = np.array([1, 0])              # Shape func 1
        self.N2 = np.array([0, 1])              # Shape func 2
        self.k  = np.zeros((2, 2))               # Local stiffness

        # Calculate the global matrices:
        #   Initialise:
        self.K = np.zeros((self.n1, self.n1))   # Global stiffness
        self.M = np.zeros((self.n1, self.n1))   # Global mass
        self.F = np.zeros(self.n1)              # Global force
        #   Assemble:
        self.construct_global_matrices()

        # Global variables:
        self.d = self.set_IC()


    # CLASS FUNCTIONS:
    def set_IC(self):
        # Calculate T initial:
        T0 = 1 + np.cos(self.x)
        d = T0  # Since the shape funcs are only defined as non-zero at single nodes, D = T0 for each node

        # However for final element we need to do:
        try:
            d[-1] = T0[-1] / self.T1
        except:
            print("IC: T1 is zero, using d[-1]= 1")
            d[-1] = 1

        # Calculate v initial:
        self.v = np.zeros((2, self.n1))
        self.v[0,:-1] = np.matmul(np.linalg.inv(self.M[:-1,:-1]), (self.F[:-1] - np.matmul(self.K[:-1, :-1], d[:-1]) ) )
        #self.v[0, -1] = 0

        # Initialise analytical solution:
        self.T_an = T0

        return d


    def construct_global_matrices(self):
        # Loop through elements:
        for n in range(self.n):
            # Construct local k_ab matrix:
            for a in range(2):
                for b in range(2):
                    self.k[a][b] = self.kappa*(1 / self.h) * (-1) ** (a + b)

            rhs = (self.N1 + self.N2) * self.f[n:n + 2] * self.h / 2  # This allows f to be spatially varying
            # Assemble:
            self.K[n:n + 2, n:n + 2] += self.k
            self.F[n:n + 2] += rhs

            self.M[n:n + 2, n:n + 2] += self.rho[n] * self.cp[n] * self.h * (1 / 6) * np.array([[2, 1], [1, 2]])

        # Add BCs:
        self.F[0] += self.q0
        self.F[-2] -= self.k[0][1] * self.T1


    def step_d(self):
        # Time marching:
        # Step 1: calculate d_tilda:
        d_tilda = self.d + self.dt * (1 - self.alpha) * self.v[0, :]

        # Step 2: Invert for v(t+1)
        LHS = self.M[:-1, :-1] + self.alpha * self.dt * self.K[:-1, :-1]
        RHS = self.F[:-1] - np.matmul(self.K[:-1, :-1], d_tilda[:-1])
        self.v[1, :-1] = np.matmul(np.linalg.inv(LHS), RHS)

        # Now calculate d(t+1):
        self.d[:-1] += self.dt * (1 - self.alpha) * self.v[0, :-1] + self.alpha * self.dt * self.v[1, :-1]

        self.v = self.v[::-1]


    def animate(self, i):

        self.set_plot_stuff()
        self.step_d()
        # Analytical solution:
        self.T_an = 1 + np.exp(-(i+1)*self.dt)*np.cos(self.x)

        self.l_an.set_ydata(self.T_an)
        self.l.set_ydata(self.d)
        self.ax.set_title(f"Time = {np.round(i*self.dt,2)}")


    def init(self):

        self.set_plot_stuff()
        self.d = self.set_IC()
        self.l_an.set_ydata(self.T_an)
        self.l.set_ydata(self.d)


    def run(self, fig, ax):
        self.ax = ax
        self.l_an, = ax.plot(self.x, self.T_an, 'k')
        self.l,    = ax.plot(self.x, self.d, 'r:')

        anim = animation.FuncAnimation(fig, self.animate,
                                       init_func=self.init,
                                       frames=40,
                                       interval=2000*self.dt, blit=False)

        return anim

    def set_plot_stuff(self):
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("Temperature")
        self.ax.legend(["Analytical", "Finite Element"])

















class Ocean_diffusion():
    # CONSTRUCTOR
    def __init__(self, nelem, T0, q0, alpha, rho=1, cp=1, kappa=1, L=20, dt=0.01):
        # Domain params:
        self.n = nelem
        self.n1 = self.n + 1
        self.L = L
        self.y = np.linspace(0, self.L, self.n1)
        self.h = (self.L) / self.n

        # Integration params:
        self.dt = dt
        self.frames = 1000
        self.fps =int(self.frames*self.dt)
        if self.fps < 1:
            self.fps = 1

        self.alpha = alpha

        # Modal params:
        self.rho = rho
        self.cp = cp
        self.kappa = kappa
        self.T0 = T0
        self.q0 = q0
        self.f = np.zeros(self.n) + 0

        # Local variables
        self.N1 = np.array([1, 0])  # Shape func 1
        self.N2 = np.array([0, 1])  # Shape func 2
        self.k = np.zeros((2, 2))  # Local stiffness

        # Calculate the global matrices:
        #   Initialise:
        self.K = np.zeros((self.n1, self.n1))  # Global stiffness
        self.M = np.zeros((self.n1, self.n1))  # Global mass
        self.F = np.zeros(self.n1)  # Global force
        #   Assemble:
        self.construct_global_matrices()

        # Global variables:
        self.d = self.set_IC_ocean()


    # CLASS FUNCTIONS:
    def set_IC_ocean(self):
        # Initial condition:
        T0 = self.y * 0
        T0[:-1] = 1300
        d = T0  # Since the shape funcs are only defined as non-zero at single nodes, D = T0 for each node

        # Initialise velocity
        self.v = np.zeros((2, self.n1))
        self.v[0, :-1] = np.matmul(np.linalg.inv(self.M[:-1, :-1]), (self.F[:-1] - np.matmul(self.K[:-1, :-1], d[:-1])))

        self.calc_analytical(t=0.00000001)
        #self.T_an[-1] = 0
        return d


    def construct_global_matrices(self):
        # Loop through elements:
        for n in range(self.n):
            # Construct local k_ab matrix:
            for a in range(2):
                for b in range(2):
                    self.k[a][b] = self.kappa * (1 / self.h) * (-1) ** (a + b)

            rhs = (self.N1 + self.N2) * self.f[n:n + 2] * self.h / 2  # This allows f to be spatially varying
            # Assemble:
            self.K[n:n + 2, n:n + 2] += self.k
            self.F[n:n + 2] += rhs

            self.M[n:n + 2, n:n + 2] += self.rho * self.cp * self.h * (1 / 6) * np.array([[2, 1], [1, 2]])

        # Add BCs:
        self.F[0] += self.q0
        self.F[-1] = self.k[0][1] * self.T0 # This currently wont do anything for this version


    def step_d(self):
        # Time marching:
        # Step 1: calculate d_tilda:
        d_tilda = self.d + self.dt * (1 - self.alpha) * self.v[0, :]

        # Step 2: Invert for v(t+1)
        LHS = self.M[:-1, :-1] + self.alpha * self.dt * self.K[:-1, :-1]
        RHS = self.F[:-1] - np.matmul(self.K[:-1, :-1], d_tilda[:-1])
        self.v[1, :-1] = np.matmul(np.linalg.inv(LHS), RHS)

        # Now calculate d(t+1):
        self.d[:-1] += self.dt * (1 - self.alpha) * self.v[0, :-1] + self.alpha * self.dt * self.v[1, :-1]

        self.v = self.v[::-1]


    def animate(self, i):
        # Step the finite difference method and calculate the analytical solution:
        self.step_d()
        self.calc_analytical(t=(i+1)*self.dt)

        # Update temperature plot
        self.l.set_xdata(self.d)
        self.l_an.set_xdata(self.T_an)
        self.Tz_ax.set_title(f"Time = {np.round(i * self.dt, 2)}")

        # Plot isotherms:
        for l in [self.iso200, self.iso400, self.iso600, self.iso800]:
            l.update(x=i+1, y=self.y[self.get_isotherm_ind(l.temp)[0]])


        return self.iso200.line, self.iso400.line, self.iso600.line, self.iso800.line,




    def init(self):
        # Temp profile plotting:
        self.d = self.set_IC_ocean()
        self.l_an.set_xdata(self.T_an)
        self.l.set_xdata(self.d)

        # Plate plotting:
        self.plate_ax.clear()
        self.set_legend()

        for l in [self.iso200, self.iso400, self.iso600, self.iso800]:
            l.line = self.plate_ax.scatter(l.xx, l.yy, c=l.colour)
            l.update(x=0, y=self.y[self.get_isotherm_ind(l.temp)[0]])

        self.set_plot_labels()


        return self.iso200.line, self.iso400.line, self.iso600.line, self.iso800.line,


    def run(self, fig, ax):
        # Designate 2 axes:
        self.Tz_ax = ax[0]
        self.plate_ax = ax[1]

        self.set_plot_labels()

        # Initial plots for temp-depth profiles of analytical and FE:
        self.l_an, = self.Tz_ax.plot(self.T_an, self.y, 'k')
        self.l,    = self.Tz_ax.plot(self.d, self.y, 'r:')

        # Generate 4 isotherm objects
        self.iso200 = Isotherm(temp=200, colour='orange', axis=self.plate_ax)
        self.iso400 = Isotherm(temp=400, colour='red',    axis=self.plate_ax)
        self.iso600 = Isotherm(temp=500, colour='green',  axis=self.plate_ax)
        self.iso800 = Isotherm(temp=600, colour='purple', axis=self.plate_ax)


        # Animate the whole system
        anim = animation.FuncAnimation(fig, self.animate,
                                       init_func=self.init,
                                       frames=self.frames,
                                       interval=self.fps, blit=False)


        return anim



    def calc_analytical(self, t):
        self.T_an =  1300 - 1300*ss.erfc((self.L - self. y) / (2 * np.sqrt((self.kappa * t) / (self.rho * self.cp))))


    def get_isotherm_ind(self, T_crit):
        min_val = np.min(np.abs(self.T_an - T_crit))
        ind = np.where(self.T_an == min_val + T_crit)
        if len(ind[0]) == 0:
            ind = np.where(self.T_an == T_crit - min_val)
        elif len(ind[0]) > 1:
            ind = ind[0]

        ind = np.array(ind)
        return ind[0]


    def set_plot_labels(self):
        self.Tz_ax.set_xlabel(r"Temperature [$^o C$]")
        self.Tz_ax.set_ylabel(r"Depth (y)")
        self.Tz_ax.set_ylim([0, 22])
        self.Tz_ax.set_xlim([-100, 1400])

        self.plate_ax.set_xlabel(r"Iteration number")
        self.plate_ax.set_ylabel(r"Depth (y)")
        self.plate_ax.set_xlim([0, self.frames])
        self.plate_ax.set_ylim([14, 20])
        self.set_legend()


    def set_legend(self):
        self.Tz_ax.legend(["Analytical", "Finite Element"], loc='lower left')
        self.plate_ax.legend([r"200 $^o C$", r"400 $^o C$", r"500 $^o C$", r"600 $^o C$"], loc='lower left')







class Isotherm():
    def __init__(self, temp, colour, axis):
        self.xx = [0]
        self.yy = [1000]
        self.colour = colour
        self.temp = temp
        self.line = axis.scatter(x=self.xx, y=self.yy, c=colour)


    def update(self, x, y):
        self.xx.append(x)
        self.yy.append(y)
        self.line.set_offsets(np.transpose(np.array([self.xx, self.yy])))