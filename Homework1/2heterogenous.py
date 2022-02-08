from wave import create_wave
from model import Model
from propagate import propagate
from gen_x import gen_x
import matplotlib.pyplot as plt

# ______________________________________________________________________________________________________________________
# 2 HETEROGENOUS MEDIUM
# ______________________________________________________________________________________________________________________
# This is very similar to the other driver script but with a few additions


# Define simulation parameters:
dx       = 0.1          # Grid spacing
L        = (0, 100)     # X domain length
Nt       = 2000         # Number of timesteps - kinda arbitrary.
x = gen_x(dx=dx, L=L)   # X domain array

BC_l = "dirichlet"      # Left boundary condition
BC_r = "neumann"        # Right boundary condition

# Define initial parameters and store in model object :
rho   = x*0 + 1                                       # Homogenous density array of same length as x
kappa = x*0 + 1                                       # Material property
c     = x*0 + 1                                       # Wave speed

# change material properties for x >60:
mask = tuple([x >60])
c[mask]     = 2
kappa[mask] = 4
m = Model(x=x, c=c, K=kappa, dx=dx, Nt=Nt, rho=rho)   # Model object holding some metadata/model params


# Now create as many waves as you want - either a vel-stress or displacement formulation
# Can set the initial conditions using a function of your choice - as in initial_conditions.py
w  = create_wave(type="displacement",  model=m,    BC_left=BC_l, BC_right=BC_r, label="Displacement")
w2 = create_wave(type="vel-stress",    model=m,    BC_left=BC_l, BC_right=BC_r, plot="v", label="Velocity")
w3 = create_wave(type="vel-stress",    model=m,    BC_left=BC_l, BC_right=BC_r, plot="T", label="Stress")

# Collect all the waves into an array for propagation
waves = [w, w2, w3]

w.march()


# Alternatively you can uncomment this code block and comment out the stuff below to view the animations one at a time
ani = propagate(waves, m, fig_title="Heterogenous case")
plt.show()





