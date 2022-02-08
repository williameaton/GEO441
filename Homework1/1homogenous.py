from wave import create_wave
from model import Model
from propagate import propagate
from gen_x import gen_x
import matplotlib.pyplot as plt

# ______________________________________________________________________________________________________________________
# 1 HOMOGENOUS MEDIUM
# ______________________________________________________________________________________________________________________
# The script below will produce one simulation for each 1(a) Dirichlet boundary conditions and (b) Neumann boundary
# conditions. The displacement, velocity and stress are all plotted.

# Define simulation parameters:
dx       = 0.1          # Grid spacing
L        = (0, 100)     # X domain length
Nt       = 2000         # Number of timesteps - kinda arbitrary.
x = gen_x(dx=dx, L=L)   # X domain array
ani = []                # Save animations to this list

# Define initial parameters and store in model object :
rho   = x*0 + 1                                       # Homogenous density array of same length as x
kappa = x*0 + 1                                       # Material property
c     = x*0 + 1                                       # Wave speed
m = Model(x=x, c=c, K=kappa, dx=dx, Nt=Nt, rho=rho) # Model object holding some metadata/model params


# Since the boundary conditions are the same on each end I will make both plots using a simple loop:
for BC in ["dirichlet", "neumann"]:

    # Now create as many waves as you want - either a vel-stress or displacement formulation
    # Can set the initial conditions using a function of your choice - as in initial_conditions.py
    w  = create_wave(type="displacement",  model=m,    BC_left=BC, BC_right=BC, label="Displacement")
    w2 = create_wave(type="vel-stress",    model=m,    BC_left=BC, BC_right=BC, plot="v", label="Velocity")
    w3 = create_wave(type="vel-stress",    model=m,    BC_left=BC, BC_right=BC, plot="T", label="Stress")

    # Collect all the waves into an array for propagation
    waves = [w, w2, w3]


    # Alternatively you can uncomment this code block and comment out the stuff below to view the animations one at a time
    ani.append(propagate(waves, m,  fig_title=BC))
plt.show()

    ### Uncomment this:
    #propagate(waves, m)
    #plt.show()




