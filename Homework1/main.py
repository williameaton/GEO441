from wave import create_wave
from model import Model
from propagate import propagate
import matplotlib.pyplot as plt
from gen_x import gen_x

# Define simulation parameters:
BC_left = "dirichlet"
BC_right = "dirichlet"
dx = 0.1
L  = (0, 100)


# Create x array
x = gen_x(dx=dx, L=L, BC_left=BC_left, BC_right=BC_right)

# Specify material properties along x
rho = x*0 + 2
kappa = x*0 + 1
c = x*0 + 1


# Define initial parameters and store in model object :
m = Model(x=x, c=c, K=kappa, dx=dx, Nt=4000, rho=rho)

# Now create as many waves as you want - either a vel-stress or displacement formulation
# Can set the initial conditions using a function of your choice - as in initial_conditions.py
w = create_wave(type="vel-stress",  model=m,  BC_left=BC_left, BC_right=BC_right, plot="T")
w2 = create_wave(type="displacement",  model=m,  BC_left=BC_left, BC_right=BC_right)

# Collect all the waves into an array for propagation
waves = [w, w2]

for n in range(1):
    w.march()


# Propagate and show:
animation = propagate(waves, m)
plt.show()

