from wave import create_wave
from model import Model
from propagate import propagate
import matplotlib.pyplot as plt
from initial_conditions import *


x = np.linspace(0, 70, 1005)
dx = x[1] - x[0]
print(len(x))

rho = x*0 + 2
kappa = x*0 + 1
c = x*0 + 1


# Define initial parameters and store in model object :
m = Model(c = c, L = 70, K =kappa, dx = dx, N = 4000, rho=rho)


# Now create as many waves as you want - either a vel-stress or displacement formulation
# Can set the initial conditions using a function of your choice - as in initial_conditions.py
w = create_wave(type="vel-stress",  model=m,  BC_left="neumann", BC_right="neumann")

# Collect all the waves into an array for propagation
waves = [w]

#for n in range(100):
#    w.march()


# Propagate and show:
animation = propagate(waves, m)
plt.show()

