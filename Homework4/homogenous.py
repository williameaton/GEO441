import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/home/we3822/Documents/princeton/GEO441/GEO441/Homework1')
from model import Model
from waveWE import create_wave
from propagate import propagate
from general_funcs import homo_ofsize


# Create initial homogenous model:
dx = 0.1
x = np.arange(0, 100+dx, dx)
rho = homo_ofsize(1, x)
K = homo_ofsize(1,x)
c = K/rho
Nt = 100000

m = Model(x, dx, c, K, Nt, rho, dt=0.03)
BC = "neumann"


w  = create_wave(type="vel-stress",  model=m,  BC_left=BC, BC_right=BC, solver="pseudospectral", label="Pseudospectral")
w2  = create_wave(type="vel-stress",  model=m,  BC_left=BC, BC_right=BC, solver="finite_difference", label="Finite Difference")

# Collect all the waves into an array for propagation
waves = [w, w2]

ani = propagate(waves, m,  fig_title=BC)
plt.show()