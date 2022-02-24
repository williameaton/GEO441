import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#sys.path.append('/home/we3822/Documents/princeton/GEO441/GEO441/Homework1')
sys.path.append('/Users/eaton/Documents/Princeton/GEO441/Homework1')
from model import Model
from waveWE import create_wave
from propagate import propagate
from general_funcs import homo_ofsize, save_anim_mp4, hetero_ofsize

# Simulation:
simulation = "heterogenous"


# Create initial homogenous model:
dx = 0.2
x = np.arange(0, 100+dx, dx)
rho = homo_ofsize(1, x)

if simulation=="homogenous":
    K = homo_ofsize(1,x)
elif simulation=="heterogenous":
    K = hetero_ofsize(bounds=[0, 60, 100], vals=[1,4], x=x)
else:
    raise ValueError("Simulation type not correct.")

c = K/rho
Nt = 1000

m = Model(x, dx, c, K, Nt, rho, dt=0.03)
BC = "neumann"


w  = create_wave(type="vel-stress",  model=m,  BC_left=BC, BC_right=BC, solver="pseudospectral", label="Pseudospectral")
w2  = create_wave(type="vel-stress",  model=m,  BC_left=BC, BC_right=BC, solver="finite_difference", label="Finite Difference")

# Collect all the waves into an array for propagation
waves = [w2, w]

ani = propagate(waves, m,  fig_title=BC, figsize=(12,8), ylims=[-0.3, 0.3])

#plt.show()
save_anim_mp4(ani, f"./videos/{simulation}.mp4", fps=50)