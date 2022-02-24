import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#sys.path.append('/home/we3822/Documents/princeton/GEO441/GEO441/Homework1')
sys.path.append('/Users/eaton/Documents/Princeton/GEO441/Homework1')
from model import Model
from waveWE import create_wave
from propagate import propagate
from general_funcs import homo_ofsize, save_anim_mp4


# Create initial homogenous model:
dx = 0.1
x = np.arange(30, 70+dx, dx)
rho = homo_ofsize(1, x)
K = homo_ofsize(1,x)
c = K/rho
Nt = 1000

m = Model(x, dx, c, K, Nt, rho, dt=0.03)
BC = "neumann"


w  = create_wave(type="vel-stress",  model=m,  BC_left=BC, BC_right=BC, solver="pseudospectral", label="Pseudospectral")
w2  = create_wave(type="vel-stress",  model=m,  BC_left=BC, BC_right="dirichlet", solver="finite_difference", label="Finite Difference")

# Collect all the waves into an array for propagation
waves = [w, w2]

ani = propagate(waves, m,  fig_title=BC, figsize=(12,8), ylims=[-0.3, 0.3])

save_anim_mp4(ani, "./videos/homogenous.mp4", fps=50)