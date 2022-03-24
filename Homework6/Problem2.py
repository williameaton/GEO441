import matplotlib.pyplot as plt
from classes import Ocean_diffusion
import matplotlib.animation as animation

# Define number of elements and timestep:
nelem = 10
dt    = 0.01

diff = Ocean_diffusion(nelem=nelem, T0=0, q0=0, alpha=0.5, rho=1, cp=1, kappa=1, L=20, dt=dt)


fig, ax = plt.subplots(2, figsize=(7,7))
fig.set_tight_layout(True)

anim = diff.run(fig=fig, ax=ax)

f = f"videos/problem2_n{nelem}_dt_{dt}.mp4"
writervideo = animation.FFMpegWriter(fps=100)
anim.save(f, writer=writervideo)
