# %matplotlib notebook
from model import model
from diffusion import diffusion
from animate import diff_animate
import matplotlib.animation as animation
import matplotlib.pyplot as plt

# DEFINE MODEL PARAMETERS - these are somewhat arbitrary
dx = 1                          # Grid spacing
L = (0, 100)                    # Domain limits
cp = 1                          # Specific heat capacity at constant pressure
rho = 1                         # Density
k = 0.7                         # Conductivity
dtc = [0.4, 0.45, 0.55, 0.6]    # Timesteps to be solved where each is multiplied by

# SOLVER METHOD:
method="forward"

# ______________________________________________________________________________________________________________________

# Initialise some variables
num = len(dtc)   # Number of different timesteps to simulate
ds = []          # Holds 'diffusion' objects
ls = []          # Holds matplotlib line objects for animation

# Initialise output figure
fig, ax = plt.subplots(num, figsize=(6,6), sharex=True)
fig.set_tight_layout(True)


# Loop for each different dt - create 'diffusion' object using object of 'model' class and append this diffusion
# object to a list called ds for sequential solving (marching) and animation
for i in range(num):
    ds.append(diffusion(model=model(dx=dx, L=L, cp=cp, rho=rho, k=k, dt=None, dtc=dtc[i]), method=method))

    # Generate matplotlib line objects to be animated (initial setup for each simulation) and add them to a list
    l, = ax[i].plot(ds[i].m.x, ds[i].T[1, :], 'k', linewidth=2)  # Plot initial setup
    ls.append(l)


# Generate animation object (this will solve at each timestep and plot on the figure and then return an animation):
ani = diff_animate(lines=ls, diff_obj=ds, fig=fig, axes=ax, interval=300, frames=100)


plt.show()

"""
# Animate and save:
f = f"./videos/{method}_diffusion.mp4"
writervideo = animation.FFMpegWriter(fps=5)
print("Saving video:")
ani.save(f, writer=writervideo)
print(f"Written to {f}")
"""

#print("For convenience I have already run and saved these movies. Alternatively if you want to view them interactively you can just run plt.show(). You will also need to uncomment '%matplotlib notebook' with the imports to view the animation in Jupyter, and the plt.close() that I use to supress the plot.")
