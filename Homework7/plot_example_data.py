# Quick python script that is used to plot the resulting output from the fortran codes: 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from plot_funcs import plotter, secs_to_a
import numpy as np

# Compile a list of 'plotter' objects that basically allow for animation
plotter_list = []
names = ['homo_k_homo_sp', 'homo_k_hetero_sp','hetero_k_homo_sp', 'hetero_k_hetero_sp']
for dir in names:
    plotter_list.append(plotter(f'example_data/{dir}'))


# Create figure and set initial plot
fig, ax = plt.subplots(figsize=(8,7))
for p in plotter_list:
    p.plot_initial_setup(ax)


# Animating functions:
def init():
    for p in plotter_list:
        p.l.set_ydata(p.init)

def animate(i):
    for p in plotter_list:
        p.l.set_ydata(p.data[i,:,1])
    ax.set_title(f"Time: {np.around(secs_to_a(i*p.timestep), 2)} years")


anim = animation.FuncAnimation(fig, animate,
                            init_func = init,
                            frames = p.no_snapshots,
                            interval = 50,
                            blit = False)
ax.legend(names)

plt.show()