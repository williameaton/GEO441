import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys
from plot_funcs import plotter, secs_to_a

directory = sys.argv[1]

p = plotter(directory)

fig, ax = plt.subplots()
p.plot_initial_setup(ax)

def init():
    p.l.set_ydata(p.init)

def animate(i):
    p.l.set_ydata(p.data[i,:,1])
    ax.set_title(f"Time: {np.around(secs_to_a(i*p.timestep), 2)} years")


anim = animation.FuncAnimation(fig, animate,
                            init_func = init,
                            frames = p.no_snapshots,
                            interval = 50,
                            blit = False)

plt.show()


