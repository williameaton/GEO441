import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys
from plot_funcs import plotter, secs_to_ka

directory = sys.argv[1]
video_name = sys.argv[2]

p = plotter(directory)

fig, ax = plt.subplots()
ax.set_ylim([-1,1])
p.plot_initial_setup(ax)

def init():
    p.l.set_ydata(p.init)

def animate(i):
    p.l.set_ydata(p.data[i,:,1])
    ax.set_title(f"Time: {np.around(i*p.timestep*p.timestep_interval, 2)} Secs")


anim = animation.FuncAnimation(fig, animate,
                            init_func = init,
                            frames = p.no_snapshots,
                            interval = 2,
                            blit = False)


if video_name != "N0_VIDE0":  
    vid_file = f"{video_name}.mp4"
    print(f"Saving animation as: {vid_file}")
    writervideo = animation.FFMpegWriter(fps=150) 
    anim.save(vid_file, writer=writervideo)   
    print("Saved.") 

else: 
    plt.show()


