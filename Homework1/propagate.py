import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
def propagate(waves, m):
    # Assuming Waves is a list of different wave-type objects:

    # Create figure
    fig, ax = plt.subplots()
    ax.set_ylim(-1, 1)
    lines = []

    # Produce initial plots to initialise each 2D Line object:
    for w in waves:
        l, = ax.plot(m.x[1:-1], w.v[1,1:-1])
        lines.append(l)

    def init():
        for i in range(len(waves)):
            # Update the wave data to its initial values:
            waves[i].set_initial_conditions()
            lines[i].set_ydata(waves[i].v[1,1:-1])
            ax.set_title(f"Time: {m.dt*0}")

        return lines,

    def animate(i):
        for j in range(len(waves)):
            # Update the wave data to its initial values:
            waves[j].march()
            lines[j].set_ydata(waves[j].v[1,1:-1])
            ax.set_title(f"Time: {np.around(m.dt*i,2)} ; dt = {m.dt}")



    # Run the animation
    anim = animation.FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=m.N,
                                   interval=1,
                                   blit=False)
    return anim