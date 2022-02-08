import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
def propagate(waves, m, fig_title=""):
    # This whole function basically animates the plot by marching through timesteps

    # Create figure
    fig, ax = plt.subplots()
    fig.suptitle(fig_title)
    ax.set_ylim(-1, 1)
    lines = []

    # Produce initial plots to initialise each 2D Line object:
    leg_list = []
    for w in waves:
        l, = ax.plot(m.x, w.plot[1, :])
        lines.append(l)
        leg_list.append(w.label)

    # Set figure legend:
    ax.legend(lines, leg_list)

    def init():
        # Function for re-setting animation when at end of animation frame loop
        for i in range(len(waves)):
            # Update the wave data to its initial values:
            waves[i].set_initial_conditions()
            lines[i].set_ydata(waves[i].plot[1,:])

            ax.set_title(f"Time: {m.dt*0}")

        return lines,


    def animate(i):
        for j in range(len(waves)):
            # Update the wave data to its initial values:
            waves[j].march()
            lines[j].set_ydata(waves[j].plot[1,:])
            ax.set_title(f"Time: {np.around(m.dt*i,2)} ; dt = {m.dt}")

        return lines,


    # Run the animation
    anim = animation.FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=m.Nt,
                                   interval=1,
                                   blit=False)
    return anim