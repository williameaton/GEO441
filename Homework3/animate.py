import numpy as np
import matplotlib.animation as animation


def diff_animate(lines, diff_obj, fig, axes, interval, frames):
    num = len(lines)                                     # Number of lines to animate
    fs = 16                                              # Font size

    # Set plot metadata - hardcoded:
    axes[-1].set_xlabel("X (m)", fontsize=fs)
    axes[-1].tick_params(axis='x', labelsize=fs)
    for k in range(num):
        axes[k].set_xlim([0,100])
        axes[k].set_ylabel("Temperature", fontsize=fs)
        axes[k].tick_params(axis='y', labelsize=fs)


    # Internal animation function
    def animate(k):
        k_0 = 5 # Wait 5 frames before starting
        if k < k_0:
            for i in range(num):
                # Just update the axes titles
                axes[i].set_title(f"dt: {np.around(diff_obj[i].m.dt, 2)} - Timestep: {k - k_0}", fontsize=16)
        else:
            for i in range(num):
                lines[i].set_ydata(diff_obj[i].T[0, :]) # Update the line object with correct T data
                axes[i].set_title(f"dt: {np.around(diff_obj[i].m.dt,2)} - Timestep: {k-k_0}", fontsize=16)
                diff_obj[i].march() # Calculate next T timestep
        return lines

    def init():
        for i in range(num):
            diff_obj[i].set_IC()                    # re-initialise T conditions
            lines[i].set_ydata(diff_obj[i].T[1, :]) # Update T lines on animation
            axes[i].set_title(f"dt: {np.around(diff_obj[i].m.dt,2)} - Timestep: -10", fontsize=16)
        return lines

    # calling the animation function
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=interval, blit=False)
    return anim