import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from classes import Diffusion

# Global matrices:
diff = Diffusion(nelem=10, T1=1, q0=0, alpha=0.5, kappa=1)


fig, ax = plt.subplots()
anim = diff.run(fig=fig, ax=ax)

f = r"videos/problem1.mp4"
writervideo = animation.FFMpegWriter(fps=5)
anim.save(f, writer=writervideo)

