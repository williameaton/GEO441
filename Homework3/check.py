# %matplotlib notebook
from main import model, diffusion, diff_animate
import matplotlib.pyplot as plt
import matplotlib.animation as animation



dx = 1  # Grid spacing
L = (0, 100)  # Domain limits
cp = 1  # Specific heat capacity at constant pressure
rho = 1  # Density
k = 0.7  # Conductivity
method="crank-nicolson"

dtc = [0.4, 0.45, 0.55, 0.6]
num = len(dtc)
ds = []
ls = []

# Create overall figure
fig, ax = plt.subplots(num, figsize=(6,6), sharex=True)
fig.set_tight_layout(True)


# Loop for each different dt
for i in range(num):
    ds.append(diffusion(model=model(dx=dx, L=L, cp=cp, rho=rho, k=k, dt=None, dtc=dtc[i]), method=method))  # Create diffusion obj

    l, = ax[i].plot(ds[i].m.x, ds[i].T[1, :], 'k', linewidth=2)  # Plot initial setup
    ls.append(l)




ani = diff_animate(lines=ls, diff_obj=ds, fig=fig, axes=ax, interval=300, frames=100)

#plt.close()
# Animate and save:
f = f"./{method}_diffusion.mp4"
writervideo = animation.FFMpegWriter(fps=5)
print("Saving video:")
ani.save(f, writer=writervideo)
print(f"Written to {f}")
print("For convenience I have already run and saved these movies. Alternatively if you want to view them interactively you can just run plt.show(). You will also need to uncomment '%matplotlib notebook' with the imports to view the animation in Jupyter, and the plt.close() that I use to supress the plot.")
