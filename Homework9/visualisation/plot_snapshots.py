import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata as gd
import matplotlib.animation as animation
import matplotlib as mpl
from general_funcs import normalise

def slice_by_time():
    pass

def interp_imshow(i, data, chl, grid, coords, xarr, zarr):
    return np.reshape(gd(points=coords, values=data[i, :, chl], xi =grid, method='linear'), (len(xarr), len(zarr)))

def init():
    im.set_data(interp_imshow(0, snapshots, chl=chl, grid=grid_2d, coords=coords, xarr=x_arr, zarr=z_arr))
    set_axtitle(ax[0], 0)

    # Update traces:
    for j in range(3):
        time_slice = [0,0]
        disp = np.array([0,0])+(2*j)
        lines[j].set_xdata(time_slice)
        lines[j].set_ydata(disp)

    return [im]

def animate(i):
    im.set_array(interp_imshow(i, snapshots, chl=chl, grid=grid_2d, coords=coords, xarr=x_arr, zarr=z_arr))
    set_axtitle(ax[0], i)
    mask = time <= i * dt * sample_DT
    time_slice = time[mask]

    # Update traces:
    for j in range(3):
        disp = traces[j, :][mask]
        lines[j].set_xdata(time_slice)
        lines[j].set_ydata(disp + 2*j)


    return [im]

def set_axtitle(ax, i):
    ax.set_title(f"Channel: {chls[chl]}; time = {i*sample_DT*5.0e-3} s")


# Data directories:
dir = "../OUTPUT_FILES/"
nx = 500
nz = 500
chl = 1
frames = 400
sample_DT = 25
dt = 0.005

chls = ["X", "Y", "Z"]

# Get coordinate data
coords = np.loadtxt(fname=f"{dir}/reduced/coords")
x_min = np.min(coords[:,0])
x_max = np.max(coords[:,0])
z_min = np.min(coords[:,1])
z_max = np.max(coords[:,1])

# Generate even grid for outputs:
x_arr = np.linspace(x_min, x_max, nx)
z_arr = np.linspace(z_min, z_max, nz)
X, Z = np.meshgrid(x_arr, z_arr)
grid_2d = (X.flatten(), Z.flatten())

# _____________________________________________________________________________________________________________________
# LOADING SNAPSHOTS:
snapshots = []
for i in range(frames):
    # Ensure snapshot number has sufficient 0's in front:
    snapshot_no = str((i+1)*sample_DT)
    while len(snapshot_no) <6:
        snapshot_no  = "0" + snapshot_no
    # Append data to array:
    snapshots.append(np.loadtxt(f"{dir}/reduced/snapshot_{snapshot_no}"))
    print(f"Loading snapshot: {i}")
# Convert to numpy array:
snapshots = np.array(snapshots)
# _____________________________________________________________________________________________________________________
# PLOTTING
fig, ax = plt.subplots(2, figsize=(13, 7))
fig.set_tight_layout(True)

traces = []
lines  = []
# Plot seismograms:
for j in range(3):
    d = np.loadtxt(f"{dir}/seismograms_01_{j+1}")
    time = d[:,0]
    traces.append(normalise((d[:,1])))
    ax[1].plot(d[:,0], normalise(d[:,1])+(2*j), 'k')
    if np.amax(np.abs(d[:,1])) > 0:
        ax[1].annotate(f"Channel: {chls[j]}", xy=(130, 205-(60*j)), xycoords='figure points')

    l, = ax[1].plot(np.array([0, 0]), np.array([0, 0]) + 2 * j, 'r')
    lines.append(l)


# Store traces for animation
dt = time[1]-time[0]
traces = np.array(traces)

ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].get_yaxis().set_visible(False)
ax[1].set_xlabel("Time [s]")
ax[1].set_xlim([d[0,0], d[-1,0]])

ax[0].set_xlabel("X [m]")
ax[0].set_ylabel("Z [m]")
set_axtitle(ax[0], 0)

ax[0].plot([0.25],[0.40], '*r') # Source
ax[0].plot([0.750000000],[0.40], 'r', marker=mpl.markers.CARETDOWNBASE) # Receiver

im  = ax[0].imshow(interp_imshow(0, snapshots, chl=chl, grid=grid_2d, coords=coords, xarr=x_arr, zarr=z_arr),
                origin='lower', extent=(x_min, x_max, z_min, z_max), cmap='RdBu', vmin=-0.2, vmax=0.2)

cbar = fig.colorbar(im, ax=ax[0])
cbar.minorticks_on()

anim = animation.FuncAnimation(fig, animate, init_func = init, frames = frames, interval = 1, blit = True)


f = r"SH.mp4"
writermp4 = animation.FFMpegWriter(fps=40)
anim.save(f, writer=writermp4)