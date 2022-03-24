# Quick python script that is used to plot the resulting output from the fortran codes: 
from time import time
import numpy as np 
import matplotlib.pyplot as plt 

# Load a list of all of the times outputted: 
times = np.loadtxt(f"snapshots/time")
no_snapshots = len(times)
data = []
for j in range(no_snapshots):
    time_name = str(int(times[j]))
    while len(time_name) <5: 
        time_name = "0" + time_name # F90 puts 0's in front of times so need to append any 
    data.append(np.loadtxt(f"snapshots/snapshot{time_name}"))

# Data now has dimensions [no. snapshots  x  len(x dim)  x  2] where the [, :, 0] is the x array and [, :, 1] is Temp
data = np.array(data)


# Create plot:
fig, ax = plt.subplots(figsize=(8,7))
for j in range(no_snapshots):
    ax.plot(data[j,:,0], data[j,:,1], alpha = 1 - 0.7*j/no_snapshots)


# Titles, labels etc...
fs = 16 # Label font size
ax.set_xlabel("X [m]", fontsize=fs)
ax.set_ylabel(r"Temperature [$^o C$]", fontsize = fs)

plt.show()