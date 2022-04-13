import numpy as np
import matplotlib.pyplot as plt

def load_trace_SV(fname):
    data = np.loadtxt(f"{fname}/seismograms_01_2")
    return data[:,0], data[:,1]


# In order: NEX, NEZ, NGLLX, NGLLZ, NO.GRID.PTS (P), NO.GRID.PTS (S)
sims = [[100, 60, 7, 7, 13, 7], [50, 30, 7, 7, 6, 3], [25, 15, 7, 7, 3, 1], [25, 15, 14, 14, 7,3]]


fig, ax = plt.subplots()
legend_str = []
for i in [0, 3]:
    s = sims[i]
    time, trace = load_trace_SV(fname=f"./compare_seismograms/SH/NEX{s[0]}_NEZ{s[1]}_NGLLX{s[2]}_NGLLZ{s[3]}")
    a = 1
    if i < 2:
        a = 1

    ax.plot(time,trace, alpha=a)
    legend_str.append(f"NEX: {s[0]}, NEZ: {s[1]}, NGLLX: {s[2]}, NGLLZ: {s[3]}")

ax.set_xlabel("Time [s]")
ax.set_xlim([25, 50])
ax.set_ylabel("Displacement")
ax.legend(legend_str)
plt.show()