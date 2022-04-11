# All of the snapshots have the first same two columns (x and z coordinates)
# So plan is to have a script that strips out those columns:
import numpy as np

# Data directories:
dir = "../OUTPUT_FILES/"
new_dir = f"{dir}/reduced/"

# _____________________________________________________________________________________________________________________
# LOADING SNAPSHOTS:
for i in range(40):
    # Ensure snapshot number has sufficient 0's in front:
    snapshot_no = str((i+1)*250)
    while len(snapshot_no) <6:
        snapshot_no  = "0" + snapshot_no

    # Load snapshot:
    data = np.loadtxt(f"{dir}/snapshot_{snapshot_no}")
    out_name = f"{new_dir}/snapshot_{snapshot_no}"
    np.savetxt(fname=out_name , X=data[:,2:])
    print(f"Saved {out_name}")

# Finally, output x and z array:
np.savetxt(fname=f"{new_dir}/coords" , X=data[:,:2])
