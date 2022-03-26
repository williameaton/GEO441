import numpy as np
class plotter():

    def __init__(self, directory):
        self.dir = directory
        # Get meta data
        self.get_meta()
        # Load data
        self.load_data()

    def get_meta(self):
        Lines = open(f'{self.dir}/meta', 'r').readlines()
        self.name = Lines[0].strip()
        self.no_snapshots = int(Lines[1])
        self.timestep = float(Lines[2])


    def load_data(self):
        # Returns data with dimensions [no. snapshots  x  len(x dim)  x  2] where the [, :, 0] is the x array and [, :, 1] is Temp

        # Load a list of all of the times outputted:
        times = np.loadtxt(f"{self.dir}/time")
        assert (self.no_snapshots == len(times))
        self.times = times

        data = []

        self.slicename = f"{self.dir}/{self.name}"
        print(f"Loading data from {self.slicename} slices...")
        for j in range(self.no_snapshots):
            time_name = str(int(times[j]))
            while len(time_name) <8:
                time_name = "0" + time_name # F90 puts 0's in front of times so need to append any
            data.append(np.loadtxt(f"{self.slicename}{time_name}"))

        self.data = np.array(data)


    def plot_initial_setup(self, ax):
        fs = 16  # Label font size
        ax.set_xlabel("X [m]", fontsize=fs)
        ax.set_ylabel(r"Temperature [$^o C$]", fontsize=fs)

        # Get initial setup for plot
        self.init = self.data[0, :, 1]
        self.init[1:] = 0

        self.l, = ax.plot(self.data[0, :, 0], self.init, 'o-', alpha=1)



def secs_to_a(secs):
    # Convert seconds to yrs:
    return secs/(60*60*24*365.25)