import numpy as np

def gen_x(dx, L):

    # Get number of grid elements required - odd if boundary types are same, even if not
    N  = int((L[1] - L[0])/dx + 1)

    # Now creating the arrays:
    x = np.linspace(L[0], L[1], N)
    return x