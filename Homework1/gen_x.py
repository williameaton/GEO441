import numpy as np

bc_bool = {"neumann":   0,
           "dirichlet": 1}

def gen_x(dx, L, BC_left, BC_right):

    # Get number of grid elements required - odd if boundary types are same, even if not
    N  = int((L[1] - L[0])/dx + 1)

    if (bc_bool[BC_left] + bc_bool[BC_right]) == 1: N += 1

    # Now creating the arrays:
    x = np.linspace(0, L[1], N)
    return x