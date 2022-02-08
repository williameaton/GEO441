import numpy as np

def gen_x(dx, L):
    # Generate x-domain array
    return np.linspace(L[0], L[1], int((L[1] - L[0])/dx + 1))