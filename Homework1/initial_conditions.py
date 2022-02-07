import numpy as np

def set_IC(x):
    # Impose initial condition where u_old and u are the same:
    u     =  np.exp(-0.1*((x-50)**2) )
    zeros  = np.zeros(len(x))
    v = np.exp(-0.1*((x-50)**2) )
    T = zeros
    return u, v, T, zeros




def set_IC2(x):
    # Impose initial condition:
    u     =  np.exp(-0.1*((x-20)**2) )*np.sin(x)
    zeros  = np.zeros(len(x))
    v = zeros
    T = zeros
    return u, v, T, zeros