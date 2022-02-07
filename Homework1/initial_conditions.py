"""import numpy as np

def set_IC(x):
    # Impose initial condition where u_old and u are the same:
    u     = np.exp(-0.1*((x-50)**2)
    v     = np.zeros(len(x))

    # Calculate T (apart from boundary values):
    T     = np.zeros(len(x))
    for i in range(1, T-1):
        T[i] = (kappa[i]/(2*dx))(u[i+1] - u[i-1])


    return u, v, T, zeros


"""