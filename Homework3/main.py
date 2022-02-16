import matplotlib.pyplot as plt
import numpy as np
import timeit
import matplotlib.animation as animation
import matplotlib.pyplot as plt


def create_A_matrix(k, prefac, dim):
    # Create A matrix:
    A = np.zeros((dim, dim))
    for i in range(1, dim - 1):
        A[i, i - 1] = k[i] * prefac[i]
        A[i, i] = 1 - (k[i + 1] + k[i]) * prefac[i]
        A[i, i + 1] = (k[i + 1]) * prefac[i]

    # boundary matrix vals:
    A[0, 0] = (1 - k[1] - k[0]) * prefac[0]
    A[0, 1] = k[1] * prefac[0]
    A[-1, -1] = (1 - k[dim - 1] - k[dim - 2]) * prefac[-1]
    A[-1, -2] = (k[dim - 1]) * prefac[-1]

    return A

def forward_matrix_march(T, A):
    T[0, :] = np.matmul(A, np.transpose(T[1, :]))
    return T

#def march_loop(T, prefac, k):
#    for i in range(1, len(x)-1):
#        T[0, i] = T[1, i] + prefac[i]*( (k[i+1] - k[i])*(T[1, i+1] - T[1, i]) + k[i]*(T[1, i+1] - 2*T[1,i] + T[1,i-1]) )
#    return T


dx  = 1                      # Grid spacing
x   = np.arange(0, 101, dx)  # X dimension
cp  = x*0 + 1                # Specific heat capacity at constant pressure
rho = x*0 + 1                # Density
k = x*0 + 0.7                  # Density

dim = len(x)

dt = 0.01

prefac = dt/(rho*cp*dx*dx)
T   = np.zeros((2, dim))

# set initial condition:
T[:, np.where(x==50)] = 1

A = create_A_matrix(k, prefac, dim)


fig, ax = plt.subplots()

l, = ax.plot(x, T[1,:], 'k')



def animate(T):
    T = forward_matrix_march(T, A)
    l.set_ydata(T[0,:])
    T = T[::-1, :]
    return T,


def init():
    T[:, :] = 0
    T[:, np.where(x == 50)] = 1
    l.set_ydata(T[1,:])
    return l,

# calling the animation function
anim = animation.FuncAnimation(fig, animate(T),
                            init_func = init,
                            frames = 500,
                            interval = 20,
                            blit = True)


plt.show()