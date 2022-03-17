import numpy as np
import matplotlib.pyplot as plt

# We want to solve for Kd = f but in this local form we want to solve for local elements:
# Define domain parameters:
n = 5      # Number of elements
h = 1 / n  # Element width
T1 = 0     # Boundary condition (right)
q0 = 0     # Boundary condition (left)
f = np.zeros(n) + 1


rho = np.zeros(n) + 1
cp = np.zeros(n) + 2

# Global matrices:
K = np.zeros((n+1,n+1))
F = np.zeros(n+1)
M = np.zeros((n+1,n+1))


# Local shape funcs/arrays:
N1 = np.array([1, 0])
N2 = np.array([0, 1])
k = np.zeros((2,2))

for n in range(n):

    # Construct local k_ab matrix:
    for a in range(2):
        for b in range(2):
            k[a][b] = (1/h)* (-1)**(a+b)

    rhs = (N1+N2)*f[n:n+2] * h/2 # This allows f to be spatially varying
    # Assemble:
    K[n:n+2, n:n+2] += k
    F[n:n+2] += rhs

    M[n:n+2, n:n+2] += rho[n]*cp[n]*h*(1/6) * np.array([[2, 1], [1, 2]])

print(F)
# Add BCs:
F[0] += q0
F[-2] -= k[0][1]*T1



# NEED TO QUERY THE ELEMENTAL FORCING MATRIX WITH LUCAS + ADD IN TIME MARCHING

D = np.matmul(np.linalg.inv(K[:-1,:-1]), F[:-1])

print(D)

fig, ax = plt.subplots()

ax.plot(D)
plt.show()