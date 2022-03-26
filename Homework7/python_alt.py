import matplotlib.pyplot as plt

from krichner_codes import *

# Constants:
NSPEC = 12
NGLL = 7
NGLOB = (NGLL-1)*NSPEC + 1

# Simulation params:
NSTEPS = 3000
time_step = 1000000000

# Model parameters:
length        = 3e3
density       = 2.5e3
conductivity  = 10e-1
heat_capacity = 0.3e3


# Get gll points, weights and lagrange derivatives:
gllxi, wgll = gll(NGLL-1)
hprime = lagrange1st(NGLL-1)




x1 = np.arange(NSPEC)*(length/NSPEC)
x2 = (1+ np.arange(NSPEC))*(length/NSPEC)

# Initialise the model value for the homogenous domain:
rho = np.zeros((NGLL, NSPEC)) + density
kappa = np.zeros((NGLL, NSPEC)) + conductivity
cp = np.zeros((NGLL, NSPEC)) + heat_capacity
jacobian = np.zeros((NGLL, NSPEC))
jac_inv = np.zeros((NGLL, NSPEC))

for i in range(NSPEC):
    for j in range(NGLL):
        jac_inv[j,i] = 2/(x2[i] - x1[i])
        jacobian[j,i] = (x2[i] - x1[i])/2

# Set up global numbering system:
ibool = np.zeros((NGLL, NSPEC), dtype=int)
x = np.zeros(NGLOB)
iglob = 0
for ispec in range(NSPEC):
    for i in range(NGLL):
        if i>0:
            iglob += 1
        ibool[i, ispec] = int(iglob)
        x[iglob] = 0.5*(1 - gllxi[i])*x1[ispec] + 0.5*(1 + gllxi[i])*x2[ispec]

mass_global = np.zeros(NGLOB)
# Construct global mass matrix:
for ispec in range(NSPEC):
    for a in range(NGLL):
        for b in range(NGLL):
            k = ibool[a,ispec]
            mass_global[k] = mass_global[k] + wgll[b]*rho[b,ispec]*cp[b,ispec]*jacobian[b,ispec]


M = np.diag(mass_global)

temperature = np.zeros(NGLOB)
dTdt = np.zeros(NGLOB)

# Set initial conditions:
temperature[0] = 10.0

# Loop in time:
K = np.zeros((NGLOB,NGLOB))

fig,ax = plt.subplots()
ax.plot(x,temperature)

rhs = np.zeros(NGLOB)

for itime in range(NSTEPS):
    print(itime)
    temperature_tilde = temperature + 0.5*time_step*dTdt

    # Enforce boundary conditions:
    temperature_tilde[0] = 10
    temperature_tilde[-1] = 0


    rhs[:] = 0
    for ispec in range(NSPEC):
        for a in range(NGLL):
            for b in range(NGLL):
                kab = 0
                for g in range(NGLL):
                    kab += wgll[g]*kappa[g,ispec]*jac_inv[g,ispec]*hprime[a,g]*hprime[b,g]

                rhs[ibool[a,ispec]] += -kab*temperature_tilde[ibool[b,ispec]]

    for v in range(NGLOB):
        dTdt[v] = rhs[v]/M[v,v]

    temperature = temperature_tilde + 0.5*time_step*dTdt

    ax.plot(x, temperature)

plt.show()