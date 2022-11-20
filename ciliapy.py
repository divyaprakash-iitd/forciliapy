import time
import numpy as np
import matplotlib.pyplot as plt
import meshfunctions as mf
import solvefem

t1 = time.time()

L = 0.2
B = L/4
Nx = 4
Ny = 10

M = mf.rectanglecor(L, B, Nx, Ny)
x = M.points[:, 0]
y = M.points[:, 1]
# Get the top and bottom boundaries
top = y == max(y)
bottom = y == min(y)

a, b, aelem = mf.shapefunctioncoefficients(M)

# Call fortran functions
# FN = np.zeros((aelem.shape[0], 2))
co = 5000
K = 10000
# Increment the indices by 1 to comply with fortran


# Run time iterations
dt = 0.001
mu = 1000

# Tip force
Ftip = 100

plt.triplot(x, y, M.simplices)
boundary = np.vstack((top, bottom)).T
FN = np.zeros((M.npoints, 2), order='F')
P = np.zeros((M.npoints, 2), order='F')
P[:, 0] = x
P[:, 1] = y

for niter in range(1, 10000, 2000):
    solvefem.solvefem.solve(M.simplices + 1, P, b, co, K, aelem,
                            boundary, FN, niter)
    plt.triplot(P[:, 0], P[:, 1], M.simplices)

elasped = time.time() - t1
plt.axis('equal')
plt.show()
