from curses.textpad import rectangle
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
niter = 10000
dt = 0.001
mu = 1000

# Tip force
Ftip = 100

plt.triplot(x, y, M.simplices)
# plt.plot(x, y, 'go')

for i in range(10000):
    FN = np.zeros((M.npoints, 2), order='F')
    FN[top, 0] = Ftip/(Nx - 1)
    # FN[36, 0] = Ftip/(Nx - 1)/2
    # FN[39, 0] = Ftip/(Nx - 1)/2

    solvefem.solvefem.computeforces(M.simplices + 1, np.array([x, y]).T, b, co, K, aelem, FN)

    # Fix the bottom nodes (Make force zero)
    FN[bottom, 1] = 0.0

    # Calculate the deformation
    for inode in range(M.npoints):
        x[inode] = x[inode] + dt/mu * FN[inode, 0]
        y[inode] = y[inode] + dt/mu * FN[inode, 1]

    # Fix the bottom nodes (Make displacement zero)
    y[bottom] = 0.0


elasped = time.time() - t1
plt.triplot(x, y, M.simplices)
# plt.plot(x, y, 'go')
plt.axis('equal')
plt.show()
