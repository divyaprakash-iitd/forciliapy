from curses.textpad import rectangle
import numpy as np
import matplotlib.pyplot as plt
import meshfunctions as mf
import solvefem

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
Ftip = 10

plt.triplot(x, y, M.simplices)
plt.plot(x, y, 'o')
plt.axis('equal')
plt.show()

FN = np.zeros((M.npoints, 2), order='F')
for i in range(10):
    solvefem.solvefem.computeforces(M.simplices + 1, np.array([x, y]).T, b, co, K, aelem, FN)
    
 