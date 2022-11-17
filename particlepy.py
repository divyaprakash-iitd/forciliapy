import time
import numpy as np
import matplotlib.pyplot as plt
import meshfunctions as mf
import solvefem

t1 = time.time()

radius = 1
N = 10

M, bpoints = mf.particlecor(radius, N)
x = M.points[:, 0]
y = M.points[:, 1]

xorg = x
yorg = y

boundary = np.isin(M.points, bpoints)[:, 0]

a, b, aelem = mf.shapefunctioncoefficients(M)

# Call fortran functions
co = 5000
K = 10000

# Run time iterations
niter = 10000
dt = 0.0001
mu = 10000

# Tip force
Ftip = 10000

plt.triplot(x, y, M.simplices)
# plt.plot(x, y, 'go')

for i in range(niter):
    FN = np.zeros((M.npoints, 2), order='F')
    if i < 2*niter/3:
        FN[boundary, 0] = Ftip * 2 * bpoints[:, 0]/radius
        FN[boundary, 1] = Ftip * bpoints[:, 1]/radius
    # FN[36, 0] = Ftip/(Nx - 1)/2
    # FN[39, 0] = Ftip/(Nx - 1)/2

    solvefem.solvefem.computeforces(M.simplices + 1, np.array([x, y]).T, b, co, K, aelem, FN)

    # Calculate the deformation
    for inode in range(M.npoints):
        x[inode] = x[inode] + dt/mu * FN[inode, 0]
        y[inode] = y[inode] + dt/mu * FN[inode, 1]

    if i % 100 == 0:
        plt.figure(figsize=(6, 6))
        # plt.triplot(M.points[:, 0], M.points[:, 1], M.simplices)
        plt.triplot(xorg, yorg, M.simplices)
        # plt.triplot(x, y, M.simplices)
        plt.xlim(-2, 2)
        plt.ylim(-1.5, 1.5)
        plt.savefig('images/plot'+str(i)+'.png')
        plt.clf()
        plt.close() 
elasped = time.time() - t1
plt.triplot(xorg, yorg, M.simplices)
# plt.triplot(x, y, M.simplices)
# plt.plot(x, y, 'go')
# plt.axis('equal')
plt.show()
