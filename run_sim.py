import time
import numpy as np
import matplotlib.pyplot as plt
import meshfunctions as mf
# import solvefem
import ibmc

DL = 0.01
DH = 0.002
L = DH/2
B = DL/10
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
co = 1/10000.0
K = 1/50000.0
E = co*(2*co+3*K)/(co+K)
nu = K/(2*(co+K))
print(E)
print(nu)

boundary = np.vstack((top, bottom)).T
FN = np.zeros((M.npoints, 2), order='F')
UN = FN
P = np.zeros((M.npoints, 2), order='F')
P[:, 0] = x
P[:, 1] = y

# plt.triplot(x,y,M.simplices)

ibmc.ibmc(M.simplices + 1, P, b, boundary, co, K, aelem, FN, UN)
