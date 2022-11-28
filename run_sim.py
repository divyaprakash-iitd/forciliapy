import time
import numpy as np
import matplotlib.pyplot as plt
import meshfunctions as mf
# import solvefem
import ibmc

DL = 0.01
DH = 0.002
L = 0.7e-3
B = 0.1e-3
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

# Generate a particle
MP, _ = mf.particlecor(L/4.0, 12)
x = MP.points[:, 0]
y = MP.points[:, 1]
PP = np.zeros((MP.npoints, 2), order='F')
PP[:, 0] = x
PP[:, 1] = y
boundary[:, :] = False
_, pb, paelem = mf.shapefunctioncoefficients(MP)

# subroutine ibmc(MM,PTM,PP,PTP,b,PTB,boundary,co,kval,aelem,ptaelem,FN,UN)
# ibmc.ibmc(M.simplices + 1, P, b, boundary, co, K, aelem, FN, UN)
ibmc.ibmc(M.simplices + 1, MP.simplices + 1, P, PP, b, pb, boundary, co, K,
          aelem, paelem, FN, UN)
