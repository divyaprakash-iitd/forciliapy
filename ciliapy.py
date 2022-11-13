from curses.textpad import rectangle
import numpy as np
import meshfunctions as mf
import solvefem

L = 0.2
B = L/4
Nx = 4
Ny = 10

M = mf.rectanglecor(L, B, Nx, Ny)

a, b, aelem = mf.shapefunctioncoefficients(M)

# Call fortran functions
# FN = np.zeros((aelem.shape[0], 2))
co = 5000
K = 10000
# Increment the indices by 1 to comply with fortran
for i in range(2):
    fn, d = solvefem.solvefem.computeforces(M.simplices + 1, M.points, b, co, K, aelem)
    print(i)
 