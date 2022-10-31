from ast import Del
from matplotlib.pyplot import axis
import numpy as np
from scipy.spatial import Delaunay

# Generate points for a triangle


def rectanglecor(L, B, Nx, Ny):
    """Genrates the delaunay triangulation for a rectangle

    Args:
        L (float): Rectangle's Length
        B (float): Rectangle's Breadth
        Nx (int): Number of nodes along x-direction
        Ny (int): Number of nodes along y-direction

    Returns:
        tuple: A tuple containing the Delaunay triangulation and the
        points defining the triangle
    """
    X = np.linspace(0, B, Nx)
    Y = np.linspace(0, L, Nx)
    X, Y = np.meshgrid(X, Y)

    X = X.reshape((X.size, -1))
    Y = Y.reshape((Y.size, -1))

    M = Delaunay(np.concatenate(X, Y), axis=1)

    return M, X, Y


def cofactor(A):
    return np.linalg.inv(A).T * np.linalg.det(A)

# Calculate shape function coefficients
def shapefunctioncoefficients(M, x, y):
