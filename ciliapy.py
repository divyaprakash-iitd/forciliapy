from curses.textpad import rectangle
import meshfunctions as mf
import ciliafor as cf

M = mf.rectanglecor(1, 0.1, 2, 10)

a, b, aelem = mf.shapefunctioncoefficients(M)

# Call fortran functions

 