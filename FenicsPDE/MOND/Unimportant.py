#Initial attempt at solving the MOND PDE in a 3D spherical domain
#The code is a modification of the example on 'Nonlinear Poisson Equation' 
#from chapter 3.2 of the "FEniCS Tutorial Vol.1" book


# Warning: from fenics import * will import both ‘sym‘ and
# ‘q‘ from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from fenics import *

# Use SymPy to compute f from the manufactured solution u
import sympy as sym

#defining the x,y,z coordinates from the coordinate array in sympy
x, y, z = sym.symbols('x[0], x[1], x[2]')

#Defining the non-linear term inside of the Divergence on the LHS of the PDE
def q(u):
    "Return nonlinear coefficient"
    return 

u = 1 + x + 2*y
f = - sym.diff(q(u)*sym.diff(u, x), x) - sym.diff(q(u)*sym.diff(u, y), y)
f = sym.simplify(f)
u_code = sym.printing.ccode(u)
f_code = sym.printing.ccode(f)
print('u =', u_code)
print('f =', f_code)


