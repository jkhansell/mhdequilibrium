from dolfin import *

mesh = UnitSquare(6, 4)
V = FunctionSpace(mesh, "Lagrange", 1)