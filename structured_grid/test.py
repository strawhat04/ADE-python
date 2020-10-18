from Solver import solver
from Meshing import struct_mesh
from BoundaryCondi import FinalMat
import numpy as np

if __name__=='__main__':
    dimensions=[2,2,2]
    mesh_elem=[5,5,5]
    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    ap,aw,ae,ad,au,ab,af,ap0,b,c=FinalMat(meshg)
    solver=solver(ap,aw,ae,ad,au,ab,af,ap0,b,c)
    solver.solve()