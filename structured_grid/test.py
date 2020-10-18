from Solver import solver
from Meshing import struct_mesh
from BoundaryCondi import FinalMat
import numpy as np
import time

if __name__=='__main__':
    dimensions=[5,5,5]
    mesh_elem=[20,20,20]

    print("for shape",mesh_elem)
    tp=time.perf_counter()
    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    ap,aw,ae,ad,au,ab,af,ap0,b,c=FinalMat(meshg)
    solver=solver(ap,aw,ae,ad,au,ab,af,ap0,b,c)
    solver.solve()
    print("Time Taken:",time.perf_counter()-tp)