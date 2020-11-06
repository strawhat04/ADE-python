from Solver import solver
from Meshing import struct_mesh
from StructuralDiscre import LinFlux
import numpy as np
import time

if __name__=='__main__':
    dimensions=[1,1,1]
    mesh_elem=[6,6,6]
    # arr=np.arange(64).reshape((4,4,4))
    # print(arr[0,0,:])

    print("for shape",mesh_elem)
    # tp=time.perf_counter()
    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    ap,aw,ae,ad,au,ab,af,ap0,b,c=LinFlux(meshg)
    solver=solver(ap,aw,ae,ad,au,ab,af,ap0,b,c)
    phi=solver.solve1()
    phi2=solver.solve2()
    print(phi[1,1,1:,1:])
    print(" #####dif#####")
    print(phi2[1,1,1:-1,1:-1])

    # print("Time Taken:",time.perf_counter()-tp)
