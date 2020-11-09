from Solver import Solver
from Meshing import struct_mesh
import numpy as np
import time

if __name__=='__main__':
    dimensions=[1,1,1]
    mesh_elem=[30,30,30]
    # arr=np.arange(64).reshape((4,4,4))
    # print(arr[0,0,:])

    print("for shape",mesh_elem)

    meshg=struct_mesh(dimensions,mesh_elem)
    meshg.generate_mesh()
    solver=Solver(meshg)
    solver.set_conditions(meshg)
    tp=time.perf_counter()
    p1=solver.solve2(meshg)
    print("slow tdma taken time:",time.perf_counter()-tp)
    tp=time.perf_counter()
    p2=solver.solve(meshg)
    print("fast tdma taken time:",time.perf_counter()-tp)
    print(np.linalg.norm(p1-p2))

    # print(phi[2:,3,:-2,:-2])
    # phi=solver.solve1()
    # phi2=solver.solve2()

    # print(" #####dif#####")
    # print(phi2[1,1,1:-1,1:-1])

    # print("Time Taken:",time.perf_counter()-tp)
