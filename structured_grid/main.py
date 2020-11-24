
#import module 
from Solver import solverCLS
from Meshing import struct_mesh
from parameters import runtime, dt
import dill
#Import python default module
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.collections import cm
import matplotlib.animation as animation


dimensions=[2,2,2]
mesh_grid=[20,18,15]
# arr=np.arange(64).reshape((4,4,4))
# print(arr[0,0,:])
xGRID_nos=mesh_grid[0]
print("for shape",mesh_grid)

meshg=struct_mesh(dimensions,mesh_grid)
meshg.generate_mesh()

#dill.dump(meshg, file = open("mesh20-18-15.pickle", "wb"))

solver=solverCLS(meshg)
solver.set_conditions(meshg)
tp=time.perf_counter()
# phi=solver.solve(meshg)
# print("slow tdma taken time:",time.perf_counter()-tp)
# tp=time.perf_counter()
phi=solver.solveByVector(meshg)
np.save('phi',phi)

print("fast tdma taken time:",time.perf_counter()-tp)
#print(np.linalg.norm(p1-phi))
#print(phi[-1,int(zGRID_nos/2), :,:])

y,z= np.meshgrid(meshg.cvygrid,meshg.cvzgrid)
fig,ax=plt.subplots(figsize=(10,5))

#PLOT COLOR MESH
cmin=np.min(phi[:,1:-1,1:-1,int(xGRID_nos/2)])
print(cmin)     
cmax=np.max(phi[:,1:-1,1:-1,int(xGRID_nos/2)])
print(cmax)
cs=ax.pcolormesh(y,z,phi[0,1:-1,1:-1,int(xGRID_nos/2)], cmap=cm.RdPu, vmin=cmin,vmax=cmax)

cbar=fig.colorbar(cs)

ax.set_xlabel("Distnace (m)")
ax.set_ylabel("Distance (m)")
ax.set_title("2D Advection-Diffusion FVM Solver\n with Dirichlet boundary Condition at x=0 and y=0")

#Func TO ANIMATE THE PLOT
def update(t):
    for txt in ax.texts:
        txt.set_visible(False)
    ax.pcolormesh(y,z,phi[t,1:-1,1:-1,int(xGRID_nos/2)], cmap=cm.RdPu, vmin=cmin,vmax=cmax)
    ax.text(0,-1/8,'at time= %f sec\nux=2, uy=2 & gamma=3'%(t*dt), size=10)

sim= animation.FuncAnimation(fig,update, frames=range(0,len(phi[:,0,0,0])), interval=500, repeat=False)
plt.show()


# print(phi[2:,3,:-2,:-2])
# phi=solver.solve1()
# phi2=solver.solve2()

# print(" #####dif#####")
# print(phi2[1,1,1:-1,1:-1])

# print("Time Taken:",time.perf_counter()-tp)
