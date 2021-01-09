
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
mesh_grid=[50,50,50]
# arr=np.arange(64).reshape((4,4,4))
# print(arr[0,0,:])
xGRID_nos=mesh_grid[0]
#print("for shape",mesh_grid)

meshg=struct_mesh(dimensions,mesh_grid)
meshg.generate_mesh()

#dill.dump(meshg, file = open("mesh16-20-15.pickle", "wb"))

solver=solverCLS(meshg)
solver.set_conditions(meshg)
tp=time.perf_counter()
# phi=solver.solve(meshg)
# print("slow tdma taken time:",time.perf_counter()-tp)
# tp=time.perf_counter()
phi=solver.solveByVector(meshg)
# print("old TDMA: ", time.perf_counter()-tp)
# tp=time.perf_counter()
# phi=solver.ADIsolver(meshg)
print("new TDMA: ", time.perf_counter()-tp)

np.save('phipure',phi)

#print("fast tdma taken time:",time.perf_counter()-tp)
#print(np.linalg.norm(p1-phi))


z,y=np.meshgrid(meshg.cvzgrid,meshg.cvygrid)
fig,ax=plt.subplots(figsize=(10,5))
print("hererere" , z.shape, y.shape, phi[0,1:-1,1:-1,int(xGRID_nos/2)].shape)
#PLOT COLOR MESH
cmin=np.min(phi[:,1:-1,1:-1,int(xGRID_nos/2)])
#print(cmin)     
cmax=np.max(phi[:,1:-1,1:-1,int(xGRID_nos/2)])
#print(cmax)
cs=ax.pcolormesh(z,y,np.transpose(phi[0,1:-1,1:-1,int(xGRID_nos/2)]), cmap=cm.RdPu, vmin=cmin,vmax=cmax)

cbar=fig.colorbar(cs)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
ax.set_xlabel("Y axis\nDistance(m)\n",fontsize=16)
ax.set_ylabel("Z axis\nDistance(m)\n",fontsize=16)
ax.set_title("YZ plane view at centroid of domain for Pure Advection\n with gamma=0, uy=1, ux=uz=0",fontsize=16)
ax.text(0,-1/7,'at time= 0 sec\n',fontsize=16)

#plt.show()
#Func TO ANIMATE THE PLOT
def update(t):
    for txt in ax.texts:
        txt.set_visible(False)
    ax.pcolormesh(z,y,np.transpose(phi[t,1:-1,1:-1,int(xGRID_nos/2)]), cmap=cm.RdPu, vmin=cmin,vmax=cmax)
    ax.text(0,-1/6,'at time= %f sec\n'%(t*dt), size=10,fontsize=16)

sim= animation.FuncAnimation(fig,update, frames=range(0,len(phi[:,0,0,0])), interval=500, repeat=False)
#sim.save('uycouse.gif', writer='imagemagick')

plt.show()


# print(phi[2:,3,:-2,:-2])
# phi=solver.solve1()
# phi2=solver.solve2()

# print(" #####dif#####")
# print(phi2[1,1,1:-1,1:-1])

# print("Time Taken:",time.perf_counter()-tp)
