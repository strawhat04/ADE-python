from Meshing import struct_mesh
from StructuralDiscre import LinFlux
#from Discretization import LinMat
import numpy as np
from matplotlib.collections import cm
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from math import exp, sin, cos, atan, sqrt

def tr(a):
	return np.transpose(a)

runtime=1
dt=0.1

xGRID_nos=20
yGRID_nos=20
zGRID_nos=20

mesh=struct_mesh([1,1,1], [xGRID_nos,yGRID_nos,zGRID_nos])
mesh.generate_mesh()
print("grid: ",mesh.xgrid)
print("cv: ", mesh.cvxgrid)

iphi=0
ap,aw,ae,ad,au,ab,af,ap0,b,c=LinFlux(mesh)

# print("aw",aw)
# print("ae",ae)
# print("ad",ad)
# print("au",au)
# print("ap0",ap0)
# print("ap",ap)

phi=iphi*np.ones((int(runtime/dt)+1,zGRID_nos+2, yGRID_nos+2, xGRID_nos+2))	 #include 2 ghost nodes for including source terms outside domain 

#initial phi
phi[0,int(zGRID_nos/2-zGRID_nos/5):int(zGRID_nos/2+zGRID_nos/5),int(yGRID_nos/2-yGRID_nos/5):int(yGRID_nos/2+yGRID_nos/5),int(xGRID_nos/2-xGRID_nos/5):int(xGRID_nos/2+xGRID_nos/5)]=10

#BOUNDARY CONDITION 

#array storing face area of each control volume element in face normal to x axis
Sx=np.dot((mesh.cvygrid[1:]-mesh.cvygrid[:-1]).reshape(-1,1), (mesh.cvzgrid[1:]-mesh.cvzgrid[:-1]).reshape(1,-1))
#array storing face area of each control volume element in face normal to y axis
Sy=np.dot((mesh.cvxgrid[1:]-mesh.cvxgrid[:-1]).reshape(-1,1), (mesh.cvzgrid[1:]-mesh.cvzgrid[:-1]).reshape(1,-1))
#array storing face area of each control volume element in face normal to z axis
Sz=np.dot((mesh.cvygrid[1:]-mesh.cvygrid[:-1]).reshape(-1,1), (mesh.cvxgrid[1:]-mesh.cvxgrid[:-1]).reshape(1,-1))

#print("Imput Boundary Conditions as: D/N x0 y0 z0 xn yn zn")
bc,x0,y0,z0,xn,yn,zn=['d',0,0,0,0,0,0]

#Computing flux coefficient for dirichlet condition
if bc=='D' or bc=='d':
	phi[:,:,:,1]=x0
	phi[:,:,1,:]=y0
	phi[:,1,:,:]=z0
	phi[:,:,:,xGRID_nos]=xn
	phi[:,:,yGRID_nos,:]=yn
	phi[:,zGRID_nos,:,:]=zn

#Computing flux coefficient for Neumann condition
if bc=='N' or bc=='n':
	aw[:,:,0]=ae[:,:,nx]=0
	b[:,:,0]+=float(x0)*Sx
	b[:,:,nx]+=float(xn)*Sx

	ad[:,0,:]=au[:,ny,:]=0
	b[:,0,:]+=float(y0)*Sy
	b[0,ny,:]+=float(yn)*Sy

	ab[0,:,:]=af[nx,:,:]=0
	b[0,:,:]+=float(z0)*Sz
	b[nx,:,:]+=float(zn)*Sz
print("phi", phi[0,int(zGRID_nos/2),:,:])
d=np.zeros((zGRID_nos,yGRID_nos,xGRID_nos))

R=np.zeros(xGRID_nos+2)		#include 2 ghost nodes for including source terms outside domain 
Q=np.zeros(xGRID_nos+2)

rdue=np.ones((zGRID_nos-2,yGRID_nos-2,xGRID_nos-2))

Q[1]=x0
Q[xGRID_nos]=xn
for t in range(1, len(phi[:,0,0])):
	rdue.fill(1)
	itera=0
	phi[t,:,:,:]=phi[t-1,:,:,:]
	while np.linalg.norm(rdue)>1e-4 : #   
		for k in range(0, zGRID_nos):	
			for j in range(1, yGRID_nos-1):
				for i in range(1, xGRID_nos-1):
					#print("in while")
					d[k,j,i]=b[k,j,i]+ap0[k,j,i]*phi[t-1,k+1,j+1,i+1]+au[k,j,i]*phi[t,k+1,j+1+1,i+1] +ad[k,j,i]*phi[t,k+1,j-1+1,i+1] + af[k,j,i]*phi[t,k+1+1,j+1,i+1] +ab[k,j,i]*phi[t,k+1-1,j+1,i+1] 
					
					R[i+1]=ae[k,j,i]/(ap[k,j,i]-aw[k,j,i]*R[i])

					Q[i+1]=(d[k,j,i]+aw[k,j,i]*Q[i])/(ap[k,j,i]-aw[k,j,i]*R[i])
					
				for i in range(xGRID_nos-1,-1,-1):	
					#print("afd", Q[i+1], R[i+1])
					phi[t,k+1,j+1,i+1]=R[i+1]*phi[t,k+1,j+1,i+2]+Q[i+1]

		#print(phi[t,1:-3,2:-2,2:-2].shape)
		rdue=phi[t,2:-2,2:-2,2:-2]-np.divide((np.multiply(aw[1:-1,1:-1,1:-1],phi[t,2:-2,2:-2,1:-3])+np.multiply(ae[1:-1,1:-1,1:-1],phi[t,2:-2,2:-2,3:-1])+np.multiply(au[1:-1,1:-1,1:-1],phi[t,2:-2,3:-1,2:-2])+np.multiply(ad[1:-1,1:-1,1:-1], phi[t,2:-2,1:-3,2:-2])+np.multiply(af[1:-1,1:-1,1:-1],phi[t,3:-1,2:-2,2:-2])+np.multiply(ab[1:-1,1:-1,1:-1], phi[t,1:-3,2:-2,2:-2])+b[1:-1,1:-1,1:-1]+ np.multiply(ap0[1:-1,1:-1,1:-1],phi[t-1,2:-2,2:-2,2:-2])),ap[1:-1,1:-1,1:-1])
		
		itera=itera+1
		#print(itera)
	print("error analysis: ", itera,np.linalg.norm(rdue))#,
	
x,y= np.meshgrid(mesh.cvxgrid,mesh.cvygrid)
fig,ax=plt.subplots(figsize=(10,5))

#PLOT COLOR MESH
cmin=np.min(phi[:,1:-1,int(yGRID_nos/2),1:-1])
print(cmin)		
cmax=np.max(phi[:,1:-1,int(yGRID_nos/2),1:-1])
print(cmax)
cs=ax.pcolormesh(x,y,phi[0,1:-1,int(yGRID_nos/2),1:-1], cmap=cm.RdPu, vmin=cmin,vmax=cmax)

cbar=fig.colorbar(cs)

ax.set_xlabel("Distnace (m)")
ax.set_ylabel("Distance (m)")
ax.set_title("2D Advection-Diffusion FVM Solver\n with Dirichlet boundary Condition at x=0 and y=0")

#Func TO ANIMATE THE PLOT
def update(t):
	for txt in ax.texts:
		txt.set_visible(False)
	ax.pcolormesh(x,y,np.transpose(phi[t,1:-1,int(yGRID_nos/2),1:-1]), cmap=cm.RdPu, vmin=cmin,vmax=cmax)
	ax.text(0,-1/8,'at time= %f sec\nux=2, uy=2 & gamma=3'%(t*dt), size=10)

sim= animation.FuncAnimation(fig,update, frames=range(0,len(phi[:,0,0])), interval=500, repeat=False)
plt.show()
#+#np.multiply(aw[1:-1,1:-1,1:-1],phi[t,1:-3,2:-2,2:-2])+np.multiply(ae[1:-1,1:-1,1:-1],phi[t,3:-1,2:-2,2:-2])+np.multiply(au[1:-1,1:-1,1:-1],phi[t,2:-2,3:-1,2:-2])+np.multiply(ad[1:-1,1:-1,1:-1], phi[t,2:-2,1:-3,2:-2])+
