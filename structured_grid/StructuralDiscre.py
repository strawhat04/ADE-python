#This module discretize the ADE equation and return the flux matrix 

from Meshing import struct_mesh
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.inf)

#advection flux coeff
def F(rho, u, faceArea):
	return rho*u*faceArea

#diffusion flux coeff
def D(gamma, dx,faceArea):
	if gamma==0 or dx==0:
		return 0
	else:
		return gamma*faceArea/dx

#Power law advection scheme
def funA(F,D):
	if D==0:
		return 0
	else:
		return max(0,(1-0.1*abs(F/(D)))**5 )
		
#control volume element length  in x direction
def dx(i):
	return mesh.cvxgrid[i+1]-mesh.cvxgrid[i]
#Grid element length in x direction 
def dxg(i):
	if i==xGRID_nos or i==0:
		return 0

	else:
		return mesh.xgrid[i]-mesh.xgrid[i-1]
#control volume element length  in y direction
def dy(j):
	return mesh.cvygrid[j+1]-mesh.cvygrid[j]
#Grid element width in y direction
def dyg(j):
	if j==yGRID_nos or j==0:
		return 0
	else:
		return mesh.ygrid[j]-mesh.ygrid[j-1]

#control volume element width in z direction
def dz(j):
	return mesh.cvzgrid[j+1]-mesh.cvzgrid[j]

#Grid element width om z direction
def dzg(j):
	if j==yGRID_nos or j==0:
		return 0
	else:
		return mesh.zgrid[j]-mesh.zgrid[j-1]


def LinFlux(struct_mesh):

	global mesh
	mesh=struct_mesh

	global xGRID_nos, yGRID_nos, zGrid_nos

	#total number of mesh grid points
	xGRID_nos=len(mesh.xgrid)
	yGRID_nos=len(mesh.ygrid)
	zGrid_nos=len(mesh.zgrid)

	#total number of CV grid points
	CV_xGRID_nos=len(mesh.cvxgrid)
	CV_yGRID_nos=len(mesh.cvygrid)
	CV_zGRID_nos=len(mesh.cvzgrid)

	#*************************************************************************
	#the variables accessed from this star box will be imported from module later
	dt=.1
	sp=0
	sc=0

	flow_u=1
	flow_v=1
	flow_w=1
	G0=1
	irho=1

	#these values to be defined at the CV grid points
	u=flow_u*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	v=flow_v*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	w=flow_w*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	gamma=G0*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	rho=irho*np.ones((CV_xGRID_nos, CV_yGRID_nos,CV_zGRID_nos))

	
	#*******************************************************************************

	#Function to inititate the flux coefficient array 
	iniFunc=lambda m: [np.zeros((xGRID_nos,yGRID_nos, zGrid_nos)) for _ in range(m)]

	#these value are defined at the mesh grid points
	ae,aw,au,ad,af,ab,ap,ap0,b,d,c=iniFunc(11)


	#loop throught every mesh elements to compute the value of flux coefficeints
	for i in range(0,xGRID_nos):
		for j in range(0,yGRID_nos):
			for k in range(0,zGrid_nos):	
				aw[i,j,k]=D(gamma[i,j,k], dxg(i),dy(j)*dz(k))*funA(F(rho[i,j,k],u[i,j,k], dy(j)*dz(k)),D(gamma[i,j,k], dxg(i), dy(j)*dz(k))) + max(F(rho[i,j,k],u[i,j,k], dy(j)*dz(k)),0)
				
				ae[i,j,k]=D(gamma[i+1,j,k], dxg(i+1),dy(j)*dz(k))*funA(F(rho[i+1,j,k],u[i+1,j,k],dy(j)*dz(k)),D(gamma[i+1,j,k], dxg(i+1), dy(j)*dz(k) )) + max(-F(rho[i+1,j,k],u[i+1,j,k],dy(j)*dz(k)),0)	
				
				ad[i,j,k]=D(gamma[i,j,k], dyg(j),dx(i)*dz(k))*funA(F(rho[i,j,k],v[i,j,k], dx(i)*dz(k)),D(gamma[i,j,k], dyg(j), dx(i)*dz(k))) + max(F(rho[i,j,k],v[i,j,k], dx(i)*dz(k)),0)

				au[i,j]=D(gamma[i,j+1,k], dyg(j+1), dx(i)*dz(k))*funA(F(rho[i,j+1,k],v[i,j+1,k],dx(i)*dz(k)),D(gamma[i,j+1,k], dyg(j+1), dx(i)*dz(k))) + max(-F(rho[i,j+1,k],v[i,j+1,k],dx(i)*dz(k)),0)	
				
				ab[i,j,k]=D(gamma[i,j,k], dzg(k),dx(i)*dy(j))*funA(F(rho[i,j,k],w[i,j,k], dx(i)*dy(j)),D(gamma[i,j,k], dzg(k), dx(i)*dy(j))) + max(F(rho[i,j,k],w[i,j,k], dx(i)*dy(j)),0)

				af[i,j,k]=D(gamma[i,j,k+1], dzg(k+1),dx(i)*dy(j))*funA(F(rho[i,j,k+1],w[i,j,k+1],dx(i)*dy(j)),D(gamma[i,j,k+1], dzg(k+1), dx(i)*dy(j))) + max(-F(rho[i,j,k+1],w[i,j,k+1],dx(i)*dy(j)),0)	

				c[i,j,k]= sp*dx(i)*dy(j)
				ap0[i,j,k]=rho[i,j,k]*dx(i)*dy(j)/dt	
				b[i,j,k]=sc*dx(i)*dy(j) # +ap0[i]*phi[t-1,i]

	
	ap=aw+ae+ad+au+ab+af+ap0-c


	#NOW COMPUTING BOUNDARY FLUX COEFFICIENT 
	
	#last element index of each array
	nx=xGRID_nos-1
	ny=yGRID_nos-1
	nz=zGrid_nos-1

	#array storing face area of each control volume element in face normal to x axis
	Sx=np.dot((mesh.cvygrid[1:]-mesh.cvygrid[:-1]).reshape(-1,1), (mesh.cvzgrid[1:]-mesh.cvzgrid[:-1]).reshape(1,-1))
	#array storing face area of each control volume element in face normal to y axis
	Sy=np.dot((mesh.cvxgrid[1:]-mesh.cvxgrid[:-1]).reshape(-1,1), (mesh.cvzgrid[1:]-mesh.cvzgrid[:-1]).reshape(1,-1))
	#array storing face area of each control volume element in face normal to z axis
	Sz=np.dot((mesh.cvygrid[1:]-mesh.cvygrid[:-1]).reshape(-1,1), (mesh.cvxgrid[1:]-mesh.cvxgrid[:-1]).reshape(1,-1))

	print("Imput Boundary Conditions as: D/N x0 y0 z0 xn yn zn")
	bc,x0,y0,z0,xn,yn,zn=input().split(' ')

	#Computing flux coefficient for dirichlet condition
	if bc=='D' or bc=='d':
		aw[0,:,:]=ae[nx,:,:]=0
		b[0,:,:]=float(x0)*np.multiply(rho[0,1:,1:],u[0,1:,1:],Sx)
		b[nx,:,:]=np.multiply(rho[0,1:,1:],u[0,1:,1:],Sx)*float(xn)

		ad[:,0,:]=au[:,ny,:]=0
		b[:,0,:]=np.multiply(rho[1:,0,1:],v[1:,0,1:],Sy)*float(y0)
		b[:,ny,:]=np.multiply(rho[1:,ny,1:],v[1:,ny,1:],Sy)*float(yn)

		ab[:,:,0]=af[:,:,nz]=0
		b[:,:,0]=np.multiply(rho[1:,1:,0],w[1:,1:,0],Sz)*float(z0)
		b[:,:,nz]=np.multiply(rho[1:,1:,nz],w[1:,1:,nz],Sz)*float(zn)

	#Computing flux coefficient for Neumann condition
	if bc=='N' or bc=='n':
		aw[0,:,:]=ae[nx,:,:]=0
		b[0,:,:]=float(x0)*Sx
		b[nx,:,:]=float(xn)*Sx

		ad[:,0,:]=au[:,ny,:]=0
		b[:,0,:]=float(y0)*Sy
		b[0,ny,:]=float(yn)*Sy

		ab[:,:,0]=af[:,:,nz]=0
		b[:,:,0]=float(z0)*Sz
		b[:,:,nz]=float(zn)*Sz

	return ap,aw,ae,ad,au,ab,af,ap0,b,c
	

