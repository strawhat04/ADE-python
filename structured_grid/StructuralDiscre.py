from Meshing import struct_mesh
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.inf)

def tr(a):
	return np.transpose(a)

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
		
#control volume element length 
def dx(i):
	return mesh.cvxgrid[i+1]-mesh.cvxgrid[i]
#Grid element length 
def dxg(i):
	if i==xGRID_nos or i==0:
		return 0
	else:
		return mesh.xgrid[i]-mesh.xgrid[i-1]

def dy(j):
	return mesh.cvygrid[j+1]-mesh.cvygrid[j]
#Grid element width
def dyg(j):
	if j==yGRID_nos or j==0:
		return 0
	else:
		return mesh.ygrid[j]-mesh.ygrid[j-1]

#control volume element width 
def dz(j):
	return mesh.cvzgrid[j+1]-mesh.cvzgrid[j]
#Grid element width
def dzg(j):
	if j==yGRID_nos or j==0:
		return 0
	else:
		return mesh.zgrid[j]-mesh.zgrid[j-1]


def LinFlux(struct_mesh):

	global mesh
	mesh=struct_mesh

	global xGRID_nos, yGRID_nos, zGrid_nos

	xGRID_nos=len(mesh.xgrid)
	yGRID_nos=len(mesh.ygrid)
	zGrid_nos=len(mesh.zgrid)

	CV_xGRID_nos=len(mesh.cvxgrid)
	CV_yGRID_nos=len(mesh.cvygrid)
	CV_zGRID_nos=len(mesh.cvzgrid)

	#*************************************************************************
	dt=.1
	sp=0
	sc=0

	flow_u=1
	flow_v=1
	flow_w=1
	G0=1
	irho=1

	u=flow_u*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	v=flow_v*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	w=flow_w*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	gamma=G0*np.ones((CV_xGRID_nos,CV_yGRID_nos,CV_zGRID_nos))
	rho=irho*np.ones((CV_xGRID_nos, CV_yGRID_nos,CV_zGRID_nos))

	#*******************************************************************************

	iniFunc=lambda m: [np.zeros((xGRID_nos,yGRID_nos, zGrid_nos)) for _ in range(m)]

	ae,aw,au,ad,af,ab,ap,ap0,b,d,c=iniFunc(11)

	matA=np.zeros((xGRID_nos,yGRID_nos))

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

	return ap,aw,ae,ad,au,ab,af,ap0,b,c,rho,u,v,w
	

