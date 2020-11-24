#This module discretize the ADE equation and return the flux matrix

#import module

#from Meshing import struct_mesh
from parameters import *

#import numpy modules
import numpy as np
import matplotlib.pyplot as plt


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
def dz(k):
	return mesh.cvzgrid[k+1]-mesh.cvzgrid[k]

#Grid element width om z direction
def dzg(k):
	if k==zGRID_nos or k==0:
		return 0
	else:
		return mesh.zgrid[k]-mesh.zgrid[k-1]


#this function will discretize the ADE equation by given input parameters
#and return the coefficient flux of system of linear equation 
def LinFlux(struct_mesh):

	#define mesh grid object globally so all other function easily access without passion again as arguments
	global mesh
	mesh=struct_mesh

	#globally define totaly number of grid points in each direction
	global xGRID_nos, yGRID_nos, zGRID_nos

	#total number of mesh grid points
	xGRID_nos=len(mesh.xgrid)
	yGRID_nos=len(mesh.ygrid)
	zGRID_nos=len(mesh.zgrid)

	#total number of CV grid points
	CV_xGRID_nos=len(mesh.cvxgrid)
	CV_yGRID_nos=len(mesh.cvygrid)
	CV_zGRID_nos=len(mesh.cvzgrid)

	#*************************************************************************
	
	#these values to be imported from the parameter module
	u,v,w=velocity_model(CV_zGRID_nos,CV_yGRID_nos,CV_xGRID_nos)
	
	gamma=diffusion_model(CV_zGRID_nos,CV_yGRID_nos,CV_xGRID_nos)
	
	rho=density_model(CV_zGRID_nos,CV_yGRID_nos,CV_xGRID_nos)

	sp,sc = source_model(zGRID_nos,yGRID_nos,xGRID_nos)


	#*******************************************************************************

	#Function to inititate the flux coefficient array
	iniFunc=lambda m: [np.zeros((zGRID_nos,yGRID_nos,xGRID_nos)) for _ in range(m)]

	#these value are defined at the mesh grid points
	ae,aw,au,ad,af,ab,ap,ap0,b,d,c=iniFunc(11)


	#loop throught every mesh elements to compute the value of flux coefficeints
	#Refer Readme to know more 
	for i in range(0,xGRID_nos):
		for j in range(0,yGRID_nos):
			for k in range(0,zGRID_nos):
				#aw=
				aw[k,j,i]=D(gamma[k,j,i], dxg(i),dy(j)*dz(k))*funA(F(rho[k,j,i],u[k,j,i], dy(j)*dz(k)),D(gamma[k,j,i], dxg(i), dy(j)*dz(k))) + max(F(rho[k,j,i],u[k,j,i], dy(j)*dz(k)),0)

				ae[k,j,i]=D(gamma[k,j,i+1], dxg(i+1),dy(j)*dz(k))*funA(F(rho[k,j,i+1],u[k,j,i+1],dy(j)*dz(k)),D(gamma[k,j,i+1], dxg(i+1), dy(j)*dz(k) )) + max(-F(rho[k,j,i+1],u[k,j,i+1],dy(j)*dz(k)),0)

				ad[k,j,i]=D(gamma[k,j,i], dyg(j),dx(i)*dz(k))*funA(F(rho[k,j,i],v[k,j,i], dx(i)*dz(k)),D(gamma[k,j,i], dyg(j), dx(i)*dz(k))) + max(F(rho[k,j,i],v[k,j,i], dx(i)*dz(k)),0)

				au[k,j,i]=D(gamma[k,j+1,i], dyg(j+1), dx(i)*dz(k))*funA(F(rho[k,j+1,i],v[k,j+1,i],dx(i)*dz(k)),D(gamma[k,j+1,i], dyg(j+1), dx(i)*dz(k))) + max(-F(rho[k,j+1,i],v[k,j+1,i],dx(i)*dz(k)),0)

				ab[k,j,i]=D(gamma[k,j,i], dzg(k),dx(i)*dy(j))*funA(F(rho[k,j,i],w[k,j,i], dx(i)*dy(j)),D(gamma[k,j,i], dzg(k), dx(i)*dy(j))) + max(F(rho[k,j,i],w[k,j,i], dx(i)*dy(j)),0)

				af[k,j,i]=D(gamma[k+1,j,i], dzg(k+1),dx(i)*dy(j))*funA(F(rho[k+1,j,i],w[k+1,j,i],dx(i)*dy(j)),D(gamma[k+1,j,i], dzg(k+1), dx(i)*dy(j))) + max(-F(rho[k+1,j,i],w[k+1,j,i],dx(i)*dy(j)),0)

				c[k,j,i]= sp[k,j,i]*dx(i)*dy(j)*dz(k)
				ap0[k,j,i]=rho[k,j,i]*dx(i)*dy(j)*dz(k)/dt
				b[k,j,i]=sc[k,j,i]*dx(i)*dy(j)*dz(k) # +ap0[i]*phi[t-1,i]

	ap=aw+ae+ad+au+ab+af+ap0-c

	return ap,aw,ae,ad,au,ab,af,ap0,b

# class conditions():
# 	def  __init__(self,struct_mesh):
#         ap,aw,ae,ad,au,ab,af,ap0,b=LinFlux(struct_mesh)
#         self.ap=ap
#         self.aw=aw
#         self.ae=ae
#         self.ad=ad
#         self.au=au
#         self.ab=ab
#         self.af=af
#         self.ap0=ap0
#         self.b=b
