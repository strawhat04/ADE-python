import numpy as np

from StructuralDiscre import LinFlux


def FinalMat(m):
	ap,aw,ae,ad,au,ab,af,ap0,b,c,rho,u,v,w= LinFlux(m)

	d=1
	n=0

	x0phi=5
	xnphi=2

	y0phi=3
	ynphi=1

	z0phi=5
	znphi=1

	fx0=3
	fxn=1

	fy0=2
	fyn=2

	fz0=3
	fzn=2

	nx=len(m.xgrid)-1
	ny=len(m.ygrid)-1
	nz=len(m.zgrid)-1


	Sx=np.dot((m.cvygrid[1:]-m.cvygrid[:-1]).reshape(-1,1), (m.cvzgrid[1:]-m.cvzgrid[:-1]).reshape(1,-1))
	Sy=np.dot((m.cvxgrid[1:]-m.cvxgrid[:-1]).reshape(-1,1), (m.cvzgrid[1:]-m.cvzgrid[:-1]).reshape(1,-1))
	Sz=np.dot((m.cvygrid[1:]-m.cvygrid[:-1]).reshape(-1,1), (m.cvxgrid[1:]-m.cvxgrid[:-1]).reshape(1,-1))
# 	print(Sx)   

	if d==1:
		aw[0,:,:]=ae[nx,:,:]=0
		b[0,:,:]=np.multiply(rho[0,1:,1:],u[0,1:,1:],Sx)*x0phi
		b[nx,:,:]=np.multiply(rho[nx,1:,1:],u[nx,1:,1:],Sx)*xnphi

		ad[:,0,:]=au[:,ny,:]=0
		b[:,0,:]=np.multiply(rho[1:,0,1:],v[1:,0,1:],Sy)*y0phi
		b[:,ny,:]=np.multiply(rho[1:,ny,1:],v[1:,ny,1:],Sy)*ynphi

		ab[:,:,0]=af[:,:,nz]=0
		b[:,:,0]=np.multiply(rho[1:,1:,0],w[1:,1:,0],Sz)*z0phi
		b[:,:,nz]=np.multiply(rho[1:,1:,nz],w[1:,1:,nz],Sz)*znphi

	if n==1:
		aw[0,:,:]=ae[nx,:,:]=0
		b[0,:,:]=fx0*Sx
		b[nx,:,:]=fxn*Sx

		ad[:,0,:]=au[:,ny,:]=0
		b[:,0,:]=fy0*Sy
		b[0,ny,:]=fyn*Sy

		ab[:,:,0]=af[:,:,nz]=0
		b[:,:,0]=fz0*Sz
		b[:,:,nz]=fzn*Sz

	return ap,aw,ae,ad,au,ab,af,ap0,b,c



