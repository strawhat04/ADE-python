import numpy as np

class struct_mesh:

	def __init__(self, dim, ele_nos):
		self.lx=dim[0]
		self.ly=dim[1]
		self.lz=dim[2]

		self.nx=ele_nos[0]
		self.ny=ele_nos[1]
		self.nz=ele_nos[2]


	def generate_mesh(self):
		self.xgrid=np.linspace(0.0,self.lx,num=self.nx, dtype=float)
		self.ygrid=np.linspace(0.0,self.ly,num=self.ny, dtype=float)
		self.zgrid=np.linspace(0.0,self.lz,num=self.nz, dtype=float)

		self.cvxgrid=np.insert((self.xgrid[1:]+self.xgrid[:-1])/2,[0,self.nx-1],[0,self.lx])
		self.cvygrid=np.insert((self.ygrid[1:]+self.ygrid[:-1])/2,[0,self.ny-1],[0,self.ly])
		self.cvzgrid=np.insert((self.zgrid[1:]+self.zgrid[:-1])/2,[0,self.nz-1],[0,self.lz])		
		



	
