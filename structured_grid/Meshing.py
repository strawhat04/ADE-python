import numpy as np

#Class for creating 3D structured mess
class struct_mesh:

	def __init__(self, dim, grid_nos):
		#instanciate the length, height and width
		self.lx=dim[0]
		self.ly=dim[1]
		self.lz=dim[2]

		#instanciate the number of mesh elements in respective direction
		self.nx=grid_nos[0]
		self.ny=grid_nos[1]
		self.nz=grid_nos[2]


	#we gonna store each axes dimension grid topology in respective grid array
	def generate_mesh(self):
		#place n grid points between 0 and length in respecive grid array to
		self.xgrid=np.linspace(0.0,self.lx,num=self.nx, dtype=float)
		self.ygrid=np.linspace(0.0,self.ly,num=self.ny, dtype=float)
		self.zgrid=np.linspace(0.0,self.lz,num=self.nz, dtype=float)

		#to create cv grid, first place points in middle of two adjecant mesh grid points and include 0 and length point at boundary
		self.cvxgrid=np.insert((self.xgrid[1:]+self.xgrid[:-1])/2,[0,self.nx-1],[0,self.lx])
		self.cvygrid=np.insert((self.ygrid[1:]+self.ygrid[:-1])/2,[0,self.ny-1],[0,self.ly])
		self.cvzgrid=np.insert((self.zgrid[1:]+self.zgrid[:-1])/2,[0,self.nz-1],[0,self.lz])
