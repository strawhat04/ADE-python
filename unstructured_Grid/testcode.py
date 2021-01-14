import meshplex
import numpy as np
import pygmsh
from MeshProcessing import MeshTopoogy
from EqDiscretization import linMatrix
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#this function is to give unique ID to every face of triangular elemental mesh
def meshcellFaceID(mesh):
	cells=mesh.cells[1].data
	#inititalizing list to store cell Face
	cellFaceID = []

	#first loop at every cell element
	for i in range(len(cells)):
		#loop over all 3 sides of triangular element
		for j in range(3):
			#keep adding those list of two points making triangular elemental side and avoid repeatation
			if [ cells[i,j-1], cells[i,j]] and [cells[i,j], cells[i,j-1]] not in cellFaceID:
				#print("in if")
				cellFaceID.append([cells[i,j-1], cells[i,j]])

	#return after converting list of face ID into numpy array
	return np.array(cellFaceID)


with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [0, 1],
            [1.0, 1],
            [1, .0],
        ],
        mesh_size=0.1
    )
    pymesh= geom.generate_mesh()


print("\n---------------\n")
print(pymesh.points)

pymesh.cellFaceID=meshcellFaceID(pymesh)

print(dir(pymesh))
# print(pymesh.cells)
# print("\n--------sorting cell:1-------------\n")

mtopo = MeshTopoogy(pymesh)
mtopo.get_neighbourCELL()

# print("neighbour element",mtopo.neighCellID)
# print("centroid",mtopo.centroid)


# print("boundary face", mtopo.BoundaryFace)
lin=linMatrix()

runtime=1.2
dt=0.05

phi=0*np.ones((int(runtime/dt)+1,len(mtopo.centroid)))

#phi[0,6]=phi[0,25]=phi[0,21]=phi[0,1]=2

fluxMat, bMat=lin.get_linMatrix(pymesh,mtopo,dt, np.zeros((len(mtopo.centroid))))
# print("flux mat", fluxMat)
# print("bMat", bMat)

for t in range(len(phi[:,0])-1):
	fluxMat, bMat=lin.get_linMatrix(pymesh,mtopo,dt, phi[t,:])
	phi[t+1,:]=np.transpose(np.linalg.solve(fluxMat,bMat))
# 	print("bmat", bMat)

	
# print(phi[-1,:])
# print(np.max(phi))

# print(phi.shape)
pmax=np.max(phi)
pmin=np.min(phi)


pymesh.field_data=phi[-1,:]
print(pymesh.field_data)
pymesh.write("out.vtk")
fig = plt.figure() 
axis = plt.axes(xlim=(0, 1), ylim=(0,1)) 

plot1=axis.tricontour(mtopo.centroid[:,0],mtopo.centroid[:,1], phi[1,:], levels=16, linewidths=0.05, colors='k')
cntr = axis.tricontourf(mtopo.centroid[:,0],mtopo.centroid[:,1], phi[1,:], levels=16, cmap ='OrRd', vmin=pmin,vmax=pmax)
fig.colorbar(cntr, ax=axis)
axis.plot(mtopo.centroid[:,0],mtopo.centroid[:,1], 'ko', ms=3)


print("Animating")
def animate(t):
    axis.clear()
    plot1=axis.tricontour(mtopo.centroid[:,0],mtopo.centroid[:,1], phi[t,:], levels=16, linewidths=0.05, colors='k')
    cntr = axis.tricontourf(mtopo.centroid[:,0],mtopo.centroid[:,1], phi[t,:], levels=16,cmap ='OrRd')
#     fig.colorbar(cntr, ax=axis)
    axis.set_title("Pure advection, Source at x=0amd vx=1, vy=0\nat t:(%f) seconds"%(dt*t))
vis=FuncAnimation(fig,animate,frames=range(0,len(phi)), interval=200,repeat=False)
vis.save('mtp1.gif', writer='imagemagick')
mesh3=meshplex.MeshTri(pymesh.points, pymesh.cells[1].data)
#mesh3.show()
plt.show()


