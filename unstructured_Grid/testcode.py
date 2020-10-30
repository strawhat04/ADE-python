import meshplex
import numpy as np
import pygmsh
from MeshProcessing import MeshTopoogy
from EqDiscretization import linMatrix
from scipy import linalg
import matplotlib.pyplot as plt

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
        mesh_size=0.4,
    )
    pymesh= geom.generate_mesh()


print("\n---------------\n")
print(pymesh.points)

pymesh.cellFaceID=meshcellFaceID(pymesh)


print(pymesh.cells)
# print("\n--------sorting cell:1-------------\n")

mtopo = MeshTopoogy(pymesh)
mtopo.get_neighbourCELL()

print("neighbour element",mtopo.neighCellID)
print("centroid",mtopo.centroid)


# print("boundary face", mtopo.BoundaryFace)
lin=linMatrix()

runtime=0.2
dt=0.1

phi=-np.ones((int(runtime/dt)+1,len(mtopo.centroid)))

# fluxMat, bMat=lin.get_linMatrix(pymesh,mtopo,dt, np.zeros((len(mtopo.centroid))))
# print("flux mat", fluxMat)
# print("bMat", bMat)

# for t in range(len(phi[:,0])-1):
# 	fluxMat, bMat=lin.get_linMatrix(pymesh,mtopo,dt, phi[t,:])
# 	phi[t+1,:]=np.transpose(np.linalg.solve(fluxMat,bMat))

	
# print(phi)
# print(np.max(phi))

# print(mtopo.centroid.shape, ans.shape)
# fig,ax=plt.subplots(figsize=(10,5))
# cs=ax.pcolormesh(mtopo.centroid, ans)
# plt.show()

#print(fluxMat.shape(),ans.shape
mesh3=meshplex.MeshTri(pymesh.points, pymesh.cells[1].data)
mesh3.show()

