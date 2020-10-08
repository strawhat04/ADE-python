import meshplex
import numpy as np
import pygmsh
from MeshProcessing import MeshTopoogy
from EqDiscretization import linMatrix

#this function is to give unique ID to every face of triangular elemental mesh
def meshcellFace(mesh):
	cells=mesh.cells[1].data
	#inititalizing list to store cell Face
	cellFace = []

	#first loop at every cell element
	for i in range(len(cells)):
		#loop over all 3 sides of triangular element
		for j in range(3):
			#keep adding those list of two points making triangular elemental side and avoid repeatation
			if [ cells[i,j-1], cells[i,j]] and [cells[i,j], cells[i,j-1]] not in cellFace:
				#print("in if")
				cellFace.append([cells[i,j-1], cells[i,j]])

	#return after converting list of face ID into numpy array
	return np.array(cellFace)



with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [0, 1],
            [1.0, 1],
            [1, .0],
        ],
        mesh_size=0.5,
    )
    mesh1= geom.generate_mesh()
print("start")
print(mesh1.points[:,0:2])

print("\n---------------\n")
print(mesh1.cells)

# c=np.sort(mesh1.cells[1].data, axis=1)
# c=c[c[:,0].argsort(kind='mergesort')]
# print(c)

cellFace=meshcellFace(mesh1)
print("cellFace:",cellFace)
#print(mesh1.cells)
# print(mesh1.cells)
# print("\n--------sorting cell:1-------------\n")

mesh = MeshTopoogy(mesh1,cellFace)
mesh.get_neighbourCELL()
mesh.get_meshDATA()
print(mesh.neighCellID)
print(mesh.commonFaceID)
mesh3=meshplex.MeshTri(mesh1.points, mesh1.cells[1].data)
print("centroid", mesh.centroid[1,0])
lin=linMatrix()
lin.get_linMatrix(mesh1,mesh,cellFace)
#print(mesh.cells)

#mesh3.show()